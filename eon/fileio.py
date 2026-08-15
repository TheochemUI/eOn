
'''
Con(figuration) i/o library.

``.con`` frames are read and written with the **readcon** package
(:class:`readcon.ConFrame`). The in-memory working type is
:class:`eon.structure.Structure` (alias ``Atoms``).
'''
import configparser
import contextlib
#from io import BytesIO as StringIO
from io import StringIO
import logging
import numpy
import os

import pickle as pickle
import readcon
import stat
import tempfile

from eon.geometry.cell import box_to_length_angle, length_angle_to_box
from eon.structure import Structure

logger = logging.getLogger('io')

def prng_state_path(config):
    '''Path of the persisted numpy PRNG pickle.

    ConfigClass.init restores from this file under path_root, so every
    writer and reset must use the same location rather than the process CWD.
    '''
    return os.path.join(config.path_root, 'prng.pkl')


def save_prng_state(path):
    state = numpy.random.get_state()
    with open(path, 'wb') as fh:
        pickle.dump(state, fh, pickle.HIGHEST_PROTOCOL)

def get_prng_state(path):
    with open(path, 'rb') as fh:
        state = pickle.load(fh)
    numpy.random.set_state(state)

# Re-export cell helpers (used by callers / POSCAR path)
__all_cell__ = ("length_angle_to_box", "box_to_length_angle")


def _process_umask():
    mask = os.umask(0)
    os.umask(mask)
    return mask


# What open() and os.mkdir() give a new file or directory. mkstemp and
# mkdtemp ignore the umask and use 0600 / 0700, so anything staged through
# them is widened back to these before being moved into place.
_UMASK = _process_umask()
DEFAULT_FILE_MODE = 0o666 & ~_UMASK
DEFAULT_DIR_MODE = 0o777 & ~_UMASK


@contextlib.contextmanager
def atomic_write(path, mode='w'):
    '''
    Rewrite a file in one step, for callers that truncate and rewrite whole.
        path: destination file
        mode: 'w' or 'wb'
    The body writes to a temporary file in the destination's own directory,
    which is renamed over the destination on a clean exit. A reader sees
    either the old contents or the new ones, and a write that fails partway
    (a full disk, a killed process) leaves the destination as it was.
    mkstemp opens at 0600, so the destination's own mode carries over, or
    the umask default for a file that does not exist yet.
    '''
    directory = os.path.dirname(os.path.abspath(path))
    fd, temp_path = tempfile.mkstemp(dir=directory,
                                     prefix='.' + os.path.basename(path) + '.',
                                     suffix='.tmp')
    try:
        with os.fdopen(fd, mode) as f:
            yield f
        try:
            os.chmod(temp_path, stat.S_IMODE(os.stat(path).st_mode))
        except OSError:
            os.chmod(temp_path, DEFAULT_FILE_MODE)
        os.replace(temp_path, path)
    except BaseException:
        try:
            os.unlink(temp_path)
        except OSError:
            pass
        raise


def _maybe_open(target, mode, attr):
    '''
    Context manager over a filename or an already open file-like object.
        target: filename or file-like object
        mode:   mode passed to open() when target is a filename
        attr:   attribute marking target as file-like ('readline' or 'write')
    A handle opened here is closed on exit; one passed in by the caller is
    left open.
    '''
    if hasattr(target, attr):
        return contextlib.nullcontext(target)
    return open(target, mode)


def _frame_to_atoms(frame):
    """Convert a readcon.ConFrame to a Structure (canonical bridge)."""
    return Structure.from_conframe(frame)


def loadcons(filename):
    frames = readcon.read_con(filename)
    return [_frame_to_atoms(f) for f in frames]


def loadxyz(filename):
    '''
    Load the first frame of an XYZ/PDB/GRO file via readcon chemfiles ingress.

    Requires a chemfiles-linked readcon (``pip install readcon-chemfiles``).
    The .con path stays on readcon-core; this is the foreign-format edge.
    '''
    if not getattr(readcon, "has_chemfiles_support", lambda: False)():
        raise RuntimeError(
            "loadxyz needs a chemfiles-linked readcon "
            "(pip install readcon-chemfiles)"
        )
    return _frame_to_atoms(readcon.read_chemfiles_first(filename))


def _tokens_are_ints(tokens):
    """True when every token parses as an int (VASP 4 counts line)."""
    if not tokens:
        return False
    for tok in tokens:
        try:
            int(tok)
        except ValueError:
            return False
    return True


def loadposcars(filename):
    '''Load every POSCAR frame from filename (VASP 4 or VASP 5).'''
    p = []
    with open(filename, 'r') as filein:
        while True:
            try:
                p.append(loadposcar(filein))
            except (ValueError, IndexError):
                # End of file: the next header line is empty.
                return p


def loadcon(filein, reset = True):
    '''
    Load a con file via readcon; return a Structure.
        filein: may be either a filename or a file-like object
    '''
    if hasattr(filein, 'readline'):
        content = filein.read()
        if reset:
            filein.seek(0)
        frames = readcon.read_con_string(content)
        if not frames:
            raise IOError("No frames found in con data")
        return _frame_to_atoms(frames[0])
    return _frame_to_atoms(readcon.read_first_frame(filein))

def _as_structure(p):
    """Copy a Structure-like object into a Structure.

    Anything carrying the five geometry attributes plus a length passes; an
    ``atom_ids`` attribute carries over, and one without gets the sequential
    numbering :class:`Structure` gives a fresh configuration.
    """
    n = len(p)
    s = Structure(n)
    s.box = numpy.asarray(p.box, dtype=float).reshape(3, 3).copy()
    s.r = numpy.asarray(p.r, dtype=float).reshape(n, 3).copy()
    s.free = numpy.asarray(p.free, dtype=float)
    s.names = list(p.names)
    s.mass = numpy.asarray(p.mass, dtype=float).reshape(n).copy()
    ids = getattr(p, 'atom_ids', None)
    if ids is not None:
        s.atom_ids = numpy.asarray(ids, dtype=numpy.uint64).reshape(-1).copy()
    return s


def _atoms_to_frame(p):
    """Convert a Structure (or Structure-like) to a readcon.ConFrame."""
    if not isinstance(p, Structure):
        p = _as_structure(p)
    return p.to_conframe()


def _path_is_compressed_con(path):
    name = os.path.basename(path).lower()
    return name.endswith(".gz") or name.endswith(".zst")


def _file_ends_with_newline(path):
    """True for an empty file too: nothing needs separating from the first frame."""
    try:
        with open(path, "rb") as fh:
            fh.seek(0, os.SEEK_END)
            if fh.tell() == 0:
                return True
            fh.seek(-1, os.SEEK_END)
            return fh.read(1) == b"\n"
    except OSError:
        return True


def savecon(fileout, p, w = 'w'):
    '''
    Save a con file via readcon.
        fileout: can be either a file name or a file-like object
        p:       Structure (or Structure-like) to save
        w:       write/append flag

    Append of an uncompressed file concatenates one serialized frame. A
    gzip or zstd target cannot be extended in place, so those still
    read the existing frames and rewrite the file.
    '''
    frame = _atoms_to_frame(p)
    if hasattr(fileout, 'write'):
        text = readcon.write_con_string([frame])
        fileout.write(text)
    elif w == 'a' and os.path.exists(fileout) and os.path.getsize(fileout) > 0:
        if _path_is_compressed_con(fileout):
            existing = readcon.read_con(fileout)
            existing.append(frame)
            readcon.write_con(fileout, existing)
        else:
            text = readcon.write_con_string([frame])
            with open(fileout, 'a') as fh:
                if not _file_ends_with_newline(fileout):
                    fh.write('\n')
                fh.write(text)
    else:
        readcon.write_con(fileout, [frame])


def load_mode(modefilein):
    '''
    Reads a mode.dat file into an N by 3 numpy array
        modefilein: may be either a file-like object of a filename
    '''
    with _maybe_open(modefilein, 'r', 'readline') as f:
        if len(f.readline().split()) == 3:
            f.seek(0)
        lines = f.readlines()
    mode = []
    for line in lines:
        l = line.split()
        if not l:
            continue
        if len(l) < 3:
            raise IOError("Malformed mode.dat line, expected three columns: %r" % line)
        for j in range(3):
            mode.append(float(l[j]))
    return numpy.array(mode).reshape(len(mode)//3, 3)

def save_mode(modefileout, displace_vector, free=None):
    '''
    Saves an Nx3 numpy array into a mode.dat file.
        modefileout:     may be either a filename or file-like object
        displace_vector: the mode (Nx3 numpy array)
        free:            optional (N,) or (N,3) mask; a 0 axis is written as 0
    17 significant digits round trip a double exactly, and match what the
    client's printf writers emit for the same value.
    '''
    vec = numpy.asarray(displace_vector, dtype=float)
    if free is not None:
        mask = numpy.asarray(free, dtype=float)
        if mask.ndim == 1:
            mask = numpy.repeat(mask.reshape(-1, 1), 3, axis=1)
        vec = vec * (mask > 0.5)
    with _maybe_open(modefileout, 'w', 'write') as f:
        for i in range(len(vec)):
            f.write("%.17g %.17g %.17g\n" % (vec[i][0], vec[i][1], vec[i][2]))


def save_results_dat(fileout, results):
    '''
    Saves a results.dat file from a dictionary
        fileout: may be either a filename or a file-like object
        results: dictionary of values, written one "<value> <key>" line each
    '''
    with _maybe_open(fileout, 'w', 'write') as f:
        for key in results:
            f.write("%s %s\n" % (results[key], key))

def modify_config(config_path, changes):
    parser = configparser.ConfigParser()
    parser.read(config_path)
    for change in changes:
        parser.set(*change)
    config_str_io = StringIO()
    parser.write(config_str_io)
    config_str_io.seek(0)
    return config_str_io

def _coerce_results_value(token):
    '''Parse one results.dat value token.

    Integers stay int. Everything float() accepts (1.0, 1e-5, nan, inf)
    becomes float. Dotted identifiers and other non-numerics stay str so
    the key is present for the caller.
    '''
    body = token.lstrip('+-')
    if body.isdigit():
        return int(token)
    try:
        return float(token)
    except ValueError:
        return token


def parse_results(filein):
    '''
    Reads a results.dat file and gives a dictionary of the values contained therein
    '''
    results = {}
    with _maybe_open(filein, 'r', 'readline') as f:
        f.seek(0)
        for line in f:
            line = line.split()
            if len(line) < 2:
                continue
            results[line[1]] = _coerce_results_value(line[0])

    return results

def loadposcar(filein):
    '''
    Load a VASP POSCAR (or one frame of a multi-frame movie).

    POSCAR is the one structure format that does not go through readcon:
    the kdb tool emits VASP 4, and movie.py writes VASP 5 for external
    viewers. The reader accepts both. After the cell, a line of integer
    counts is VASP 4 (species taken from the comment); a line of species
    names followed by the counts is VASP 5.

    filein: filename or file-like object
    '''
    with _maybe_open(filein, 'r', 'readline') as f:
        # Line 1: comment, often the species names.
        AtomTypes = f.readline().split()
        # Line 2: scaling of coordinates
        scale = float(f.readline())
        # Lines 3-5: the box
        box = numpy.zeros((3, 3))
        for i in range(3):
            line = f.readline().split()
            box[i] = numpy.array([float(line[0]), float(line[1]), float(line[2])]) * scale
        # Line 6 is either VASP 4 counts or VASP 5 species names.
        line = f.readline().split()
        if _tokens_are_ints(line):
            NumAtomsPerType = [int(tok) for tok in line]
        else:
            if line:
                AtomTypes = line
            counts_line = f.readline().split()
            NumAtomsPerType = [int(tok) for tok in counts_line]
        # Now have enough info to make the Structure object.
        num_atoms = sum(NumAtomsPerType)
        p = Structure(num_atoms)
        # Fill in the box.
        p.box = box
        # Selective Dynamics (optional) or Cartesian/Direct.
        sel = f.readline()[0]
        selective_flag = (sel == 's' or sel == 'S')
        if not selective_flag:
            car = sel
        else:
            car = f.readline()[0]
        direct_flag = not (car == 'c' or car == 'C' or car == 'k' or car == 'K')
        atom_index = 0
        for i in range(len(NumAtomsPerType)):
            for j in range(NumAtomsPerType[i]):
                p.names[atom_index] = AtomTypes[i]
                line = f.readline().split()
                if(selective_flag):
                    assert len(line) >= 6
                else:
                    assert len(line) >= 3
                pos = line[0:3]
                if selective_flag:
                    sel = line[3:6]
                    flags = []
                    for flag in sel:
                        flags.append(1.0 if flag in ('T', 't') else 0.0)
                    while len(flags) < 3:
                        flags.append(flags[0] if flags else 1.0)
                    p.free[atom_index] = flags
                p.r[atom_index] = numpy.array([float(q) for q in pos])
                if direct_flag:
                    p.r[atom_index] = numpy.dot(p.r[atom_index], p.box)
                else:
                    p.r[atom_index] *= scale
                atom_index += 1
        return p


def saveposcar(fileout, p, w='w', direct = False):
    '''
    Save a VASP 5 POSCAR (species on the comment line and again before
    the counts). movie.py uses this for visualization; loadposcar reads
    this layout and the VASP 4 layout kdb produces.
        fileout: name to save it under
        p:       Structure (or Structure-like) to save
        w:       write/append flag
    '''
    with _maybe_open(fileout, w, 'write') as poscar:
        atom_types = []
        num_each_type = {}
        rows_each_type = {}
        for i, name in enumerate(p.names):
            if not name in atom_types:
                atom_types.append(name)
                num_each_type[name] = 1
                rows_each_type[name] = [i]
            else:
                num_each_type[name] += 1
                rows_each_type[name].append(i)
        # The header names the types once and then gives one count per type,
        # so a reader walks the coordinate block type by type. A configuration
        # whose species interleave writes its atoms in a different order than
        # it holds them.
        order = [i for name in atom_types for i in rows_each_type[name]]
        poscar.write(" ".join(atom_types)+'\n') #Line 1: Atom type
        poscar.write("1.0\n") #Line 2: scaling
        for i in range(3):
            poscar.write(" ".join(['%20.14f' % s for s in p.box[i]])+'\n')  #lines 3-5: box
        poscar.write(" ".join(atom_types)+'\n') #Line 6: Atom type
        poscar.write(" ".join(['%s' % num_each_type[key] for key in atom_types])+'\n')
        poscar.write('Selective Dynamics\n') #line 7: selective dynamics
        if direct:
            poscar.write('Direct\n')  #line 8 cartesian coordinates
            ibox = numpy.linalg.inv(numpy.array(p.box))
            positions = numpy.dot(p.r, ibox)
        else:
            poscar.write('Cartesian\n') #line 8 cartesian coordinates
            positions = p.r
        free = numpy.asarray(p.free, dtype=float)
        for i in order:
                posline = " ".join(['%20.14f' % s for s in positions[i]]) + " "
                if free.ndim == 1:
                    axis_free = (free[i] > 0.5, free[i] > 0.5, free[i] > 0.5)
                else:
                    axis_free = (free[i, 0] > 0.5, free[i, 1] > 0.5, free[i, 2] > 0.5)
                for flag in axis_free:
                    posline += '   T' if flag else '   F'
                poscar.write(posline+'\n')


from configparser import ConfigParser as SCP
class ini(SCP):

    def __init__(self, filenames):
        self.loaded = False
        self.filenames = filenames
        SCP.__init__(self)

    def read(self):
        self.loaded = True
        SCP.read(self, self.filenames)

    def get(self, section, option, default="ini_no_default", **kwargs):
#    def get(self, section, option, default="ini_no_default"):
        if not self.loaded:
            self.read()
        try:
            SCP.read(self, self.filenames)
            #value = SCP.get(self, section, option, raw=True, **kwargs)
            value = SCP.get(self, section, option, **kwargs)
            #value = SCP.get(self, section, option, raw=True)
            #value = SCP.get(self, section, option)
        except:
            if default == "ini_no_default":
                raise NameError("Section or option missing, no default specified")
            return default
        try:
            return int(value)
        except ValueError:
            pass
        try:
            return float(value)
        except ValueError:
            pass
        if value.lower() == 'true':
            return True
        if value.lower() == 'false':
            return False
        return value

    def getint(self, *args):
        raise NotImplementedError("Use the get function with this ConfigParser wrapper")
    def getfloat(self, *args):
        raise NotImplementedError("Use the get function with this ConfigParser wrapper")
    def getboolean(self, *args):
        raise NotImplementedError("Use the get function with this ConfigParser wrapper")

    def set(self, section, option, value):
        if section not in self.sections():
            self.add_section(section)
        SCP.set(self, section, option, str(value))
        if type(self.filenames) == str:
            name = self.filenames
        else:
            name = self.filenames[-1]
        with atomic_write(name) as configfile:
            self.write(configfile)


class Dynamics:
    """ The Dynamics class handles I/O for the dynamics.txt file of an aKMC simulation. """

    def __init__(self, filename):
        self.filename = filename
        if not os.path.exists(filename):
            f = open(self.filename, 'w')
            header = "%12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s\n" % ('step-number', 'reactant-id', 'process-id', 'product-id', 'step-time', 'total-time', 'barrier', 'rate', 'energy')
            f.write(header)
            f.write("-" * len(header))
            f.write("\n")
            f.close()
            self.next_step = 0

        # read last lines of the file to determine iteration nr
        else:
            f = open(self.filename,'r')
            f.seek(0,2)	#seek to EOF
            fsize = f.tell()
            # seek 1024 bytes back (or to beginning of file if fsize < 1024 )
            # last line must be contained in this block
            f.seek( max( fsize - 1024 , 0 ) , 0)
            lines = f.readlines()
            self.next_step = int ( lines[-1].split()[0] ) + 1 # determine iteration nr of next step

    def append(self, reactant_id, process_id, product_id, step_time, total_time, barrier, rate, energy):
        f = open(self.filename, 'a')
        f.write("%12d  %12d  %12d  %12d  %12e  %12e  %12f  %12e  %12f\n" % (self.next_step, reactant_id, process_id, product_id, step_time, total_time, barrier, rate, energy))
        f.close()
        self.next_step += 1

    def append_sb(self, reactant_id, process_id, product_id, step_time, total_time, basin_id, rate, energy):
        f = open(self.filename, 'a')
        f.write("%12d  %12d  %12d  %12d  %12e  %12e  %12d  %12e  %12f\n" % (self.next_step, reactant_id, process_id, product_id, step_time, total_time, basin_id, rate, energy))
        f.close()
        self.next_step += 1

    def get(self):
        f = open(self.filename, 'r')
        lines = f.readlines()[2:]
        f.close()
        data = []
        for line in lines:
            split = line.split()
            data.append({"reactant":    int(split[1]),
                         "process":     int(split[2]),
                         "product":     int(split[3]),
                         "steptime":    float(split[4]),
                         "totaltime":   float(split[5]),
                         "barrier":     float(split[6]),
                         "prefactor":   float(split[7])})
        return data

def load_potfiles(pot_dir):
    ret = {}
    if os.path.isdir(pot_dir):
        for name in os.listdir(pot_dir):
            if os.path.isdir(name):
                continue
            a = open(os.path.join(pot_dir, name), 'r')
            b = StringIO("".join(a.readlines()))
            c = os.stat(os.path.join(pot_dir, name)).st_mode
            ret[name] = (b,c)
    return ret

class TableException(Exception):
    pass

class Table:
    """
    A class that provides a nice io abstraction for table like data.  The data
    is saved in a pretty printed format. Also provides nice data retrival
    methods.

    >>> t = Table("sample.tbl", ['id', 'name', 'age' ])
    >>> t.eagerwrite = False
    >>> t.add_row({'id':0,'name':"Sam","age":24})
    >>> t.add_row({'id':1,'name':"David","age":50})
    >>> t.add_row({'id':2,'name':"Anna","age":21})
    >>> t #doctest: +NORMALIZE_WHITESPACE
        id name  age
        -- ----- ---
        0  Sam   24
        1  David 50
        2  Anna  21

    Rows can be accessed directly:
    >>> t.rows[1] #doctest: +SKIP
        {'age': 50, 'id': 1, 'name': 'David'}
    >>> t.max_value('age') #doctest: +NORMALIZE_WHITESPACE
        50
    >>> t.min_row('age') #doctest: +NORMALIZE_WHITESPACE +SKIP
        {'age': 21, 'id': 2, 'name': 'Anna'}
    >>> sorted(t.min_row('id').items()) #doctest: +NORMALIZE_WHITESPACE
        [('age', 24), ('id', 0), ('name', 'Sam')]
    >>> len(t) #doctest: +NORMALIZE_WHITESPACE
        3
    >>> sum(t.getcolumn('age')) #doctest: +NORMALIZE_WHITESPACE
        95
    >>> t.write() #doctest: +SKIP

    The table can be loaded from disk without specifying columns. This is
    slightly unsafe because the columns can't be checked, but it could cut down
    on the verbosity in some places.
    >>> t2 = Table("sample.tbl") #doctest: +SKIP
"""

    #XXX: This is the number of digits that a floating point number gets
    #     serialized with. Should it be some sort of config option?
    #     Or is there just a good default?

    def __init__(self, filename, columns=None, overwrite=False):
        self.filename = filename
        self.columns = columns
        self.rows = []
        self.columntypes = {}
        self.columnwidths = {}
        self.initialized = False
        self.overwrite = overwrite

        self.floatprecision = 6
        self.eagerwrite = True

    def init(self):
        """Checks to see if self.filename exists. If it does self.rows
           will be initialized from disk."""
        self.initialized = True
        if os.path.isfile(self.filename) and not self.overwrite:
            self.read(self.filename)
        else:
            if self.columns is None:
                raise TableException("Columns are not optional for new tables")

            for c in self.columns:
                self.columnwidths[c] = len(c)

    def read(self, filename):
        self.eagerwrite = False
        f = open(self.filename, "r")
        filecolumns = f.readline().split()
        if self.columns != None:
            if filecolumns != self.columns:
                raise TableException("Column name mismatch: %s" % filename)
        else:
            self.columns = filecolumns

        for c in self.columns:
            self.columnwidths[c] = len(c)

        # skip comment line
        f.readline()

        for line in f:
            fields = line.split()
            row = {}
            coli = 0
            for field in fields:
                try:
                    field = int(field)
                except ValueError:
                    try:
                        field = float(field)
                    except ValueError:
                        field = field.strip()
                        pass
                row[self.columns[coli]] = field
                coli += 1
            self.add_row(row)
        f.close()
        self.eagerwrite = True

    def __repr__(self):
        if not self.initialized:
            self.init()
        f = StringIO()
        self.writefilehandle(f)
        return f.getvalue()

    def __len__(self):
        if not self.initialized:
            self.init()
        return len(self.rows)

    def __iter__(self):
        for row in self.rows:
            yield row

    def write(self):
        if not self.initialized:
            self.init()
        #print("into table write: ",self.filename)
        with atomic_write(self.filename) as f:
            self.writefilehandle(f)

    def writefilehandle(self, filehandle):
        f = filehandle
        line = ' '.join([ "%-*s"%(self.columnwidths[c], c) for c in self.columns ])
        f.write(line+"\n")

        line = ''
        for c in self.columns:
            line += '-'*self.columnwidths[c]+' '
        f.write(line+'\n')

        for row in self.rows:
            line = ""
            for c in self.columns:
                if self.columntypes[c] == float:
                    line += "%#-*.*G " % (self.columnwidths[c],self.floatprecision,
                                         row[c])
                else:
                    line += "%-*s " % (self.columnwidths[c],row[c])
            f.write(line+"\n")

    def add_row(self, row):
        if not self.initialized:
            self.init()
        mismatched_columns = set(self.columns).symmetric_difference(set(row.keys()))
        if len(mismatched_columns) != 0:
            raise TableException("Mismatched columns %s" % str(mismatched_columns))

        if len(self.rows) == 0:
            for c in row:
                self.columntypes[c] = type(row[c])
        else:
            for c in row:
                if type(row[c]) != self.columntypes[c]:
                    raise TableException("Type mismatch for column %s" % c)

        for c in row:
            if self.columntypes[c] == float:
                self.columnwidths[c] = max(self.columnwidths[c], self.floatprecision+5)
            else:
                self.columnwidths[c] = max(self.columnwidths[c], len(str(row[c])))

        self.rows.append(row)
        if self.eagerwrite:
            self.write()

    def delete_row(self, column, value):
        if not self.initialized:
            self.init()
        rows_to_delete = []
        for row in self.rows:
            if row[column] == value:
                rows_to_delete.append(row)
#        map(self.rows.remove, rows_to_delete)
        for row in rows_to_delete:
            self.rows.remove(row)
        if self.eagerwrite:
            self.write()
        return len(rows_to_delete)

    def delete_row_func(self, column, func):
        if not self.initialized:
            self.init()

        rows_to_delete = []
        for row in self.rows:
            if func(row[column]):
                rows_to_delete.append(row)
#        map(self.rows.remove, rows_to_delete)
        for row in rows_to_delete:
            self.rows.remove(row)
        if self.eagerwrite:
            self.write()
        return len(rows_to_delete)

    def find_value(self, column, func):
        if not self.initialized:
            self.init()
        value = None
        for row in self.rows:
            if value is None:
                value = row[column]
                continue
            value = func(value, row[column])
        return value

    def find_row(self, column, func):
        if not self.initialized:
            self.init()
        value = None
        for row in self.rows:
            if value is None:
                value = row
                continue

            if func(row[column],value[column])==row[column]:
                value = row
        return value

    def min_value(self, column):
        return self.find_value(column, min)
    def min_row(self, column):
        return self.find_row(column, min)
    def max_value(self, column):
        return self.find_value(column, max)
    def max_row(self, column):
        return self.find_row(column, max)

    def get_row(self, column, value):
        if not self.initialized:
            self.init()

        for row in self.rows:
            if row[column] == value:
                return row

        return None

    def get_rows(self, column, value):
        if not self.initialized:
            self.init()

        result = []
        for row in self.rows:
            if row[column] == value:
                result.append(row)

        return result

    def get_column(self, column):
        if not self.initialized:
            self.init()
        results = []
        for row in self.rows:
            results.append(row[column])
        return results

if __name__=='__main__':
    import doctest
    doctest.testmod()
