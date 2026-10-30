
""" The statelist module. """

import logging
import os
from pathlib import Path

from eon import atoms
from eon import fileio as io
from eon.config import ConfigClass # Typing

logger = logging.getLogger('statelist')


def _intact_prefix(data, header_lines, min_fields, exact):
    """Bytes to keep, and whether the final line was torn.

    A final line past the header is torn when it has no trailing newline
    or the wrong number of fields. ``exact`` requires that count;
    otherwise the line needs at least ``min_fields``.
    """
    if not data:
        return data, False
    if data.endswith(b"\n"):
        end = len(data) - 1
        prev = data.rfind(b"\n", 0, end)
        start = 0 if prev < 0 else prev + 1
        if data[:start].count(b"\n") < header_lines:
            return data, False
        count = len(data[start:end].split())
        ok = count == min_fields if exact else count >= min_fields
        if ok:
            return data, False
        return data[:start], True
    prev = data.rfind(b"\n")
    start = 0 if prev < 0 else prev + 1
    if data[:start].count(b"\n") < header_lines:
        return data, False
    return data[:start], True


def _last_line_start(handle, size, ends_nl):
    """Byte offset of the final line. ``ends_nl`` skips a trailing newline."""
    skip = 1 if ends_nl else 0
    found = 0
    pos = size
    while pos > 0:
        step = min(4096, pos)
        pos -= step
        handle.seek(pos)
        chunk = handle.read(step)
        idx = len(chunk)
        while True:
            idx = chunk.rfind(b"\n", 0, idx)
            if idx < 0:
                break
            found += 1
            if found > skip:
                return pos + idx + 1
    return 0


def _preceding_newlines(handle, start, need):
    if need <= 0 or start <= 0:
        return 0
    handle.seek(0)
    remaining = start
    count = 0
    while remaining > 0 and count < need:
        chunk = handle.read(min(4096, remaining))
        if not chunk:
            break
        remaining -= len(chunk)
        count += chunk.count(b"\n")
    return count


def _torn_keep(path, header_lines, min_fields, exact):
    """Cut off a torn tail. Return the length to restore if the append fails.

    Only the tail, and enough of the front to see the header, is read.
    """
    path = Path(path)
    if not path.is_file():
        return 0
    size = path.stat().st_size
    if size == 0:
        return 0
    with open(path, "rb") as handle:
        handle.seek(size - 1)
        ends_nl = handle.read(1) == b"\n"
        start = _last_line_start(handle, size, ends_nl)
        if _preceding_newlines(handle, start, header_lines) < header_lines:
            return size
        if not ends_nl:
            keep = start
        else:
            handle.seek(start)
            fields = len(handle.read(size - 1 - start).split())
            bad = fields != min_fields if exact else fields < min_fields
            keep = start if bad else size
    if keep != size:
        os.truncate(path, keep)
    return keep


def _append_line(path, line, header_lines, min_fields, exact):
    """Append one full line, after cutting off a torn tail.

    The line is a single write, followed by a flush and an fsync. If
    that raises, the file is truncated back to the length it had once
    the torn tail was gone.
    """
    restore = _torn_keep(path, header_lines, min_fields, exact)
    try:
        with open(path, "a") as handle:
            handle.write(line)
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        try:
            os.truncate(path, restore)
        except OSError:
            pass
        raise


class StateList:
    """ The StateList class.  Serves as an interface to State objects and StateList metadata. """
    def __init__(self, StateClass, initial_state = None, config: ConfigClass = None):
        ''' Check to see if state_path exists and that state zero exists.
            Initializes state zero when passed a initial_state only if state
            zero doesn't already exist. '''

        if config is None:
            raise TypeError("statelist requires a ConfigClass instance")
        self.config = config
        self.path = self.config.path_states
        self.epsilon_e = self.config.comp_eps_e
        self.epsilon_r = self.config.comp_eps_r
        self.use_identical = self.config.comp_use_identical
        self.StateClass = StateClass

        # Paths
        root = Path(self.path)
        self.state_table_path = str(root / "state_table")

        # Create the statelist directory if it does not exist.
        if not root.is_dir():
            logger.warning("State list path does not exist; Creating: %s" % self.path)
            root.mkdir(parents=True)
            Path(self.state_table_path).write_text("")

        # Create the zero state directory if it does not exist.
        if not (root / "0").is_dir():
            if initial_state is None:
                raise IOError("Missing zeroth state directory and no reactant provided")
            self.StateClass(
                statepath=str(root / "0"),
                statenumber=0,
                statelist=self,
                previous_state_num=-1,
                reactant_path=initial_state,
                config=config,
            )

        # Other class variables.
        self.states = {}
        self.state_table = None

    def get_num_states(self):
        """ Returns the number of lines in the state_table file. """
        self.load_state_table(force=True)
        return len(self.state_table)

    def get_product_state(self, state_number, process_id):
        ''' Returns a State object referenced by state_number and process_id. '''
        # If the number of states in state_table is zero
        #we need to add the zero state and energy to the state table.
        if self.get_num_states() == 0:
            zst = self.get_state(0)
            self.append_state_table(zst.get_energy())

        # Load the state object containing the process we want the product for.
        st = self.get_state(state_number)
        st.load_process_table()

        # Get the state number for the product.
        newstnr = st.procs[process_id]['product']

        # If the product id is not initialized, make sure it is not a copy of an existing state.
        # Otherwise, create it, connect it to st, and return it.
        if newstnr == -1:

            # Make a list of states for which we need to compare self.configurations.
            enew = st.procs[process_id]['product_energy']
            energetically_close = []
            for id in range(self.get_num_states()):
                if abs(self.get_state(id).get_energy() - enew) < self.epsilon_e:
                    energetically_close.append(id)

            # Perform distance checks on the energetically close self.configurations.
            if len(energetically_close) > 0:
                pnew = st.get_process_product(process_id)
                for id in energetically_close:
                    p = self.get_state(id).get_reactant()
                    if atoms.match(p, pnew, self.config.comp_eps_r, self.config.comp_neighbor_cutoff, True, check_rotation=self.config.comp_check_rotation, use_identical=self.config.comp_use_identical, remove_translation=self.config.comp_remove_translation):
                        if id == state_number:
                            logging.warning("State %i process %i leads back to initial state",
                                            state_number, process_id)
                        self.register_process(st.number, id, process_id)
                        return self.get_state(id)

            # The id for the new state is the number of states.
            newstnr = self.get_num_states()

            # Create the new state object.
            newst = self.StateClass(
                statepath=self.state_path(newstnr),
                statenumber=newstnr,
                statelist=self,
                previous_state_num=state_number,
                reactant_path=st.proc_product_path(process_id),
                config=self.config,
            )
            self.register_process(st.number, newstnr, process_id)

            # Append the new state to the state table.
            self.append_state_table(st.procs[process_id]['product_energy'])

        # The product state already exists, so get it.
        else:
            newst = self.get_state(newstnr)

        # Return the product state.
        return newst

    def get_state(self, state_number):
        ''' Returns a state object. '''
        if state_number in self.states:
            return self.states[state_number]
        st = self.StateClass(
            statepath=str(Path(self.path) / str(state_number)),
            statenumber=state_number,
            statelist=self,
            config=self.config,
        )
        self.states[state_number] = st
        return st

    def load_state_table(self, force=False):
        """Load the state table from disk.

        If already loaded and *force* is false, do nothing. When *force* is
        true (used by :meth:`get_num_states`), re-read so states appended by
        another process are seen; a stale count would mint a state number
        that is already taken.
        """
        if self.state_table is not None and not force:
            return
        data = Path(self.state_table_path).read_bytes()
        kept, torn = _intact_prefix(data, 0, 2, True)
        if torn:
            logger.warning(
                "state_table %s has a torn final line; dropping it",
                self.state_table_path,
            )
        lines = kept.decode().splitlines()
        self.state_table = []
        for l in lines:
            parts = l.split()
            if not parts:
                continue
            self.state_table.append(float(parts[1]))

    def save_state_table(self):
        if self.state_table != None:
            with io.atomic_write(self.state_table_path) as f:
                for i in range(len(self.state_table)):
                    f.write("% 7d %16.5f\n" % (i, self.state_table[i]))

    def append_state_table(self, energy):
        number = self.get_num_states()
        if energy is None:
            raise IOError(
                "State %i has no reactant energy in its info file; it cannot "
                "be added to the state table" % number
            )
        line = "% 7d %16.5f\n" % (number, energy)
        _append_line(self.state_table_path, line, 0, 2, True)
        if self.state_table != None:
            self.state_table.append(energy)

    def state_path(self, state_number):
        """ Utility function to return the compiled path of a state, whether it exists or not. """
        return str(Path(self.path) / str(state_number))

if __name__ == "__main__":
    pass
