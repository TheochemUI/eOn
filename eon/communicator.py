
from array import array
import os
import shutil
import logging
logger = logging.getLogger('communicator')

from time import sleep, time
from subprocess import Popen, PIPE
import tarfile
from io import StringIO
import pickle as pickle
import glob
import re
import numpy
import sys

from eon.config import config as EON_CONFIG
from eon.config import ConfigClass # Typing

def tryint(s):
    try:
        return int(s)
    except:
        return s

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)


def get_communicator(config: ConfigClass = EON_CONFIG):
    # This is an ugly hack to "remember" a communicator as it isn't possible to construct
    # the MPI communicator multiple times and it needs to remember its object level variables.
    if hasattr(get_communicator, 'comm'):
        return get_communicator.comm

    if config.comm_type=='cluster':
        comm = Script(config.path_scratch, config.comm_job_bundle_size,
                                   config.comm_script_name_prefix,
                                   config.comm_script_path,
                                   config.comm_script_queued_jobs_cmd,
                                   config.comm_script_cancel_job_cmd,
                                   config.comm_script_submit_job_cmd)
    elif config.comm_type=='local':
        comm = Local(config.path_scratch, config.comm_local_client,
                                  config.comm_local_ncpus, config.comm_job_bundle_size)
    elif config.comm_type=='mpi':
        comm = MPI(config.path_scratch, config.comm_job_bundle_size)
    else:
        logger.error(str(config.comm_type)+" is an unknown communicator.")
        raise ValueError()
    get_communicator.comm = comm
    return comm


class NotImplementedError(Exception):
    pass


class CommunicatorError(Exception):
    pass


class EONClientError(Exception):
    """An EON client finished without outputting results, it probably crashed."""
    pass


class Communicator:
    def __init__(self, scratchpath, bundle_size=1, config: ConfigClass = EON_CONFIG):
        self.config = config
        if not os.path.isdir(scratchpath):
            # should probably log this event
            os.makedirs(scratchpath)
        self.scratchpath = scratchpath
        self.bundle_size = bundle_size

    def submit_jobs(self, data, invariants):
        '''Throws CommunicatorError if fails.'''
        raise NotImplementedError()

    def get_results(self, results):
        '''Returns a list of dictionaries containing the results.'''
        raise NotImplementedError()

    def get_queue_size(self):
        '''Returns the number of items waiting to run in the queue.'''
        raise NotImplementedError()

    def cancel_state(self, statenumber):
        '''Returns the number of workunits that were canceled.'''
        raise NotImplementedError()

    def get_bundle_size(self, job_path):
        if not isinstance(job_path, list):
            # List files in job_path.
            fnames = os.listdir(job_path)
        else:
            # job_path is already a list of filenames.
            fnames = job_path
        # Count results*.dat files.
        pattern = re.compile(r'results(?:_\d+)?.dat$')
        size = sum(1
                   for fname in fnames
                   if pattern.match(fname))
        is_bundle = not (size == 1 and "results.dat" in fnames)
        return size, is_bundle

    def unbundle(self, resultpath, keep_result):
        '''This method unbundles multiple jobs into multiple single
           jobs so the akmc script can process them.

           If the job did not return results (probably because it
           crashed or was canceled), this method will raise
           EONClientError.

        '''
        # These are the files in the result directory that we keep.
        jobpaths = [ os.path.join(resultpath,d) for d in os.listdir(resultpath)
                    if os.path.isdir(os.path.join(resultpath,d)) ]

        regex = re.compile(r"(\w+)_(\d+)(\.\w+)")
        for jobpath in jobpaths:
            basename, dirname = os.path.split(jobpath)
            if not keep_result(dirname):
                continue
            # Need to figure out how many jobs were bundled together
            # and then create the new job directories with the split files.
            bundle_size, is_bundle = self.get_bundle_size(jobpath)

            if bundle_size == 0:
                logger.error("Client running in %s returned no results. "
                             "Check its output for errors." % jobpath)
                # GH: just log the error and continue instead of quitting
                #raise EONClientError("Client running in %s returned no results. "
                #                     "Check its output for errors." % jobpath)
                continue

            results = [{'name': dirname} for i in range(bundle_size)]

            if not is_bundle:
                # Only a single task inside this job, no need to unbundle.
                for filename in glob.glob(os.path.join(jobpath, "*.*")):
                    if not (filename.endswith(".con") or
                            filename.endswith(".dat")):
                        continue
                    rootname, fname = os.path.split(filename)
                    f = open(filename,'r')
                    filedata = StringIO(f.read())
                    f.close()

                    # add result to results
                    results[0][fname] = filedata
                    results[0]['number'] = 0
            else:
                # Several tasks bundled inside this job, we need to unbundle.
                filenames = glob.glob(os.path.join(jobpath,"*_[0-9]*.*"))
                for filename in filenames:
                    if not (filename.endswith(".con") or
                            filename.endswith(".dat")):
                        continue

                    # parse filename
                    rootname, fname = os.path.split(filename)
                    match = regex.match(fname)
                    if not match:
                        continue
                    parts = match.groups()
                    index = int(parts[1])
                    key = parts[0]+parts[2]

                    # Load data into stringIO object (should we just return filehandles?)
                    try:
                        f = open(filename,'r')
                        filedata = StringIO(f.read())
                        f.close()
                    except (IOError, OSError):
                        logger.exception("Failed to read file %s" % filename)
                        continue

                    # add result to results
                    results[index][key] = filedata
                    results[index]['number'] = index

            # XXX: UGLY: We need a way to check if there are no results.
            if not any([ filename.startswith('results') for filename in list(results[0].keys())]):
                logger.warning("Failed to find a result.dat file for %s",results[0]['name'])
                results = []
            yield results

    def make_bundles(self, data, invariants):
        '''This method is a generator that bundles together multiple jobs into a single job.
           Example usage:
               for jobpath in self.make_bundles(data, invariants):
                   do_stuff()'''

        # Split jobpaths in to lists of size self.bundle_size.
        chunks = [ data[i:i+self.bundle_size] for i in range(0, len(data), self.bundle_size) ]
        for chunk in chunks:
            # create the bundle's directory

            job_path = os.path.join(self.scratchpath, chunk[0]['id'])
            os.mkdir(job_path)

            for filename in list(invariants.keys()):
                f = open(os.path.join(job_path, filename), 'w')
                file_contents, file_permissions = invariants[filename]
#                f.write(invariants[filename].getvalue())
                f.write(file_contents.getvalue())
                f.close()
                os.chmod(os.path.join(job_path, filename), file_permissions)

            # Concatenate all of the displacement and modes together.
            n = 0
            for job in chunk:
                for basename in list(job.keys()):
                    splitname = basename.rsplit(".", 1)
                    if len(splitname)!=2:
                        continue
                    if self.bundle_size == 1:
                        filename = basename
                    else:
                        filename = "%s_%d.%s" % (splitname[0], n, splitname[1])
                    f = open(os.path.join(job_path, filename), 'w')
                    f.write(job[basename].getvalue())
                    f.close()
                n += 1

            # Returns the jobpath to the new bigger workunit.
            yield job_path


class MPI(Communicator):
    def __init__(self, scratchpath, bundle_size, config: ConfigClass = EON_CONFIG):
        Communicator.__init__(self, scratchpath, bundle_size, config = config)
        from mpi4py.MPI import COMM_WORLD
        self.comm = COMM_WORLD

        self.client_ranks = [ int(r) for r in os.environ['EON_CLIENT_RANKS'].split(":") ]
        self.config.comm_job_buffer_size = len(self.client_ranks)

        self.resume_jobs = []
        if os.path.isdir(self.scratchpath):
            self.resume_jobs = [ d for d in os.listdir(self.scratchpath) if os.path.isdir(os.path.join(self.scratchpath,d)) ]
        logger.info("Found %i jobs to resume in %s", len(self.resume_jobs), self.scratchpath)

    def submit_jobs(self, data, invariants):
        ready_ranks = self.get_ready_ranks()
        for jobpath in self.make_bundles(data, invariants):
            rank = ready_ranks.pop()
            tmp = numpy.empty(1, dtype='i')
            self.comm.Recv(tmp, source=rank, tag=1)
            #buf = array('c', jobpath+'\0')
            buf = array('b')
            bufval = jobpath+'\0'
            buf.frombytes(bufval.encode())
            self.comm.Send(buf, rank)

    def run_resume_jobs(self):
        if len(self.resume_jobs) == 0: return
        ready_ranks = self.get_ready_ranks()
        while True:
            if len(self.resume_jobs) == 0: break
            if len(ready_ranks) == 0: break

            jobdir = self.resume_jobs.pop()
            rank = ready_ranks.pop()

            jobpath = os.path.join(self.scratchpath,jobdir)
            tmp = numpy.empty(1, dtype='i')
            self.comm.Recv(tmp, source=rank, tag=1)
            #buf = array('c', jobpath+'\0')
            buf = array('b')
            bufval = jobpath+'\0'
            buf.frombytes(bufval.encode())
            self.comm.Send(buf, rank)

    def get_ready_ranks(self):
        ready_ranks = []
        for rank in self.client_ranks:
            ready = self.comm.Iprobe(rank, tag=1)
            if ready:
                # logger.info("Rank %i is ready" % rank)
                ready_ranks.append(rank)
        return ready_ranks

    def get_queue_size(self):
        self.run_resume_jobs()
        nready = len(self.get_ready_ranks())
        nclients = len(self.client_ranks)
        qs = nclients - nready

        return qs

    def get_results(self, resultspath, keep_result):

        print("into get_results")
        '''Moves work from scratchpath to results path.'''
#        from mpi4py.MPI import ANY_SOURCE, Status
        import mpi4py.MPI as MPI
        print("after import")

        status = MPI.Status()
        while self.comm.Iprobe(source=MPI.ANY_SOURCE, tag=0, status=status):
            #buf = array('c', '\0'*1024)
            buf = numpy.array(['\0']*1024, dtype="S1")
            self.comm.Recv([buf, MPI.CHARACTER], source=status.source, tag=0)
            #print("after Recv")
            #print("buf: ",buf)
            strterm = numpy.where(buf == b'\0')
            strindex = strterm[0][0]
            #jobdir = buf[:buf.index('\0')].tostring()
            jobdir = buf[:strindex].tostring()
            #print("jobdir: ",jobdir.decode())
            jobdir = os.path.split(jobdir)[1].decode()
            #print("jobdir: ",jobdir)

            if self.config.debug_keep_all_results:
                shutil.copytree(os.path.join(self.scratchpath,jobdir),
                                os.path.join(self.config.path_root, self.config.debug_results_path, jobdir))
            dest_dir = os.path.join(resultspath, jobdir)
            shutil.move(os.path.join(self.scratchpath,jobdir), dest_dir)
        for bundle in self.unbundle(resultspath, keep_result):
            for result in bundle:
                yield result


    def get_number_in_progress(self):
        return int(os.environ['EON_NUMBER_OF_CLIENTS'])

    def cancel_state(self, state):
        #XXX: how to support this...
        return 0


class Local(Communicator):
    def __init__(self, scratchpath, client, ncpus, bundle_size, config: ConfigClass = EON_CONFIG):
        Communicator.__init__(self, scratchpath, bundle_size, config = config)

        # number of cpus to use
        self.ncpus = ncpus

        # path to the client
        if '/' in client:
            self.client = os.path.abspath(client)
            if not os.path.isfile(self.client):
                logger.error("Can't find client: %s", client)
                raise CommunicatorError("Can't find client binary: %s"%client)
        else:
            # is the client in the local directory?
            if os.path.isfile(client):
                self.client = os.path.abspath(client)
            # is the client in the path?
            elif sum([ os.path.isfile(os.path.join(d, client)) for d in
                       os.environ['PATH'].split(':') ]) != 0:
                self.client = client
            else:
                logger.error("Can't find client: %s", client)
                raise CommunicatorError("Can't find client binary: %s"%client)

        self.joblist = []

        import atexit
        # don't let clients hang around if the script dies
        atexit.register(self.cleanup)

    def cleanup(self):
        '''Kills the running eonclients.'''
        import signal
        for job in self.joblist:
            p = job[0]
            try:
                os.kill(p.pid, signal.SIGKILL)
            except OSError:
                pass

    def get_results(self, resultspath, keep_result):
        '''Moves work from scratchpath to results path.'''
        jobdirs = [ d for d in os.listdir(self.scratchpath)
                    if os.path.isdir(os.path.join(self.scratchpath,d)) ]

        for jobdir in jobdirs:
            if self.config.debug_keep_all_results:
                shutil.copytree(os.path.join(self.scratchpath,jobdir), os.path.join(self.config.path_root, self.config.debug_results_path,jobdir))
            dest_dir = os.path.join(resultspath, jobdir)
            shutil.move(os.path.join(self.scratchpath,jobdir), dest_dir)
        for bundle in self.unbundle(resultspath, keep_result):
            for result in bundle:
                yield result

        # Clean out scratch directory
        for name in os.listdir(self.scratchpath):
            path_name = os.path.join(self.scratchpath, name)
            if not os.path.isdir(path_name):
                continue
            shutil.rmtree(path_name)

    def check_job(self, job):
        p, jobpath = job
        if p.returncode == 0:
            logger.info('Job finished: %s' % jobpath)
            return True
        else:
            stdout, stderr = p.communicate()
            errmsg = "job failed: %s: %s" % (jobpath, stderr)
            logger.warning(errmsg)

    def submit_jobs(self, data, invariants):
        '''Run up to ncpu number of clients to process the work in jobpaths.
           The job directories are moved to the scratch path before the calculation
           is run. This method doesn't return anything.'''

        for jobpath in self.make_bundles(data, invariants):
            # move the job directory to the scratch directory
            # update jobpath to be in the scratch directory
            fstdout = open(os.path.join(jobpath, "stdout.dat"),'w')
            p = Popen(self.client, cwd=jobpath, stdout=fstdout, stderr=PIPE)
            #commands.getoutput("renice -n 20 -p %d" % p.pid)
            self.joblist.append((p,jobpath))

            while len(self.joblist) == self.ncpus:
                for i in range(len(self.joblist)):
                    p = self.joblist[i][0]
                    retval = p.poll()
                    if retval is None:
                        continue
                    else:
                        self.check_job(self.joblist[i])
                        self.joblist.pop(i)
                        break
                sleep(0.1)

        # wait for everything to finish
        for job in self.joblist:
            p = job[0]
            p.wait()
            self.check_job(job)

    def cancel_state(self, state):
        return 0

    def get_queue_size(self):
        return 0

    def get_number_in_progress(self):
        return 0


class Script(Communicator):

    def __init__(self, scratch_path, bundle_size, name_prefix, scripts_path,
                 queued_jobs_cmd, cancel_job_cmd, submit_job_cmd, config: ConfigClass = EON_CONFIG):
        Communicator.__init__(self, scratch_path, bundle_size, config=config)

        self.queued_jobs_cmd = os.path.join(scripts_path, queued_jobs_cmd)
        self.cancel_job_cmd = os.path.join(scripts_path, cancel_job_cmd)
        self.submit_job_cmd = os.path.join(scripts_path, submit_job_cmd)
        self.job_id_path = os.path.join(scratch_path, "script_job_ids")

        self.name_prefix = name_prefix

        # read in job ids
        try:
#            f = open(self.job_id_path, "r")
            f = open(self.job_id_path, "rb")
            self.jobids = pickle.load(f)
            f.close()
        except IOError:
            self.jobids = {}
            pass

    def save_jobids(self):
#        f = open(self.job_id_path, "w")
        f = open(self.job_id_path, "wb")
        pickle.dump(self.jobids, f)
        f.close()

    def get_results(self, resultspath, keep_result):
        '''Moves work from scratchpath to results path.'''
        # queued_jobs.sh jobid1 jobid2 jobid 3
        # the inverse of the jobids returned is
        # job dirs needs to map
        queued_jobs = self.get_queued_jobs()

        finished_jobids = set(self.jobids.keys()) - set(self.get_queued_jobs())

        finished_eonids = []
        for jobid in finished_jobids:
            finished_eonids.append(int(self.jobids.pop(jobid)))

        jobdirs = [ d for d in os.listdir(self.scratchpath)
                    if os.path.isdir(os.path.join(self.scratchpath,d))
                    if int(d.rsplit('_', 1)[-1]) in finished_eonids ]

        #try to return jobs in order
        sort_nicely(jobdirs)

        for jobdir in jobdirs:
            if self.config.debug_keep_all_results:
                shutil.copytree(os.path.join(self.scratchpath,jobdir), os.path.join(self.config.path_root, self.config.debug_results_path,jobdir))
            dest_dir = os.path.join(resultspath, jobdir)
            shutil.move(os.path.join(self.scratchpath,jobdir), dest_dir)

        for bundle in self.unbundle(resultspath, keep_result):
            for result in bundle:
                yield result

    def check_command(self, status, output, cmdname):
        if status != 0:
            logger.error(output)
            raise CommunicatorError("'%s' returned a non-zero exit status"%cmdname)

    def submit_jobs(self, data, invariants):
        for jobpath in self.make_bundles(data, invariants):
            # submit_job.sh jobname jobpath
            # should return a jobid
            # need to associate this jobid with our jobid
            jobpath = os.path.realpath(jobpath)
            jobname = "%s_%s" % (self.name_prefix, os.path.basename(jobpath))
            eon_jobid = jobname.rsplit('_',1)[-1]

            cmd = "%s %s %s" % (self.submit_job_cmd, jobname, jobpath)
#            status, output = commands.getstatusoutput(cmd)
            p = Popen([self.submit_job_cmd,jobname,jobpath], stdout=PIPE, stderr=PIPE)
            output, error = p.communicate()
            output = output.decode()
            error = error.decode()
            status = p.returncode
            self.check_command(status, output, cmd)

            jobid = int(output.strip())
            self.jobids[jobid] = eon_jobid

            # XXX: It is probably slow to save after EVERY job submission,
            #      but is slow better than losing jobs?
            self.save_jobids()

    def cancel_state(self, state):
        # cancel_job.sh jobid
        if len(list(self.jobids.keys())) == 0:
            return 0
        for job_id in list(self.jobids.keys()):
            cmd = "%s %i" % (self.cancel_job_cmd, job_id)
            job_id_string = "%s" % (job_id)
#            status, output = commands.getstatusoutput(cmd)
            p = Popen([self.cancel_job_cmd, job_id_string], stdout=PIPE, stderr=PIPE)
            output, error = p.communicate()
            output = output.decode()
            error = error.decode()
            status = p.returncode
            self.check_command(status, output, cmd)

            if status != 0:
                logger.warn("Job cancel failed with error: %s" % output)
        self.jobids = {}
        self.save_jobids()
        shutil.rmtree(self.config.path_scratch)
        os.makedirs(self.config.path_scratch)
        return len(list(self.jobids.keys()))

    def get_queued_jobs(self):
#        status, output = commands.getstatusoutput(self.queued_jobs_cmd)
        p = Popen([self.queued_jobs_cmd,''], stdout=PIPE, stderr=PIPE)
        output, error = p.communicate()
        output = output.decode()
        error = error.decode()
        status = p.returncode
        self.check_command(status, output, self.queued_jobs_cmd)
        queued_job_ids = []
        for line in output.split("\n"):
            try:
                queued_job_ids.append(int(line))
            except ValueError:
                continue
        return list(set(self.jobids).intersection(queued_job_ids))

    def get_number_in_progress(self):
        return 0

    def get_queue_size(self):
        return len(self.get_queued_jobs())
