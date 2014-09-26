##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------

from array import array
import config
import os
import shutil
import logging
logger = logging.getLogger('communicator')

from time import sleep, time
import subprocess
import commands
import tarfile
from cStringIO import StringIO
import cPickle as pickle
import glob
import re
import numpy

# To ensure backward compatibility
import sys
if sys.version_info < (2, 5):
    def any(iterable):
        for element in iterable:
            if element:
                return True
        return False

import re

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


def get_communicator():
    # This is an ugly hack to "remember" a communicator as it isn't possible to construct
    # the MPI communicator multiple times and it needs to remember its object level variables.
    if hasattr(get_communicator, 'comm'):
        return get_communicator.comm

    if config.comm_type=='boinc':
        comm = BOINC(config.path_scratch, config.comm_boinc_project_dir,
                config.comm_boinc_wu_template_path, config.comm_boinc_re_template_path,
                config.comm_boinc_appname, config.comm_boinc_results_path,
                config.comm_job_bundle_size, config.comm_boinc_priority)
    elif config.comm_type=='cluster':
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
    elif config.comm_type=='arc':
        comm = ARC(config.path_scratch, config.comm_job_bundle_size, 
                                config.comm_client_path, config.comm_blacklist,
                                config.comm_num_submit)
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
    def __init__(self, scratchpath, bundle_size=1):
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
        return size

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

        regex = re.compile("([\w]+)_([\d]+)(\.[\w]+)")
        for jobpath in jobpaths:
            basename, dirname = os.path.split(jobpath)
            if not keep_result(dirname):
                continue
            # Need to figure out how many jobs were bundled together
            # and then create the new job directories with the split files.
            bundle_size = self.get_bundle_size(jobpath)

            if bundle_size == 0:
                logger.error("Client running in %s returned no results. "
                             "Check its output for errors." % jobpath)
                raise EONClientError("Client running in %s returned no results. "
                                     "Check its output for errors." % jobpath)

            results = [{'name': dirname} for i in xrange(bundle_size)]

            filenames = glob.glob(os.path.join(jobpath,"*_[0-9]*.*"))
            if len(filenames) == 0:
                for filename in glob.glob(os.path.join(jobpath, "*.*")):
                    rootname, fname = os.path.split(filename)
                    f = open(filename,'r')
                    filedata = StringIO(f.read())
                    f.close()

                    # add result to results
                    results[0][fname] = filedata
                    results[0]['number'] = 0

            for filename in filenames:
                if filename[-3:] != 'con' and filename[-3:] != 'dat':
                    continue
                try:
                    # parse filename
                    rootname, fname = os.path.split(filename)

                    match = regex.match(fname)
                    if not match:
                        continue
                    parts = match.groups()

                    index = int(parts[1])
                    key = parts[0]+parts[2]

                    # Load data into stringIO object (should we just return filehandles?)
                    f = open(filename,'r')
                    filedata = StringIO(f.read())
                    f.close()

                    # add result to results
                    results[index][key] = filedata
                    results[index]['number'] = index
                except:
                    logger.exception("Failed to handle file %s" % filename)

            # XXX: UGLY: We need a way to check if there are no results.
            if not any([ filename.startswith('results') for filename in results[0].keys() ]):
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

            for filename in invariants.keys():
                f = open(os.path.join(job_path, filename), 'w')
                f.write(invariants[filename].getvalue())
                f.close()

            # Concatenate all of the displacement and modes together.
            n = 0
            for job in chunk:
                for basename in job.keys():
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

class BOINC(Communicator):
    def __init__(self, scratchpath, boinc_project_dir, wu_template, result_template,
                 appname, boinc_results_path, bundle_size, priority):
        '''
        This constructor modifies sys.path to include the BOINC python
        modules.  It then tries to connect to the BOINC mysql database raising
        exceptions if there are problems connecting. It also creates a file
        named uniqueid in the scratchpath to identify BOINC jobs as belonging
        to this akmc run if it doesn't exist. It then reads in the uniqueid
        file and stores that as an integer in self.uniqueid.
        '''

        Communicator.__init__(self, scratchpath, bundle_size)
        self.wu_template = wu_template
        self.result_template = result_template
        self.appname = appname
        self.boinc_project_dir = boinc_project_dir
        self.boinc_results_path = boinc_results_path
        self.priority = priority
        os.environ['BOINC_PROJECT_DIR'] = self.boinc_project_dir
        import sys
        sys.path.insert(0, os.path.join(self.boinc_project_dir, 'py'))
        try:
            import Boinc.database
            import Boinc.db_base
            import Boinc.boinc_db
        except ImportError:
            raise CommunicatorError("The Boinc python module could not be imported.\n"
                "Perhaps the boinc project path is set incorrectly?")

        self.database = Boinc.database 
        self.boinc_db_constants = Boinc.boinc_db
        self.db_base = Boinc.db_base

        try:
            self.database.connect_default_config()
        except:
            # XXX: This error handling is maybe a little ugly, but provides all the information
            # that you would want to know. The exception that connect_default_config() throws
            # is not helpful. It often will just say that it couldn't parse a xml file when it
            # really means it can't find the project's config.xml file.
            import traceback
            # print the traceback from connect_default_config
            traceback.print_exc()
            # raise a nice human readable error
            raise CommunicatorError("Couldn't connect to the BOINC database.")

        self.dbconnection = self.db_base.dbconnection
        self.cursor = self.dbconnection.cursor()

        # generate our unique id if it doesn't already exist.
        uniqueid_path = os.path.join(self.scratchpath, "uniqueid")

        if not os.path.isfile(uniqueid_path):
            f = open(uniqueid_path, 'w')
            import random 
            uid = random.randint(0, 2147483647)
            f.write("%s\n" % uid)
            f.close()
            logger.debug("Wrote new unique id %i to %s" % (uid, uniqueid_path))
        try:
            f = open(uniqueid_path)
        except IOError:
            raise CommunicatorError("Unable to open the uniqueid file: %s" % uniqueid_path)

        try:
            self.uniqueid = int(f.read().strip())
            logger.debug("Read in unique id %i from %s" % (self.uniqueid, 
                uniqueid_path))
        except ValueError:
            raise CommunicatorError("Trouble converting uniqueid value in %s to integer" % uniqueid_path)

        self.average_flops = self.get_average_flops()
        logger.debug("Current average flops per wu is %.2e", self.average_flops)

    def get_average_flops(self):
        'This function might be slow with large result tables and without '
        'mysql indices on result.cpu_time, result.workunits, result.hostid, '
        'and workunit.batch.'

        # number of wus to average over
        limit = 500
        query = "select r.cpu_time*h.p_fpops " \
                "from workunit w, result r, host h "\
                "where r.workunitid=w.id and r.hostid=h.id and w.batch=%i "\
                "and cpu_time>0 limit %i" % (self.uniqueid, limit)

        self.cursor.execute(query)

        rows = self.cursor.fetchall()
        if rows:
            average_flops = 0.0
            counter = 0
            for row in rows:
                average_flops += row.values()[0]
                counter += 1
            average_flops /= counter
        else:
            # 2e11 flops is about a 100 second job (assuming 2 gigaflop cpu)
            average_flops = 2e11

        return average_flops

    def get_queue_size(self):
        server_state = self.boinc_db_constants.RESULT_SERVER_STATE_UNSENT
        query = 'select count(*) from result where batch=%i and server_state=%i'
        query = query % (self.uniqueid, server_state)
        self.cursor.execute(query)
        row = self.cursor.fetchone()
        number_unsent = row['count(*)']
        return number_unsent

    def get_number_in_progress(self):
        server_state = self.boinc_db_constants.RESULT_SERVER_STATE_IN_PROGRESS
        query = 'select count(*) from result where batch=%i and server_state=%i'
        query = query % (self.uniqueid, server_state)
        self.cursor.execute(query)
        row = self.cursor.fetchone()
        number_in_progress = row['count(*)']
        return number_in_progress

    def cancel_state(self, statenumber):
        # XXX: This function might be too expensive. Probably needs to be
        #      profiled later. It has to get all of the result rows that correspond
        #      to a unique id find which state they correspond to by parsing the
        #      wu name and then update the rows that correspond to statenumber.

        state_unsent = self.boinc_db_constants.RESULT_SERVER_STATE_UNSENT
        q1 = "select id,workunitid,name from result where batch=%i and server_state=%i"
        q1 = q1 % (self.uniqueid, state_unsent)

        self.cursor.execute(q1)

        if self.cursor.rowcount == 0:
            return 0

        result_ids = []
        workunit_ids = []
        while 1:
            row = self.cursor.fetchone()
            if row == None:
                break

            result_id = row['id']
            workunit_id = row['workunitid']
            name = row['name']
            result_statenumber = name.split('_')[1]

            if statenumber == int(result_statenumber):
                result_ids.append(str(result_id))
                workunit_ids.append(str(workunit_id))

        resultid_string = '('+','.join(result_ids) + ')'
        state_over = self.boinc_db_constants.RESULT_SERVER_STATE_OVER
        outcome_not_needed = self.boinc_db_constants.RESULT_OUTCOME_DIDNT_NEED
        error_mask_cancelled = self.boinc_db_constants.WU_ERROR_CANCELLED
        q1 = "update result set server_state=%i, outcome=%i where id in %s"
        q1 = q1 % (state_over, outcome_not_needed, resultid_string)
        self.cursor.execute(q1)

        workunitid_string = '('+','.join(workunit_ids) + ')'
        q2 = "update workunit set error_mask=%i, transition_time=%i where id in %s" 
        q2 = q2 % (error_mask_cancelled, int(time()), workunitid_string)
        self.cursor.execute(q2)
        num_cancelled_wu = self.cursor.rowcount

        self.db_base.dbconnection.commit()

        return num_cancelled_wu

    def dir_hier_path(self, filename):
        cmd = os.path.join(self.boinc_project_dir,"bin","dir_hier_path")
        path = commands.getoutput("%s %s" % (cmd, filename))
        return path

    def submit_jobs(self, jobdata, invariants):
        now = time()
        chunks = [ jobdata[i:i+self.bundle_size] for i in
                   range(0, len(jobdata), self.bundle_size) ]
        for jobs in chunks:
            wu_name = "%i_%s" % (self.uniqueid, jobs[0]['id'])
            tarname = "%s.tgz" % wu_name
            tarpath = self.dir_hier_path(tarname)
            tar = tarfile.open(tarpath, "w:gz")

            jobfiles = {}
            n=0
            for job in jobs:
                job.pop('id')
                for origname,data in job.iteritems():
                    splitname = origname.rsplit(".",1)
                    newname = "%s_%d.%s" % (splitname[0], n, splitname[1])
                    jobfiles[newname] = data
                n += 1

            # Add the files in job and in invariants to the tar file
            for filelist in (jobfiles, invariants):
                for filename, filehandle in filelist.iteritems():
                    info = tarfile.TarInfo(name=filename)
                    info.size = len(filehandle.getvalue())
                    info.mtime = now
                    filehandle.seek(0)
                    tar.addfile(info, filehandle);
            tar.close()

            self.create_work(tarpath, wu_name)

    def create_work(self, tarpath, wu_name):
        create_wu_cmd = os.path.join('bin', 'create_work')

        # XXX: make sure permissions are correct
        # this should be a config option for the boinc group
        mode = 0666
        os.chmod(tarpath, mode)

        arglist = [create_wu_cmd]
        arglist.append("-appname")
        arglist.append(self.appname)
        arglist.append("-wu_name")
        arglist.append(wu_name)
        arglist.append("-wu_template")
        arglist.append(self.wu_template)
        arglist.append("-result_template")
        arglist.append(self.result_template)
        arglist.append("-batch")
        arglist.append(str(self.uniqueid))
        arglist.append("-priority")
        arglist.append(str(self.priority))
        arglist.append("-rsc_fpops_est")
        arglist.append(str(self.average_flops))
        # last arguments are the filenames
        arglist.append("%s.tgz" % wu_name)

        p = subprocess.Popen(arglist, cwd=self.boinc_project_dir,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        retval = p.wait()

        if retval != 0:
            stdout, stderr = p.communicate()
            errstr = "Problem submitting the BOINC workunit:\nstdout: %s\nstderr %s" % \
                        (stdout, stderr)
            raise CommunicatorError(errstr)

    def get_results(self, resultspath, keep_result):
        all_results = os.listdir(self.boinc_results_path)
        all_results = [ f for f in all_results if '_' in f ]
        my_results = [ f for f in all_results
                                if f.split('_')[0] == str(self.uniqueid) and
                                'no_output_files' not in f]

        for resultfile in my_results:
            # jobname looks like state#_job#
            # jobname is everything but the first underscore records
            jobname = '_'.join(resultfile.split('_')[1:])
            resultpath = os.path.join(self.boinc_results_path, resultfile)
            if not keep_result(jobname):
                os.remove(resultpath)
                continue

            try:
                tar = tarfile.open(resultpath)
                bundle_size = self.get_bundle_size(tar.getnames())
                results = [ {'name':jobname} for i in range(bundle_size) ]
                if config.debug_keep_all_results:
                    rp = os.path.join(config.path_root,config.debug_results_path)
                    if not os.path.isdir(rp):
                        os.mkdir(rp)
                    os.mkdir(os.path.join(rp,jobname))
                for tarinfo in tar:
                    try:
                        index = int(tarinfo.name.split('_')[-1].split('.')[0])
                    except:
                        logger.exception("Failed to process file %s in tar" % tarinfo.name)
                        continue
                    splitname = tarinfo.name.rsplit(".",1) 
                    newfilename = "%s.%s" % (splitname[0].rsplit("_",1)[0],splitname[1])
                    # Read the file in the tar archive into a stringio; you cannot return the filehandle
                    # that extractfile returns as it will be closed when tar.close() is called.
                    fh = StringIO(tar.extractfile(tarinfo).read())
                    if config.debug_keep_all_results:
                        f = open(os.path.join(config.path_root, config.debug_results_path, jobname, tarinfo.name),'w')
                        f.write(fh.getvalue())
                        f.close()
                        fh.seek(0)
                    results[index][newfilename] = fh
                    results[index]["number"] = index
                tar.close()
                os.remove(resultpath)
            except:
                logger.exception(
                        "Something tar-file related went wrong with file %s" % resultpath)
                try:
                    os.remove(resultpath)
                except:
                    logger.exception("Failed to remove %s" % resultpath)
                continue

            for result in results:
                yield result

class MPI(Communicator):
    def __init__(self, scratchpath, bundle_size):
        Communicator.__init__(self, scratchpath, bundle_size)

        from mpi4py.MPI import COMM_WORLD
        self.comm = COMM_WORLD

        self.client_ranks = [ int(r) for r in os.environ['EON_CLIENT_RANKS'].split(":") ]
        config.comm_job_buffer_size = len(self.client_ranks)

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
            buf = array('c', jobpath+'\0')
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
            buf = array('c', jobpath+'\0')
            self.comm.Send(buf, rank)

    def get_ready_ranks(self):
        ready_ranks = []
        for rank in self.client_ranks:
            ready = self.comm.Iprobe(rank, tag=1)
            if ready:
                logger.info("Rank %i is ready" % rank)
                ready_ranks.append(rank)
        return ready_ranks

    def get_queue_size(self):
        self.run_resume_jobs()
        nready = len(self.get_ready_ranks())
        nclients = len(self.client_ranks)
        qs = nclients - nready

        return qs

    def get_results(self, resultspath, keep_result):
        '''Moves work from scratchpath to results path.'''
        from mpi4py.MPI import ANY_SOURCE, Status

        status = Status()
        while self.comm.Iprobe(ANY_SOURCE, tag=0, status=status):
            buf = array('c', '\0'*1024)
            self.comm.Recv(buf, source=status.source, tag=0)
            jobdir = buf[:buf.index('\0')].tostring()
            jobdir = os.path.split(jobdir)[1]

            if config.debug_keep_all_results:
                shutil.copytree(os.path.join(self.scratchpath,jobdir),
                                os.path.join(config.path_root, config.debug_results_path, jobdir))
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
    def __init__(self, scratchpath, client, ncpus, bundle_size):
        Communicator.__init__(self, scratchpath, bundle_size)

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
            if config.debug_keep_all_results:
                shutil.copytree(os.path.join(self.scratchpath,jobdir), os.path.join(config.path_root, config.debug_results_path,jobdir))
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
            p = subprocess.Popen(self.client,cwd=jobpath,
                    stdout=fstdout, stderr=subprocess.PIPE)
            #commands.getoutput("renice -n 20 -p %d" % p.pid)
            self.joblist.append((p,jobpath))

            while len(self.joblist) == self.ncpus:
                for i in range(len(self.joblist)):
                    p = self.joblist[i][0]
                    retval = p.poll()
                    if retval == None:
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
                 queued_jobs_cmd, cancel_job_cmd, submit_job_cmd):
        Communicator.__init__(self, scratch_path, bundle_size)

        self.queued_jobs_cmd = os.path.join(scripts_path, queued_jobs_cmd)
        self.cancel_job_cmd = os.path.join(scripts_path, cancel_job_cmd)
        self.submit_job_cmd = os.path.join(scripts_path, submit_job_cmd)
        self.job_id_path = os.path.join(scratch_path, "script_job_ids")

        self.name_prefix = name_prefix

        # read in job ids
        try:
            f = open(self.job_id_path, "r")
            self.jobids = pickle.load(f)
            f.close()
        except IOError:
            self.jobids = {}
            pass

    def save_jobids(self):
        f = open(self.job_id_path, "w")
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
            if config.debug_keep_all_results:
                shutil.copytree(os.path.join(self.scratchpath,jobdir), os.path.join(config.path_root, config.debug_results_path,jobdir))
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
            status, output = commands.getstatusoutput(cmd)
            self.check_command(status, output,cmd)

            jobid = int(output.strip())
            self.jobids[jobid] = eon_jobid

            # XXX: It is probably slow to save after EVERY job submission,
            #      but is slow better than losing jobs?
            self.save_jobids()

    def cancel_state(self, state):
        # cancel_job.sh jobid
        if len(self.jobids.keys()) == 0:
            return 0
        for job_id in self.jobids.keys():
            cmd = "%s %i" % (self.cancel_job_cmd, job_id)
            status, output = commands.getstatusoutput(cmd)
            if status != 0:
                logger.warn("Job cancel failed with error: %s" % output)
        self.jobids = {}
        self.save_jobids()
        shutil.rmtree(config.path_scratch) 
        os.makedirs(config.path_scratch)
        return len(self.jobids.keys())

    def get_queued_jobs(self):
        # get_queued_jobs.sh
        # should return the jobids of the jobs in the queue
        status, output = commands.getstatusoutput(self.queued_jobs_cmd)
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


class ARC(Communicator):

    def __init__(self, scratchpath, bundle_size, client_path, blacklist, num_submit):
        self.init_completed = False

        Communicator.__init__(self, scratchpath, bundle_size)

        try:
            import arclib
            self.arclib = arclib
        except ImportError:
            raise CommunicatorError("ARCLib can't be imported. Check if PYTHONPATH is set correctly")

        self.arclib.SetNotifyLevel(self.arclib.WARNING)

        self.client_path = client_path
        self.num_submit = num_submit
        self.blacklist = blacklist

        self.queue_info = None

        # Check grid certificate proxy
        try:
            c = self.arclib.Certificate(self.arclib.PROXY)
        except self.arclib.CertificateError, msg:
            raise CommunicatorError(str(msg) + ".\n\nForgot to run grid-proxy-init?\n")
        if c.IsExpired():
            raise CommunicatorError("Grid proxy has expired!")
        logger.info("Grid proxy is valid for " + c.ValidFor())

        # Get a list of jobs, and find their statuses.
        self.busy_clusters = []
        self.active_jobs = []
        self.active_jobsfilename = os.path.join(self.scratchpath, "active_jobs.txt")
        if os.path.isfile(self.active_jobsfilename):
            jobids = {}
            f = open(self.active_jobsfilename, "r")
            for line in f:
                (jid, jname) = line.split('#')
                jobids[jid] = jname[:-1]  # (Remove trailing '\n' from name).
        else:
            jobids = {}
        if jobids:
            for info in self.arclib.GetJobInfo(jobids.keys()):
                job = {"id": info.id, "name": jobids[info.id]}
                if info.status in [ "FINISHED", "FAILED" ]:
                    job["stage"] = "Done"
                    job["success"] = (info.status == "FINISHED")
                elif info.status in [ "DELETED", "KILLED", "KILLING" ]:
                    job["stage"] = "Aborted" # Supposed to disappear by itself soonish
                elif info.status in [ "ACCEPTING", "ACCEPTED", "PREPARING", "PREPARED", "SUBMITTING", "INLRMS:Q" ]:
                    job["stage"] = "Queueing"
                    if info.cluster not in self.busy_clusters:
                        self.busy_clusters.append(info.cluster)
                elif info.status == "":
                    # XXX: The most common reason for info.status == "" is that
                    # the job was submitted so recently that ARC info.sys.
                    # hasn't picked up the job yet, which is why I decided
                    # to consider it to be "Queueing". But it could also be the
                    # ARC info.sys being down, or other problems.
                    job["stage"] = "Queueing"
                    if info.cluster not in self.busy_clusters:
                        self.busy_clusters.append(info.cluster)
                else:
                    job["stage"] = "Running"

                if job["stage"] != "Aborted":
                    self.active_jobs.append(job)
        #logger.info("Job %s / %s found in state %s (%s)" % (job["name"], job["id"], job["stage"], info.status))

        # loads the server job-queue
        self.job_queue_filename = os.path.join(self.scratchpath, "job_queue.pickle")
        if os.path.isfile(self.job_queue_filename):
            f = open(self.job_queue_filename, "r")
            self.job_queue = pickle.load(f)
            f.close()
        else:
            self.job_queue = []
        self.init_completed = True


    def __del__(self):
        """
            Remember jobs for future invocations.
            """
        logger.debug("ARC.__del__ invoked!")

        if self.init_completed:
            # Save the jobs to two files for the future, but only if
            # __init__() was successful - if it wasn't we might not have
            # read all the jobs from previous runs, in which case
            # information on those jobs would be overwritten.
            # (And no jobs could have been submitted or retrieved;
            # init_completed = False means we crashed at an early stage, so
            # there's nothing new to save)
            f = open(self.active_jobsfilename, "w")
            for j in self.active_jobs:
                if j["stage"] not in ["Aborted", "Retrieved"]:
                    f.write(j["id"] + '#' + j["name"] + '\n')
            f.close()

            f = open(self.job_queue_filename, "w")
            pickle.dump(self.job_queue, f)
            f.close()

    def create_wrapper_script(self):
        '''Create a wrapper script to execute a job.  Return path to script.'''

        s = """
            #!/bin/bash

            set -e # Immediately exit if any command fail

            ls -l
            if [ -f client-bin ]; then
            # It seems we got a EON client binary as an inputfile. It does
            # (probably) not have execute bit set, and it might be a
            # sym-link to a file we don't own, so we have to make a copy
            # and change permissions of the copy instead.
            export PATH=$PATH:$PWD
            cp client-bin client
            chmod +x client
            fi
            ls -l
            tar jxvf $1.tar.bz2
            cd $1
            client
            cd $HOME
            tar jcvf $1.tar.bz2 $1
            echo $1
            """
        script_path = os.path.join(self.scratchpath, 'wrapper.sh')
        try:
            f = open(script_path, "w")
            f.write(s)
            f.close()
        except Exception, msg:
            raise CommunicatorError("Can't create wrapper script: %s" % msg)

        return script_path


    def create_tarball(self, src, dest):
        """Pack directory 'src' into tar.bz2 file 'dest', makeing sure it will unpack into
            a dir called basename(src), rather than path/to/src"""

        # Remove trailing '/'; it would cause trouble with
        # os.path.dirname() and os.path.basename()
        if src[-1] == '/':
            src = src[:-1]

        # Since we'll change directory for a while, we should make sure we
        # use absolute paths, rather than relative ones:
        src = os.path.abspath(src)
        dest = os.path.abspath(dest)

        dirname = os.path.dirname(src)
        basename = os.path.basename(src)

        cwd = os.getcwd()
        try:
            os.chdir(dirname)
            tarball = tarfile.open(dest, 'w:bz2')
            tarball.add(basename)
            tarball.close()
        finally:
            os.chdir(cwd)


    def open_tarball(self, filename, dest_dir):
        """Pack upp tar.bz2 file beneth the directory dest_dir."""
        try:
            tarball = tarfile.open(filename, 'r:bz2')

            # For security reasons, filter out filenames that might end up
            # outside of 'dest_dir':
            files = tarball.getmembers()
            good_files = [ f for f in files if f.name[0:2] != '..' and f.name[0] != '/' ]

            for f in good_files:
                tarball.extract(path=dest_dir, member=f)
            tarball.close()
        except IOError:
            logger.warning("Could not open tarball %s" % (filename))


    def create_job(self, job_path, wrapper_path):
        '''Prepare a job who's inputfiles are found in 'job_path'.
            Return a string to be used as input for xsrl and job name.'''

        # Remove trailing '/'; it would cause trouble with os.path.basename()
        if job_path[-1] == '/':
            job_path = job_path[:-1]

        basename = os.path.basename(job_path)
        tarball_path = os.path.join(self.scratchpath, basename + ".tar.bz2")
        self.create_tarball(job_path, tarball_path)
        shutil.rmtree(job_path)

        s = "&"
        s += "(executable=%s)" % os.path.basename(wrapper_path)
        s += "(arguments=%s)" % basename

        s += "(inputFiles="
        if self.client_path:
            s += "(%s %s)" % ("client-bin", self.client_path)
        s += "(%s %s)" % (os.path.basename(wrapper_path), wrapper_path)
        s += "(%s %s)" % (os.path.basename(tarball_path), tarball_path)
        s += ")"

        s += "(outputFiles="
        s += "(%s '')" % os.path.basename(tarball_path)
        s += ")"

        if config.debug_keep_all_results:
            s += "(stdout=stdout)"
            s += "(stderr=stderr)"
            s += "(gmlog=gmlog)"

        if not self.client_path:
            s += "(runTimeEnvironment=APPS/CHEM/EON2)"

        jobname = "%s" % basename
        s += "(jobName=%s)" % jobname

        # requirements to the host
        # 500 MB memory
        s += '(Memory="500")'
        # 120 minutes
        s += '(wallTime="120")'

        logger.debug("job_string: " + s)

        return s, jobname, tarball_path


    def submit_job_to_target(self, hostname):
        if(len(self.job_queue)):

            job_string, job_name, job_tarball = self.job_queue.pop(0)
            xrsl = self.arclib.Xrsl(job_string)

            # fetch the target which hostname passed was 
            target = None
            all_targets = self.arclib.ConstructTargets(self.queue_info, xrsl)
            for t in all_targets:
                if t.cluster.hostname == hostname:
                    target = self.arclib.PerformStandardBrokering([t])
                    break
            if target:
                try:
                    job_id = self.arclib.SubmitJob(xrsl, target)
                    os.remove(job_tarball)
                except self.arclib.JobSubmissionError, msg:
                    raise CommunicatorError(msg)
                except self.arclib.XrslError, msg:
                    raise CommunicatorError(msg)

                self.arclib.AddJobID(job_id, job_name)
                self.active_jobs.append({"id": job_id, "name": job_name, "stage":"Queueing"})
                logger.info("Submitted " + job_id)
            else:
                logger.info("Failed to submit to " + hostname)                


    def get_targets_hostnames(self, xrsl):
        """Get list of hostnames of clusters+queues we can submit to."""
        hostnames = []
        targets = self.get_targets(xrsl)
        for t in targets:
            hostnames.append(t.cluster.hostname)
        return hostnames


    def get_targets(self, xrsl):
        """Get list of clusters+queues we can submit to."""

        if not self.queue_info:
            self.queue_info = self.arclib.GetQueueInfo()

        targets_initial = self.arclib.ConstructTargets(self.queue_info, xrsl)

        logger.debug("List of targets: " + ', '.join([ t.cluster.hostname for t in targets_initial ]))

        if self.blacklist:
            targets_not_bl = []
            for t in targets_initial:
                if t.cluster.hostname not in self.blacklist:
                    targets_not_bl.append(t)
            logger.debug("List of targets after blacklisting: " + ', '.join([ t.cluster.hostname for t in targets_not_bl ]))
        else:
            targets_not_bl = targets_initial

        if self.busy_clusters:
            targets = []
            for t in targets_not_bl:
                if t.cluster.hostname not in self.busy_clusters:
                    targets.append(t)
            logger.debug("List of targets after queing check: " + ', '.join([ t.cluster.hostname for t in targets ]))
        else:
            targets = targets_not_bl
        return self.arclib.PerformStandardBrokering(targets)


    def submit_jobs(self, data, invariants):
        '''Throws CommunicatorError if fails.'''
        # add new jobs to queue
        wrapper_path = self.create_wrapper_script()
        for job_path in self.make_bundles(data, invariants):
            job_string, job_name, job_tarball = self.create_job(job_path, wrapper_path)
            self.job_queue.append([job_string, job_name, job_tarball])

        # are there jobs to submit
        if len(self.job_queue): 
            # uses first element in queue to get list of free resources
            xrsl = self.arclib.Xrsl(self.job_queue[0][0])
            targets_hostnames = self.get_targets_hostnames(xrsl)
        else:
            targets_hostnames = []

        for target_hostname in targets_hostnames:
            submitted_to_target = 0
            while (submitted_to_target < self.num_submit):
                self.submit_job_to_target(target_hostname)
                submitted_to_target = submitted_to_target + 1


    def get_job_output(self, jobid, resultspath):
        """Fetch the output files of a job.
            The files are put in a subdirectory of resultspath,
            and the full path of the subdirectory is returned."""

        n = jobid.split('/')[-1]
        outputdir = os.path.join(resultspath, n)
        if not os.path.isdir(outputdir):
            os.makedirs(outputdir)

        ftp = self.arclib.FTPControl()
        ftp.DownloadDirectory(jobid, outputdir)

        return outputdir


    def get_results(self, resultspath, keep_result):
        '''Returns a list of directories containing the results.'''

        result_dirs = []
        done = [ j for j in self.active_jobs if j["stage"] == "Done" ]
        for job in done:
            jid = job["id"]
            jname = job["name"]

            p = self.get_job_output(jid, self.scratchpath)
            tarball = os.path.join(p, "%s.tar.bz2" % jname)

            if job["success"]:
                self.open_tarball(tarball, resultspath)
                logger.info("Fetched %s / %s" % (jname, jid)) 
            else:
                logger.warning("Job %s / %s FAILED.\nOutput files can be found in %s" % (jname, jid, p)) 

            job["stage"] = "Retrieved"

            self.arclib.RemoveJobID(jid) # Remove from ~/.ngjobs
            self.arclib.CleanJob(jid) # Remove from ARC sever

            if not config.debug_keep_all_results:
                os.remove(tarball) # keep server scratch dir clean
                shutil.rmtree(p)

        for bundle in self.unbundle(resultspath, keep_result):
            for result in bundle:
                yield result


    def get_queue_size(self):
        '''Returns the number of items waiting to run in the server queue.'''
        return len(self.job_queue)


    def cancel_job(self, job):
        '''Removed job from cluster queue'''
        try:
            self.arclib.CleanJob(job["id"])
        except:
            logger.warning('Did not clean up the cancelled job!')
            return

        # we are sure that the job is removed from the cluster
        # remove id from the local list
        self.arclib.RemoveJobID(job["id"])
        job["stage"] = "Aborted"
        return


    def cancel_state(self, statenumber):
        '''Returns the number of workunits that were canceled.'''

        logger.debug("cancel_state called with statenumber = %i (%s)" % (int(statenumber), type(statenumber)))

        n = 0
        for j in self.active_jobs:
            sn = j["name"].split('_')[0]
            #            if int(sn) == int(statenumber) and j["stage"] not in [ "Aborted", "Retrieved" ]:
            if int(sn) <= int(statenumber) and j["stage"] not in [ "Aborted", "Retrieved" ]:
                self.cancel_job(j)
                logger.debug("Canceling job %s / %s" % (j["name"], j["id"]))
                n += 1

        self.job_queue = []
        return n

