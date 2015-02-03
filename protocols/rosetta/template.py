#!/usr/bin/env python2
# This work is licensed under the terms of the MIT license. See LICENSE for the full text.
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=240:00:00
#$ -t 1-#$#numclusterjobs#$#
#$ -l arch=linux-x64
#$ -l mem_free=#$#mem_free#$#
#$ -l netapp=1G,scratch=1G

# 1-#$#numclusterjobs#$# total (in case less is set above)

import socket
import sys

print "Python version:", sys.version
print "Hostname:", socket.gethostname()

from datetime import *
import os
import subprocess
import time
import sys
import shutil
import inspect
import gzip
import tempfile
import math

try:
    import cPickle as pickle
except:
    print 'cPickle not available, using regular pickle module'
    import pickle

script_name = '#$#scriptname#$#.py'

print "Script:", script_name

# Constants
tasks_per_process = #$#tasks_per_process#$#
total_number_tasks = #$#numjobs#$#

cluster_rosetta_bin = '#$#cluster_rosetta_bin#$#'
cluster_rosetta_db = '#$#cluster_rosetta_db#$#'

local_rosetta_bin = '#$#local_rosetta_bin#$#'
local_rosetta_db = '#$#local_rosetta_db#$#'
local_scratch_dir = '/tmp'

job_pickle_file = 'data/job_dict.pickle'

app_name = '#$#appname#$#'

zip_rosetta_output = True

generic_rosetta_args = [
    #$#rosetta_args_list#$#
]

sge_task_id = 0
run_locally = True
run_on_sge = False
if os.environ.has_key("SGE_TASK_ID"):
    sge_task_id = long(os.environ["SGE_TASK_ID"])
    run_locally = False
    run_on_sge = True

job_id = 0
if os.environ.has_key("JOB_ID"):
    job_id=long(os.environ["JOB_ID"])

def roundTime(dt=None, roundTo=1):
    """Round a datetime object to any time period (in seconds)
    dt : datetime.datetime object, default now.
    roundTo : Closest number of seconds to round to, default 1 second.
    Author: Thierry Husson 2012 - Use it as you want but don't blame me.
    http://stackoverflow.com/questions/3463930/how-to-round-the-minute-of-a-datetime-object-python/10854034#10854034
    """
    if dt == None : dt = datetime.now()
    seconds = total_seconds(dt - dt.min)
    # // is a floor division, not a comment on following line:
    rounding = (seconds+roundTo/2) // roundTo * roundTo
    return dt + timedelta(0,rounding-seconds,-dt.microsecond)

def total_seconds(td):
    '''
    Included in python 2.7 but here for backwards-compatibility for old Python versions
    '''
    return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6

class Reporter:
    '''
    Used to report and estimate remaining runtime when script is run locally
    '''
    def __init__(self, task, entries='files', print_output=True):
        self.print_output = print_output
        self.start = time.time()
        self.entries = entries
        self.lastreport = self.start
        self.task = task
        self.report_interval = 1 # Interval to print progress (seconds)
        self.n = 0
        self.completion_time = None
        if self.print_output:
            print '\nStarting ' + task
        self.total_count = None # Total tasks to be processed
    def set_total_count(self, x):
        self.total_count = x
    def report(self, n):
        self.n = n
        t = time.time()
        if self.print_output and self.lastreport < (t-self.report_interval):
            self.lastreport = t
            if self.total_count:
                percent_done = float(self.n)/float(self.total_count)
                time_now = time.time()
                est_total_time = (time_now-self.start) * (1.0/percent_done)
                time_remaining = est_total_time - (time_now-self.start)
                minutes_remaining = math.floor(time_remaining/60.0)
                seconds_remaining = int(time_remaining-(60*minutes_remaining))
                print "  Processed: "+str(n)+" "+self.entries+" (%.1f%%) %02d:%02d" % (percent_done*100.0, minutes_remaining, seconds_remaining)
            else:
                print "  Processed: "+str(n)+" "+self.entries

    def increment_report(self):
        self.report(self.n+1)
    def decrement_report(self):
        self.report(self.n-1)
    def add_to_report(self,x):
        self.report(self.n+x)
    def done(self):
        self.completion_time = time.time()
        if self.print_output:
            print 'Done %s, processed %d %s, took %.3f seconds\n' % (self.task,self.n,self.entries,self.completion_time-self.start)
    def elapsed_time(self):
        if self.completion_time:
            return self.completion_time - self.start
        else:
            return time.time() - self.start

def run_single(task_id, rosetta_bin, rosetta_db, scratch_dir=local_scratch_dir, verbosity=1, move_output_files=False):
    global job_dict

    time_start = roundTime()

    if verbosity >= 1:
        print 'Starting time:',time_start
        print 'Task id:',task_id

    if os.path.isfile(job_pickle_file):
        p = open(job_pickle_file,'r')
        job_dict = pickle.load(p)
        p.close()

    job_dirs = sorted(job_dict.keys())
    
    job_dir = job_dirs[task_id]

    if verbosity >= 1:
        print 'Job dir:', job_dir

    # Make temporary directories
    if not os.path.isdir(scratch_dir):
        os.mkdir(scratch_dir)
    tmp_data_dir = tempfile.mkdtemp(prefix='%d.%d_data_' % (job_id,task_id), dir=scratch_dir)
    tmp_output_dir = tempfile.mkdtemp(prefix='%d.%d_output_' % (job_id,task_id), dir=scratch_dir)

    if verbosity >= 1:
        print 'Temporary data dir:', tmp_data_dir
        print 'Temporary output dir:', tmp_output_dir

    args=[
        os.path.join(rosetta_bin, app_name),
    ]

    # Append specific Rosetta database path if this argument is included
    if len(rosetta_db) > 0:
        args.append('-database')
        args.append(rosetta_db)

    def copy_file_helper(original_file):
        new_file = os.path.join(tmp_data_dir, os.path.basename(original_file))
        shutil.copy(original_file, new_file)
        value = os.path.relpath(new_file, tmp_output_dir)
        if verbosity>=1:
            print 'Copied file to local scratch:', os.path.basename(original_file)
        return value

    flags_dict = job_dict[job_dir]
    for flag in flags_dict:
        if flag.startswith('FLAGLIST'):
            for split_flag in flags_dict[flag]:
                args.append(str(split_flag))
            continue

        # Process later if input file list
        if flag == 'input_file_list':
            continue

        # Check if argument is a rosetta script variable
        if flag == '-parser:script_vars':
            args.append(flag)
            script_vars = flags_dict[flag]
            assert( not isinstance(script_vars, basestring) )
            for varstring in script_vars:
                assert( '=' in varstring )
                name, value = varstring.split('=')
                if os.path.isfile(value):
                    value = copy_file_helper( os.path.abspath(value) )
                args.append( '%s=%s' % (name, value) )
            continue

        # Add only the value if NOAPPEND, otherwise also add the key here
        if not flag.startswith('NOAPPEND'):
            args.append(flag)

        # Check if argument is a file or directory
        # If so, copy to temporary data directory
        value = str(flags_dict[flag])
        if os.path.isfile(value):
            value = copy_file_helper( os.path.abspath(value) )

        elif os.path.isdir(value):
            original_dir = os.path.abspath(value)
            new_dir = os.path.join(tmp_data_dir, os.path.basename(original_dir))
            shutil.copytree(original_dir, new_dir)
            value = os.path.relpath(new_dir, tmp_output_dir)
            if verbosity >= 1:
                print 'Copied dir to local scratch:', os.path.basename(original_dir)

        if not isinstance(flags_dict[flag], basestring):
            # Is a list
            for list_item in flags_dict[flag]:
                args.append(list_item)
        else:
            args.append(value)

    if 'input_file_list' in flags_dict:
        input_list_file = os.path.join(tmp_data_dir, 'structs.txt')
        tmp_pdb_dir = os.path.join(tmp_data_dir, 'input_list_pdbs')
        os.mkdir(tmp_pdb_dir)
        args.append('-l')
        args.append( os.path.relpath(input_list_file, tmp_output_dir) )
        f = open(input_list_file, 'w')
        pdbs = sorted(flags_dict['input_file_list'])
        for i, input_pdb in enumerate(pdbs):
            inner_tmp_pdb_dir = tempfile.mkdtemp(prefix='pdb_dir_', dir=tmp_pdb_dir)
            new_input_file = os.path.join(inner_tmp_pdb_dir, os.path.basename(input_pdb))
            shutil.copy(input_pdb, new_input_file)
            if new_input_file.endswith('.gz'):
                new_input_file = unzip_file(new_input_file)
            f.write( os.path.abspath(new_input_file) )
            f.write('\n')
        f.close()

    args.extend(generic_rosetta_args)

    job_dir_path = os.path.join(os.getcwd(),str(job_dir))
    outfile_path = os.path.join(tmp_output_dir,'rosetta.out')
    
    if verbosity >= 1:
        print 'Args:'
        print args
        print ''

    # Run Rosetta
    rosetta_outfile = open(outfile_path, 'w')

    rosetta_env = os.environ.copy()

    # Add to LD path (needed for Rosetta MySQL)
    add_extra_ld_path = #$#add_extra_ld_path#$#
    if add_extra_ld_path:
        if "LD_LIBRARY_PATH" in rosetta_env:
            rosetta_env["LD_LIBRARY_PATH"] = "#$#extra_ld_path#$#:" + rosetta_env["LD_LIBRARY_PATH"]
        else:
            rosetta_env["LD_LIBRARY_PATH"] = "#$#extra_ld_path#$#"

    # Check that Rosetta binary exists
    assert( os.path.isfile( args[0] ) )
    # Output to file
    rosetta_process = subprocess.Popen(args, stdout=rosetta_outfile, stderr=subprocess.STDOUT, close_fds = True, cwd=tmp_output_dir, env=rosetta_env )

    # Output to terminal for debugging
    # rosetta_process=subprocess.Popen(args, close_fds = True, cwd=tmp_output_dir, env=rosetta_env )

    return_code = rosetta_process.wait()

    rosetta_outfile.close()

    if verbosity>=1:
        print 'Rosetta return code:', return_code, '\n'

    if zip_rosetta_output and os.path.isfile(outfile_path):
        zip_file(outfile_path)

    if not os.path.isdir(job_dir_path):
        print 'Making new job output directory: ', job_dir_path
        os.makedirs(job_dir_path)

    # Move files to job_dir from scratch dir recursively
    def move_file_helper(d, copy_to_dir):
        for x in os.listdir(d):
            if x.startswith('tmppdbnocopy'):
                # Skip temporary output PDB files
                continue
            x=os.path.abspath(os.path.join(d, x))
            if os.path.isfile(x):
                if x.endswith('.pdb'):
                    x=zip_file(x)
                shutil.copy(x, copy_to_dir)
                os.remove(x)
            elif os.path.isdir(x):
                new_copy_to_dir = os.path.join(copy_to_dir, x) 
                if not os.path.isdir( new_copy_to_dir ):
                    os.makedirs(new_copy_to_dir)
                move_file_helper( os.path.join(d, x), new_copy_to_dir )

    move_file_helper(tmp_output_dir, job_dir_path)
            
    # Delete temporary directories
    shutil.rmtree(tmp_output_dir)
    shutil.rmtree(tmp_data_dir)

    # Check if on SGE to move special output files and calculate RAM usage                                                             
    ram_usage = None
    ram_usage_type = None
    if run_on_sge and move_output_files:
        try:
            # Encase this in try block in case script_name is wrong
            shutil.move("%s.o%d.%d" % (script_name,job_id,sge_task_id), job_dir_path)
            shutil.move("%s.e%d.%d" % (script_name,job_id,sge_task_id), job_dir_path)
        except IOError:
            print 'Failed moving script files, check stored name'

        qstat_p = subprocess.Popen(['/usr/local/sge/bin/linux-x64/qstat', '-j', '%d' % job_id],
                                   stdout=subprocess.PIPE)
        out, err = qstat_p.communicate()

        # Import slow re module only if we made it this far
        import re
        for line in out.split(os.linesep):
            m = re.match('(?:usage\s+%d[:]\s+.*?)(?:maxvmem[=])(\d+[.]\d+)([a-zA-Z]+)(?:.*?)' % sge_task_id, line)
            if m:
                ram_usage = float(m.group(1))
                ram_usage_type = m.group(2)

    time_end = roundTime()
    if verbosity>=1:
        print 'Ending time:', time_end
        print "Elapsed time:", time_end-time_start
        if ram_usage:
            print 'Max virtual memory usage: %.1f%s' % (ram_usage, ram_usage_type)

    return time_end

def run_local():
    global job_dict

    # Only import multiprocessing here because it may not be available in old cluster python environments
    from multiprocessing import Pool
    import multiprocessing

    # Integer argument allows number of jobs to run to be specified
    # Jobs will be run single-threaded for debugging
    if len(sys.argv) > 1:
        jobs_to_run = int(sys.argv[1])
    else:
        jobs_to_run = None

    class MultiWorker:
        def __init__(self, task, func):
            self.reporter = Reporter(task)
            self.func = func
            self.pool = Pool()
            self.number_finished = 0
        def cb(self, time_return):
            self.number_finished += 1
            self.reporter.report(self.number_finished)
        def addJob(self,argsTuple):
            if jobs_to_run:
                # Testing
                self.cb( self.func(argsTuple[0], argsTuple[1], argsTuple[2], argsTuple[3], argsTuple[4]) )
            else:
                # Multiprocessing
                self.pool.apply_async(self.func, argsTuple, callback=self.cb)

        def finishJobs(self):
            self.pool.close()
            self.pool.join()
            self.reporter.done()

    worker = MultiWorker('running script locally', run_single)

    if os.path.isfile(job_pickle_file):
        p = open(job_pickle_file, 'r')
        job_dict = pickle.load(p)
        p.close()

    num_jobs = len(job_dict.keys())
    worker.reporter.set_total_count(num_jobs)

    for i in xrange(0, num_jobs): # Manually specify which jobs to run here, if you desire. Or pass as argument
        if jobs_to_run and i >= jobs_to_run:
            break
        worker.addJob( (i, local_rosetta_bin, local_rosetta_db, local_scratch_dir, 1) )

    worker.finishJobs()

def run_cluster():
    # Determine which tasks this job will run
    tasks_to_run = []
    starting_task_id = tasks_per_process * ( int(sge_task_id) - 1 )
    tasks_to_run = range(
        starting_task_id,
        min(starting_task_id + tasks_per_process, total_number_tasks)
    )

    if len(tasks_to_run) == 0:
        print 'ERROR: No tasks to run!!!'

    for i, task_id in enumerate(tasks_to_run):
        if len(tasks_to_run) != 1:
            print 'Running subtask %d (%d of %d)' % (task_id, i+1, len(tasks_to_run))
        if i+1 == len(tasks_to_run):
            run_single(task_id, cluster_rosetta_bin, cluster_rosetta_db, scratch_dir='/scratch', move_output_files=True)
        else:
            run_single(task_id, cluster_rosetta_bin, cluster_rosetta_db, scratch_dir='/scratch')

def zip_file(file_path):
    if os.path.isfile(file_path):
        f_in = open(file_path, 'rb')
        f_out_name = file_path + '.gz'
        f_out = gzip.open(f_out_name, 'wb')
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()
        os.remove(file_path)
        return f_out_name

def unzip_file(file_path):
    if os.path.isfile(file_path) and file_path.endswith('.gz'):
        f_in = gzip.open(file_path,'rb')
        f_out_name = file_path[:-3]
        f_out = open(f_out_name, 'wb')
        f_out.writelines(f_in)
        f_out.close()
        f_in.close()
        os.remove(file_path)
        return f_out_name

if __name__=='__main__':
    if run_locally:
        run_local()
    else:
        run_cluster()
