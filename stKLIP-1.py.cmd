#!/bin/csh -f
#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = /u/home/b/blewis34/job-outputs/stKLIP.joblog.$JOB_ID
#$ -o /u/home/b/blewis34/job-outputs/stKLIP.joblog.$JOB_ID
#  error           = Merged with joblog
#$ -j y
#  The following items pertain to the user program
#  user program    = /u/home/b/blewis34/speckle-stats/stKLIP-H2.py
#  arguments       = 
#  program input   = Specified by user program
#  program output  = Specified by user program
#  Resources requested
#
#$ -l h_data=16000M,h_rt=24:00:00,exclusive
# #
#  Name of application for log
#$ -v QQAPP=job
#  Email address to notify
#$ -M blewis34@mail
#  Notify at beginning and end of job
#$ -m bea
#  Job is not rerunable
#$ -r n
#
# Initialization for serial execution
#
  unalias *
  set qqversion = 
  set qqapp     = "job serial"
  set qqidir    = /u/home/b/blewis34/speckle-stats
  set qqjob     = stKLIP-H2.py
  set qqodir    = /u/home/b/blewis34/speckle-stats
  cd     /u/home/b/blewis34/speckle-stats
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "UGE job for stKLIP-H2.py"
  echo ""
  echo "  script directory:"
  echo "    "/u/home/b/blewis34/speckle-stats
  echo "  Submitted to UGE:"
  echo "    "$qqsubmit
  echo "  SCRATCH directory:"
  echo "    "$qqscratch
#
  echo ""
  echo "script started on:   "` hostname -s `
  echo "script started at:   "` date `
  echo ""
#
# Run the user program
#
  source /u/local/Modules/default/init/modules.csh
  module load intel/13.cs
  module load python/3.7.2
  source /u/home/b/blewis34/.bash_profile
#
  echo stKLIP-H2.py "" \>\& stKLIP.output.$JOB_ID
  echo ""

#PUT SIM PARAMS HERE
  /usr/bin/time python3 /u/home/b/blewis34/speckle-stats/stKLIP-H2.py MEDIS_1sec_Aug2020.h5 2,3,4,5,7,10 1,2,3,4,5,7,10,15,20 85,135,125,175 >& /u/home/b/blewis34/job-outputs/stKLIP.output.$JOB_ID
#
  echo ""
  echo "script finished at:  "` date `
#
# Cleanup after serial execution
#
  source /u/local/bin/qq.sge/qr.runtime
#
  echo "-------- /u/home/b/blewis34/job-outputs/stKLIP.joblog.$JOB_ID --------" >> /u/local/apps/queue.logs/job.log.serial
  if (`wc -l /u/home/b/blewis34/job-outputs/stKLIP.joblog.$JOB_ID  | awk '{print $1}'` >= 1000) then
	head -50 /u/home/b/blewis34/job-outputs/stKLIP.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
	echo " "  >> /u/local/apps/queue.logs/job.log.serial
	tail -10 /u/home/b/blewis34/job-outputs/stKLIP.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  else
	cat /u/home/b/blewis34/job-outputs/stKLIP.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  endif
  exit (0)
