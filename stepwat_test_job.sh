#!/bin/bash

#---Example JANUS job, for use with qsub ie: "qsub -q janus-debug testjob.sh" wherein janus-debug is the queue to run it on
#---Make sure to call "use Torque" & "use Moab" beforehand
 
#---Set job parameters first... lines starting with #PBS are treated by bash as comments, but are interpreted by qsub as arguments.

#---Set the name of the job
#PBS -N job_stepwat_test

#---Set a walltime for the job. The time format is HH:MM:SS.
#PBS -l walltime=0:05:00

#---Select 1 node, 1 processor per node.  1 processor total...
#PBS -l nodes=1:ppn=1

#---Join the output and error files
#PBS -j oe

#---The following commands will be executed when the script is run
. /curc/tools/utils/dkinit
use NCAR-Parallel-Intel
use R-2.13

# gets us to the correct working directory...
cd /projects/donovanm/STEPPEWAT_test/testing
pwd

./stepwat -f files.in -s -q -e
