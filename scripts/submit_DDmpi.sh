#!/bin/bash

# change directory to work directory
#$ -cwd

# define job log file
#$ -o parallel_job_penode.joblog.$JOB_ID

# Merged job error file error   with job log
#$ -j y

# If number of nodes is not passed to qsub, uncomment following line and define number of nodes as needed (e.g. 5)
# #$ -pe node* 5

# Request resources per node:
# h_data = minimum memory per core - use at least 4Gb
# h_rt = time in HHH:MM:SS
# exclusive = reserves the entire node (needed for openmp)
# highp = add highp to run on sponsor nodes
#$ -l h_data=4096M,h_rt=8:00:00,exclusive

#  Email address to notify
#$ -M $USER@mail
#  Notify at beginning and end of job
#$ -m bea
#  Job is not rerunable
#$ -r n

# echo info on joblog:
echo "Job $JOB_ID started on:   "` hostname -s `
echo "Job $JOB_ID started on:   "` date `
echo " "
echo "PE_HOSTFILE: "
cat $PE_HOSTFILE | awk '{print $1" "$2}'


# set the environment
. /u/local/Modules/default/init/modules.sh
module load gcc/4.7.2  > /dev/null 2>&1
module load openmpi > /dev/null 2>&1 

# run the program
$MPI_BIN/mpiexec -pernode ./DDmpi >& parallel_job_penode.output.$JOB_ID