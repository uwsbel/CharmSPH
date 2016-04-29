#!/bin/bash

#PBS -l nodes=1:ppn=32
#PBS -l walltime=05:00:00
#PBS -N testjob

#PBS -e $PBS_JOBID.err
#PBS -o $PBS_JOBID.out

. /opt/modules/default/init/bash
module swap PrgEnv-cray PrgEnv-gnu
module add craype-hugepages2M

cd $PBS_O_WORKDIR

aprun -n 32 ./charmrun +p32 ./charmsph -x 1.5 -y 1.5 -z 1.5 -t 20000 -h 0.025 -dt 0.0002 -csm 4 -wp 200 -lbp 200 -wb 0 -mv 1 > out.$PBS_JOBID
