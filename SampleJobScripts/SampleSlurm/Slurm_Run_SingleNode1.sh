#!/bin/bash

#SBATCH -N 1 -n 16 -t 0-4:00:00
#SBATCH -o charm-run-log.o%j

cd /srv/home/felipegb94/repositories/CharmSPH
make clean
make


./charmrun +p16 ./charmsph -x 1.5 -y 1.5 -z 1.5 -t 20000 -h 0.025 -dt 0.0002 -csm 4 -wp 200 -lbp 200 -wb 0 -mv 1 -w 0

