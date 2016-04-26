#!/bin/bash

#SBATCH -N 1 -n 8 -t 0-5:00:00
#SBATCH -o charm-run-log.o%j

cd /srv/home/felipegb94/repositories/CharmSPH
make clean
make


./charmrun +p8 ./charmsph -x 2 -y 2 -z 2 -t 20000 -h 0.025 -dt 0.0002 -wb 1 -csm 4 -wp 250 -lbp 200 -wb 0
