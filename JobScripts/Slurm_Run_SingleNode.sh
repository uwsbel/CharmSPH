#!/bin/bash

#SBATCH -N 1 -n 32 -t 0-20:00:00
#SBATCH -o charm-run-log.o%j

cd /srv/home/felipegb94/repositories/CharmSPH
make clean
make


./charmrun +p32 ./charmsph -x 2 -y 2 -z 2 -t 20000 -h 0.02 -dt 0.0002 -wb 1 -csm 4 -wp 250 -lbp 200
./charmrun +p16 ./charmsph -x 2 -y 2 -z 2 -t 20000 -h 0.02 -dt 0.0002 -wb 1 -csm 4 -wp 250 -lbp 200
./charmrun +p8 ./charmsph -x 2 -y 2 -z 2 -t 20000 -h 0.02 -dt 0.0002 -wb 1 -csm 4 -wp 250 -lbp 200

