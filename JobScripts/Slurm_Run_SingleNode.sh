#!/bin/bash

#SBATCH -N 1 -n 12 -t 0-00:15:00
#SBATCH -o charm-run-log.o%j

cd /home/felipegb94/repositories/CharmSPH
make clean
make


./charmrun +p12 ./charmsph 1 1 1 1000