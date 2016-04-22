#!/bin/bash

#SBATCH -N 4 --ntasks-per-node=32 --gres=cputype:amd:1 -t 0-01:00:00
#SBATCH -o charm-run-log.o%j

cd /home/felipegb94/repositories/CharmSPH
make clean
make


./charmrun +p128 ./charmsph 2 2 2 20000