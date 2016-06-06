#!/bin/bash

#SBATCH --reservation=felipegb94_4
#SBATCH --ntasks=8 --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --gres=cputype:amd:1 
#SBATCH -t 0-01:00:00
#SBATCH -o charm-run-log.o%j

cd /srv/home/felipegb94/CharmSPH


./charmrun +p504 ./charmsph +ppn 63 +commap 0 +pemap 1-63 -x 1.6 -y 1.6 -z 1.6 -h 0.025 -w 0

