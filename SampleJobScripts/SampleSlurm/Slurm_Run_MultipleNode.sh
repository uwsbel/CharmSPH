#!/bin/bash

#SBATCH --reservation=felipegb94_4
#SBATCH --ntasks=6 --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --gres=cputype:amd:1 
#SBATCH -t 0-01:00:00
#SBATCH -o charm-run-log.o%j

cd /srv/home/felipegb94/CharmSPH


./charmrun +p378 ./charmsph +ppn 63 +commap 0 +pemap 1-63 -x 1.6 -y 1.6 -z 1.6 -w 0

