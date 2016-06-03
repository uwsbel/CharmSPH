#!/bin/bash

# SimParams
csm=4
H=0.025
X=1.6
Y=1.6
Z=1.6
W=0
WP=20
T=101

#SBATCH --reservation=felipegb94_4
#SBATCH --ntasks=8 --ntasks-per-node=2
#SBATCH --cpus-per-task=32
#SBATCH --gres=cputype:amd:1 
#SBATCH -t 0-00:30:00

numTasks=8
cpusPerNode=64
cpusPerTask=32
numCommThreads=1
numNodes=4 # Number of physical nodes
let numProcs="$cpusPerNode * $numNodes" # Total number of processors available
let ppn="$cpusPerTask - $numCommThreads"
let N="$cpusPerNode / $cpusPerTask"

SimID="n"$numTasks"_N"$N"_csm"$csm"_h"$H
OutputDir="output/"$SimID"/"

cd /srv/home/felipegb94/CharmSPH
mkdir "output"
mkdir $OutputDir

printf "{\
		\"numCharmNodes\":  %d,\n\
	    \"numCharmNodesPerNode\": %d,\n\
	    \"numPEPerCharmNode\": %d,\n \
	    \"numCommThreads\":  %d,\n\
	    \"ppnPerCharmNode\":  %d,\n\
	    \"numNodes\":  %d,\n\
	    \"maxPPN\": %d\n\
	    }\n"\
	    $numTasks $N $cpusPerTask $numCommThreads $ppn $numNodes $cpusPerNode > "output/"$SimID"/EnvParams.json"


./charmrun +p$numProcs ./charmsph +ppn $cpusPerTaks +commap 0,32 +pemap 1-31:33-63 \
	-x $X -y $Y -z $Z -h $H -w $W -id $SimID -wp $WP -t $T


