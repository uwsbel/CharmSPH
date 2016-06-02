#!/bin/bash


# Slurm Command: 
#./charmrun +p378 ./charmsph +ppn 63 +commap 0 +pemap 1-63 -x 1.6 -y 1.6 -z 1.6 -w 0
 

#!/bin/bash


maxPPN=32
numNodes=4 # Number of physical nodes
N=4 # Number of Charm++ nodes per physical node
let n="$numNodes * $N" # Total number of charm++ nodes
let d="$maxPPN / $N" # Total number of PE's per charm++ node
let ppn="$d - 1"
commThreads=1 # Number of communication threads (usually 1)

# SimParams
csm=4
X=1.6
Y=1.6
Z=1.6
H=0.025
W=0
SimID="n"$n"_N"$N"_csm"$csm"_h"$H
OutputDir="output/"$SimID"/"


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
	    $n $N $d $commThreads $ppn $numNodes $maxPPN > "output/"$SimID"/EnvParams.json"


#PBS -l nodes=1:ppn=32:xe
#PBS -l walltime=00:01:00
#PBS -N $SimID

#PBS -e $OutputDir$PBS_JOBID.err
#PBS -o $OutputDir$PBS_JOBID.out

## aprun -n 16 -N 4 -d 8 ./charmsph +ppn 7 +commap 0,8,16,24 +pemap 1-7:9-15:17-23:25-31
## -n = Total number of Charm++ Nodes
## n/N = Total number of Physical Nodes
## -N = Number of Charm++ nodes per physicalnode
## -d = Total number of processing elements
## +ppn = Number computation processing elements, 
## d - ppn = Number of communication elements
## +commap = Where to map communication elements
## +pemat = Where to map processing elements


# aprun -n $n -N $N -d $d ./charmsph +ppn $ppn +commap 0 +pemap 1-31 -x $X -y $Y -z $Z -h $H -w $W -id $SimID