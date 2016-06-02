#!/bin/bash

echo "Creating Simulation ID.."


# Slurm Command: 
#./charmrun +p378 ./charmsph +ppn 63 +commap 0 +pemap 1-63 -x 1.6 -y 1.6 -z 1.6 -w 0
 



# Torque command: 
#PBS -l nodes=4:ppn=32:xe
#PBS -l walltime=01:00:00
#PBS -N 4_N2_4h
# aprun -n 16 -N 4 -d 8 ./charmsph +ppn 7 +commap 0,8,16,24 +pemap 1-7:9-15:17-23:25-31
# -n = Total number of Charm++ Nodes
# n/N = Total number of Physical Nodes
# -N = Number of Charm++ nodes per physicalnode
# -d = Total number of processing elements
# +ppn = Number computation processing elements, 
# d - ppn = Number of communication elements
# +commap = Where to map communication elements
# +pemat = Where to map processing elements

# +p = numNodes * ppn
maxPPN=32
numNodes=4 # Number of physical nodes
N=4 # Number of Charm++ nodes per physical node
let n="$numNodes * $N" # Total number of charm++ nodes
let d="$maxPPN / $N" # Total number of PE's per charm++ node
commThreads=1 # Number of communication threads (usually 1)

# SimParams
csm=4
x=1.6
y=1.6
z=1.6
h=0.05
w=0

SimID="n"$n"_N"$N"_csm"$csm

mkdir "output"
mkdir "output/"$SimID

printf "{\
		\"numCharmNodes\":  %d,\n\
	    \"numCharmNodesPerNode\": %d,\n\
	    \"numPEPerCharmNode\": %d,\n \
	    \"numCommThreads\":  %d,\n\
	    \"numNodes\":  %d,\n\
	    \"maxPPN\": %d\n\
	    }\n"\
	    $n $N $d $commThreads $numNodes $maxPPN > "output/"$SimID"/EnvParams.json"





# #PBS -e $PBS_JOBID.err
# #PBS -o $PBS_JOBID.out

# . /opt/modules/default/init/bash

# cd $PBS_O_WORKDIR


# aprun -n 8 -N 2 -d 16  ./charmsph +ppn 15 +commap 0,16 +pemap 1-15:17-31




