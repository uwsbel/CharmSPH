
maxPPN=32
numNodes=1 # Number of physical nodes
N=1 # Number of Charm++ nodes per physical node
n=1 # Total charm++ nodes (numNodes * N)
d=16 # total PE's per charm++ node
ppn=16 # PE's for computation
commThreads=0
let p="$d * $n"

csm=4
H=0.05
X=0.8
Y=0.8
Z=0.8
W=1
WB=1
WP=250
T=10000
#SimID="p"$p"_n"$n"_N"$N"_csm"$csm"_h"$H
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
         $numTasks $N $cpusPerTask $numCommThreads $ppn $numNodes $cpusPerNode > "output/"$SimID"/EnvParams.json"

./charmrun +p${p} ./charmsph -x ${X} -y ${Y} -z ${Z} -h ${H} -w ${W} -wb ${WB} -wp ${WP} -t ${T} -csm ${csm}

