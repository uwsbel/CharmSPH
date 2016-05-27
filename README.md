# CharmSPH
Charm++ implementation of Smoothed Particle Hydrodynamics for distributed memory architectures

## Table of Contents
1. [Setting up CharmSPH & Charm++](#charmsetup)
	* Setup in Euler
	* [Setup in Blue Waters](#bwsetup)   
2. [Running CharmSPH](#runningcharmsph)
3. [Development Tips](#developmenttips)
4. [Parameters](#parameters)
5. [Outputs](#outputs)
6. [Algorithm Details](#algorithmdetails)
7. [Adding a feature to CharmSPH](#addfeature)
8. [Development Tips](#developmenttips)

<a name="charmsetup">
## Setting up CharmSPH and Charm++
</a>

### In Euler

```
	module load openmpi
```

```
	cd charm
 ./build charm++ mpi-linux-x86_64   mpicxx smp  -j16  --with-production --enable-tracing	
```

<a name="bwsetup">
### In Blue Waters
</a>

1. Download Charm++. CharmSPH was developed and tested using Charm++ 6.6.1, so we recommend going with that version. Run  the following commands on your home directory:

	```
	wget http://charm.cs.illinois.edu/distrib/charm-6.6.1.tar.gz
	tar -zxvf charm-6.6.1.tar.gz
	cd charm
	```
	
2. Select the gnu compiler environment and load a hugepages module. Apparently there are some issues in craype-hugepages8M, so the Blue Waters staff recommended me to compile charm++ with craype-hugepages2M:

	```
	module swap PrgEnv-cray PrgEnv-gnu
	module add craype-hugepages2M
	```

3. Blue Charm++. Refer to [Blue Waters Charm++ Reference Page](https://bluewaters.ncsa.illinois.edu/charm) for more information about the following commands:

	```
	./build charm++ gni-crayxe smp --with-production --enable-tracing -j8
	```
This step takes a few minutes.

4. Create an environment variable `CHARMDIR` that points to Charm++, charmc compiler. This is needed so no changes to the `Makefile` are needed. To do this open your `.bashrc` file, which should be on your home directory

	```
	export CHARMDIR=~/charm/
	```
	
5. Try out one of their examples and see if its working. But at this point you should be set to run an Charm++ program in BW.


<a name="runningcharmsph">
## Running CharmSPH
</a>
### Running on a Slurm Cluster
1. Create an interactive job. In slurm this would be: 
```
srun -u bash -i
```
This gives us a single job with a single core. To get AMD only nodes use the following command instead: 
```
srun -u --gres=cputype:amd:1 bash -i
```
2. Now we want to allocate the resources needed to actually run the Charm program. These would be multiple nodes with an specific number of cores per node. To do this
```
salloc -N 4 --ntasks-per-node=64 --gres=cputype:amd:1 bash -i
```
Here we allocated 4 nodes with 64 threads per node (AMD nodes have 32 physical cores, and 64 logical cores).


### Running CharmSPH on Blue Waters

1. Login into your blue waters account. 
2. 

<a name="parameters">
## Input Parameters
</a>

```
  /**
   *  Input Parameters
   *    * -x = x domain dimension
   *    * -y = y domain dimension
   *    * -z = z domain dimension
   *    * -t = t is the total number of time steps
   *    * -h = h is the particle interaction radius (cutoff radius)
   *    * -dt = dt delta t at every time step
   *    * -mv = mv is the estimate of the maximum velocity of the particles in the model.
   *    * -wp = Write period. After every wp steps we write output
   *    * -wb = Write boundary. 1 if you want to write the boundary, 0 if not.
   *    * -csm = Cell Size Multiplier. How much should we multiply
   */
```

<a name="outputs">
## Outputs
</a>

At every run certain parameters (look at parameters section) are set, which determine the output format. The following steps happen through the program.

1. output folder is created if it does not exist.
2. Inside the output folder a file with the following format is created:

```
    charmsph_h_CellSizeMult_numCores_dt_t.json
```

This file will contain the simulation parameters used.


<a name="algorithmdetails">
## Algorithm Details
</a>

### Hybrid decomposition
Cells are a dense 3D chare array that represent a spatial decomposition of the 3D simulation space. They are responsible for sending positions and collecting the forces for their atoms. Computes, on the other hand, form a sparse 6-dimensional array of chares. Such a representation makes it convenient for a pair of cells with coordinates (x1, y1, z1) and (x2, y2, z2) to use a compute with coordinates (x1, y1, z1, x2, y2, z2) to calculate forces for their atoms.

<a name="addfeature">
## Add feature to CharmSPH instruction
</a>
Refer to this [tutorial on branches and merging] (https://github.com/Kunena/Kunena-Forum/wiki/Create-a-new-branch-with-git-and-manage-branches) for more info about these instructions. To add a feature to CharmSPH you will have to clone the repository, create a new branch for the feature, complete and test the feature, merge to master, and finally remove the feature branch. To this end you have to execute the following commands:

```
git clone https://github.com/uwsbel/CharmSPH.git
git checkout -b [name_of_your_new_branch]
git push origin [name_of_your_new_branch]
git remote add [name_of_your_remote] 
```

The previous commands will create a new branch from master where you can work on you feature. If you execute `git branch` you will see a list of all the branches you have locally available. If you want to get all branches to your local copy execute `git fetch`. Finally once you are done working and testing the feature you can follow these commands:

```
git merge [name_of_your_remote]/develop
git branch -d [name_of_your_new_branch]
git push origin :[name_of_your_new_branch]
```

The last two command delete the local and remote branch you were working on.

<a name="developmenttips">
## Development Tips
</a>

### Running serially: Single Core, Single Cell & Compute Chare
When developing the SPH side of CharmSPH you might want to make the code serial so it is easier to debug and track. In order to do this go to `Cell.cc` and look for the `createComputes()` function. Here you will have to uncomment and comment a few sections of the code. The comments in the file indicate what to comment and uncomment. Furthermore make sure that the cell size is the same as the domain size. For example if you want a 1x1x1 domain then make sure to hard code these values in `cellSize` at `Main.cc`.

### Development Workflows

Charm++ is an unusual framework to develop for. Setting it up so a developer can code in an IDE is complicated, Windows development is not fully supported and if you are learning Charm++ running a program in a 2-4 cores machine (your laptop) is not the same as running your program in a development node (20-64 cores) or a few nodes. This section of the documentation includes different development setups various lab members use. 

### Felipe's development - OSX + Sublime Text 2(Edit) and Build/Run in Development Node

If you are not proficient in vim, nano, or emacs (like me) or simply cannot live without Sublime Text (like me) and you develop in a Unix-based system then you might want to use my development setup. 

Note: This setup might also work in Windows if you have Cygwin and you download the rsync tool for Cygwin. If you want to set it up let me know and I might be able to help.

My setup relies on the terminal command `rsync`. This basically syncs two directories a local and remote directory, so you can quickly reflect any changes made to your code in the remote copy of your code. Under `DevelopmentUtils/` there is a bash script called `remote-sync.sh`. The script is pretty self explanatory. Here is how my script looks: 

```
#!/bin/sh

# If you have ssh without password setup
USER= # USERNAME TO REMOTE HOST GOES HERE
HOST= # REMOTE HOST ALIAS GOES HERElagrange

# Path to repository in your local machine
LOCAL_DIR_CharmSPH=/Users/felipegb94/sbel/repos/CharmSPH
# Path to repository in your euler account
REMOTE_DIR_CharmSPH=/home/felipegb94/repositories/CharmSPH

# Note 1: If you have --exclude rsync will ignore those folder
rsync --verbose --recursive --times \
    --exclude ".*" --exclude "Debug/" --exclude "build/" --exclude "data/" --exclude "armadillo_bits" --exclude "armadillo"\
        $LOCAL_DIR_CharmSPH/ $USER@$HOST:$REMOTE_DIR_CharmSPH/

```

Once you clone the script under `DevelopmentUtils`, change the variables `USER`, `HOST`, `LOCAL_DIR_CharmSPH` and `REMOTE_DIR_CharmSPH` to the appropriate values and finally add the following line to your `.bashrc` file:

```
alias rs='/PATH/TO/remote-sync.sh'
```

Finally if you don't have ssh without a password you should set it up because it makes the syncing process much nicer. The following steps are the ones I am constantly doing:

1. Edit code.
2. In terminal session number 1 (local session in any directory) type: `rs`
3. In terminal session number 2 (remote session in `CharmSPH` directory) type: `make`
4. Run charm++ code.

The `rs` command only take a long time the first time you run it. After that every time you run it only take a few seconds.


