# CharmSPH
Charm++ implementation of Smoothed Particle Hydrodynamics for distributed memory architectures

## Hybrid decomposition
Cells are a dense 3D chare array that represent a spatial decomposition of the 3D simulation space. They are responsible for sending positions and collecting the forces for their atoms. Computes, on the other hand, form a sparse 6-dimensional array of chares. Such a representation makes it convenient for a pair of cells with coordinates (x1, y1, z1) and (x2, y2, z2) to use a compute with coordinates (x1, y1, z1, x2, y2, z2) to calculate forces for their atoms.

## CharmSPH Development tips
1. When developing the SPH side of CharmSPH you might want to make the code serial so it is easier to debug and track. In order to do this go to `Cell.cc` and look for the `createComputes()` function. Here you will have to uncomment and comment a few sections of the code. The comments in the file indicate what to comment and uncomment.

## Charm++ Development Workflows
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
