#!/bin/sh

# When you have no way to run a program in your local machine (due to hardware restrictions, 
# amount of data output, etc), but you still wish to develop locally you can use the 'rsync' 
# command to keep a local directory and a remote directory in sync. This script will be useful 
# when you modify a file locally and want the changes to be reflected in the remote copy, you just 
# have to run this script.
# For ease of usage I would add the following line to you .bashrc or .bash_profile files, so you 
# can call this script from any folder by just typing `rs`:
#   alias rs='/PATH/TO/remote-sync.sh'
# Example:
#   alias rs='/Users/felipegb94/sbel/repos/DevelopmentUtils/remote-sync.sh 


# If you have ssh without password setup
USER= # USERNAME TO REMOTE HOST GOES HERE
HOST= # REMOTE HOST ALIAS GOES HERElagrange

# If you don't have ssh without password setup. I recommend setup the ssh without password so the 
# sync process means only typing rs and not typing rs and then your password.
USER2= # USERNAME TO REMOTE HOST GOES HERE
HOST2= # URL TO HOST GOES HERE


# Path to repository in your local machine
LOCAL_DIR_CharmSPH=/Users/felipegb94/sbel/repos/CharmSPH
# Path to repository in your euler account
REMOTE_DIR_CharmSPH=/home/felipegb94/repositories/CharmSPH

# Note 1: If you have --exclude rsync will ignore those folder
rsync --verbose --recursive --times \
    --exclude ".*" --exclude "Debug/" --exclude "build/" --exclude "data/" --exclude "armadillo_bits" --exclude "armadillo"\
        $LOCAL_DIR_CharmSPH/ $USER@$HOST:$REMOTE_DIR_CharmSPH/
