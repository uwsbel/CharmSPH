#!/bin/sh

# When you have output data in a remote directory and you which to visualize that data (using 
# paraview GUI) or use it locally you might want to keep your data directories in sync. The
# following script allows you to do that.
# For ease of usage I would add the following line to you .bashrc or .bash_profile files, so you 
# can call this script from any folder by just typing `sd`:
#   alias sd='/PATH/TO/remote-syncdata.sh'
# 


# If you have ssh without password setup
USER= # USERNAME TO REMOTE HOST GOES HERE
HOST= # REMOTE HOST ALIAS GOES HERElagrange

# If you don't have ssh without password setup. I recommend setup the ssh without password so the 
# sync process means only typing rs and not typing rs and then your password.
USER2= # USERNAME TO REMOTE HOST GOES HERE
HOST2= # URL TO HOST GOES HERE


# Path to repository in your local machine
LOCAL_DIR_CharmSPH=/Users/felipegb94/sbel/repos/CharmSPH/data
# Path to repository in your euler account
REMOTE_DIR_CharmSPH=/home/felipegb94/repositories/CharmSPH/data


rsync --verbose --recursive --times \
	--exclude ".*" --exclude "Debug/" \
	$USER@$HOST:$REMOTE_DIR_CharmSPH $LOCAL_DIR_CharmSPH/



