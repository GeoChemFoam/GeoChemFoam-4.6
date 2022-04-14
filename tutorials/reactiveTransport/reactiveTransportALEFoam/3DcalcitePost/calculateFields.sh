#!/bin/bash

# Load user environment variables
source $HOME/.bashrc

#export $GCFOAM_DIR/lib

NP=24

MPIRUN=mpirun

cd ../temp

mapFields ../3DcalcitePost -case ../temp -sourceTime latestTime
mv 0/pointMotionU* 0/pointMotionU

./runCase0.sh
