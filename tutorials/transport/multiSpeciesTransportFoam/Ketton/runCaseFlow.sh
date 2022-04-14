#!/bin/bash

# Load user environment variables 
source ~/.bashrc

NP=8

MPIRUN=mpirun 

rm -f *.out

cp -r 0_orig 0

# Decompose
echo -e "DecomposePar"
decomposePar > decomposePar.out

# Run simpleFoam in parallel
echo -e "Run simpleFoam in parallel"
$MPIRUN -np $NP simpleFoam -parallel  > simpleFoam.out

# ReconstructPar
echo -e "reconstructPar"
reconstructPar -latestTime > reconstructPar.out

rm -rf 0
rm -rf processor*
echo -e "processPoroPerm"
processPoroPerm > poroPerm.out

