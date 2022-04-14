#!/bin/bash

# Load user environment variables 
source ~/.bashrc

NP=8

MPIRUN=mpirun 

rm -rf *.out

mv [1-9]* 0
cp 0_orig/A 0/.
cp 0_orig/T 0/.
cp 0_orig/AT 0/.
rm -rf 0/uniform

# Decompose
echo -e "DecomposePar"
decomposePar > decomposePar.out

# Run reactiveTransportFoam in parallel
echo -e "Run reactiveTransportFoam in parallel"
$MPIRUN -np $NP reactiveTransportFoam -parallel  > reactiveTransportFoam.out

# ReconstructPar
echo -e "reconstructPar"
reconstructPar > reconstructPar.out

rm -rf processor*

echo "process concentration"
processConcentration > processConcentration.out
