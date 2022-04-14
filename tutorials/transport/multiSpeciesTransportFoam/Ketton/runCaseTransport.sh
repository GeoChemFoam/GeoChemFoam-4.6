#!/bin/bash

# Load user environment variables 
source ~/.bashrc

NP=8

MPIRUN=mpirun 

rm -rf *.out

mv [1-9]* 0
cp 0_orig/Na+ 0/.
cp 0_orig/Cl- 0/.
rm -rf 0/uniform

# Decompose
echo -e "DecomposePar"
decomposePar > decomposePar.out

# Run multiSpeciesTransportFoam in parallel
echo -e "Run multiSpeciesTransportFoam in parallel"
$MPIRUN -np $NP multiSpeciesTransportFoam -parallel  > multiSpeciesTransportFoam.out

# ReconstructPar
echo -e "reconstructPar"
reconstructPar > reconstructPar.out

rm -rf processor*

echo "process concentration"
processConcentration > processConcentration.out
