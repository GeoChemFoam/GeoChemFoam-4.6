#!/bin/bash

# ###### DO NOT MAKE CHANGES FROM HERE ###################################

set -e

NP=8

MPIRUN=mpirun 

# Decompose
echo -e "DecomposePar"
decomposePar > decomposePar.out

# Run simpleFoam in parallel
echo -e "Run simpleFoam in parallel"
$MPIRUN -np $NP simpleDBSFoam -parallel  > simpleFoam.out

for i in processor*; do mv $i/[1-9]*/U $i/0/.; done
for i in processor*; do mv $i/[1-9]*/p $i/0/.; done
for i in processor*; do mv $i/[1-9]*/phi $i/0/.; done

rm -rf processor*/[1-9]*

echo -e "Run heatTransportSimpleFoam in parallel"
$MPIRUN -np $NP heatTransportSimpleFoam -parallel  > heatTransportSimpleFoam.out 

echo -e "reconstructPar"
reconstructPar -latestTime > reconstructPar.out

cp [1-9]*/T 0/.

rm -rf [1-9]*

rm -rf processor*

echo -e "processHeatTransfer"
processHeatTransfer > processHeatTransfer.out




