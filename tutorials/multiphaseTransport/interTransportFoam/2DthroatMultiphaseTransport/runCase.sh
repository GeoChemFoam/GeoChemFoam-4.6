#!/bin/bash

set -e

cp 0_org/alpha1 0/alpha1
cp 0_org/T 0/T
echo -e "blockMesh"
blockMesh > blockMesh.out
echo -e "setFields"
setFields > setFields.out
echo "decomposePar"
decomposePar > decomposePar.out
echo "Run interTransportFoam in parallel"
mpiexec -np 4 interTransportFoam -parallel > interTransportFoam.out
echo "reconstructPar"
reconstructPar > reconstructPar.out
rm -rf process*

echo "process phase concentration"
processPhaseConcentration > processPhaseConcentration.out
echo "process interface transfer"
processInterfaceTransfer > processInterfaceTransfer.out 

