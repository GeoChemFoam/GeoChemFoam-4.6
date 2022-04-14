#!/bin/bash

set -e

echo "blockMesh"
blockMesh > blockMesh.out
cp 0/alpha1.org 0/alpha1
echo "setFields"
setFields > setFields.out
echo "decomposePar"
decomposePar > decomposePar.out
echo "run interGCFoam in parallel"
mpiexec -np 4 interGCFoam  -parallel > interGCFoam.out
echo "reconstructPar"
reconstructPar > reconstructPar.out
rm -rf proc*
