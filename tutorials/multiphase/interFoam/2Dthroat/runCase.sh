#!/bin/bash

set -e

echo "blockMesh"
blockMesh > blockMesh.out
cp 0/alpha1.org 0/alpha1
echo "setFields"
setFields > setFields.out
echo "decomposePar"
decomposePar > decomposePar.out
echo "run interFoam in parallel"
mpiexec -np 4 interFoam  -parallel > interFoam.out
echo "reconstructPar"
reconstructPar > reconstructPar.out
rm -rf proc*
