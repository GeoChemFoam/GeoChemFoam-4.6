#!/bin/bash


set -e

echo -e "decomposePar"
decomposePar > decomposePar.out
echo -e "Run interTransferFoam in parallel"
mpiexec -np 8 interTransferFoam -parallel > interTransferFoam
echo -e "reconstructPar"
reconstructPar > reconstructPar.out
rm -rf processor*
