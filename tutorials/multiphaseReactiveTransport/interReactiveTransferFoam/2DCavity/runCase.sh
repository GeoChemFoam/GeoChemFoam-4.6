#!/bin/bash

set -e

echo -e "decomposePar"
decomposePar > decomposePar.out
echo -e "Run interReactiveTranferFoam in parallel"
mpiexec -np 8 interReactiveTransferFoam -parallel
echo -e "reconstructPar"
reconstructPar > reconstructPar.out
rm -rf processor*
