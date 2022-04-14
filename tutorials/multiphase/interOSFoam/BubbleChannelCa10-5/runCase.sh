#!/bin/bash

set -e

echo "decomposePar"
decomposePar > decomposePar.out
echo "Run interOSFoam in parallel"
mpiexec -np 8 interOSFoam  -parallel > interOSFoam.out
echo "reconstructPar"
reconstructPar > reconstructPar.out
rm -rf proc*
echo "write capillary relaxation steps"
./writeStep.sh
