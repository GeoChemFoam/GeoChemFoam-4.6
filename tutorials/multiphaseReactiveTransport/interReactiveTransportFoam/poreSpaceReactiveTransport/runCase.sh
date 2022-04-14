#!/bin/bash

set -e


cp system/controlDictRun system/controlDict
echo "decomposePar"
decomposePar > decomposePar.out
echo "Run interReactiveTransportFoam in parallel"
mpiexec -np 8 interReactiveTransportFoam -parallel > interReactiveTransportFoam.out
echo "reconstructPar"
reconstructPar > reconstructPar.out
rm -rf processor*
