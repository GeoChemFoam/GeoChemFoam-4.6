#!/bin/bash

set -e


cp system/controlDictRun system/controlDict
cp system/fvSolutionRun system/fvSolution
mpiexec -np 24 reactiveTransportDBSFoam -parallel
reconstructPar
rm -rf 0/polyMesh
rm -rf processor*
processSolidArea
