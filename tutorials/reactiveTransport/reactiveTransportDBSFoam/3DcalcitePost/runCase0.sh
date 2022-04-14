#!/bin/bash



set -e

cp system/controlDict0 system/controlDict
cp system/fvSolution0 system/fvSolution

mpiexec -np 24 reactiveTransportDBSFoam -parallel
reconstructPar -latestTime

rm -rf processor*/0
rm -rf processor*/0.5/uniform
for i in processor*; do mv "$i/0.5" "$i/0"; done

rm -rf 0
rm -rf 0.5/uniform
mv 0.5 0

