#!/bin/bash

set -e

cp constant/dynamicMeshDictAMR constant/dynamicMeshDict

python << END
import os
nRef=200
maxRef=1
lowRef=0.01
upRef=0.99

os.system('sed -i "s/nRef/'+str(nRef)+'/g" constant/dynamicMeshDict')
os.system('sed -i "s/lowRef/'+str(lowRef)+'/g" constant/dynamicMeshDict')
os.system('sed -i "s/upRef/'+str(upRef)+'/g" constant/dynamicMeshDict')
os.system('sed -i "s/refLevel/'+str(maxRef)+'/g" constant/dynamicMeshDict')
END

cp system/controlDictRun system/controlDict
cp system/fvSolutionRun system/fvSolution
mpiexec -np 24 reactiveTransportDBSFoam -parallel
reconstructParMesh -noZero
rm -rf processor*
processSolidArea


