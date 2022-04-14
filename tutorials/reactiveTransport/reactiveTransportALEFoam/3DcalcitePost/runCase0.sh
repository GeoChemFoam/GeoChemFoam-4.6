#!/bin/bash




# ###### DO NOT MAKE CHANGES FROM HERE ###################################



set -e

cp system/controlDict0 system/controlDict
cp system/fvSolution0 system/fvSolution

decomposePar
mpiexec -np 24 reactiveTransportALEFoam -parallel
reconstructPar
rm -rf processor*

cp 0.5/H+ 0/.
cp 0.5/U 0/.
cp 0.5/p 0/.
cp 0.5/phi 0/.
cp 0.5/pointMotionU 0/.
cp 0.5/cellMotionU 0/.

rm -rf 0.5 

