#!/bin/bash

# Load user environment variables 
source ~/.bashrc

NP=8

source $OF4X_DIR/OpenFOAM-4.x/etc/bashrc WM_LABEL_SIZE=64 FOAMY_HEX_MESH=yes

MPIRUN=mpirun 

echo -e "make stl"
cd constant/triSurface
tar -xf Ketton.raw.tar.gz
python raw2stl.py
rm Ketton.raw
cd ../..


# Create background mesh
echo -e "Create background mesh"
blockMesh  > blockMesh.out

# Decompose background mesh
echo -e "Decompose background mesh"
decomposePar > decomposeBlockMesh.out
rm -rf processor*/0/*

# Run snappyHexMesh in parallel
echo -e "Run snappyHexMesh in parallel"
$MPIRUN -np $NP snappyHexMesh -overwrite -parallel  > snappyHexMesh.out

echo -e "splitMeshRegions and decompose fields"
reconstructParMesh -constant > reconstructParMesh.out
$MPIRUN -np $NP splitMeshRegions -largestOnly -overwrite -parallel > splitMeshRegions.out

# renumberMesh
echo -e "renumberMesh"
$MPIRUN -np $NP renumberMesh -overwrite -parallel > renumberMesh.out


echo -e "checkMesh"
$MPIRUN -np $NP checkMesh -parallel > checkMesh.out

echo -e "transformPoints" 
transformPoints -scale '(5.3e-6 5.3e-6 5.3e-6)' > transformPoints.out

rm -rf processor*
