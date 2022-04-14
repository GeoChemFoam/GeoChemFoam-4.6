#!/bin/bash

# Load user environment variables 
source $HOME/.bashrc

NP=24

source $OF4X_DIR/OpenFOAM-4.x/etc/bashrc WM_LABEL_SIZE=64 WM_COMPILER_TYPE=ThirdParty FOAMY_HEX_MESH=yes

MPIRUN=mpirun 

Nx=134
Ny=75
Nz=10

Nlevel=1

cp system/snappyHexMeshDictLMR system/snappyHexMeshDict 
sed -i "s/nlev/$Nlevel/g" system/snappyHexMeshDict

cp constant/triSurface/calcitePost_org.stl constant/triSurface/calcitePost.stl

rm -rf 0
cp system/blockMeshDictRun system/blockMeshDict
sed -i "s/nx/$Nx/g" system/blockMeshDict
sed -i "s/ny/$Ny/g" system/blockMeshDict
sed -i "s/nz/$Nz/g" system/blockMeshDict

# Create background mesh
echo -e "Create background mesh"
blockMesh  > blockMesh.out

# Decompose background mesh
echo -e "Decompose background mesh"
decomposePar > decomposeBlockMesh.out

# Remove fields on this stage
rm -rf ./processor*/0/*

# Run snappyHexMesh in parallel
echo -e "Run snappyHexMesh in parallel"
$MPIRUN -np $NP snappyHexMesh -overwrite -parallel  > snappyHexMesh.out

# reconstruct mesh to fields decomposition
echo -e "splitMeshRegions and decompose fields"
reconstructParMesh -constant > reconstructParMesh.out
#decomposePar -fields  > decomposeFields.out
$MPIRUN -np $NP splitMeshRegions -largestOnly -overwrite -parallel > splitMeshRegions.out

# renumberMesh
echo -e "renumberMesh"
$MPIRUN -np $NP renumberMesh -overwrite -parallel > renumberMesh.out


echo -e "checkMesh"
$MPIRUN -np $NP checkMesh -parallel > checkMesh.out


rm -rf *.out processor*

cp -r 0_org 0
