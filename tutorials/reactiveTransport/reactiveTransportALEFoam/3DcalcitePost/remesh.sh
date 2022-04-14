#!/bin/bash

# Load user environment variables
source $HOME/.bashrc

#export $GCFOAM_DIR/lib

NP=24

MPIRUN=mpirun

cd ../temp

./removeInternalFaces.sh

surfaceMeshTriangulate -patches '(movingWalls)' constant/triSurface/calcitePost.stl
cp constant/triSurface/calcitePost.stl ../3DcalcitePost/constant/triSurface/calcitePost.stl

source $OF4X_DIR/OpenFOAM-4.x/etc/bashrc WM_LABEL_SIZE=64 WM_COMPILER_TYPE=ThirdParty FOAMY_HEX_MESH=yes


rm -rf 0 0.* [1-9]*

./delete.sh

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

cp -r constant/polyMesh 0/.

