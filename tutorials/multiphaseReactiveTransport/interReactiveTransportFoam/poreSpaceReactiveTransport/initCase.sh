#!/bin/bash

set -e

source $OF4X_DIR/OpenFOAM-4.x/etc/bashrc WM_LABEL_SIZE=64 FOAMY_HEX_MESH=yes


cp system/controlDictInit system/controlDict

echo "blockMesh"
blockMesh > blockMesh.out
echo "Run snappyHexMesh in parallel"
snappyHexMesh -overwrite > snappyHexMesh.out
echo "transformPoints"
transformPoints -scale '(0.01 0.01 0.01)' > transformPoints.out
echo "checkMesh"
checkMesh > checkMesh.out

