#!/bin/bash

set -e

rm -f *.csv
rm -rf 0.* *e-*
rm -rf processor*
rm -f constant/polyMesh/boundary
rm -f constant/polyMesh/faces
rm -f constant/polyMesh/neighbour
rm -f constant/polyMesh/owner
rm -f constant/polyMesh/points
rm -f *.out
