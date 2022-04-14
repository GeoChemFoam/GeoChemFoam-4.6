#!/bin/bash




# ###### DO NOT MAKE CHANGES FROM HERE ###################################


set -e

#cp -rf constant/polyMesh 0/.
python << END
import os
import h5py
import numpy as np

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

ndir = 0

print("create voxelized mesh")
for dir in os.listdir(os.getcwd()): 
  if (is_number(dir)):
    ndir +=1
    print(dir)
    os.system('cp -r '+dir+' ../temp/.')
    os.system('cp ../temp/system/snappyHexMeshDictNoSnap ../temp/system/snappyHexMeshDict')
    os.system('./remesh.sh ') 
    os.system('mapFields ../3DcalcitePost -case ../temp -sourceTime '+dir+' > logMapField.out')
    print('rm -rf '+dir)
    os.system('rm -rf '+dir)
    os.system('mv ../temp/0 ./'+dir)


print("processMeshCellCenters")
os.system('processMeshCellCenters > log.processCellCenters')


END


