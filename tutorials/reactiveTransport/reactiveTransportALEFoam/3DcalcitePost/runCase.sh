#!/bin/bash




# ###### DO NOT MAKE CHANGES FROM HERE ###################################



set -e

cp system/fvSolutionRun system/fvSolution

rm -rf polyMesh_old
cp -r constant/polyMesh polyMesh_old

cp -r ../3DcalcitePost ../temp

python << END
import os


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


a=0;
s=str(a)
removeOld = 0

for n in range(0,4000):
  os.system('cp system/controlDictRun system/controlDict')
  os.system('sed -i "s/var/'+str(n)+'/g" system/controlDict') 
  while a<n:
    os.system('decomposePar')
    os.system('mpiexec -np 24 reactiveTransportALEFoam -parallel')
    os.system('reconstructPar -latestTime' )
    os.system('rm -rf processor*')
    if removeOld == 1:
        os.system('rm -rf '+s)
    for directories in os.listdir(os.getcwd()): 
      if (is_number(directories)):
        if (float(directories)>a):
          a=float(directories)
          s=directories
    os.system('rm polyMesh_old/points')
    os.system('cp '+s+'/polyMesh/points polyMesh_old/.')
    os.system('cp polyMesh_old/* '+s+'/polyMesh/.')
    os.system('cp -r '+s+' ../temp/.')
    os.system( './remesh.sh') 
    os.system('./calculateFields.sh')
    os.system('cd ../3DcalcitePost')
    os.system('rm -rf '+s)
    os.system('mv ../temp/0 ./'+s)
    os.system('cp '+s+'/polyMesh/* polyMesh_old/.')
    if a<n:
      removeOld=1
    else:
      removeOld=0
END

processPoroSurf

rm -rf ../temp
rm -rf polyMesh_old
