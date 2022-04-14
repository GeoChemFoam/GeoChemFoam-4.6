#!/bin/bash

set -e

cd constant/triSurface

echo "make hdf5 image"
python raw2hdf5.py
mv Bentheimer400-5mum.hdf5 ../../.

cd ../..

rm -rf 0
cp -r 0_orig 0
mkdir constant/polyMesh
cp system/blockMeshDict constant/polyMesh/.


python << END
import numpy as np
import h5py
import array
import os
f=h5py.File("Bentheimer400-5mum.hdf5",'r')
my_array=f['image'][()]

nx=206
ny=200
nz=200

nlevel=2
nx1 = nx/nlevel
ny1 = ny/nlevel
nz1 = nz/nlevel

ncell = nx1*ny1*nz1 
eps = np.zeros((nx1,ny1,nz1),dtype=float)

print('calculate eps')
for k in range (0, nz1):
    for kk in range (0,nlevel):
        for j in range (0, ny1):
            for jj in range (0,nlevel):
                for i in range (0, nx1):
                    for ii in range (0,nlevel):
                        eps[i,j,k] += (0.001*my_array[i*nlevel+ii,k*nlevel+kk,j*nlevel+jj]/255+(1-my_array[i*nlevel+ii,k*nlevel+kk,j*nlevel+jj]/255))/nlevel/nlevel/nlevel

print('create eps')
f=open('0/eps','a')
f.seek(0) #get to the first position
f.write("FoamFile"+'\n')
f.write("{"+'\n')
f.write("    version     2.0;"+'\n')
f.write("    format      ascii;"+'\n')
f.write("    class       volScalarField;"+'\n')
f.write("    object      eps;"+'\n')
f.write("}"+'\n')
f.write(""+'\n')
f.write("dimensions      [0 0 0 0 0 0 0];"+'\n')
f.write("internalField   nonuniform List<scalar>"+'\n')
f.write(str(ncell)+'\n')
f.write("("+'\n')
for k in range (0,nz1):
    for j in range (0, ny1):
        for i in range (0, nx1):
            f.write(str(eps[i,j,k])+'\n')
f.write(")"+'\n')
f.write(";"+'\n')
f.write(""+'\n')
f.write("boundaryField"+'\n')
f.write("{"+'\n')
f.write("    front"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("    back"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("    left"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("    right"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("    top"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("    bottom"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("}"+'\n')
f.close()


os.system('sed -i "s/nx/'+str(nx1)+'/g" constant/polyMesh/blockMeshDict')
os.system('sed -i "s/ny/'+str(ny1)+'/g" constant/polyMesh/blockMeshDict')
os.system('sed -i "s/nz/'+str(nz1)+'/g" constant/polyMesh/blockMeshDict')

os.system("blockMesh")

os.system('smoothSolidSurface')
os.system('cellSet')
os.system('refineHexMesh refinementRegion -overwrite')


os.system('processMeshCellCenters')

ix = np.zeros(nx*ny*nz)
iy = np.zeros(nx*ny*nz)
iz = np.zeros(nx*ny*nz)
file = open("0/cellCenters","r")
Lines = file.readlines()
count =0
wbool=0
for line in Lines:
  ls = line.strip()
  if (ls==")"):
      break
  if (wbool==1):
      x=float(ls.split("(")[1].split(")")[0].split()[0])
      ix[count] = np.floor((x-2.5e-6)/5e-6);
      y=float(ls.split("(")[1].split(")")[0].split()[1])
      iy[count] = np.floor((y-2.5e-6)/5e-6);
      z=float(ls.split("(")[1].split(")")[0].split()[2])
      iz[count] = np.floor((z-2.5e-6)/5e-6);
      count +=1
  if (ls=="("):
      wbool=1

ncell = count

newEps = np.zeros(ncell)

print("eps")
file = open("0/eps","r")
Lines = file.readlines()
count =0
wbool=0
for line in Lines:
  ls = line.strip()
  if (ls==")"):
      break
  if (wbool==1):
      epsVal = float(ls)
      if (epsVal<0.99) and (epsVal>0.01):
          val= my_array[ix[count].astype(int),iz[count].astype(int),iy[count].astype(int)]
          newEps[count] = 0.001*val/255.0 + (1-val/255.0)
      else:
          newEps[count] = epsVal
      count +=1
  if (ls=="("):
      wbool=1

os.system('rm 0/eps')

f=open('0/eps','a')
f.seek(0) #get to the first position
f.write("FoamFile"+'\n')
f.write("{"+'\n')
f.write("    version     2.0;"+'\n')
f.write("    format      ascii;"+'\n')
f.write("    class       volScalarField;"+'\n')
f.write("    object      eps;"+'\n')
f.write("}"+'\n')
f.write(""+'\n')
f.write("dimensions      [0 0 0 0 0 0 0];"+'\n')
f.write("internalField   nonuniform List<scalar>"+'\n')
f.write(str(ncell)+'\n')
f.write("("+'\n')

for k in range(0,ncell):
            f.write(str(newEps[k])+'\n')
f.write(")"+'\n')
f.write(";"+'\n')
f.write(""+'\n')
f.write("boundaryField"+'\n')
f.write("{"+'\n')
f.write("    front"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("    back"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("    left"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("    right"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("    top"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("    bottom"+'\n')
f.write("    {"+'\n')
f.write("        type zeroGradient;"+'\n')
f.write("    }"+'\n')
f.write("}"+'\n')
f.close()

os.system('smoothSolidSurface')

END

rm 0/cellCenters
rm Bentheimer*
