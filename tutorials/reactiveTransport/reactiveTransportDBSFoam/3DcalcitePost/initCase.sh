#!/bin/bash


set -e



rm -rf 0


cp -r 0_org 0
rm 0/eps
cp constant/polyMesh/blockMeshDictRun constant/polyMesh/blockMeshDict
cp constant/dynamicMeshDictRun constant/dynamicMeshDict

python << END
import numpy as np
import h5py
import array
import os
f=h5py.File("calcitePost.h5",'r')
my_array=f['/t0/channel0'][()]

nz=10;
ny=75;
nx=134;

ncell = nx*ny*nz 
p=int(900/ny)
q=int(1608/nx)
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
eps=np.zeros((nx,ny),dtype=float)
for j in range (0,ny):
    for jj in range (0,q):
        for i in range (0,nx):
            for ii in range (0,p):	
                eps[i,j]+=(my_array[i*p+ii,j*q+jj]/255.0+ (1-my_array[i*p+ii,j*q+jj]/255.0)*0.001)/p/q
for k in range (0,nz):
    for j in range (0, ny):
        for i in range (0, nx):
            f.write(str(eps[i,j])+'\n')
f.write(")"+'\n')
f.write(";"+'\n')
f.write(""+'\n')
f.write("boundaryField"+'\n')
f.write("{"+'\n')
f.write("    frontAndBack"+'\n')
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

os.system('sed -i "s/nx/'+str(nx)+'/g" constant/polyMesh/blockMeshDict')
os.system('sed -i "s/ny/'+str(ny)+'/g" constant/polyMesh/blockMeshDict')
os.system('sed -i "s/nz/'+str(nz)+'/g" constant/polyMesh/blockMeshDict')

Kf = 1e12 
os.system('cp constant/transportPropertiesRun constant/transportProperties')
os.system('sed -i "s/val/'+str(Kf)+'/g" constant/transportProperties')



END

blockMesh
#smoothSolidSurface
decomposePar
