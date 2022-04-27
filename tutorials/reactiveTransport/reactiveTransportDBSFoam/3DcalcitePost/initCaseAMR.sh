#!/bin/bash


set -e

cp constant/dynamicMeshDictAMR constant/dynamicMeshDict
cp system/controlDictAMRInit system/controlDict
cp system/fvSolutionAMRInit system/fvSolution

rm -rf 0
cp -r 0_org 0
rm 0/eps
cp constant/polyMesh/blockMeshDictRun constant/polyMesh/blockMeshDict

python << END
import numpy as np
import h5py
import array
import os
f=h5py.File("calcitePost.h5",'r')
my_array=f['/t0/channel0'][()]


nz=10
ny=75
nx=134

maxRef=1
lowRef=0.01

dt=0.5/maxRef
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
                eps[i,j]+=(my_array[i*p+ii,j*q+jj]/255.0+ (1-my_array[i*p+ii,j*q+jj]/255.0)*1e-4)/p/q
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

os.system('sed -i "s/nRef/1/g" constant/dynamicMeshDict')
os.system('sed -i "s/lowRef/'+str(lowRef)+'/g" constant/dynamicMeshDict')
os.system('sed -i "s/upRef/1/g" constant/dynamicMeshDict')
os.system('sed -i "s/refLevel/'+str(maxRef)+'/g" constant/dynamicMeshDict')

os.system('sed -i "s/dt/'+str(dt)+'/g" system/controlDict')

os.system("blockMesh")

os.system("decomposePar")
os.system("mpiexec -np 24 reactiveTransportDBSFoam -parallel")

os.system("reconstructParMesh -latestTime")

os.system("rm -rf 0")

os.system("mv 0.5 0")

os.system("rm 0/phi")
os.system("cp 0_org/U 0/.")
os.system("cp 0_org/p 0/.")
os.system("cp 0_org/H+ 0/.")

os.system("processMeshCellCenters")
 
f=h5py.File("calcitePost.h5",'r')
my_array=f['/t0/channel0'][()]

nz=20*maxRef
ny=150*maxRef
nx=268*maxRef


p=int(900/ny)
q=int(1608/nx)

dx=266.7*1e-5/nx
dy=149.6*1e-5/ny

eps=np.zeros((nx,ny),dtype=float)
for j in range (0,ny):
    for jj in range (0,q):
        for i in range (0,nx):
            for ii in range (0,p):	
                eps[i,j]+=(my_array[i*p+ii,j*q+jj]/255.0+ (1-my_array[i*p+ii,j*q+jj]/255.0)*1e-3)/p/q
                
ix = np.zeros(nx*ny*nz)
iy = np.zeros(nx*ny*nz)
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
      ix[count] = np.floor((x+108.35*1e-5)/dx);
      y=float(ls.split("(")[1].split(")")[0].split()[1])
      iy[count] = np.floor((y+49.8*1e-5)/dy);
      count +=1
  if (ls=="("):
      wbool=1

ncell = count

newEps = np.zeros(ncell)

file = open("0/eps","r")
Lines = file.readlines()
count =0
wbool=0
for line in Lines:
  ls = line.strip()
  if (ls==")"):
      break
  if (wbool==1):
      epsVal=float(ls)
      if (epsVal<1) and (epsVal>0.01):
        newEps[count] = eps[ix[count].astype(int),iy[count].astype(int)]
      else:
        newEps[count] = epsVal
      count +=1
  if (ls=="("):
      wbool=1
           



f=open('0/eps2','a')
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

Kf = 1e12
os.system('cp constant/transportPropertiesRun constant/transportProperties')
os.system('sed -i "s/val/'+str(Kf)+'/g" constant/transportProperties')


END

rm -f 0/cellCenters
rm -f 0/ddt*
rm -rf 0/uniform
mv 0/eps2 0/eps
rm -rf processor*/0

for i in processor*; do mv "$i/0.5" "$i/0"; done

rm processor*/0/phi
rm -f processor*/0/ddt*
rm -rf processor*/0/uniform

rm -rf 0/R
rm -rf processor*/0/R

decomposePar -fields

