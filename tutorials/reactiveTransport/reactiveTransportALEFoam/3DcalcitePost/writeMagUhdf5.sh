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

dimX = 268  # dimension of image in X
dimY = 150  # dimension of image in Y
dimZ = 20

print("count dir")
ndir = 0

for dir in os.listdir(os.getcwd()): 
  if (is_number(dir)):
    ndir +=1

os.system('rm -f 3DcalictePostALE.hdf5')

time = np.zeros(ndir)
C = np.zeros((dimX, dimY, dimZ, ndir))

tcount =0
for dir in os.listdir(os.getcwd()): 
  if (is_number(dir)):
    print(dir)
    time[tcount]=float(dir)
    x = np.zeros(dimX*dimY*dimZ)
    y = np.zeros(dimX*dimY*dimZ)
    z = np.zeros(dimX*dimY*dimZ)
    file = open(dir+"/cellCenters","r")
    Lines = file.readlines()
    count =0
    wbool=0
    for line in Lines:
      ls = line.strip()
      if (ls==")"):
        break
      if (wbool==1):
        x[count]=float(ls.split("(")[1].split(")")[0].split()[0]) 
        y[count]=float(ls.split("(")[1].split(")")[0].split()[1])
        z[count]=float(ls.split("(")[1].split(")")[0].split()[2])
        count +=1
      if (ls=="("):
        wbool=1

    file = open(dir+"/H+","r")
    Lines = file.readlines()
    count =0
    wbool=0
    for line in Lines:
      ls = line.strip()
      if (ls==")"):
        break
      if (wbool==1):
        a = np.round((x[count]+0.0010885)/1e-5)-1
        b = np.round((y[count]+0.000503)/1e-5)-1
        c = np.round((z[count]+0.000105)/1e-5)-1
        C[a.astype(int), b.astype(int),c.astype(int),tcount] = float(ls)
        count +=1
      if (ls=="("):
        wbool=1
    tcount += 1



f = h5py.File("3DcalcitePostALE.hdf5","w")

v = np.argsort(time)
time = time[v]
C=C[:,:,:,v]

f.create_dataset('time', data=time, dtype="float", compression="gzip")
f.create_dataset('C', data=C, dtype="float", compression="gzip")



f.close()




END


