import os
import h5py
import numpy as np

# USER INPUT
image = "Bentheimer400-5mum"
xdim = 400
ydim = 400 
zdim = 412
voxelsize_um =5 
# END USER INPUT

wd = os.getcwd()
target_direc = wd

Dim_size = np.array((zdim, ydim, xdim), dtype=np.int)

f = open(target_direc + "/" + image + ".raw", 'rb')
img_arr = np.fromfile(f, dtype=np.uint8)
img = img_arr.reshape(Dim_size[0], Dim_size[1], Dim_size[2])

g = h5py.File(target_direc + "/" + image + ".hdf5", 'w')
g.create_dataset('image', data=img, dtype="uint16", compression="gzip")
g.attrs['voxelsize_microns'] = voxelsize_um
f.close()
g.close()
