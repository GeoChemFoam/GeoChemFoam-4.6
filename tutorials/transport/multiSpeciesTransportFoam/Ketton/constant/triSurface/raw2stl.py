#script to raw image and create an stl mesh for DNS
import matplotlib.pyplot as plt
import numpy as np
import skimage
from skimage import data
from skimage.measure import marching_cubes_lewiner
from stl import mesh

# User Input
target_direc = ""
image_name = 'Ketton'
output_direc = ""

x_dim=512
y_dim=512
z_dim=512 

# END User Input

img = np.fromfile(target_direc+image_name+'.raw', dtype='uint8', sep="")
img = np.reshape(img, (z_dim, y_dim, x_dim))
img_padded = np.pad(img, pad_width=1,mode='constant', constant_values=0)

verts, faces, normals, values = skimage.measure.marching_cubes_lewiner(img_padded, 55, step_size=2, allow_degenerate=1)

# Create the mesh
imagemesh = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        imagemesh.vectors[i][j] = verts[f[j],:]

# Write the mesh to file "cube.stl"
imagemesh.save(output_direc+image_name+'_meshed.stl')
