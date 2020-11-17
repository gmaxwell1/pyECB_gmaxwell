""" 
filename: spherical_plot_triangulation.py

This script can be used to generate a 3d, spherical plot of all fitted directions from 
a set of experiments, which have already been analyzed using the extract_linear_slopes.py script.
The slopes are displayed as color coded scatter plot at the desired directions on the unit sphere. 


Author: Nicholas Meinhardt (QZabre)
        nmeinhar@student.ethz.ch
        
Date: 13.11.2020
"""

#%%
# standard library imports
import numpy as np
import os
import matplotlib.pyplot as plt

# local imports
from modules.analysis_tools import *
from modules.interpolation_tools import *


#%%
# read in fit parameters from file, which originates from a series of measurements
data_directory = './test_data/config_tests_20_11_10'
data_filename = 'linear_fits.csv'
data_filepath = os.path.join(data_directory, data_filename)

data = extract_fit_parameters_from_file(data_filepath)
filenames, slopes, slopes_std, offsets, offsets_std, B_directions, phis, phis_std, thetas, thetas_std, B_max_1, B_max_2 = data

complemented = False

#%%
# add reverse directions 
complemented_data = add_reverse_directions(phis, phis_std, thetas, thetas_std, slopes, slopes_std, offsets, offsets_std)
phis, phis_std, thetas, thetas_std, slopes, slopes_std, offsets, offsets_std = complemented_data

complemented = True

# %%
# apply Delaunay triangulation on the surface of the sphere, which corresponds to finding convex Hull
vertices_simplices, inidces_simplices, points = delaunay_triangulation_spherical_surface(phis, thetas)

# %%
#  create 3d-plot
component = 1


# transform angles to Carthesian coordinates
xx = np.sin(np.radians(thetas))*np.cos(np.radians(phis))
yy = np.sin(np.radians(thetas))*np.sin(np.radians(phis))
zz = np.cos(np.radians(thetas))

# generate figure with 3d-axis
fig = plt.figure(figsize = 1.5*plt.figaspect(1.))
ax = fig.add_subplot(111, projection='3d')

# add triangles as surfaces
# ax.plot_trisurf(points[:,0], points[:,1], points[:,2], triangles = inidces_simplices, alpha=0.8)


# plot arrows for x, y, z axis
length_axes = 2.8
ax.quiver(length_axes/2, 0, 0, length_axes, 0, 0, color='k',
          arrow_length_ratio=0.08, pivot='tip', linewidth=1.1)
ax.quiver(0, length_axes/2, 0, 0, length_axes, 0, color='k',
          arrow_length_ratio=0.08, pivot='tip', linewidth=1.1)
ax.quiver(0, 0, length_axes/2, 0, 0, length_axes, color='k',
          arrow_length_ratio=0.08, pivot='tip', linewidth=1.1)
ax.text(1.6, 0, 0, 'x')
ax.text(0, 1.65, 0, 'y')
ax.text(0, 0, 1.6, 'z')

# create a sphere
u, v = np.mgrid[0:2*np.pi:16j, 0:np.pi:40j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)
ax.plot_surface(x, y, z, color='k', rstride=1, cstride=1,
                    alpha=0.05, antialiased=False, vmax=2)  

# plot equator
u, v = np.mgrid[0:2*np.pi:40j, np.pi/2:np.pi/2:1j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)
ax.plot_wireframe(x, y, z, color='k', linewidth=0.5)

# plot data
mask = np.array([not isnan(slopes[0, component, coil]) for coil in range(3)])
valid_coil = np.arange(3)[mask][0]  
c = ax.scatter(xx,yy,zz, c = slopes[:,component, valid_coil], cmap='inferno', s=15)

cbar_ax = fig.add_axes([1.0, 0.15, 0.05, 0.7])
cbar = fig.colorbar(c, cax=cbar_ax)
cbar.set_label('linear slope $d B_{} / d I_{}$'.format(['x','x','z'][component], valid_coil+1))

# add triangles as lines
add_triangles_to_3dplot(ax, points, inidces_simplices, spherical=True, colored_triangles=False)

# axis settings
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_zlim([-1,1])
ax.set_xlabel('$B_x$ [mT]')
ax.set_ylabel('$B_y$ [mT]')
ax.set_zlabel('$B_z$ [mT]')

# switch off axes and planes in background, rotate to nice position.
ax.set_axis_off()
ax.view_init(10, 30)

# # test interpolation with a test vector
# t = np.pi /4
# p = np.pi / 3
# test_vector = np.array([np.sin(t)*np.cos(p), np.sin(t)*np.sin(p), np.cos(t)]) # for testing interpolation scheme

# # plot test vector 

# test_vector = 1 * test_vector / np.linalg.norm(test_vector)
# ax.quiver(test_vector[0], test_vector[1], test_vector[2], 
#             1.1*test_vector[0], 1.1*test_vector[1], 1.1*test_vector[2],
#             arrow_length_ratio=0.2, color='b')

# # find the associtated triangle, edge or vertex that contains the test vector
# indices_interpol = find_associated_triangle(test_vector, points, inidces_simplices)
# highlight_interpolation_area_in_3dplot(ax, points, indices_interpol, spherical = True, color='g')

plt.tight_layout()

# if complemented:
#     image_path = os.path.join(data_directory, 
#                     '3d_slopes_{}{}_polar_complemented.png'.format(['x','x','z'][component], valid_coil+1))
# else:
#     image_path = os.path.join(data_directory, 
#                     '3d_slopes_{}{}_polar_measured.png'.format(['x','x','z'][component], valid_coil+1))
# fig.savefig(image_path, dpi=300)

plt.show()

# %%
