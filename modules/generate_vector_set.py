""" 
filename: generate_set_of_vec.py

This script can be used to generate configuration files containing magnetic field
vectors for various directions and the same magnetic field strength and test the interpolation scheme.

In Part 1, the vectors can be generated and stored both in Carthesian and spherical coordinates.

In Part 2, a test vector can be plotted including the associated triangle/edge/vertex.
Moreover,the interpolation scheme can be tested for both single and multiple inputs.


Author: Nicholas Meinhardt (QZabre)
        nmeinhar@student.ethz.ch
        
Date: 11.11.2020
"""
#%%
# import modules
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt

try:
    from modules.generate_configurations import generate_configs_half_sphere, plot_vectors
except ModuleNotFoundError:
    import sys
    sys.path.insert(1, os.path.join(sys.path[0], '..'))
    from modules.generate_configurations import generate_configs_half_sphere, plot_vectors, generate_test_points_whole_sphere
    from modules.interpolation_tools import find_associated_triangle, highlight_interpolation_area_in_3dplot, apply_linear_interpolation
    

# Part 1: --------------------------------------------------------------
#%%
# generate configurations

# Settings:
n_sectors = 16    # >= 4 required, number of points along the equator. 
elevation_factor_equator = 0    # elevation useful to prevent redundancy for vectors on equator
magnitude = 20      # desired magnitude 
save_format = 's'   # 's' for spherical and 'c' for Carthesian coordinates
test_vector = np.array([0, 1, 1]) # for testing interpolation scheme

# estimate vectors and respective angles in upper and lower hemisphere,
# make sure that equator is treated only once
_, vectors_upper, thetas_upper, phis_upper = generate_configs_half_sphere(n_sectors, magnitude=magnitude, 
                                    upper=True, elevation_factor_equator = elevation_factor_equator,
                                    include_equator=True)
_, vectors_lower, thetas_lower, phis_lower = generate_configs_half_sphere(n_sectors, magnitude=magnitude, 
                                    upper=False, elevation_factor_equator = elevation_factor_equator,
                                    include_equator=False)

# combine both hemispheres
vectors = np.append(vectors_upper, vectors_lower, axis=0)
thetas = np.append(thetas_upper, thetas_lower, axis=0)
phis = np.append(phis_upper, phis_lower, axis=0)

n_vectors = len(vectors)
print('number of vectors: {}'.format(n_vectors))


# plot all considered vectors on a sphere 
fig, ax, inidces_simplices, points = plot_vectors(vectors, magnitude=magnitude, phis=phis, thetas=thetas, 
                                                add_tiangulation= True)

plt.show()


# save the combinations to csv file
if save_format == 's':
    df = pd.DataFrame({ 'convention': np.full(n_vectors, 's' ), 
                        'magnitude [mT]': np.full(n_vectors, magnitude), 
                        'theta [°]': np.degrees(thetas),
                        'phi [°]': np.degrees(phis)})
else:
    df = pd.DataFrame({ 'convention': np.full(n_vectors, 'c' ), 
                        'Bx [mT]': vectors[:,0], 
                        'By [mT]': vectors[:,1],
                        'Bz [mT]': vectors[:,2]})

 
# try:
#     # start with parent directory
#     directory = '../config_files'
#     output_file_name = 'configs_wholeSphere_mag{}_size{}.csv'.format(magnitude, n_vectors)
#     data_filepath = os.path.join(directory, output_file_name)

#     df.to_csv(data_filepath, index=False, header=True)
# except:
#     # remain in current workind directory
#     directory = './config_files'
#     output_file_name = 'configs_wholeSphere_mag{}_size{}.csv'.format(magnitude, n_vectors)
#     data_filepath = os.path.join(directory, output_file_name)

#     df.to_csv(data_filepath, index=False, header=True)


# Part 2: --------------------------------------------------------------
#%%
# test interpolation with a test vector
t = np.pi /16
p = np.pi / 2.1
test_vector = np.array([np.sin(t)*np.cos(p), np.sin(t)*np.sin(p), np.cos(t)]) # for testing interpolation scheme


# plot all considered vectors on a sphere 
fig, ax, inidces_simplices, points = plot_vectors(vectors, magnitude=magnitude, phis=phis, thetas=thetas, 
                                                add_tiangulation= True)

# plot test vector 
test_vector = magnitude * test_vector / np.linalg.norm(test_vector)
ax.quiver(test_vector[0], test_vector[1], test_vector[2], 
            1.1*test_vector[0], 1.1*test_vector[1], 1.1*test_vector[2],
            arrow_length_ratio=0.2, color='b')

# find the associtated triangle, edge or vertex that contains the test vector
indices_interpol = find_associated_triangle(test_vector, points, inidces_simplices)
highlight_interpolation_area_in_3dplot(ax, points, indices_interpol, spherical = True, color='g')

image_path = os.path.join('../config_files', 'generated_points_{}_test_interpolation.png'.format(n_sectors))
fig.savefig(image_path, dpi=300)


plt.show()

#%%
# define relevant quantity that should be interpolated, length must be equal to length(points) as only restriction 
values = np.tile(points,3).reshape(-1,3,3)


y = apply_linear_interpolation(test_vector, points, values, inidces_simplices)
print(y)


multiple_vecs = generate_test_points_whole_sphere(64, magnitude)
print(multiple_vecs.shape)
y_multiple = np.array([apply_linear_interpolation(test_vector, points, values, inidces_simplices) for x in multiple_vecs])
# print(y_multiple)



# %%
