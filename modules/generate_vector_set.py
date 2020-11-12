""" 
filename: generate_set_of_vec.py

This script can be used to generate configuration files containing magnetic field
vectors for various directions and the same magnetic field strength.
The vectors can be generated and stored both in Carthesian and spherical coordinates.


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
    from modules.generate_configurations import generate_configs_half_sphere, plot_vectors
    

#%%
# generate configurations

# Settings:
n_sectors = 32      # >= 4 required, number of points along the equator. 
elevation_factor_equator = 0    # elevation useful to prevent redundancy for vectors on equator
magnitude = 20      # desired magnitude 
save_format = 's'   # 's' for spherical and 'c' for Carthesian coordinates

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
fig, ax = plot_vectors(vectors, magnitude=magnitude, phis=phis, thetas=thetas, add_tiangulation= True)

# image_path = os.path.join('../config_files', 'generated_points_{}.png'.format(n_sectors))
# fig.savefig(image_path, dpi=300)
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


try:
    # start with parent directory
    directory = '../config_files'
    output_file_name = 'configs_wholeSphere_mag{}_size{}.csv'.format(magnitude, n_vectors)
    data_filepath = os.path.join(directory, output_file_name)

    df.to_csv(data_filepath, index=False, header=True)
except:
    # remain in current workind directory
    directory = './config_files'
    output_file_name = 'configs_wholeSphere_mag{}_size{}.csv'.format(magnitude, n_vectors)
    data_filepath = os.path.join(directory, output_file_name)

    df.to_csv(data_filepath, index=False, header=True)

# %%
