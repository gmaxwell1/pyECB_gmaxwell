""" 
filename: data_analysis_accuracy_characterization.py

This file contains functions that are used for extracting data from measurements from the Metrolab sensor 
and analysing the accuracy of the linear actuation matrix that we have used so far<<<<<< .

Author: Maxwell Guerne, Nicholas Meinhardt (Qzabre)
        gmaxwell@student.ethz.ch
        nmeinhar@student.ethz.ch
        
Date: 09.10.2020
"""
#%%
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# local imports 
try:
    from modules.analysis_tools import *
except ModuleNotFoundError:
    import sys
    sys.path.insert(1, os.path.join(sys.path[0], '..'))
finally:
    from modules.analysis_tools import *
    from modules.interpolation_tools import find_start_of_saturation




#%%
directory = r'data_sets\config_tests_20_11_10'
threshold = 10
currents, B_measured, B_expected = collectAndExtract(directory, threshold, remove_saturation=True)


#%%
# find best rotation of expected fields that represents a misalignment between magnet and sensor
p, pcov, B_expected_rotated = fit_rotation_matrix(B_measured, B_expected)
print('fitted angles: alpha = {:.2f}°, beta = {:.2f}°, gamma = {:.2f}°'.format(p[0], p[1], p[2]))

# illustrate the effect of the found rotation by printing its 
_ = rotation_on_basis_vectors(*p, verbose=True)


#%%
# B_reference = B_expected
B_reference = B_expected_rotated


# errors in theta and phi
delta_theta = np.abs(get_theta(B_reference) - get_theta(B_measured))
delta_phis = np.abs(get_phi(B_reference) - get_phi(B_measured))

print('mean error theta: {:.2f}°, std: {:.2f}°'.format(np.mean(delta_theta), np.std(delta_theta)))
print('mean error phi: {:.2f}°, std: {:.2f}°'.format(np.mean(delta_phis), np.std(delta_phis)))

# angles between expectations and measurements
dot = np.array([np.dot(B_measured[i], B_reference[i]) for i in range(len(B_measured))])
norms_measured = np.linalg.norm(B_measured, axis=1)
norms_expected = np.linalg.norm(B_reference, axis=1)
alphas = np.degrees(np.arccos(dot / (norms_measured * norms_expected)))

print('mean angular error: {:.2f}°, std: {:.2f}°'.format(np.mean(alphas), np.std(alphas)))
print('min / max angular error: {:.2f}° / {:.2f}°'.format(np.min(alphas), np.max(alphas)))
print('median angular error: {:.2f}°'.format(np.median(alphas)))



# %%
# generate histogram
# the histogram of the data
n, bins, patches = plt.hist(alphas, 50)

plt.xlabel('Angular Error, $\\alpha$ [°]')
plt.ylabel('Counts')
plt.text(15, 200, '$\\mu$ = {:.2f}°,\n'.format(np.mean(alphas))+ 
                '$\\sigma$ = {:.2f}°,\n'.format(np.std(alphas))+ 
                'min = {:.2f}°,\n'.format(np.min(alphas))+ 
                'max = {:.2f}°,\n'.format(np.max(alphas))+ 
                'median = {:.2f}°,\n'.format(np.median(alphas)))
plt.savefig('angular_error_0mT_excl_saturation.png', dpi=300)
plt.show()

# %%
average_error = np.mean(np.abs(norms_expected - norms_measured))
std_error = np.std(np.abs(norms_expected - norms_measured))
max_error = np.amax(np.abs(norms_expected - norms_measured))
min_error = np.amin(np.abs(norms_expected - norms_measured))
median_error = np.median(np.abs(norms_expected - norms_measured))
print('average absolute difference between predicted \nand measured magnetic field strength:'
      f' {average_error:.2f} mT')
print(f'standard deviation: {std_error:.2f} mT')
print(f'Median deviation: {median_error:.2f} mT')
print(f'maximum/minimum deviation: {max_error:.2f} mT, {min_error:.2f} mT')

# %%
fig, ax = plt.subplots(2,sharex=True)
ax[0].plot(norms_measured, alphas, linestyle='', marker='.', label='angle error', color='C2')
ax[1].set_xlabel('measured field magnitude [mT]')
ax[1].set_xticks([10,20,40,60,80,100,120])
ax[1].set_xticks([30,50,70,90,110],minor=True)
ax[1].set_xticklabels([10,20,40,60,80,100,120])
ax[0].set_ylabel('$\\Delta$angle [°]', color='C2')
ax[1].plot(norms_measured, np.abs(norms_expected - norms_measured), linestyle='', marker='.', label='magnitude error', color='C1')
ax[1].set_ylabel('$\\Delta$B [mT]', color='C1')
ax[0].legend()
ax[1].legend()
plt.savefig('angular_error_0mT_excl_saturation.png', dpi=300)


# %%
