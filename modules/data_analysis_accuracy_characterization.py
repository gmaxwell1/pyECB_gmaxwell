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

########## local imports ##########
try:
    from modules.analysis_tools import *
except ModuleNotFoundError:
    import sys
    sys.path.insert(1, os.path.join(sys.path[0], '..'))
finally:
    from modules.analysis_tools import *
    from modules.interpolation_tools import find_start_of_saturation


#%%
def collectAndExtract(directory, B_min, remove_saturation = True,
                        verbose=False, fraction_cutoff = 0.02,
                        affine_fct = lambda x, a, b: a*x + b):

    # collect all csv-files in this directory
    filenames = []
    [filenames.append(file) for file in os.listdir(directory) if file.endswith(".csv")]

    # if in list, remove the linear_fits file from previous fits
    try:
        filenames.remove('linear_fits.csv')
    except ValueError:
        pass

    # initialize lists to store the relevant data
    # B_measured = []
    # B_expected = []
    # list_phi = []
    # list_theta = []
    # list_slopes = []
    # list_offsets = []
    # list_B_max_1 = []
    # list_B_max_2 = []

    # loop through all csv files in a dictionary and fit the data
    for i in range(len(filenames)):
        if verbose:
            print(filenames[i])

        # read in raw measurment data
        data_filepath = os.path.join(directory, filenames[i])
        I, mean_data, std_data, expected_fields = extract_raw_data_from_file(data_filepath)

        # -> this could be used to exclude data above saturation, but could be left out
        # estimate minimum and maximum indices of region within which the linear relation between current and field holds  
        # even though find_start_of_saturation offers the possibility to specify the considered component, 
        # keep the default stting, which detects the field component that has the greatest absolute field values. 
        # This should work fine for situations, where one component is dominating. 
        i_min, i_max = find_start_of_saturation(I, mean_data, std_data, fraction_cutoff=fraction_cutoff,
                                fitting_fct = affine_fct)

        # estimate field magnitudes
        magnitudes = np.linalg.norm(mean_data, axis=1)

        # set up mask to only keep data with magnitudes larger than B_min
        mask_keep = magnitudes >= B_min

        # optionally: also remove potentially saturated part
        if remove_saturation:
            mask_keep[:i_min+1] = False
            mask_keep[i_max:] = False

        # collect all relevant data
        if i == 0:
            B_measured = mean_data[mask_keep]
            B_expected = expected_fields[mask_keep]
        else:
            B_measured = np.append(B_measured, mean_data[mask_keep], axis=0)
            B_expected = np.append(B_expected, expected_fields[mask_keep], axis=0)

    print(B_measured.shape)
    # raw_data = pd.read_csv(filepath).to_numpy()
    return B_measured, B_expected


#%%
directory = '../test_data/config_tests_20_11_03'
B_measured, B_expected = collectAndExtract(directory, 10)


#%%
# errors in theta and phi
delta_theta = get_theta(B_expected) - get_theta(B_measured)
phis_theta = get_phi(B_expected) - get_phi(B_measured)


# angles between expectations and measurements
dot = np.array([np.dot(B_measured[i], B_expected[i]) for i in range(len(B_measured))])
norms_measured = np.linalg.norm(B_measured, axis=1)
norms_expected = np.linalg.norm(B_expected, axis=1)
alphas = np.degrees(np.arccos(dot / (norms_measured * norms_expected)))

print('mean angular error: {:.2f}°, std: {:.2f}°'.format(np.mean(alphas), np.std(alphas)))