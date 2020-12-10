""" 
filename: data_analysis_accuracy_characterization.py

This file contains functions that are used for extracting data from measurements from the Metrolab sensor 
and analysing the accuracy of the linear actuation matrix that we have used so far<<<<<< .

Author: Maxwell Guerne, Nicholas Meinhardt (Qzabre)
        gmaxwell@student.ethz.ch
        nmeinhar@student.ethz.ch
        
Date: 09.10.2020
"""

import os
import numpy as np
import pandas as pd



def collectAndExtract(directory):
    # collect all csv-files in this directory
    filenames = []
    [filenames.append(file) for file in os.listdir(directory) if file.endswith(".csv")]

    # if in list, remove the linear_fits file from previous fits
    try:
        filenames.remove('linear_fits.csv')
    except ValueError:
        pass

    # initialize lists to store the relevant data
    B_measured = []
    B_expected = []
    list_phi = []
    list_theta = []
    list_slopes = []
    list_offsets = []
    list_B_max_1 = []
    list_B_max_2 = []

    # loop through all csv files in a dictionary and fit the data
    for data_filename in filenames:
        if verbose:
            print(data_filename)

       
        # collect fit parameters
        list_filenames.append(data_filename)
        list_B_directions.append(B_direction)
        list_phi.append(phi)
        list_phi_std.append(phi_std)
        list_theta.append(theta)
        list_theta_std.append(theta_std)
        list_slopes.append(slopes)
        list_slopes_std.append(slopes_std)
        list_offsets.append(offsets)
        list_offsets_std.append(offsets_std)
        list_B_max_1.append(B_max_1)
        list_B_max_2.append(B_max_2)
        
        
    raw_data = pd.read_csv(filepath).to_numpy()
