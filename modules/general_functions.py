""" 
filename: general_functions.py

The following functions are of general use for file management and simple calculations. 
They also contain transformations between the different coordinate systems in use,
namely stage coordinates, (Metrolab) sensor coordinates and magnet coordinates. 

Author: Nicholas Meinhardt (QZabre)
        nmeinhar@student.ethz.ch
        
Date: 20.10.2020
"""

########## Standard library imports ##########
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import serial
from numpy.linalg import norm
import os
import sys
from datetime import datetime


def angle_wrt_z(vec):
    """
    Return angle (radian) of vector with respect to z axis.
    """
    mag = norm(vec)
    return np.arccos(vec[2]/mag)

def inplane_angle_wrt_x(vec):
    """
    Return angle (radian) of vector with respect to z axis.
    """
    return np.arctan2(vec[1], vec[0])

def sensor_to_magnet_coordinates(data):
    """
    Transform from Metrolab sensor coordinates to magnet coordinates, using the following transformation:

    x -> -y
    y -> z
    z -> -x

    Arg: data (ndarray) can be 1d or multidimensional array, where the last dimension must be of length 3 and contain 
    x,y,z data. 
    """
    data_magnet_coords = np.zeros_like(data)
    
    # treat 1d arrays and multi-dimensional arrays differently
    if len(data.shape)==1:
        data_magnet_coords[0] = -data[2]
        data_magnet_coords[1] = -data[0]
        data_magnet_coords[2] = data[1]
    else:
        data_magnet_coords[...,0] = -data[...,2]
        data_magnet_coords[...,1] = -data[...,0]
        data_magnet_coords[...,2] = data[...,1]

    return data_magnet_coords

def magnet_to_sensor_coordinates(data):
    """
    Transform from magnet coordinates to Metrolab sensor coordinates, using the following transformation:

    x -> -z
    y -> -x
    z -> y

    Arg: data (ndarray) can be 1d or multidimensional array, where the last dimension must be of length 3 and contain 
    x,y,z data. 
    """
    data_magnet_coords = np.zeros_like(data)
    
    # treat 1d arrays and multi-dimensional arrays differently
    if len(data.shape)==1:
        data_magnet_coords[0] = -data[1]
        data_magnet_coords[1] = data[2]
        data_magnet_coords[2] = -data[0]
    else:
        data_magnet_coords[...,0] = -data[...,1]
        data_magnet_coords[...,1] = data[...,2]
        data_magnet_coords[...,2] = -data[...,0]

    return data_magnet_coords

def transform_between_sensor_stage_coordinates(data):
    """
    Transform from magnet coordinates to Metrolab sensor coordinates, using the following transformation:

    x -> -y
    y -> -x

    Arg: data (ndarray) can be 1d or multidimensional array, where the last dimension must be of length 3 and contain 
    x,y,z data. 
    """
    data_magnet_coords = np.zeros_like(data)
    
    # treat 1d arrays and multi-dimensional arrays differently
    if len(data.shape)==1:
        data_magnet_coords[0] = -data[1]
        data_magnet_coords[1] = -data[0]
        data_magnet_coords[2] = data[2]
    else:
        data_magnet_coords[...,0] = -data[...,1]
        data_magnet_coords[...,1] = -data[...,0]
        data_magnet_coords[...,2] = data[...,2]

    return data_magnet_coords

def save_in_dir(means, directory, label, stds=None, coords=False, now=False):
    """
    Write the provided means and (optional) standard deviations to the file 'means_'+label+'.csv' 
    or 'yyy_mm_dd_hh-mm-ss_means_'+label+'.csv' in provided directory.

    The coords flag distinguishes between magnetic field strength (False) and spatial coordinates (True) as source of the data.

    Args:
    - means (array or list): measured mean values of B-field or coordinates 
    - directory (string): valid path of the directory in which the file should be stored
    - label (string or float): used to label the csv file  
    - stds (array or list): measured standard deviations of B-field or coordinates. 
    Should have at least the same size as means.
    - coords (bool): Flag to switch between B-field (False) and spatial coordinates (True)
    - verbose (bool): switching on/off print-statements for displaying progress
    - now (bool): if True, the current date time is added to the filename, 
    such that it reads 'yyy_mm_dd_hh-mm-ss_means_'+label+'.csv'
    """
    # Under Linux, user rights can be set with the scheme below,
    # where 755 means read+write+execute for owner and read+execute for group and others.
    # However, note that you can only set the file’s read-only flag with it under Windows.
    access_rights = 0o755
    os.chmod(directory, access_rights)

    if now:
        time_stamp = datetime.now().strftime("%y_%m_%d_%H-%M-%S") 
        output_file_name = "{}_means_{}.csv".format(time_stamp, label)
    else:
        output_file_name = "means_{}.csv".format(label)
    data_filepath = os.path.join(directory, output_file_name)

    if stds is not None and not coords:
        df = pd.DataFrame({ 'mean Bx [mT]':  means[:, 0], 
                            'std Bx [mT] ': stds[:, 0], 
                            'mean By [mT]':  means[:, 1], 
                            'std By [mT] ': stds[:, 1],
                            'mean Bz [mT]':  means[:, 2], 
                            'std Bz [mT] ': stds[:, 2]})
        df.to_csv(data_filepath, index=False, header=True)
    
    elif not coords:
        df = pd.DataFrame({ 'mean Bx [mT]':  means[:, 0], 
                            'mean By [mT]':  means[:, 1], 
                            'mean Bz [mT]':  means[:, 2]})
        df.to_csv(data_filepath, index=False, header=True)

    elif coords and stds is None:
        df = pd.DataFrame({ 'Index':  np.arange(len(means[:, 0])) + 1, 
                            'x [mm]':  means[:, 0], 
                            'y [mm]':  means[:, 1], 
                            'z [mm]':  means[:, 2]})
        df.to_csv(data_filepath, index=False, header=True)
        
    elif not coords:
        df = pd.DataFrame({ 'Index':  np.arange(len(means[:, 0])) + 1, 
                            'x [mm]':  means[:, 0], 
                            'std x [mm]':  stds[:, 0], 
                            'y [mm]':  means[:, 1], 
                            'std y [mm]':  stds[:, 1], 
                            'z [mm]':  means[:, 2],
                            'std z [mm]':  stds[:, 2]})
        df.to_csv(data_filepath, index=False, header=True)


def ensure_dir_exists(directory, access_rights=0o755, purpose_text='', verbose=False):
    """
    Ensure that the directory exists and create the respective folders if it doesn't.

    Prints message that either a folder has been created or that it already exists. 
    Note that only read and write options can be set under Windows, rest ignored.

    Args:
    - directory (string) is a path which should exist after calling this function
    - access_rights: set rights for reading and writing.
    - purpose_text (str): add more information what the dictionary was created for

    Return:
    - 0 if directory needed to be created
    - 1 if directory already exists
    - -1 if there was an exception
    """
    try:
        os.mkdir(directory, access_rights)
        os.chmod(directory, access_rights)
        if verbose:
            print('Created directory {}: {}'.format(
                purpose_text, os.path.split(directory)[1]))
        return 0
    except FileExistsError:
        if verbose:
            print('Folder already exists, no new folder created.')
        return 1
    except Exception as e:
        # if verbose:
        # this is important and should always be printed - gmaxwell, 8.10.2020
        print('Failed to create new directory due to {} error'.format(e))
        return -1