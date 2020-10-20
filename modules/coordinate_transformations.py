""" 
filename: coordinate_transformations.py

The following functions are to transform between the three coordinate systems in use,
namely stage coordinates, (Metrolab) sensor coordinates and magnet coordinates. 

Author: Nicholas Meinhardt (QZabre)
        nmeinhar@student.ethz.ch

Date: 20.10.2020
"""

import numpy as np



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
        data_magnet_coords[:,0] = -data[:,2]
        data_magnet_coords[:,1] = -data[:,0]
        data_magnet_coords[:,2] = data[:,1]

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
        data_magnet_coords[:,0] = -data[:,1]
        data_magnet_coords[:,1] = data[:,2]
        data_magnet_coords[:,2] = -data[:,0]

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
        data_magnet_coords[:,0] = -data[:,1]
        data_magnet_coords[:,1] = -data[:,0]
        data_magnet_coords[:,2] = data[:,2]

    return data_magnet_coords
