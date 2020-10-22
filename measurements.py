"""
filename: measurements.py

This script is meant to be used to measure magnetic field values with the Hall
sensor cube/Metrolab THM1176. It is adapted to interface with the ECB, e.g. to set current values and 
measure the generated magnetic field of the vector magnet.

Author: Nicholas Meinhardt, Maxwell Guerne-Kieferndorf (QZabre)
        nmeinhar@student.ethz.ch, gmaxwell@student.ethz.ch

Date: 20.10.2020
"""
########## Standard library imports ##########
import numpy as np
import serial
from time import sleep, time
import matplotlib.pyplot as plt
import os
import pandas as pd
from datetime import datetime

########## local imports ##########
from conexcc.conexcc_class import *
# from modules.calibrate_cube import get_new_mean_data_set, find_center_axis, angle_calib
# from modules.plot_hall_cube import plot_many_sets, plot_stage_positions, plot_set, plot_sensor_positions
from modules.general_functions import ensure_dir_exists
from MetrolabTHM1176.thm1176 import MetrolabTHM1176Node
from modules.MetrolabMeasurements import readoutMetrolabSensor, get_mean_dataset_MetrolabSensor

__all__ = [
    'calibration',
    'newMeasurementFolder',
    'measure',
    'saveDataPoints'
]

########## sensor cube/Conexcc ports ##########
# port_sensor = 'COM3'
z_COM_port = 'COM6' # z-coordinate controller
y_COM_port = 'COM5'

def newMeasurementFolder(defaultDataDir='data_sets', sub_dir_base='z_field_meas', verbose=False):
    """
    This function creates a new directory to store data from a measurement run.

    Args:
        defaultDataDir (str, optional): The directory where you want to store the subdirectories containing the actual data files. Defaults to 'data_sets'.
        sub_dir_base (str, optional): The specific subdirectory base name that will be suffixed with a number. Defaults to 'z_field_meas'.
        verbose (bool, optional): Whether it should tell you everything that's going on. Defaults to False.

    Returns:
        sub_dirname, dataDir (str): name of subdirectory where data is stored and the absolute path to it.
    """
    index = 1
    cwd = os.getcwd()
    if verbose:
        print(cwd)
    sub_dirname = sub_dir_base + '_' + str(index)
    dataDir = os.path.join(cwd, defaultDataDir, sub_dirname)
    if verbose:
        print(dataDir)
    # iterate though postfixes until a new directory name is found
    while ensure_dir_exists(dataDir, verbose=verbose):
        index = index + 1
        sub_dirname = sub_dir_base + '_' + str(index)
        dataDir = os.path.join(cwd, defaultDataDir, sub_dirname)
        if verbose:
            print(dataDir)

    return sub_dirname, dataDir


def calibration(node: MetrolabTHM1176Node, meas_height=1.5):
    # initialize actuators
    CC_Z = ConexCC(com_port=z_COM_port, velocity=0.4, set_axis='z', verbose=False)
    CC_Y = ConexCC(com_port=y_COM_port, velocity=0.4, set_axis='y', verbose=False)
    CC_Z.wait_for_ready()
    CC_Y.wait_for_ready()
    
    cal_pos_z = 20
    start_pos_z = CC_Z.read_cur_pos()
    total_distance_z = cal_pos_z-start_pos_z
    
    cal_pos_y = 0
    start_pos_y = CC_Y.read_cur_pos()
    total_distance_y = abs(cal_pos_y-start_pos_y)
    
    print('Moving to calibration position...')
    CC_Z.move_absolute(new_pos=cal_pos_z)
    CC_Y.move_absolute(new_pos=cal_pos_y)
    if total_distance_y > total_distance_z:
        progressBar(CC_Y, start_pos_y, total_distance_y)
    else:
        progressBar(CC_Z, start_pos_z, total_distance_z)

    char = input('\nPress enter to start (zero-gauss chamber!) calibration (any other key to skip): ')
    if char == '':
        node.calibrate()
    input('Press enter to continue measuring')
    
    # change meas_offset_z according to the following rule: 
    # height above sensor = meas_offset_z + 0.9 - 7.66524 (last number is "zero")
    meas_offset_z = meas_height - 0.9 + 7.66524
    start_pos_z = CC_Z.read_cur_pos()
    total_distance_z = abs(meas_offset_z-start_pos_z)
    
    meas_offset_y = 14.9
    start_pos_y = CC_Y.read_cur_pos()
    total_distance_y = abs(meas_offset_y-start_pos_y)
    
    CC_Z.move_absolute(new_pos=meas_offset_z)
    CC_Y.move_absolute(new_pos=meas_offset_y)
    
    if total_distance_y > total_distance_z:
        progressBar(CC_Y, start_pos_y, total_distance_y)
    else:
        progressBar(CC_Z, start_pos_z, total_distance_z)


def progressBar(CC: ConexCC, start_pos, total_distance):
    while not CC.is_ready():
            sleep(0.02)
            pos = CC.read_cur_pos()
            ratio = (abs(pos-start_pos) / total_distance)
            left = int(ratio * 30)
            right = 30-left
            print('\r[' + '#' * left + ' ' * right + ']', f' {ratio * 100:.0f}%', sep='', end='', flush=True)
    
    print('\n')


def measure(node: MetrolabTHM1176Node, N=50, average=False):
    """
    starts communication with hall sensor cube, measures the magnetic field with the specified sensor (change specific_sensor variable if necessary)
    
    Args:
    - node (MetrolabTHM1176Node): represents the Metrolab THM 1176 sensor
    - N: number of data points collected for each average
    - specific_sensor: sensor from which data will be fetched, only with hall sensor cube

    Returns:
    - mean_data: mean measurement data of 'specific_sensor' (averaged over N measurements)
    - std_data: standard deviation in each averaged measurement
    (where mean, std are returned as ndarrays of shape (1, 3) for the 3 field directions)
    """
    # if dataDir is not None:
    #     ensure_dir_exists(dataDir, verbose=False)
               
    if average:
        # measure average field at one point with the specific
        mean_data, std_data = get_mean_dataset_MetrolabSensor(node, sampling_size=N)
        ret1, ret2 = mean_data, std_data
    else:
        # This option is more for getting time-field measurements.
        mean_data, std_data = get_mean_dataset_MetrolabSensor(node, sampling_size=1)
        ret1, ret2 = mean_data, std_data

    return ret1, ret2

def timeResolvedMeasurement(period=0.001, averaging=1, block_size=1, duration=10, return_temp_data=False):
    """
    Measure magnetic flux density over time.

    Args:
        period (float, optional): Trigger period in seconds. Defaults to 0.001 (1 ms).
        averaging (int, optional): The arithmetic mean of this number of measurements is taken before they are fetched. 
                                   Results in smoother measurements. Defaults to 1.
        block_size (int, optional): How many measurements should be fetched at once. Defaults to 1.
        duration (int, optional): Total duration of measurement. Defaults to 10.
        return_temp_data (bool, optional): Temperature measured by sensor will or will not be returned. This is a dimensionless value between 0 and 64k. 
                                           Not calibrated, mainly useful for detecting problems due to temperature fluctuations. Defaults to False.

    Raises:
        ValueError: If for some reason the number of time values and B field values is different.

    Returns:
        lists of floats: Bx, By, Bz, timeline (x, y and z components of B field, times of measurements)
        list of floats: temp (temperature values as explained above)
    """
    with MetrolabTHM1176Node(period=period, block_size=block_size,
                             range='0.3 T', average=averaging, unit='MT') as node:
        node.start_acquisition()
        sleep(duration)
        node.stop = True
        
        flux_temp_data = [node.data_stack[key] for key in item_name]
        timeline = node.data_stack['Timestamp']
        
        t_offset = timeline[0]
        for ind in range(len(timeline)):
            timeline[ind] = round(timeline[ind]-t_offset, 3)
        
        Bx = flux_temp_data[0]
        By = flux_temp_data[1]
        Bz = flux_temp_data[2]
        temp = flux_temp_data[3]
        
        try:
            if (len(Bx) != timeline or len(By) != timeline or len(Bz) != timeline):
                raise ValueError("length of Bx, By, Bz do not match that of the timeline")  
            else:
                if return_temp_data:
                    return Bx, By, Bz, temp, timeline
                
                return Bx, By, Bz, timeline
        except:
            return [0], [0], [0], [0]
        

def saveDataPoints(I, mean_data, std_data, expected_fields, directory, data_filename_postfix='B_field_vs_I'):
    """
    Saves input data points to a .csv file

    Args:
    - I, mean_values, std_values, expected_values are ndarrays of shape (#measurements, 3), 
    containing applied current, experimentally estimated/expected mean values and standard deviations 
    for x,y,z-directions.
    - directory: valid path of directory where the image should be saved
    - data_filename_postfix: The image is saved as '%y_%m_%d_%H-%M-%S_'+ data_filename_postfix +'.png'

    """
    
    if directory is not None:
        ensure_dir_exists(directory, verbose=False)
        
    try:
        if len(I[0]) == 3:
            # depending on which function in main_menu.py was used to measure
            df = pd.DataFrame({'channel 1 [A]': I[:, 0],
                            'channel 2 [A]': I[:, 1],
                            'channel 3 [A]': I[:, 2],
                            'mean Bx [mT]': mean_data[:, 0],
                            'mean By [mT]': mean_data[:, 1],
                            'mean Bz [mT]': mean_data[:, 2],
                            'std Bx [mT]': std_data[:, 0],
                            'std By [mT]': std_data[:, 1],
                            'std Bz [mT]': std_data[:, 2],
                            'expected Bx [mT]': expected_fields[:, 0],
                            'expected By [mT]': expected_fields[:, 1],
                            'expected Bz [mT]': expected_fields[:, 2]})

    except:
        df = pd.DataFrame({'I (all Channels) [A]': I,
                           'mean Bx [mT]': mean_data[:, 0],
                           'mean By [mT]': mean_data[:, 1],
                           'mean Bz [mT]': mean_data[:, 2],
                           'std Bx [mT]': std_data[:, 0],
                           'std By [mT]': std_data[:, 1],
                           'std Bz [mT]': std_data[:, 2],
                           'expected Bx [mT]': expected_fields[:, 0],
                           'expected By [mT]': expected_fields[:, 1],
                           'expected Bz [mT]': expected_fields[:, 2]})

    now = datetime.now().strftime('%y_%m_%d_%H-%M-%S')
    output_file_name = '{}_{}.csv'.format(now, data_filename_postfix)
    file_path = os.path.join(directory, output_file_name)
    df.to_csv(file_path, index=False, header=True)



if __name__ == '__main__':
    
    CC_Z = ConexCC(com_port=z_COM_port, velocity=0.4, set_axis='z', verbose=False)
    CC_Z.wait_for_ready()
    
    start_pos_z = CC_Z.read_cur_pos()
    print(start_pos_z)
    CC_Z.move_absolute(new_pos=8.25)
    CC_Z.wait_for_ready()
    print(CC_Z.read_cur_pos())

    #     print(curr_pos_z)
    #     c = input('raise again?')
    
    