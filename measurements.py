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


def calibration(node: MetrolabTHM1176Node):
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
        
    meas_offset_z = 8.55
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
            sleep(0.2)
            pos = CC.read_cur_pos()
            ratio = (abs(pos-start_pos) / total_distance)
            left = int(ratio * 30)
            right = 30-left
            print('\r[' + '#' * left + ' ' * right + ']', f' {ratio * 100:.0f}%', sep='', end='', flush=True)
    
    print('\n')


def measure(node: MetrolabTHM1176Node, dataDir=None, N=50, average=False):
    """
    starts communication with hall sensor cube, measures the magnetic field with the specified sensor (change specific_sensor variable if necessary)
    
    Args:
    - node (MetrolabTHM1176Node): represents the Metrolab THM 1176 sensor
    - dataDir: directory where measurements will be stored (entire path). Make sure that you know that the directory exists!
    - N: number of data points collected for each average
    - specific_sensor: sensor from which data will be fetched, only with hall sensor cube

    Returns: 
    - meas_time: measurement time (s)
    - meas_data: measured fields (x, y, z componenents)
    XOR
    - mean_data: mean measurement data of 'specific_sensor' (averaged over N measurements)
    - std_data: standard deviation in each averaged measurement
    (where mean, std are returned as ndarrays of shape (1, 3) for the 3 field directions)
    """
    if dataDir is not None:
        ensure_dir_exists(dataDir, verbose=False)
               
    if average:
        # measure average field at one point with the specific
        mean_data, std_data = get_mean_dataset_MetrolabSensor(node, sampling_size=N, directory=dataDir)
        ret1, ret2 = mean_data, std_data
    else:
        # This option is more for getting time-field measurements.
        meas_time, meas_data = readoutMetrolabSensor(node, measure_runs=N, directory=dataDir)
        # see .\modules\MetrolabMeasurements.py for more details on these function
        ret1, ret2 = meas_time, meas_data

    return ret1, ret2


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
    
    calibration(MetrolabTHM1176Node(unit="MT", sense_range_upper="0.1 T"))
    
