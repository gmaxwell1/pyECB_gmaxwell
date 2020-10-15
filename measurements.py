"""
filename: measurements.py

This script is meant to be used to measure magnetic field values with the Hall
sensor cube. It is adapted to inteface with the ECB, e.g. to set current values and 
measure the generated magnetic field of the vector magnet.

Author: Nicholas Meinhardt, Maxwell Guerne-Kieferndorf (QZabre)
        nmeinhar@student.ethz.ch, gmaxwell@student.ethz.ch

Date: 08.10.2020
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
from modules.conexcc_control import *
from modules.calibrate_cube import get_new_mean_data_set, find_center_axis, angle_calib
from modules.plot_hall_cube import plot_many_sets, plot_stage_positions, plot_set, plot_sensor_positions
from modules.serial_reader import get_new_data_set, ensure_dir_exists


__all__ = [
    'newMeasurementFolder',
    'measure',
    'makePlots',
    'saveDataPoints'
]

########## sensor cube port ##########
port_sensor = 'COM3'


def newMeasurementFolder(defaultDataDir='data_sets', sub_dir_base='z_field_meas', verbose=False):
    """
    This function is supposed to create a new (unique) subdirectory to store data for a measurement run.
    """
    index = 1
    cwd = os.getcwd()
    sub_dirname = sub_dir_base + '_' + str(index)
    dataDir = os.path.join(cwd, defaultDataDir, sub_dirname)
    # iterate though postfixes until a new directory name is found
    while ensure_dir_exists(dataDir, verbose=verbose):
        index = index + 1
        sub_dirname = sub_dir_base + '_' + str(index)
        dataDir = os.path.join(cwd, defaultDataDir, sub_dirname)

    return sub_dirname


def measure(sub_dirname='z_field_meas_set',on_stage=True, N=50, specific_sensor=55):
    """
    starts communication with hall sensor cube, measures the magnetic field with the specified sensor (change specific_sensor variable if necessary)

    Args:
    -subdirname: folder where measurements will be stored
    -specific_sensor: sensor from which data will be fetched
    -N: number of data points collected for each average

    Returns: 
    -mean measurement data of 'specific_sensor' (averaged over N measurements)
    -standard deviation in each averaged measurment
    (where mean, std and abs(std/mean) are returned as ndarrays of shape (1, 3) for the 3 field directions)
    -directory where the data will be saved to
    """

    # establish temporary connection to calibration cube: open serial port; baud rate = 256000
    with serial.Serial(port_sensor, 256000, timeout=2) as cube:
        # measure field with all sensors
        mean_data, std_data, _, directory = get_new_mean_data_set(N, sub_dirname, cube, no_enter=Ture, on_stage=True)
        # if on_stage:
        #     resp = 1
        #     path = ''
        #     while resp == 1:
        #         # N measurements of all 64 sensors are made and saved to a csv file.
        #         resp, directory, csvfile = get_new_data_set(measure_runs=N, sub_dirname=sub_dirname, cube=cube, no_enter=True,verbose=False, on_stage=on_stage)
        #         path = os.path.join(directory, csvfile)
        # else:
        #     path = get_new_data_set(measure_runs=N, sub_dirname=sub_dirname, cube=cube, no_enter=True, verbose=False, on_stage=on_stage)
    
    # return path
    # see .\modules\calibrate_cube.py for more details on this function
    return mean_data[specific_sensor-1, :], std_data[specific_sensor-1, :], directory


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


# make plots of the data collected (most likely from the function 'measure')
def makePlots(I, mean_data, std_data, expected_fields):

    # %%
    # create a simple plot
    fig, axs = plt.subplots(3, 1, sharex=True)
    fig.set_size_inches(6, 6)

    # set axis labels
    axs[-1].set_xlabel('$I$ [A]')
    ylabels = ['$B_z$ [mT]', '$|B|$ [mT]', '$\\theta$ [°]']

    # calculate magnitudes of measured and expected magnetic field
    mean_magnitudes = np.linalg.norm(mean_data, axis=1)
    expected_magnitudes = np.linalg.norm(expected_fields, axis=1)

    # collect plot data.
    # Note: errorbars display std, estimate errors for magnitudes (and angle) using propagation of uncertainty,
    # assuming that the measured fields in x,y,z direction are independent variables
    plot_mean_data = [mean_data[:, 2],
                      mean_magnitudes,
                      np.arccos(mean_data[:, 2]/mean_magnitudes)]

    plot_std_data = [std_data[:, 2],
                     np.sqrt(np.einsum('ij,ij->i', mean_data **
                                       2, std_data**2))/mean_magnitudes,
                     np.zeros_like(mean_magnitudes)]

    plot_expected_data = [expected_fields[:, 2],
                          expected_magnitudes,
                          np.zeros_like(expected_magnitudes)]

    # acutal plotting of Bz, |B|, theta-angle wrt z-axis
    for i in range(len(axs)):
        axs[i].set_ylabel(ylabels[i])

        axs[i].errorbar(I, plot_mean_data[i], yerr=plot_std_data[i],
                        linestyle='', marker='.', capsize=2, label='measured @ ~ 3.5 mm')
        axs[i].plot(I, plot_expected_data[i], linestyle='--',
                    marker='.', label='expected @ ~ 3 mm')
        axs[i].legend()

    plt.tight_layout()

    # save the figure
    output_file_name = '{}_{}.png'.format(now, data_filename_postfix)
    file_path = os.path.join(directory, output_file_name)
    fig.savefig(file_path, dpi=300)


if __name__ == '__main__':
    
    subDir = newMeasurementFolder(sub_dir_base='z_field_meas')
    t1 = time()
    where = measure(subDir, N=5)
    t2 = time()
    print('time: ', t2-t1)
    print(where)

