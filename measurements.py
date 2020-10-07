"""
filename: measurements.py

This script is meant to be used to measure magnetic field values with the Hall
sensor cube. It is adapted to inteface with the ECB, e.g. to set current values and 
measure the generated magnetic field of the vector magnet.

Author: Nicholas Meinhardt, Maxwell Guerne-Kieferndorf (QZabre)
        nmeinhar@student.ethz.ch, gmaxwell@student.ethz.ch

Date: 07.10.2020
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
from modules.serial_reader import get_new_data_set


__all__ = [
        'measure',
        'makePlots',
        'saveDataPoints'
        ]

########## set measurement parameters and folder name ##########

N = 50  # number of measurements per sensor for averaging 
specific_sensor = 55

########## initialize sensor ##########
port_sensor = 'COM3'

# # initialize actuators
# init_pos = np.array([8.544, 4.256, 4.0])
# COM_ports = ['COM7','COM6','COM5']
# CC_X, CC_Y, CC_Z = setup(init_pos, COM_ports = COM_ports)

# # manually adjust z 
# z_offset = 4.0
# CC_Z.move_absolute(z_offset)


########## starts communication with hall sensor cube, measures the ##########
#          magnetic field with the specified sensor (change specific_sensor 
#          variable if necessary) 
#          returns: mean measurement data (averaged over N measurements)
#                   standard deviation in each averaged measurment 
#                   directory where the data will be saved to

def measure(folder_name='first_characterization_prototype_along_z'):

    # establish temporary connection to calibration cube: open serial port; baud rate = 256000
    with serial.Serial(port_sensor, 256000, timeout=2)  as cube: 
        # measure field with all sensors
        mean_data, std_data, _, directory = get_new_mean_data_set(N, filename=folder_name, cube=cube, 
                                                                        no_enter=True, on_stage=True)

    # see .\modules\calibrate_cube.py for more details on this function
    return mean_data[specific_sensor-1,:], std_data[specific_sensor-1,:], directory


def saveDataPoints(I, mean_data, std_data, expected_fields, directory, data_filename_postfix='B_field_vs_I'):
        # 'B_vs_I_in_plane'
        # save the results
        df = pd.DataFrame({ 'I [A]': I, 
                            'mean Bx [mT]': mean_data[:,0],
                            'mean By [mT]': mean_data[:,1],
                            'mean Bz [mT]': mean_data[:,2],
                            'std Bx [mT]': std_data[:,0],
                            'std By [mT]': std_data[:,1],
                            'std Bz [mT]': std_data[:,2],
                            'expected Bx [mT]': expected_fields[:,0],
                            'expected By [mT]': expected_fields[:,1],
                            'expected Bz [mT]': expected_fields[:,2]})

        now = datetime.now().strftime('%y_%m_%d_%H-%M-%S')
        output_file_name = '{}_{}.csv'.format(now, data_filename_postfix) 
        file_path = os.path.join(directory, output_file_name)
        df.to_csv(file_path, index=False, header=True)



# make plots of the data collected (most likely from the function 'measure')
def makePlots(I, mean_data, std_data, expected_fields):
    
    #%%
    # create a simple plot
    fig, axs = plt.subplots(3,1, sharex=True)
    fig.set_size_inches(6,6)

    # set axis labels
    axs[-1].set_xlabel('$I$ [A]')
    ylabels = ['$B_z$ [mT]', '$|B|$ [mT]', '$\\theta$ [°]']

    # calculate magnitudes of measured and expected magnetic field
    mean_magnitudes = np.linalg.norm(mean_data, axis=1)
    expected_magnitudes = np.linalg.norm(expected_fields, axis=1)

    # collect plot data. 
    # Note: errorbars display std, estimate errors for magnitudes (and angle) using propagation of uncertainty,
    # assuming that the measured fields in x,y,z direction are independent variables
    plot_mean_data = [mean_data[:,2],
                    mean_magnitudes,
                    np.arccos(mean_data[:,2]/mean_magnitudes)]

    plot_std_data = [std_data[:,2],
                    np.sqrt(np.einsum('ij,ij->i', mean_data**2, std_data**2))/mean_magnitudes,
                    np.zeros_like(mean_magnitudes)]

    plot_expected_data = [expected_fields[:,2],
                    expected_magnitudes,
                    np.zeros_like(expected_magnitudes)]

    # acutal plotting of Bz, |B|, theta-angle wrt z-axis
    for i in range(len(axs)):
        axs[i].set_ylabel(ylabels[i])

        axs[i].errorbar(I, plot_mean_data[i], yerr=plot_std_data[i],
                                linestyle='', marker='.', capsize = 2, label = 'measured @ ~ 3.5 mm')
        axs[i].plot(I, plot_expected_data[i], linestyle='--', marker='.', label = 'expected @ ~ 3 mm')
        axs[i].legend()

    plt.tight_layout()

    # save the figure
    output_file_name = '{}_{}.png'.format(now, data_filename_postfix) 
    file_path = os.path.join(directory, output_file_name)
    fig.savefig(file_path, dpi = 300)

