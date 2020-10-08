#%%
import numpy as np
import serial
from time import sleep, time
import matplotlib.pyplot as plt
import os
import pandas as pd
from datetime import datetime

from modules.conexcc_control import *
from modules.calibrate_cube import get_new_mean_data_set, find_center_axis, angle_calib
from modules.plot_hall_cube import plot_many_sets, plot_stage_positions, plot_set, plot_sensor_positions
from modules.serial_reader import get_new_data_set
from modules.analysis_tools import estimate_std_magnitude, estimate_std_theta


#%%
# set measurement parameters and folder name 
N = 50  # number of measurements per sensor for averaging 
specific_sensor = 55
folder_name = 'first_characterization_prototype_along_z' #'first_characterization_prototype_plane'
data_filename_postfix = 'B_vs_I_along_z'    # 'B_vs_I_in_plane'

#%%
# initialize sensor 
port_sensor = 'COM4'

# # initialize actuators
# init_pos = np.array([8.544, 4.256, 4.0])
# COM_ports = ['COM7','COM6','COM5']
# CC_X, CC_Y, CC_Z = setup(init_pos, COM_ports = COM_ports)

# # manually adjust z 
# z_offset = 4.0
# CC_Z.move_absolute(z_offset)


#%%
# initialize arrays for data 
I = np.linspace(-5,10, num = 16) # insert actual current array here 

mean_data_specific_sensor = np.zeros((len(I),3))
std_data_specific_sensor = np.zeros((len(I),3))
expected_fields = np.zeros((len(I),3))


# establish permanent connection to calibration cube: open serial port; baud rate = 256000
with serial.Serial(port_sensor, 256000, timeout=2)  as cube: 

    for i in range(len(I)):
        # apply current I[i]
        # ...
        B_expected = [1,1,1] 


        # measure field with all sensors and create an image
        mean_data, std_data, _, directory = get_new_mean_data_set(N, filename=folder_name, cube=cube, 
                                                                    no_enter=True, on_stage=True)

        # collect measured and expected magnetic field          
        mean_data_specific_sensor[i] = mean_data[specific_sensor-1,:]
        std_data_specific_sensor[i] = std_data[specific_sensor-1,:]
        expected_fields[i] = B_expected


# save the results
df = pd.DataFrame({ 'I [A]': I, 
                    'mean Bx [mT]': mean_data_specific_sensor[:,0],
                    'mean By [mT]': mean_data_specific_sensor[:,1],
                    'mean Bz [mT]': mean_data_specific_sensor[:,2],
                    'std Bx [mT]': std_data_specific_sensor[:,0],
                    'std By [mT]': std_data_specific_sensor[:,1],
                    'std Bz [mT]': std_data_specific_sensor[:,2],
                    'expected Bx [mT]': expected_fields[:,0],
                    'expected By [mT]': expected_fields[:,1],
                    'expected Bz [mT]': expected_fields[:,2]})

now = datetime.now().strftime('%y_%m_%d_%H-%M-%S')
output_file_name = '{}_{}.csv'.format(now, data_filename_postfix) 
file_path = os.path.join(directory, output_file_name)
df.to_csv(file_path, index=False, header=True)


#%%
expected_fields= np.zeros_like(mean_data_specific_sensor)
expected_fields[:,2] = 0.017*I

def generate_plots(I, mean_values, std_values, expected_values, data_filename_postfix = 'B_vs_I', 
                        save_plot_values=False, save_image = True, ):
    """
    Generate plots of B vs I containing errorbars with mean values and standard deviations. 
    By default, the 

    Args: 
    - I, mean_values, std_values, expected_values are ndarrays of shape (#measurements, 3), 
    containing applied current, experimentally estimated/expected mean values and standard deviations 
    for x,y,z-directions.
    - data_filename_postfix: The image is saved as '%y_%m_%d_%H-%M-%S_'+ data_filename_postfix +'.png'
    - save_image: flag to save or not save the image
    - save_plot_values: flag to save or not save the data that are plotted to an additional file, 
    including the estimated errors  

    Return:
    
    """
    # create a simple plot
    fig, axs = plt.subplots(3,1, sharex=True)
    fig.set_size_inches(6,6)

    # set axis labels
    axs[-1].set_xlabel('$I$ [A]')
    ylabels = ['$B_z$ [mT]', '$|B|$ [mT]', '$\\theta$ [°]']


    # calculate magnitudes of measured and expected magnetic field
    mean_magnitudes = np.linalg.norm(mean_data_specific_sensor, axis=1)
    expected_magnitudes = np.linalg.norm(expected_fields, axis=1)

    # collect plot data. 
    # Note: errorbars display std, estimate errors for magnitudes (and angle) using propagation of uncertainty,
    # assuming that the measured fields in x,y,z direction are independent variables
    plot_mean_data = [mean_data_specific_sensor[:,2],
                    mean_magnitudes,
                    np.degrees(np.arccos(mean_data_specific_sensor[:,2]/mean_magnitudes))]

    plot_std_data = [std_data_specific_sensor[:,2],
                    estimate_std_magnitude(mean_data_specific_sensor, std_data_specific_sensor),
                    np.degrees(estimate_std_theta(mean_data_specific_sensor, std_data_specific_sensor))]

    plot_expected_data = [expected_fields[:,2],
                    expected_magnitudes,
                    np.degrees(np.arccos(expected_fields[:,2]/expected_magnitudes))]

    # actual plotting 
    for i in range(len(axs)):
        axs[i].set_ylabel(ylabels[i])

        axs[i].errorbar(I, plot_mean_data[i], yerr=plot_std_data[i],
                                linestyle='', marker='.', capsize = 2, label = 'measured @ ~ 3.5 mm')
        axs[i].plot(I, plot_expected_data[i], linestyle='--', marker='.', label = 'expected @ ~ 3 mm')
        axs[i].legend()

    # further adjustments
    axs[2].set_ylim((0,180))
    plt.tight_layout()

    # save the figure
    now = datetime.now().strftime('%y_%m_%d_%H-%M-%S')
    output_file_name = '{}_{}.png'.format(now, data_filename_postfix) 
    file_path = os.path.join(directory, output_file_name)
    fig.savefig(file_path, dpi = 300)

    fig.show()


