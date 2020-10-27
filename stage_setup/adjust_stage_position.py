"""
filename: adjust_stage_position.py

This script is meant to calibrate the measurement stage, such that the sensor in use is centered
and a defined distance above the pole tips. Specifically for use with the Hall Sensor cube.

Author: Nicholas Meinhardt (QZabre)
        nmeinhar@student.ethz.ch


Date: 09.10.2020
"""

# %%
########## Standard library imports ##########
import numpy as np
import serial
from time import sleep, time
import matplotlib.pyplot as plt
import os

########## local imports ##########
try:
    from modules.conexcc_control import setup, reset_to
except ModuleNotFoundError:
    import sys
    sys.path.insert(1, os.path.join(sys.path[0], '..'))
finally:
    from modules.conexcc_control import setup, reset_to, get_coords
    from conexcc.conexcc_class import *
    from modules.calibrate_cube import angle_calib_cube, av_single_sens, find_center_axis_with_cube
    from modules.calibration import find_center_axis
    from modules.plot_hall_cube import plot_angle_spherical, plot_angle
    from modules.serial_reader import get_new_data_set
    from MetrolabTHM1176.thm1176 import MetrolabTHM1176Node
    from modules.MetrolabMeasurements import get_mean_dataset_MetrolabSensor


# %%
# set measurement parameters and folder name
sampling_size = 25  # number of measurements per sensor for averaging

# %%
# initialize actuators
init_pos = np.array([5.0, 15.9, 8.25])
COM_ports = ['COM7', 'COM6', 'COM5']
CC_X, CC_Y, CC_Z = setup(init_pos, COM_ports=COM_ports)

#%%
# optionally only initialize actuator of z-axis
CC_Z = ConexCC(com_port='COM5',
                velocity=0.4, set_axis='z', verbose=False)
CC_Z.wait_for_ready()
CC_Z.move_absolute(8.0)

# %%
# manually adjust stage position
z_offset = 7.6
new_pos = [5.0, 15.9, z_offset]
_ = reset_to(new_pos, CC_X, CC2=CC_Y, CC3=CC_Z)

# %%
# set the bounds for x and y that are used during the calibration process
limits_x = [0, 15]
limits_y = [10, 23]

# set the bounds for x and y that are used during the calibration process, relative to mid position
# mid = [6.217,3.022]
# distance = 2
# limits_x = [mid[0] - distance, mid[0] + distance]
# limits_y = [mid[1] - distance, mid[1] + distance]

# set minimum set size, i.e. the precision of the calibration process
min_step_size = 1e-3
grid_number = 15
update_factor = 3

# choose between quick and less quick way 
extended = False

# %%
with MetrolabTHM1176Node() as node:

    # find center axis
    center_pos, fx = find_center_axis(CC_X, CC_Y, node, extended=extended, 
                                      limits_x=limits_x, limits_y=limits_y, grid_number=grid_number,
                                      min_step_size=min_step_size, update_factor=update_factor)

    print('\nEstimated center position: ({:.4f}, {:.4f})'.format(center_pos[0], center_pos[1]))
    start = np.append(center_pos, z_offset)
    reset_to(start, CC_X, CC2=CC_Y, CC3=CC_Z)

    field = get_mean_dataset_MetrolabSensor(node, sampling_size)[0]
    print('final field: {} mT'.format(np.round(field, 3)))
    print('in-plane field: {:.3f} mT'.format(np.sqrt(field[0]**2 + field[1]**2)))


#%%
# plot the current angles theta and phi in a spherical plot
with MetrolabTHM1176Node() as node:
    for _ in range(10):
        B = get_mean_dataset_MetrolabSensor(node, sampling_size)[0]
        plot_angle_spherical(B)




#%% 
# this part is only relevant when using Calibration Cube
#  -----------------------------------------------------------------


# # settings when using calibration cube
# specific_sensor = 55

# # initialize sensor
# port_sensor = 'COM4'

# establish permanent connection to calibration cube: open serial port; baud rate = 256000
# with serial.Serial(port_sensor, 256000, timeout=2) as cube:

#     # find center axis
#     center_pos, fx = find_center_axis_with_cube(CC_X, CC_Y, cube, specific_sensor=specific_sensor, extended=extended, 
#                                       limits_x=limits_x, limits_y=limits_y, grid_number=grid_number,
#                                       min_step_size=min_step_size, update_factor=update_factor)

#     print('\nEstimated center position: ({:.4f}, {:.4f})'.format(center_pos[0], center_pos[1]))
#     start = np.append(center_pos, z_offset)
#     reset_to(start, CC_X, CC2=CC_Y, CC3=CC_Z)

#     field = av_single_sens(cube, specific_sensor, 50)
#     print('final field: {} mT'.format(np.round(field, 3)))
#     print('in-plane field: {:.3f} mT'.format(np.sqrt(field[0]**2 + field[1]**2)))

#%%
# # plot the current angles theta and phi in a spherical plot
# with serial.Serial(port_sensor, 256000, timeout=2) as cube:
#     for _ in range(10):
#         B = get_new_data_set(cube=cube, specific_sensor=specific_sensor, no_enter=True)
#         plot_angle_spherical(B)

# # %%
# with serial.Serial(port_sensor, 256000, timeout=2) as cube:
#     B = av_single_sens(cube, specific_sensor, N)
# print(np.round(B,4))