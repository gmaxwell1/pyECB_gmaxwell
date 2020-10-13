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
    from modules.conexcc_control import setup, reset_to
    from modules.calibrate_cube import find_center_axis, angle_calib, av_single_sens
    from modules.plot_hall_cube import plot_angle_spherical, plot_angle
    from modules.serial_reader import get_new_data_set


# %%
# set measurement parameters and folder name
N = 50  # number of measurements per sensor for averaging
specific_sensor = 55

# %%
# initialize actuators
init_pos = np.array([3.0, 0.1, 13.2])
COM_ports = ['COM7', 'COM6', 'COM5']
CC_X, CC_Y, CC_Z = setup(init_pos, COM_ports=COM_ports)

# initialize sensor
port_sensor = 'COM4'

# %%
# manually adjust z
z_offset = 11.5
CC_Z.move_absolute(z_offset)

# %%
# set the bounds for x and y that are used during the calibration process
limits_x = [5.5, 9.5]
limits_y = [0.2, 2.2]

# set minimum set size, i.e. the precision of the calibration process
min_step_size = 1e-3
grid_number = 10

# choose between quick and less quick way 
extended = True

# %%
# establish permanent connection to calibration cube: open serial port; baud rate = 256000
with serial.Serial(port_sensor, 256000, timeout=2) as cube:

    # find center axis
    center_pos, fx = find_center_axis(CC_X, CC_Y, cube, specific_sensor=specific_sensor, extended=extended, 
                                      limits_x=limits_x, limits_y=limits_y, grid_number=grid_number,
                                      min_step_size=min_step_size)

    print('\nEstimated center position: ({:.4f}, {:.4f})'.format(center_pos[0], center_pos[1]))
    start = np.append(center_pos, z_offset)
    reset_to(start, CC_X, CC2=CC_Y, CC3=CC_Z)

    field = av_single_sens(cube, specific_sensor, 50)
    print('final field: {} mT'.format(np.round(field, 3)))
    print('in-plane field: {:.3f} mT'.format(np.sqrt(field[0]**2 + field[1]**2)))


# %%
with serial.Serial(port_sensor, 256000, timeout=2) as cube:
    while True:
        B = get_new_data_set(cube=cube, specific_sensor=15, no_enter=True)
        plot_angle_spherical(B)

# %%
