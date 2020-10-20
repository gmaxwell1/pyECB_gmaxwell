"""
filename: 2d_scan.py

This script is meant to perform a 2d scan of the magnetic field using a single sensor,
specifically for use with the Hall Sensor cube.

Author: Nicholas Meinhardt (QZabre)
        nmeinhar@student.ethz.ch


Date: 13.10.2020
"""

# %%
########## Standard library imports ##########
import numpy as np
import serial
from time import time, sleep
import os

########## local imports ##########
try:
    from modules.conexcc_control import setup, reset_to
except ModuleNotFoundError:
    import sys
    sys.path.insert(1, os.path.join(sys.path[0], '..'))
finally:
    from modules.conexcc_control import setup, reset_to
    from modules.calibrate_cube import grid_2D_cube
    from modules.calibration import grid_2D
    from MetrolabTHM1176.thm1176 import MetrolabTHM1176Node

# %%
# set measurement parameters and folder name
sampling_size = 20 # number of measurements per sensor for averaging

directory = './test_data/2d_scans'

# number of grid points per dimension
grid_number = 25

# %%
# initialize actuators
init_pos = np.array([0., 0, 8])
COM_ports = ['COM7', 'COM6', 'COM5']
CC_X, CC_Y, CC_Z = setup(init_pos, COM_ports=COM_ports)


# %%
# manually adjust stage position
z_offset = 8
new_pos = [0, 0, z_offset]
_ = reset_to(new_pos, CC_X, CC2=CC_Y, CC3=CC_Z)

# %%
# set the bounds for x and y that are used for the scan
limits_x = [0.0, 25]
limits_y = [0.0, 25]

# set the bounds for x and y that are used for the scan, relative to mid position
# mid = [7.8866, 0.0166]
# distance = 2
# limits_x = [mid[0] - distance, mid[0] + distance]
# limits_y = [0, 4]


#%%
# perform actual 2d scan
with MetrolabTHM1176Node() as node:

    positions_corrected, B_field, filepath = grid_2D(CC_X, CC_Y, node, z_offset, 
                                      xlim=limits_x, ylim=limits_y, grid_number=grid_number,
                                      sampling_size=sampling_size, save_data=True, directory=directory)

#%%
# this part uses the Calibration Cube as Sensor
# --------------------------------------------------------------------

# # initialize sensor
# specific_sensor = 55
# port_sensor = 'COM4'

# # establish permanent connection to calibration cube: open serial port; baud rate = 256000
# with serial.Serial(port_sensor, 256000, timeout=2) as cube:

#     positions_corrected, B_field, filepath = grid_2D_cube(CC_X, CC_Y, cube, specific_sensor, z_offset, 
#                                       xlim=limits_x, ylim=limits_y, grid_number=grid_number,
#                                       sampling_size=sampling_size, save_data=True, directory=directory)

