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
from time import time
import os

########## local imports ##########
try:
    from modules.conexcc_control import setup, reset_to
except ModuleNotFoundError:
    import sys
    sys.path.insert(1, os.path.join(sys.path[0], '..'))
finally:
    from modules.conexcc_control import setup, reset_to
    from modules.calibrate_cube import grid_2D

# %%
# set measurement parameters and folder name
sampling_size = 10 # number of measurements per sensor for averaging
specific_sensor = 55

directory = './test_data/2d_scans'


# %%
# initialize actuators
init_pos = np.array([7.8866, 0.0166, 15.8])
COM_ports = ['COM7', 'COM6', 'COM5']
CC_X, CC_Y, CC_Z = setup(init_pos, COM_ports=COM_ports)

# initialize sensor
port_sensor = 'COM4'

# %%
# manually adjust z
z_offset = 15.8
CC_Z.move_absolute(z_offset)

# %%
# set the bounds for x and y that are used during the calibration process
# limits_x = [5.0, 9.0]
# limits_y = [0.1, 4.0]

# set the bounds for x and y that are used during the calibration process, relative to mid position
mid = [7.8866, 0.0166]
distance = 2
limits_x = [mid[0] - distance, mid[0] + distance]
limits_y = [0, 4]

# set minimum set size, i.e. the precision of the calibration process
grid_number = 50


#%%
# establish permanent connection to calibration cube: open serial port; baud rate = 256000
with serial.Serial(port_sensor, 256000, timeout=2) as cube:

    # find center axis
    positions_corrected, B_field, filepath = grid_2D(CC_X, CC_Y, cube, specific_sensor, z_offset, 
                                      xlim=limits_x, ylim=limits_y, grid_number=grid_number,
                                      sampling_size=sampling_size, save_data=True, directory=directory)



# %%