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
    from main_comm import *

# %%
# set measurement parameters and folder name
sampling_size = 15 # number of measurements per sensor for averaging

directory = './data_sets/2d_scans_different_fields/set3'

# number of grid points per dimension
grid_number = 20

# %%
# initialize actuators
init_pos = np.array([3, 13, 8.3])
# ports for Magnebotix PC
COM_ports = ['COM4', 'COM5', 'COM6']
CC_X, CC_Y, CC_Z = setup(init_pos, COM_ports=COM_ports)


# %%
# manually adjust stage position
# z_offset = 8.3
# new_pos = 
# _ = reset_to(new_pos, CC_X, CC2=CC_Y, CC3=CC_Z)


# %%
# set the bounds for x and y that are used for the scan
limits_x = [3.5, 6.1]
limits_y = [1.2, 3.8]

# set the bounds for x and y that are used for the scan, relative to mid position
# mid = [7.8866, 0.0166]
# distance = 2
# limits_x = [mid[0] - distance, mid[0] + distance]
# limits_y = [0, 4]
desCurrents = [0] * 8
# set currents in coils
openConnection()
currentConfig = [1,1,1]
# integer value, mA
currentStrength = 2000
for i in range(len(currentConfig)):
        desCurrents[i] = currentStrength * currentConfig[i]
enableCurrents()
demagnetizeCoils()
sleep(0.3)
setCurrents(desCurrents, b'1')

#%%
# perform actual 2d scan
with MetrolabTHM1176Node(block_size=30, period=0.01, range='0.3 T', average=1) as node:
# node = MetrolabTHM1176Node(block_size=30, period=0.01, range='0.3 T', average=1)
    # CC_Y.move_absolute(0.0)
    # CC_Z.move_absolute(21.0)
    # state = False
    # while not state:
    #     state = CC_Z.is_ready() and CC_Y.is_ready()

    # enter = input('press enter to calibrate')
    # if enter == '':
    #     node.calibrate()

    # input('press any key to continue')

    # CC_Y.move_absolute(13.0)
    # CC_Z.move_absolute(z_offset)
    # state = False
    # while not state:
    #     state = CC_Z.is_ready() and CC_Y.is_ready()
    filename_suffix = f'2d_scan_({currentConfig[0]}_{currentConfig[1]}_{currentConfig[2]})'
    positions_corrected, B_field, filepath = grid_2D(CC_X, CC_Y, node, 8.3, xlim=limits_x, ylim=limits_y, grid_number=grid_number,
                                                     sampling_size=sampling_size, save_data=True,suffix=filename_suffix, directory=directory)
disableCurrents()

closeConnection()

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

