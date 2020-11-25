""" 
filename: calibration.py

The following functions can be used for preparing precise measurements with the Metrolab sensor.

Author: Nicholas Meinhardt (QZabre)
        nmeinhar@student.ethz.ch
        
Date: 20.10.2020
"""
# %%
# standard library imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import serial
from numpy.linalg import norm
import os
import sys
from datetime import datetime

# local imports
from modules.conexcc_control import all_ready, setup, check_no_motion, get_coords, check_validity
from conexcc.conexcc_class import *
from modules.plot_hall_cube import plot_angle, plot_angle_spherical
from modules.general_functions import save_in_dir, transform_between_sensor_stage_coordinates, ensure_dir_exists
from MetrolabTHM1176.thm1176 import MetrolabTHM1176Node
from modules.MetrolabMeasurements import get_mean_dataset_MetrolabSensor
from measurements import measure


# %%

def grid_2D(CC_X: ConexCC, CC_Y: ConexCC, node: MetrolabTHM1176Node, height, xlim=None, ylim=None,
            grid_number=50, sampling_size=10, verbose=False, save_data=False, suffix='2d_scan', directory=None):
    """
    Sweep over xy-plane and measure the field on each grid point, return the field values and coordinates.

    Note: 
    - The coordinates of sensor and actuators are not the same, this function accounts for this fact
    and returns all data in SENSOR COORDINATE SYSTEM.
    - It can happen that an invalid position that is out of the actuator's bounds should be reached.
    In this case, an according message is printed to the terminal and the actuator does not move. 

    Args:
    - CC_X, CC_Y (conexcc_class instances): represent x,y actuators
    - node (MetrolabTHM1176Node): represents the Metrolab THM 1176 sensor
    - height (float): z-component of the 2d-plane, this value will be carried over to all positions
    - xlim, ylim (None or tuple/list of floats): Either None or (min, max) as limits for x and y direction.
    If no limits are provided, the min and max values of the respective actuator is used.
    - grid_number (int): number of points +1 to sweep over per axis per iteration
    - sampling_size (int): Number of samples considered for the mean value of the field estimated 
    with the chosen sensor at a fixed position
    - verbose (bool): switching on/off print-statements for displaying progress
    - save_data (bool): flag to switch on/off saving of data
    - suffix (str): added at the end of the filename, by default '2d_scan'
    - directory (string): valid path of folder in which the data are stored

    Returns: positions_corrected, B_field
    - positions_corrected is an ndarray of shape ((grid_number+1)**2, 3) containing the stage positions 
    in sensor coordinates
    - B_field: is an ndarray of shape ((grid_number+1)**2, 3) containing the measured fields
    - file_path: path of the data file 
    which is the final distance between gird points
    """
    # set min and max limits for x,y axis
    if xlim == None:
        xmin = CC_X.min_limit
        xmax = CC_X.max_limit
    else:
        xmin = xlim[0]
        xmax = xlim[1]
    if ylim == None:
        ymin = CC_Y.min_limit
        ymax = CC_Y.max_limit
    else:
        ymin = ylim[0]
        ymax = ylim[1]

    # set step size accordingly
    x_step = (xmax-xmin) / grid_number
    y_step = (ymax-ymin) / grid_number

    # move actuator(or cube, respectively) to the minimum position and wait until controllers are ready again.
    CC_X.move_absolute(xmin)
    CC_Y.move_absolute(ymin)
    all_ready(CC_X, CC2=CC_Y, timeout=60)

    # set step size accordingly
    x_step = (xmax-xmin) / grid_number
    y_step = (ymax-ymin) / grid_number

    # collect some grid
    if verbose:
        print('initial position:({:.4f}, {:.4f})'.format(
            CC_X.read_cur_pos(), CC_Y.read_cur_pos()))
        print('\nstep size (x): {:.5f}\nstep size (y): {:.5f}'.format(
            x_step, y_step))
        print("\nNow sweeping along xy- plane: ")

    # initialize array to store the measured magnetic field components and the positions
    B_field = np.zeros((grid_number+1, grid_number+1, 3))
    positions = np.zeros((grid_number+1, grid_number+1, 3))

    # sweep along xy-plane and estimate the inplane field component
    for i in range(grid_number + 1):    # along x
        print('{} / {} rows'.format(i, grid_number))

        if i != 0:
            success = CC_X.move_relative(x_step)
            # if the actuator would have to move out of bounds, use absolute movement
            if not success:
                if CC_X.read_cur_pos() + x_step <= CC_X.min_limit:
                    CC_X.move_absolute(CC_X.min_limit)
                elif CC_X.read_cur_pos() + x_step >= CC_X.max_limit:
                    CC_X.move_absolute(CC_X.max_limit)

        for j in range(grid_number + 1):    # along y
            if j != 0:
                success = CC_Y.move_relative((-1)**i * y_step)
                # if the actuator would have to move out of bounds, use absolute movement
                if not success:
                    if CC_Y.read_cur_pos() + (-1)**i * y_step <= CC_Y.min_limit:
                        CC_Y.move_absolute(CC_Y.min_limit)
                    elif CC_Y.read_cur_pos() + (-1)**i * y_step >= CC_Y.max_limit:
                        CC_Y.move_absolute(CC_Y.max_limit)

            all_ready(CC_X, CC2=CC_Y, verbose=verbose)
            check_no_motion(CC_X, CC2=CC_Y, verbose=verbose)

            # save position and measured field value
            positions[i, j, 0:2] = get_coords(CC_X, CC2=CC_Y)[0]
            B_field[i, j, :] = measure(node, N=sampling_size, average=True)[0]

            if verbose:
                if i % 2 == 0:
                    print('Position ({},{}): ({:.4f}, {:.4f})'.format(i, j,
                                                                      (xmin + x_step * i), (ymin + y_step * j)))
                else:
                    print('Position ({},{}): ({:.4f}, {:.4f})'.format(i, j,
                                                                      (xmin + x_step * i), (ymax - y_step * j)))

    # coordinate systems of sensor and axes differ, account for this and reshuffle the position array
    # x and y must be changed and get an additional (-1)-factor each
    positions_corrected = transform_between_sensor_stage_coordinates(positions)

    # eventually reshape the arrays
    positions_corrected = positions_corrected.reshape((grid_number + 1)**2, 3)
    B_field = B_field.reshape((grid_number + 1)**2, 3)

    # add the z-component as a position
    positions_corrected[:, 2] = height

    # save data if desired
    if save_data:
        df = pd.DataFrame({'x [mm]': positions_corrected[:, 0],
                           'y [mm]': positions_corrected[:, 1],
                           'z [mm]': positions_corrected[:, 2],
                           'mean Bx [mT]': B_field[:, 0],
                           'mean By [mT]': B_field[:, 1],
                           'mean Bz [mT]': B_field[:, 2]})

        if directory is None:
            directory = os.getcwd()
        ensure_dir_exists(directory)

        now = datetime.now().strftime("%y_%m_%d_%H-%M-%S")
        output_file_name = "{}_{}.csv".format(now, suffix)
        data_filepath = os.path.join(directory, output_file_name)

        try:
            df.to_csv(data_filepath, index=False, header=True)
        except FileNotFoundError:
            data_filepath = os.path.join(os.getcwd(), output_file_name)
            df.to_csv(data_filepath, index=False, header=True)

    return positions_corrected, B_field, data_filepath


def search_with_Metrolab(CC1: ConexCC, CC2: ConexCC, node: MetrolabTHM1176Node, min_step_size=5*1e-3,
                         xlim=None, ylim=None, grid_number=10, sampling_size=10,
                         verbose=False, update_factor=1):
    """
    Search the (x,y)-positon with minimum magnetic field along the xy-plane for constant hight z
    using the Metrolab sensor. 

    Note: It can happen that an invalid position that is out of the actuator's bounds should be reached.
    In this case, an according message is printed to the terminal and the actuator does not move. 
    This might yield an incorrect calibration! 

    Args:
    - CC1, CC2 (conexcc_class instances): represent x,y actuators
    - node (MetrolabTHM1176Node): represents the Metrolab THM 1176 sensor
    - min_step_size (float): stop searching for minimum once the step size is smaller than this value.
    - xlim, ylim (None or tuple/list of floats): Either None or (min, max) as limits for x and y direction.
        If no limits are provided, the min and max values of the respective actuator is used.
    - grid_number (int): controls the stepsize when sweeping along x- and y-axis. The larger the number, 
      the smaller the stepsize. 
      NOTE: grid_number >= 2 is required to not raise errors, > 2 to not get trapped in while-loop!
    - sampling_size (int): Number of samples considered for the mean value of the field estimated 
      with the chosen sensor at a fixed position
    - verbose (bool): switching on/off print-statements for displaying progress
    - update_factor (float): after one iteration, the bounds xmin and xmax are updated
    as xmin = xpos - update_factor* xstep and xmax = xpos + update_factor* xstep (same for y).
    Thus, the smaller the update_factor, the faster one reaches the min_step_size. However,
    it becomes more likely to miss the actual minimum if update_factor is too small. 

    Returns: xpos, ypos, xmax/grid_number
    - xpos, ypos (float): final center positions along x and y at which the in-plane field is 
    supposed to be minimal 
    - precision (float): Final precision, which is the maximum of the final range of values along x and y 
    within which the center position is assumed devided by grid_number, 
    """
    # set min and max limits for x,y axis
    if xlim == None:
        xmin = CC1.min_limit
        xmax = CC1.max_limit
    else:
        xmin = xlim[0]
        xmax = xlim[1]
    if ylim == None:
        ymin = CC2.min_limit
        ymax = CC2.max_limit
    else:
        ymin = ylim[0]
        ymax = ylim[1]

    # move actuator(or cube, respectively) to the minimum position and wait until controllers are ready again.
    CC1.move_absolute(xmin)
    CC2.move_absolute(ymin)
    all_ready(CC1, CC2=CC2)

    # Sweep over plane by considering each a min and max x-value (y-value), find position with minimum and
    # redefine min/max values as this position +- one step size. As long as grid_number > 2 and
    # if min/max remain between allowed limits, the procedudure should converge (?)
    xrang = xmax-xmin
    yrang = ymax-ymin
    while ((xrang/grid_number) >= min_step_size) or ((yrang/grid_number) >= min_step_size):
        if verbose:
            print('\nstep size (x): {:.5f}\nstep size (y): {:.5f}'.format(
                xrang/grid_number, yrang/grid_number))

        # sweep along x-axis and measure B_y field component
        if verbose:
            print("\nNow sweeping x: ")
        fsx = []
        for i in range(int(grid_number + 1)):
            CC1.move_absolute(xrang/grid_number * i + xmin)
            all_ready(CC1, CC2=CC2)
            fsx.append(get_mean_dataset_MetrolabSensor(
                node, sampling_size)[0][1])

            if verbose:
                print('Position {}: {:.5f}'.format(
                    i, (xrang/grid_number * i + xmin)))
                print(
                    'in-plane field at position {}: {:.3f} mT'.format(i, fsx[-1]))

        # move to x-position with smallest (absolute value of) B_y
        fminx = np.argmin(abs(np.asarray(fsx)))
        xpos = (xmin+(xrang/grid_number)*fminx)
        CC1.move_absolute(xpos)
        all_ready(CC1, CC2=CC2)
        if verbose:
            print('Now moving to position {}: {:.5f}'.format(fminx, xpos))

        # update xmin, xmax as the x-position with smallest B_y +- one stepsize
        xmin = xpos - update_factor*(xrang/grid_number)
        xmax = xpos + update_factor*(xrang/grid_number)
        xmin, xmax = check_validity(xmin, xmax, CC1)
        xrang = xmax-xmin
        if verbose:
            print('New xmin: {:.5f} new xmax: {:.5f} new xrang: {:.5f}'.format(
                xmin, xmax, xrang))

        # sweep along y-axis and measure B_x field component
        if verbose:
            print("\nNow sweeping y: ")
        fsy = []
        for i in range(int(grid_number + 1)):
            CC2.move_absolute(yrang/grid_number * i + ymin)
            all_ready(CC1, CC2=CC2)
            fsy.append(get_mean_dataset_MetrolabSensor(
                node, sampling_size)[0][0])
            if verbose:
                print('Position {}: {:.5f}'.format(
                    i, (yrang/grid_number * i + ymin)))
                print(
                    'in-plane field at position {}: {:.3f} mT'.format(i, fsy[-1]))

        # move to x-position with smallest (absolute value of) B_x
        fminy = np.argmin(abs(np.asarray(fsy)))
        ypos = (ymin+(yrang/grid_number)*fminy)
        CC2.move_absolute(ypos)
        all_ready(CC1, CC2=CC2)
        if verbose:
            print('Now moving to position {}: {:.5f}'.format(fminy, ypos))

        # update ymin, ymax as the y-position with smallest B_x +- one stepsize
        ymin = ypos - update_factor*(yrang/grid_number)
        ymax = ypos + update_factor*(yrang/grid_number)
        ymin, ymax = check_validity(ymin, ymax, CC2)
        yrang = ymax-ymin
        if verbose:
            print('New ymin: {:.5f}  new ymax: {:.5f} new yrang: {:.5f}'.format(
                ymin, ymax, yrang))

    return xpos, ypos, np.max([xrang/grid_number, yrang/grid_number])


def find_center_axis(CC1, CC2, node: MetrolabTHM1176Node, sampling_size=10, min_step_size=5*1e-3,
                     limits_x=[0, 10], limits_y=[0, 10], grid_number=10, verbose=True,
                     extended=False, update_factor=1):
    """
    Find the xy-position of minimum in-plane magnetic field and estimate the field vector at this position,
    both using the Metrolab sensor.

    Args:
    - CC1, CC2 (conexcc_class instances): represent x,y actuators
    - node (MetrolabTHM1176Node): represents the Metrolab THM 1176 sensor
    - sampling_size (int): number of estimates of B-field used for average
    - min_step_size (float): stop searching for minimum of B-field along xy-plane 
    once the step size is smaller than this value.
    - limits_x, limits_x (array of length 2): minimum and maximum values of x (y) that are considered 
    when searching for the center position.
    - verbose (bool): switching on/off print-statements for displaying progress

    Return:
    - center_position = [xpos, ypos] (ndarray of floats): contains the position of the minimum in-plane field
    - field (float): measured B-field vector at x0
    """
    if extended:
        # xpos, ypos, precision = search_extended(CC1, CC2, cube, specific_sensor, min_step_size=min_step_size,
        #                            xlim=limits_x, ylim=limits_y, verbose=verbose, grid_number=grid_number,
        #                            update_factor=update_factor)
        raise NotImplementedError
    else:
        xpos, ypos, precision = search_with_Metrolab(CC1, CC2, node, min_step_size=min_step_size,
                                                     xlim=limits_x, ylim=limits_y, verbose=verbose,
                                                     grid_number=grid_number, update_factor=update_factor)
    center_position = np.array([xpos, ypos])
    field = get_mean_dataset_MetrolabSensor(node, sampling_size)[0]

    return center_position, field


def estimate_theta_error(node: MetrolabTHM1176Node, sampling_size=20):
    """
    Return the error in theta [degrees], which is the angle between the measured magnetic field 
    and the (magnet's) z axis. 

    Note: generate z-field with (1,1,1) configuration first

    Args:
    - node (MetrolabTHM1176Node): represents the Metrolab THM 1176 sensor
    - sampling_size (int): number of estimates of B-field used for average
    """
    field_vector = get_mean_dataset_MetrolabSensor(node, sampling_size)[0]
    mag = np.linalg.norm(field_vector)
    return 90 - np.degrees(np.arccos(field_vector[2]/mag))
