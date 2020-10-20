""" 
filename: calibrate_cube.py

The following functions can be used for preparing precise measurements with the Hall Sensor Cube.

Author: Jona Buehler 2020

Documentation and Updates by Nicholas Meinhardt (Qzabre)
                             nmeinhar@student.ethz.ch
        
Date: 09.10.2020
"""
#%%
########## Standard library imports ##########
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import serial
from numpy.linalg import norm
import os
import sys
from datetime import datetime

########## local imports ##########
from modules.serial_reader import get_new_data_set, direct_readout, MeasurementError
from modules.conexcc_control import all_ready, setup, check_no_motion, get_coords, check_validity
from conexcc.conexcc_class import *
from modules.plot_hall_cube import plot_angle, plot_angle_spherical
from modules.general_functions import angle_wrt_z, transform_between_sensor_stage_coordinates, save_in_dir, ensure_dir_exists


# %%
def search_with_cube(CC1: ConexCC, CC2: ConexCC, cube, specific_sensor, min_step_size=5*1e-3, xlim=None, ylim=None,
           grid_number=10, sampling_size=10, verbose=False, update_factor=1):
    """
    Search the (x,y)-positon with minimum magnetic field along the xy-plane for constant hight z. 

    Note: It can happen that an invalid position that is out of the actuator's bounds should be reached.
    In this case, an according message is printed to the terminal and the actuator does not move. 
    This might yield an incorrect calibration! 

    Args:
    - CC1, CC2 (conexcc_class instances): represent x,y actuators
    - cube (serial.Serial): represents the magnetic field sensor.
    - specific_sensor (int in [1,64]): ID of a specific Hall sensor of the whole cube.
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
            fsx.append(av_single_sens(cube, specific_sensor, sampling_size)[1])

            if verbose:
                print('Position {}: {:.5f}'.format(
                    i, (xrang/grid_number * i + xmin)))
                print('Field at position {}: {:.3f}'.format(i, fsx[-1]))

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
            fsy.append(av_single_sens(cube, specific_sensor, sampling_size)[0])
            if verbose:
                print('Position {}: {:.5f}'.format(
                    i, (yrang/grid_number * i + ymin)))
                print('Field at position {}: {:.3f}'.format(i, fsy[-1]))

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

def search_extended(CC_X: ConexCC, CC_Y: ConexCC, cube, specific_sensor, min_step_size=5*1e-3, xlim=None, ylim=None,
           grid_number=10, sampling_size=10, verbose=False, update_factor=1.0):
    """
    Search the (x,y)-positon with minimum magnetic field along the xy-plane for constant hight z.
    In contrast to the search-function, all points on the grid are measured to find the minimum, 
    which takes longer than only searching along one axis. 

    Note: It can happen that an invalid position that is out of the actuator's bounds should be reached.
    In this case, an according message is printed to the terminal and the actuator does not move. 
    This might yield an incorrect calibration! 

    Args:
    - CC_X, CC_Y (conexcc_class instances): represent x,y actuators
    - cube (serial.Serial): represents the magnetic field sensor.
    - specific_sensor (int in [1,64]): ID of a specific Hall sensor of the whole cube.
    - min_step_size (float): stop searching for minimum once the step size is smaller than this value.
    - xlim, ylim (None or tuple/list of floats): Either None or (min, max) as limits for x and y direction.
    If no limits are provided, the min and max values of the respective actuator is used.
    - grid_number (int): number of points to sweep over per axis per iteration, thereby controlling the step size of sweeps.
    NOTE: > 2 to not get trapped in while-loop!
    - sampling_size (int): Number of samples considered for the mean value of the field estimated 
      with the chosen sensor at a fixed position
    - verbose (bool): switching on/off print-statements for displaying progress
    - update_factor (float): after one iteration, the bounds xmin and xmax are updated
    as xmin = xpos - update_factor* xstep and xmax = xpos + update_factor* xstep (same for y).
    Thus, the smaller the update_factor, the faster one reaches the min_step_size. However,
    it becomes more likely to miss the actual minimum if update_factor is too small. 

    Returns: xpos, ypos, precision
    - xpos, ypos: 
    - precision: maximum of xrang/grid_number and yrang/grid_number, 
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

    # move actuator(or cube, respectively) to the minimum position and wait until controllers are ready again.
    CC_X.move_absolute(xmin)
    CC_Y.move_absolute(ymin)
    all_ready(CC_X, CC2=CC_Y, timeout=60)

    print('initial position:({:.4f}, {:.4f})'.format(CC_X.read_cur_pos(),CC_Y.read_cur_pos()))


    # Sweep over plane by considering each a min and max x-value (y-value), find position with minimum and
    # redefine min/max values as this position +- one step size. As long as num_prec > 2 and
    # if min/max remain between allowed limits, the procedudure should converge (?)
    xrang = xmax-xmin
    yrang = ymax-ymin
    while ((xrang/grid_number) >= min_step_size) or ((yrang/grid_number) >= min_step_size):

        # define new step size
        x_step = xrang/grid_number
        y_step = yrang/grid_number
        if verbose:
            print('\nstep size (x): {:.5f}\nstep size (y): {:.5f}'.format(x_step, y_step))
        
        # initialize array to store the in-plane magnitudes of magnetic field
        inplane_magnitude = np.zeros((grid_number+1, grid_number+1))

        # sweep along xy-plane and estimate the inplane field component
        if verbose:
            print("\nNow sweeping along xy- plane: ")
        for i in range(grid_number + 1):    # along x
            if i != 0:
                CC_X.move_relative(x_step)
            for j in range(grid_number + 1):    # along y     
                if j != 0:
                    CC_Y.move_relative((-1)**i * y_step)
                all_ready(CC_X, CC2=CC_Y)

                field = av_single_sens(cube, specific_sensor, sampling_size)
                inplane_magnitude[i,j] = np.sqrt(field[0]**2 + field[1]**2)

                if verbose:
                    if i % 2 ==0: 
                        print('Position ({},{}): ({:.4f}, {:.4f})'.format(i,j, 
                                                (xmin + x_step * i), (ymin + y_step * j)))
                    else:
                        print('Position ({},{}): ({:.4f}, {:.4f})'.format(i,j, 
                                                (xmin + x_step * i), (ymax - y_step * j)))
                    print('measured field: {:.3f} mT'.format(np.sqrt(field[0]**2 + field[1]**2)))

        # find position of minimum magnitude of in-plane field 
        i_min, j_min = np.unravel_index(np.argmin(inplane_magnitude), inplane_magnitude.shape)

        # move to minimum position
        xpos = xmin+ i_min * x_step
        if i % 2 == 0:
            ypos = ymin + j_min * y_step
        else:
            ypos = ymax - j_min * y_step
        CC_X.move_absolute(xpos)
        CC_Y.move_absolute(ypos)
        if verbose:
            print('Now moving to position ({},{}): ({:.4f}, {:.4f})'.format(i_min, j_min, xpos, ypos))

        # update boundaries, such that the neighboring grid points become minimum/maximum values
        xmin = xpos - update_factor*x_step
        xmax = xpos + update_factor*x_step
        xmin, xmax = check_validity(xmin, xmax, CC_X)
        xrang = xmax - xmin
        ymin = ypos - update_factor*y_step
        ymax = ypos + update_factor*y_step
        ymin, ymax = check_validity(ymin, ymax, CC_Y)
        yrang = ymax - ymin

        if verbose:
            print('New xmin: {:.4f} new xmax: {:.4f} new xrang: {:.4f}'.format(xmin, xmax, xrang))
            print('New ymin: {:.4f} new ymax: {:.4f} new yrang: {:.4f}'.format(ymin, ymax, yrang))

    return xpos, ypos, np.max([xrang, yrang])/grid_number


def av_single_sens(cube, specific_sensor, N, max_number_attempts=10):
    """
    Measure magnetic field vector N times using a specific sensor and return mean values.

    It may happen that no valid results are returned from the sensor. In this case, still try to gather N
    valid samples by repeating the faulty measurement up to a maximum of max_number_attempts repetitions.

    Note: This function measures the magnetic field without specifying the position. It is recommended
    to ensure that the sensor is not moving beforehand.

    Args:
    - cube (serial.Serial): represents the magnetic field sensor
    - specific_sensor (int): number of a desired sensor, specific_sensor in [1, 64]  
    - N (int): number of B-field estimates that are averaged
    - max_number_attempts (int): maximum number of attempts to replace faulty measurements by valid outcomes.

    Returns 1d-ndarray with 3 entries, corresponding to the mean magnetic field vector at this position.
    """
    field_samples = []
    for _ in range(N):
        field = get_new_data_set(cube=cube, specific_sensor=specific_sensor)
        # if something does not work, get_new_data might return 1 instead of vector
        if not isinstance(field, int):
            field_samples.append(field)
    # if there are too few data because of a failure, get more
    if len(field_samples) < N:
        number_attempts = max_number_attempts
        while len(field_samples) < N and number_attempts <= max_number_attempts:
            field = get_new_data_set(
                cube=cube, specific_sensor=specific_sensor)
            # if something does not work, get_new_data might return 1 instead of vector
            if not isinstance(field, int):
                field_samples.append(field)

    # return the mean field vector
    return np.mean(field_samples, axis=0)


def find_center_axis_with_cube(CC1, CC2, cube, N=10, min_step_size=5*1e-3, specific_sensor=54, limits_x=[0, 10],
                     limits_y=[0, 10], grid_number = 10, verbose=True, extended=True, update_factor=1):
    """
    Find the xy-position of minimum in-plane magnetic field and estimate the field vector at this position.

    Args:
    - CC1, CC2 (conexcc_class instances): represent x,y actuators
    - cube (serial.Serial): represents the magnetic field sensor.
    - N (int): number of estimates of B-field used for average
    - min_step_size (float): stop searching for minimum of B-field along xy-plane 
    once the step size is smaller than this value.
    - specific_sensor (int in [1,64]): number of sensor that is used when searching the center position. 
    - limits_x, limits_x (array of length 2): minimum and maximum values of x (y) that are considered 
    when searching for the center position.
    - verbose (bool): switching on/off print-statements for displaying progress

    Return:
    - x0 = [xpos, ypos] (ndarray of floats): contains the position of the minimum in-plane field
    - f0 (float): measured B-field vector at x0
    """
    if verbose:
        print("\nTrying to find center axis using sensor # {} ...\n".format(
            specific_sensor))
    if extended:
        xpos, ypos, precision = search_extended(CC1, CC2, cube, specific_sensor, min_step_size=min_step_size,
                                   xlim=limits_x, ylim=limits_y, verbose=verbose, grid_number=grid_number,
                                   update_factor=update_factor)
    else:
        xpos, ypos, precision = search_with_cube(CC1, CC2, cube, specific_sensor, min_step_size=min_step_size, 
                                   xlim=limits_x, ylim=limits_y, verbose=verbose, grid_number=grid_number,
                                   update_factor=update_factor)
    x0 = np.array([xpos, ypos])
    f0 = av_single_sens(cube, specific_sensor, N)

    return x0, f0

def angle_calib_cube(desired, cube, specific_sensor=54, N=10, visual_feedback=True, eps=0.5,
                max_number_trials=100, spherical=True, verbose=True):
    """
    Compare measured to desired angle in a while-loop. 
    Leave the loop when the difference between both angles is less than eps degrees 
    or when the maximum number of trials is reached

    Note: Pass all angles in degrees!

    Args: 
    - desired (float): desired angle of magnetic field with respect to z-axis
    - cube (serial.Serial): represents the magnetic field sensor.
    - specific_sensor (int): id in [1,64] of sensor that should be used for calibration
    - N (int): sampling size of magnetic field estimates used for average
    - visual_feedback (bool): flag to switch on/off a plotted visualization of measured vector,
    showing an normalized vector and the measured angle with z-axis
    - eps (float): acceptable difference in angle
    - max_number_trials (int): to avoid an infinite while-loop, a maximum number of trials can be set.
    - spherical (bool): flag to switch on spherical plot, else 'normal' 3d plot is shown.
    - verbose (bool): switching on/off print-statements for displaying progress

    Return 0 
    """
    diff = 5*eps
    i = 0
    while abs(diff) > eps and i < max_number_trials:
        vec = av_single_sens(cube, specific_sensor, N)
        ang = angle_wrt_z(vec)
        diff = desired-np.degrees(ang)
        if visual_feedback and spherical:
            plot_angle_spherical(vec)
        elif visual_feedback:
            plot_angle(vec)
        i += 1
        if verbose:
            print("\r Angle to z axis: {:.2f} °; Difference to desired angle ({:.2f} °) is {:.2f} °".format(
                np.degrees(ang), desired, diff))
    if verbose:
        print("Calibration of angle successfull!")
    return 0

def get_new_mean_data_set(measure_runs, cube, specific_sensor=None, omit_64=False, verbose=False, max_num_retrials=5,
                            save_raw_data= False, save_mean_data=False, directory=None):
    """
    Estimate field vectors with all sensors measure_runs-times and return the mean and std 
    for the specific_sensor only or for all sensors (default). 

    Note: Sensor #64 has produced incorrect results in past, thus it can be omitted using the omit_True flag

    Args: 
    - measure_runs (int): sampling size to estimate mean magnetic field vector and std, 
    i.e. number of times all sensor are read out in series before averaging 
    - cube (serial.Serial): represents the magnetic field sensor.
    - specific_sensor (None or int in [1,64]): If None, the mean values of all sensors are returned. 
    If a number is provided, only the data from the sensor with this ID will be returned. 
    - directory (string): valid path of the folder where data files can be stored (provided save_mean_data=True)
    - max_num_retrials (int): when errors occur during the serial measurement of all sensors,
    the measurement is repeated up to max_number_retrials times before raising a MeasurementError exception.
    - save_raw_data (bool): if True, the raw data estimated during each measurement round for each sensor are saved as a csv-file.
    - save_mean_data (bool): if True, the mean data averaged over all measurement rounds for each sensor are saved as a csv-file.
    - omit_64 (bool): if True, sensor #64 will be omitted, else all 64 sensors are considered.
    - verbose (bool): switching on/off print-statements for displaying progress

    Return: 
    - mean_data, std_data (ndarrays): Mean magnetic field and its standard deviation as a vector, 
    either for all sensors or for a single sensor only. If specific_sensor=None, both arrays are of shape (number sensors, 3),
    else they are 1d arrays of length 3.

    Exceptions:
    - MeasurementError, if no raw data could be aquired after max_num_retrials repetitions. 
    """
    # perform measurement and collect the raw data 
    for _ in range(max_num_retrials):
        try:
            _, meas_data = direct_readout(cube, measure_runs=measure_runs, save_data=save_raw_data, directory=directory,
                                    fname_postfix='raw_data', verbose=verbose, omit_64=omit_64)
        except MeasurementError:
            pass
        else:
            break
    
    # estimate the mean and std from raw data for each sensor
    try:
        mean_data = np.mean(meas_data, axis=0)
        std_data = np.std(meas_data, axis=0)
    # if it was not possible to obtain valid measurement results after max_num_retrials, raise MeasurementError, too
    except UnboundLocalError:
        raise MeasurementError
    
    # save mean data if desired
    if save_mean_data:
        # if no directory is provided, just take current working directory
        if directory is None:
            directory = os.getcwd()
        save_in_dir(mean_data, directory, 'data', stds=std_data, now=True)

    # depending on the specific_sensor flag, return mean and std fields either for only this sensor or for all sensors
    if specific_sensor is not None:
        return mean_data[specific_sensor, :], std_data[specific_sensor, :]
    else:
        return mean_data, std_data


def get_new_mean_data_set_old(N, sub_dirname=None, cube=None, no_enter=False, on_stage=False, omit_64=False,
                          verbose=False):
    """
    Estimate field vectors N-times with all sensors and calculate mean, std and abs(std/mean) as vectors for each sensor.

    Note: Sensor #64 has produced incorrect results in past, thus it can be omitted using the omit_True flag

    Args: 
    - N (int): number of times all 64 (or 63) sensor are read out in series
    - cube (serial.Serial): represents the magnetic field sensor.
    - sub_dirname (string): name of folder where data files are stored
    - no_enter (bool): if True, measurement starts automatically, else the user is asked to press enter to start.
    - on_stage (bool): flag used to set the action upon occurence an error when reading a measurement outcome 
          from the sensor. If False, continue measuring and write a "'Read Failure', 0,0,0,0"-line to file. 
          If True, the measurement is stopped entirely and the output file is deleted. 
    - omit_64 (bool): if True, sensor #64 will be omitted, else all 64 sensors are considered.
    - verbose (bool): switching on/off print-statements for displaying progress

    Return:
    - if on_stage=True: mean_data, std_data, perc_data, directory 
    (where mean, std and abs(std/mean) are returned as ndarrays of shape (number_sensors, 3) for the 3 field directions)
    - else: mean_data, std_data, perc_data
    """
    # measure field using all sensors N times and save results to csv file
    if on_stage:
        resp = 1
        path = ''
        while resp == 1:
            resp, directory, csvfile = get_new_data_set(measure_runs=N, sub_dirname=sub_dirname, cube=cube, verbose=verbose,
                                                        no_enter=no_enter, on_stage=on_stage, omit_64=omit_64)
            path = os.path.join(directory, csvfile)
    else:
        path = get_new_data_set(measure_runs=N, sub_dirname=sub_dirname, cube=cube, no_enter=no_enter,
                                on_stage=on_stage, omit_64=omit_64, verbose=verbose)

    # import the measurement data from csv file
    dataD = pd.read_csv(path)
    if sys.version_info[0] == 3:
        data = dataD.to_numpy()
    else:
        data = dataD.values

    # initialize arrays for measurement outcomes, adapt length of arrays depending on omit_64 flag
    if omit_64:
        number_sensors = 63
    else:
        number_sensors = 64
    x_mean = np.zeros(number_sensors)
    x_std = np.zeros(number_sensors)
    x_perc = np.zeros(number_sensors)
    y_mean = np.zeros(number_sensors)
    y_std = np.zeros(number_sensors)
    y_perc = np.zeros(number_sensors)
    z_mean = np.zeros(number_sensors)
    z_std = np.zeros(number_sensors)
    z_perc = np.zeros(number_sensors)

    # collect results for each sensor and save mean, std and abs(std/mean)
    for i in range(number_sensors):
        sensor_i_data = np.zeros((N, 5))
        for k in range(N):
            if (data[i+k*number_sensors, :].dtype == 'float64'):
                sensor_i_data[k, :] = data[i+k*number_sensors, :]
            else:
                # print this message in any case!
                print("could not convert data properly! wrong data type: ",
                      data[i+k*number_sensors, :].dtype)
                sensor_i_data[k, :] = 0
        x_mean[i] = np.mean(sensor_i_data[:, 2])
        x_std[i] = np.std(sensor_i_data[:, 2])
        x_perc[i] = np.divide(x_std[i], abs(x_mean[i]))
        y_mean[i] = np.mean(sensor_i_data[:, 3])
        y_std[i] = np.std(sensor_i_data[:, 3])
        y_perc[i] = np.divide(y_std[i], abs(y_mean[i]))
        z_mean[i] = np.mean(sensor_i_data[:, 4])
        z_std[i] = np.std(sensor_i_data[:, 4])
        z_perc[i] = np.divide(z_std[i], abs(z_mean[i]))

    mean_data = np.zeros((number_sensors, 3))
    mean_data[:, 0] = x_mean
    mean_data[:, 1] = y_mean
    mean_data[:, 2] = z_mean

    std_data = np.zeros((number_sensors, 3))
    std_data[:, 0] = x_std
    std_data[:, 1] = y_std
    std_data[:, 2] = z_std

    perc_data = np.zeros((number_sensors, 3))
    perc_data[:, 0] = x_perc
    perc_data[:, 1] = y_perc
    perc_data[:, 2] = z_perc

    if on_stage:
        if verbose:
            print("Aquired Mean Data Set: Done!")
        return mean_data, std_data, perc_data, directory
    else:
        if verbose:
            print("Saving...")
        cwd = os.getcwd()
        directory = cwd + '\\' + str(N) + "_means\\"
        save_in_dir(mean_data, directory, "autosave")
        return mean_data, std_data, perc_data


def grid_2D_cube(CC_X: ConexCC, CC_Y: ConexCC, cube, specific_sensor, height, xlim=None, ylim=None,
           grid_number=50, sampling_size=10, verbose=False, save_data=False, directory=None):
    """
    Sweep over xy-plane and measure the field on each grid point, return the field values and coordinates.
    
    Note: 
    - The coordinates of sensor and actuators are not the same, this function accounts for this fact
    and returns all data in SENSOR COORDINATE SYSTEM.
    - It can happen that an invalid position that is out of the actuator's bounds should be reached.
    In this case, an according message is printed to the terminal and the actuator does not move. 

    Args:
    - CC_X, CC_Y (conexcc_class instances): represent x,y actuators
    - cube (serial.Serial): represents the magnetic field sensor.
    - specific_sensor (int in [1,64]): ID of a specific Hall sensor of the whole cube.
    - height (float): z-component of the 2d-plane, this value will be carried over to all positions
    - xlim, ylim (None or tuple/list of floats): Either None or (min, max) as limits for x and y direction.
    If no limits are provided, the min and max values of the respective actuator is used.
    - grid_number (int): number of points +1 to sweep over per axis per iteration
    - sampling_size (int): Number of samples considered for the mean value of the field estimated 
    with the chosen sensor at a fixed position
    - verbose (bool): switching on/off print-statements for displaying progress
    - save_data (bool): flag to switch on/off saving of data
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
        print('initial position:({:.4f}, {:.4f})'.format(CC_X.read_cur_pos(),CC_Y.read_cur_pos()))
        print('\nstep size (x): {:.5f}\nstep size (y): {:.5f}'.format(x_step, y_step))
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
                if CC_X.read_cur_pos() + x_step < CC_X.min_limit:
                    CC_X.move_absolute(CC_X.min_limit)
                elif CC_X.read_cur_pos() + x_step > CC_X.max_limit:
                    CC_X.move_absolute(CC_X.max_limit)

        for j in range(grid_number + 1):    # along y    
            if j != 0:
                success = CC_Y.move_relative((-1)**i * y_step)
                # if the actuator would have to move out of bounds, use absolute movement
                if not success:
                    if CC_Y.read_cur_pos() + (-1)**i * y_step < CC_Y.min_limit:
                        CC_Y.move_absolute(CC_Y.min_limit)
                    elif CC_Y.read_cur_pos() + (-1)**i * y_step > CC_Y.max_limit:
                        CC_Y.move_absolute(CC_Y.max_limit)
                
            all_ready(CC_X, CC2=CC_Y, verbose=verbose)
            check_no_motion(CC_X, CC2=CC_Y, verbose=verbose)

            # save position and measured field value 
            positions[i,j,0:2] = get_coords(CC_X, CC2=CC_Y)[0]
            B_field[i,j,:] = av_single_sens(cube, specific_sensor, sampling_size)

            if verbose:
                if i % 2 ==0: 
                    print('Position ({},{}): ({:.4f}, {:.4f})'.format(i,j, 
                                            (xmin + x_step * i), (ymin + y_step * j)))
                else:
                    print('Position ({},{}): ({:.4f}, {:.4f})'.format(i,j, 
                                            (xmin + x_step * i), (ymax - y_step * j)))

    # coordinate systems of sensor and axes differ, account for this and reshuffle the position array
    # x and y must be changed and get an additional (-1)-factor each
    positions_corrected = transform_between_sensor_stage_coordinates(positions)

    # eventually reshape the arrays
    positions_corrected = positions_corrected.reshape((grid_number+1)**2,3)
    B_field = B_field.reshape((grid_number+1)**2,3)

    # add the z-component as a position
    positions_corrected[:,2] = height

    # save data if desired
    if save_data:
        df = pd.DataFrame({ 'x [mm]': positions_corrected[:,0], 
                            'y [mm]': positions_corrected[:,1], 
                            'z [mm]': positions_corrected[:,2], 
                            'mean Bx [mT]': B_field[:,0],
                            'mean By [mT]': B_field[:,1],
                            'mean Bz [mT]': B_field[:,2]})

        if directory is None:
            directory = os.getcwd()
        ensure_dir_exists(directory)

        now = datetime.now().strftime("%y_%m_%d_%H-%M-%S") 
        output_file_name = "{}_2d_scan.csv".format(now)
        data_filepath = os.path.join(directory, output_file_name)

        try:
            df.to_csv(data_filepath, index=False, header=True)
        except FileNotFoundError:
            data_filepath = os.path.join(os.getcwd(), output_file_name)
            df.to_csv(data_filepath, index=False, header=True)
    
    return positions_corrected, B_field, data_filepath



# %%
if __name__ == "__main__":
    
    # %%
    # set up actuators
    reset = np.array([3.0, 3.0, 0])  # np.array([0,0,0])
    COM_ports = ['COM7', 'COM6', 'COM5']
    CC_X, CC_Y, CC_Z = setup(reset, COM_ports=COM_ports)

    # %%
    # set up the sensor
    port = 'COM4'
    sensor_id = 54
    N = 10

    with serial.Serial(port, 256000, timeout=2) as cube:
        # for sensor_id in range(1,65):
        #     measured_field = av_single_sens(cube, sensor_id, N)
        #     print('id{}: {}'.format(sensor_id, measured_field))

        # xpos, ypos, something = search(CC_X, CC_Y, cube, specific_sensor = 63, sampling_size = N)

        # angle_calib(0, cube, visual_feedback=True)

        mean_data, std_data, perc_data, _ = get_new_mean_data_set(
            N, cube=cube, on_stage=True, no_enter=True)

    # print('(xpos, ypos) = ({:.2f}, {:.2f})'.format(xpos, ypos))
    # print('xmax/num_prec = {:.4f}'.format(something))
    print(np.shape(mean_data))

    [print('sensor {}: {} +- {} ({})'.format(i+1, np.around(mean_data[i, :], 2),
                                             np.around(std_data[i, :], 2),
                                             np.around(perc_data[i, :], 2))) for i in range(64)]

#%%
