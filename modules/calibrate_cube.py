# Author Jona Buehler 2020
# Documentation and Updates by Nicholas Meinhardt

#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import serial
from numpy.linalg import norm
import os
import sys
# from math import copysign
# from time import sleep
# from scipy.optimize import root

from modules.serial_reader import get_new_data_set
from modules.conexcc_control import all_ready, save_in_dir, setup #, get_coords, correct_reset 
from conexcc.conexcc_class import *
from modules.plot_hall_cube import plot_set, plot_angle, plot_angle_spherical


#%%
def check_validity(min_val, max_val, CC1:ConexCC):
    """
    Compare min_val and max_val with hard limits of the actuator. If values are within these bounds, 
    return them. Else replace the invalid boundary with the corresponding limit set by the actuator.

    Args: 
    - min_val < max_val (float)
    - CC1, CC2 are instances of conexcc_class to control an actuator

    Returns: valid_min_val, valid_max_val
    """
    #  all correct
    if min_val >= CC1.min_limit and max_val <= CC1.max_limit:
        return min_val, max_val
    # min_val is not valid
    elif min_val < CC1.min_limit and max_val <= CC1.max_limit:
        return CC1.min_limit, max_val
    # max_val is not valid
    elif min_val >= CC1.min_limit and max_val > CC1.max_limit:
        return min_val, CC1.max_limit
    # neither boundary value is valid
    elif min_val < CC1.min_limit and max_val > CC1.max_limit:
        return CC1.min_limit, CC1.max_limit



def search(CC1:ConexCC, CC2:ConexCC, cube, specific_sensor, min_step_size=5*1e-3, xlim=None, ylim=None, 
                            num_prec = 10.0, sampling_size = 10, verbose=False):
    """
    Search the (x,y)-positon with minimum magnetic field along the xy-plane for constant hight z. 

    Note: It can happen that an invalid position that is out of the actuator's bounds should be reached.
    In this case, an according message is printed to the terminal and the actuator does not move. 
    This might yield an incorrect calibration! 

    Args:
    - CC1, CC2 are instances of conexcc_class for x and y axis
    - cube: instance of serial.Serial class, representing the magnetic field sensor.
    - specific_sensor (int in [1,64]): ID of a specific Hall sensor of the whole cube.
    - min_step_size: stop searching for minimum once the step size is smaller than this value.
    - xlim, ylim: Either None or (min, max) tuple or list as limits for x and y direction.
        If no limits are provided, the min and max values of the respective actuator is used.
    - num_prec (float): controls the stepsize when sweeping along x- and y-axis. The larger the number, 
      the smaller the stepsize. 
      NOTE: num_prec >= 2 is required to not raise errors, > 2 to not get trapped in while-loop!
    - sampling_size: Number of samples considered for the mean value of the field estimated 
      with the chosen sensor at a fixed position
    - verbose: switching on/off print-statements for displaying progress

    Returns: xpos, ypos, xmax/num_prec
    - xpos, ypos: 
    - xmax/num_prec:  why?
    """
    # set min and max limits for x,y axis
    if xlim==None:
        xmin = CC1.min_limit
        xmax=CC1.max_limit
    else:
        xmin=xlim[0]
        xmax=xlim[1]
    if ylim==None:
        ymin=CC2.min_limit
        ymax=CC2.max_limit
    else:
        ymin=ylim[0]
        ymax=ylim[1]

    # move actuator(or cube, respectively) to the minimum position and wait until controllers are ready again.
    CC1.move_absolute(xmin)
    CC2.move_absolute(ymin)
    all_ready(CC1, CC2=CC2)

    # Sweep over plane by considering each a min and max x-value (y-value), find position with minimum and 
    # redefine min/max values as this position +- one step size. As long as num_prec > 2 and 
    # if min/max remain between allowed limits, the procedudure should converge (?)
    xrang=xmax-xmin
    yrang=ymax-ymin
    while ((xrang/num_prec) >= min_step_size) or ((yrang/num_prec) >= min_step_size): 
        if verbose:
            print('\nstep size (x): {:.5f}\nstep size (y): {:.5f}'.format(xrang/num_prec, yrang/num_prec))

        # sweep along x-axis and measure B_y field component
        if verbose:
            print("\nNow sweeping x: ")
        fsx=[]
        for i in range(int(num_prec + 1)): 
            CC1.move_absolute(xrang/num_prec * i + xmin)
            all_ready(CC1, CC2=CC2)
            fsx.append(av_single_sens(cube, specific_sensor, sampling_size)[1])

            if verbose:
                print('Position {}: {:.5f}'.format(i, (xrang/num_prec * i + xmin)))
                print('Field at position {}: {:.3f}'.format(i, fsx[-1]))

        # move to x-position with smallest (absolute value of) B_y
        fminx=np.argmin(abs(np.asarray(fsx)))
        xpos=(xmin+(xrang/num_prec)*fminx)
        CC1.move_absolute(xpos)
        all_ready(CC1, CC2=CC2)
        if verbose:
            print('Now moving to position {}: {:.5f}'.format(fminx, xpos))

        # update xmin, xmax as the x-position with smallest B_y +- one stepsize
        xmin=xpos-(xrang/num_prec)
        xmax=xpos+(xrang/num_prec)
        xmin, xmax = check_validity(xmin, xmax, CC1)
        xrang=xmax-xmin 
        if verbose:
            print('New xmin: {:.5f} new xmax: {:.5f} new xrang: {:.5f}'.format(xmin, xmax, xrang))
        

        # sweep along y-axis and measure B_x field component
        if verbose:
            print("\nNow sweeping y: ")
        fsy=[]
        for i in range(int(num_prec + 1)):
            CC2.move_absolute(yrang/num_prec * i + ymin)
            all_ready(CC1, CC2=CC2)
            fsy.append(av_single_sens(cube, specific_sensor, sampling_size)[0])
            if verbose:
                print('Position {}: {:.5f}'.format(i, (yrang/num_prec * i + ymin)))
                print('Field at position {}: {:.3f}'.format(i, fsy[-1]))

        # move to x-position with smallest (absolute value of) B_x
        fminy = np.argmin(abs(np.asarray(fsy)))
        ypos = (ymin+(yrang/num_prec)*fminy)
        CC2.move_absolute(ypos)
        all_ready(CC1, CC2=CC2)
        if verbose:
            print('Now moving to position {}: {:.5f}'.format(fminy, ypos))

        # update ymin, ymax as the y-position with smallest B_x +- one stepsize
        ymin = ypos-(yrang/num_prec)
        ymax = ypos+(yrang/num_prec)
        ymin, ymax = check_validity(ymin, ymax, CC2)
        yrang = ymax-ymin
        if verbose:
            print('New ymin: {:.5f}  new ymax: {:.5f} new yrang: {:.5f}'.format(ymin, ymax, yrang))

    return xpos, ypos, xrang/num_prec

def av_single_sens(cube, specific_sensor, N, max_number_attempts = 10):
    """
    Measure magnetic field vector N times using a specific sensor and return mean values.

    It may happen that no valid results are returned from the sensor. In this case, still try to gather N
    valid samples by repeating the faulty measurement up to a maximum of max_number_attempts repetitions.

    Note: This function measures the magnetic field without specifying the position. It is recommended
    to ensure that the sensor is not moving beforehand.

    Args:
    - cube: instance of serial.Serial class, representing the magnetic field sensor
    - specific_sensor (int): number of a desired sensor, specific_sensor in [1, 64]  
    - N (int): number of B-field estimates that are averaged
    - max_number_attempts (int) is the maximum number of attempts to replace faulty measurements by valid outcomes.

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
            field = get_new_data_set(cube=cube, specific_sensor=specific_sensor)
            # if something does not work, get_new_data might return 1 instead of vector
            if not isinstance(field, int):
                field_samples.append(field)

    # return the mean field vector
    return np.mean(field_samples, axis=0)


def find_center_axis(CC1, CC2, cube, N = 10, min_step_size=5*1e-3, specific_sensor=54, limits_x = [0,10], 
                        limits_y = [0,10], verbose=True):
    """
    Find the xy-position of minimum in-plane magnetic field and estimate the field vector at this position.

    Args:
    - CC1, CC2 are instances of conexcc_class for x and y axis
    - cube: instance of serial.Serial class, representing the magnetic field sensor.
    - N: number of estimates of B-field used for average
    - min_step_size: stop searching for minimum of B-field along xy-plane 
    once the step size is smaller than this value.
    - specific_sensor: number of sensor that is used when searching the center position. 
    - limits_x, limits_x: array of length 2 with minimum and maximum values of x (y) that are considered 
    when searching for the center position.
    - verbose: switching on/off print-statements for displaying progress

    Return:
    - x0 = [xpos, ypos]: ndarray containing the position of the minimum in-plane field
    - f0: measured B-field vector at x0
    """
    if verbose:
        print("\nTrying to find center axis using sensor # {} ...\n".format(specific_sensor))

    xpos, ypos, precision = search(CC1, CC2, cube, specific_sensor, min_step_size=min_step_size, 
                                        xlim = limits_x, ylim = limits_y,
                                        verbose=verbose)
    x0 = np.array([xpos, ypos])
    f0 = av_single_sens(cube, specific_sensor, N)

    return x0, f0

def angle_wrt_z(vec):
    """
    Return angle (radian) of vector with respect to z axis.
    """
    mag = norm(vec)
    return np.arccos(vec[2]/mag)

def angle_calib(desired, cube, specific_sensor=54, N=10, visual_feedback=True, eps=0.5, 
                    max_number_trials = 100, spherical = True, verbose=True):
    """
    Compare measured to desired angle in a while-loop. 
    Leave the loop when the difference between both angles is less than eps degrees 
    or when the maximum number of trials is reached

    Note: Pass all angles in degrees!

    Args: 
    - desired: desired angle of magnetic field with respect to z-axis
    - cube: cube: instance of serial.Serial class, representing the magnetic field sensor.
    - specific_sensor (int): id in [1,64] of sensor that should be used for calibration
    - N: sampling size of magnetic field estimates used for average
    - visual_feedback (bool): flag to switch on/off a plotted visualization of measured vector,
    showing an normalized vector and the measured angle with z-axis
    - eps: acceptable difference in angle
    - max_number_trials: to avoid an infinite while-loop, a maximum number of trials can be set.
    - spherical: flag to switch on spherical plot, else 'normal' 3d plot is shown.
    - verbose: switching on/off print-statements for displaying progress
    Return 0 
    """
    diff=5*eps
    i = 0
    while abs(diff)>eps and i < max_number_trials:
        vec = av_single_sens(cube, specific_sensor, N)
        ang = angle_wrt_z(vec)
        diff=desired-np.degrees(ang)
        if visual_feedback and spherical:
            plot_angle_spherical(vec)
        elif visual_feedback:
            plot_angle(vec) 
        i += 1
        if verbose:
            print("\r Angle to z axis: {:.2f} °; Difference to desired angle ({:.2f} °) is {:.2f} °".format(np.degrees(ang), desired, diff))
    if verbose:
        print("Calibration of angle successfull!")
    return 0

def get_new_mean_data_set(N, filename=None, cube=None, no_enter=False, on_stage=False, omit_64 = False, 
                            verbose=False):
    """
    Estimate field vectors N-times with all sensors and calculate mean, std and abs(std/mean) as vectors for each sensor.
    
    Note: Sensor #64 has produced incorrect results in past, thus it can be omitted using the omit_True flag

    Args: 
    - N: number of times all 64 (or 63) sensor are read out in series
    - cube: instance of serial.Serial class, representing the magnetic field sensor.
    - filename: name of folder where data files are stored
    - no_enter (bool): if True, measurement starts automatically, else the user is asked to press enter to start.
    - on_stage (bool): flag used to set the action upon occurence an error when reading a measurement outcome 
	  from the sensor. If False, continue measuring and write a "'Read Failure', 0,0,0,0"-line to file. 
	  If True, the measurement is stopped entirely and the output file is deleted. 
    - omit_64: if True, sensor #64 will be omitted, else all 64 sensors are considered.
    - verbose: switching on/off print-statements for displaying progress

    Return:
    - if on_stage=True: mean_data, std_data, perc_data, directory 
    (where mean, std and abs(std/mean) are returned as ndarrays of shape (number_sensors, 3) for the 3 field directions)
    - else: mean_data, std_data, perc_data
    """
    # measure field using all sensors N times and save results to csv file
    if on_stage:
        resp=1
        path=''
        while resp==1:
            resp, directory, csvfile = get_new_data_set(measure_runs = N, filename=filename, cube=cube, verbose=verbose,
                                                    no_enter=no_enter, on_stage=on_stage, omit_64=omit_64)
            path = os.path.join(directory,csvfile)
    else:
        path = get_new_data_set(measure_runs = N, filename=filename, cube=cube, no_enter=no_enter, 
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
            if (data[i+k*number_sensors,:].dtype == 'float64'):
                sensor_i_data[k,:] = data[i+k*number_sensors,:]
            else:
                # print this message in any case!
                print("could not convert data properly! wrong data type: ", data[i+k*number_sensors,:].dtype)
                sensor_i_data[k,:] = 0
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
    mean_data[:,0] = x_mean
    mean_data[:,1] = y_mean
    mean_data[:,2] = z_mean

    std_data = np.zeros((number_sensors, 3))
    std_data[:,0] = x_std
    std_data[:,1] = y_std
    std_data[:,2] = z_std

    perc_data = np.zeros((number_sensors, 3))
    perc_data[:,0] = x_perc
    perc_data[:,1] = y_perc
    perc_data[:,2] = z_perc

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

#%%
if __name__ == "__main__":
    #N statistical values
    N = int(100)
    mean_data, std_data, perc_data = get_new_mean_data_set(N, no_enter=True)

    x_av_perc = np.mean(perc_data[:,0])
    y_av_perc = np.mean(perc_data[:,1])
    z_av_perc = np.mean(perc_data[:,2])

    #print(x_mean, "+-", x_std)
    #print(y_mean, "+-", y_std)
    #print(z_mean, "+-", z_std)
    print()
    print("Average percentual offset in")
    print("x direction: ", x_av_perc, " %")
    print("y direction: ", y_av_perc, " %")
    print("z direction: ", z_av_perc, " %")
    print(np.max(norm(mean_data, axis=1)))

    #plot_set(mean_data, Pos=np.array([0,0,0]), Vecs=False, Cont=True, single_comp='x')
    #plot_set(mean_data, Pos=np.array([0,0,0]), Vecs=False, Cont=True, single_comp='y')
    #plot_set(mean_data, Pos=np.array([0,0,0]), Vecs=False, Mag_label=True, Cont=True, single_comp='xy')
    #plot_set(mean_data, Pos=np.array([0,0,0]), Vecs=False, Cont=True, single_comp='z')
    plot_set(mean_data, Pos=np.array([0,0,0]))



    # ------------------------Testing area---------------------------------------------------
    #%%
    # set up actuators
    reset= np.array([3.0, 3.0, 0])#np.array([0,0,0])
    COM_ports = ['COM7', 'COM6', 'COM5']
    CC_X, CC_Y, CC_Z = setup(reset, COM_ports = COM_ports)

    #%%
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

        mean_data, std_data, perc_data, _ = get_new_mean_data_set(N, cube= cube, on_stage=True, no_enter=True)

    # print('(xpos, ypos) = ({:.2f}, {:.2f})'.format(xpos, ypos))
    # print('xmax/num_prec = {:.4f}'.format(something))
    print(np.shape(mean_data))

    [print('sensor {}: {} +- {} ({})'.format(i+1, np.around(mean_data[i,:], 2), 
                                                    np.around(std_data[i,:], 2), 
                                                    np.around(perc_data[i,:], 2))) for i in range(64)]




# %%
