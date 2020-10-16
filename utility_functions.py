"""
filename: utility_functions.py

This collection of functions has functions that can be used to manipulate the currents on the ECB channels. 
The user can choose from various methods to set currents on the different channels ('coils' in Pantec's terminology)
and thus generate a magnetic field with the vector magnet. For example, we can sweep through current values on the 
different channels and simultaneously measure the actual magnetic field to compare the theory with the actual results.

Author: Maxwell Guerne-Kieferndorf (QZabre)
        gmaxwell@student.ethz.ch

Date: 12.10.2020
"""
########## Standard library imports ##########
import numpy as np
import math
from time import time, sleep

########## local imports ##########
import transformations as tr
from main_comm import *
from measurements import *
from modules.analysis_tools import generate_plots

##########  Current parameters ##########

desCurrents = [0, 0, 0, 0, 0, 0, 0, 0]  # in milliamps
currDirectParam = b'1'

##########  Vector magnet properties ##########

windings = 508  # windings per coil
resistance = 0.47  # resistance per coil


def sweepCurrents(config_list=None, config='z', start_val=0, end_val=1, steps=5):
    """
    sweeps all currents in the (1,1,1) or (1,0,-1) configuration, meaning we have either a z field or an x-y-plane field, measures magnetic field 
    and stores the measured values in various arrays

    Args:
        config_list (numpy array of size 3, optional): current configuration entered by user, Defaults to np.array([0,0,1])
        config (str, optional): Choose from multiple possible 'configurations' (coil1, coil2, coil3) of currents. Possible values:
            - 'r': randomly generated configuration.
            - 'z': (1,1,1)
            - 'xy0': (1,0,-1)
        Defaults to 'z'.
        start_val (int, optional): [description]. Defaults to 0.
        end_val (int, optional): [description]. Defaults to 1.
        steps (int, optional): [description]. Defaults to 5.
    """
    global currDirectParam
    global desCurrents

    # initialization of all arrays
    all_curr_steps = np.linspace(start_val, end_val, steps)
    mean_values = np.zeros((steps, 3))
    stdd_values = np.zeros((steps, 3))
    expected_fields = np.zeros((steps, 3))
    all_curr_vals = np.zeros((steps, 3))
    directory = ''
    # some possible configurations
    current_direction = np.ndarray(3)
    if config_list is not None:
        max_val = abs(np.amax(config_list))
        current_direction[0] = config_list[0] / max_val
        current_direction[1] = config_list[1] / max_val
        current_direction[2] = config_list[2] / max_val
    elif config == 'z':
        current_direction[0] = 1
        current_direction[1] = 1
        current_direction[2] = 1
    elif config == 'xy0':
        current_direction[0] = 1
        current_direction[1] = 0
        current_direction[2] = -1
    # r for randomized
    elif config == 'r':
        arg = np.random.uniform(low=-1.0, high=1.0, size=3)
        max_val = abs(np.amax(arg))
        current_direction[0] = arg[0] / max_val
        current_direction[1] = arg[1] / max_val
        current_direction[2] = arg[2] / max_val
    else:
        print('invalid input!')
        return

    # create subdirectory to save measurements
    subDirBase = '({}_{}_{})_field_meas'.format(str(int(10*current_direction[0])), str(
        int(10*current_direction[1])), str(int(10*current_direction[0])))
    folder, filePath = newMeasurementFolder(sub_dir_base=subDirBase)

    enableCurrents()
    # iterate through all possible steps
    for i in range(steps):
        # set the current on each channel
        for k in range(3):
            desCurrents[k] = current_direction[k]*all_curr_steps[i]
        all_curr_vals[i] = current_direction*all_curr_steps[i]

        # tentative estimation of resulting B field
        B_expected = tr.computeMagField(
            current_direction*all_curr_steps[i], windings)

        setCurrents(desCurrents, currDirectParam)
        # Let the field stabilize
        sleep(0.1)
        # collect measured and expected magnetic field (of the specified sensor in measurements)
        print('measurement nr. ', i+1)
        # see measurements.py for more details
        mean_data, std_data = measure(filePath, N=15)
        mean_values[i] = mean_data
        stdd_values[i] = std_data
        expected_fields[i] = B_expected

    # end of measurements
    disableCurrents()
    # saving data section (prepared for plotting)
    saveDataPoints((all_curr_vals / 1000), mean_values,
                   stdd_values, expected_fields, filePath)


def rampVectorField(theta, phi, start_mag, finish_mag, steps):
    """
    Ramps magnetic field from start_magn to finish_magn in a specified number of steps and over a specified duration (sort of analogous to sweepCurrent)
    Measure the magnetic field values, save them to a file.


    Args:
        theta (float): Give the direction of the magnetic field. Angle between field vector and z axis
        phi (float): angle between field vector projection on xy plane and x axis
        start_mag (float): Magnitude in mT
        finish_mag (float): Magnitude in mT, finish_mag > start_mag must hold!!!
        steps (int): number of measurement steps
    """
    global currDirectParam
    global desCurrents

    magnitudes = np.linspace(start_mag, finish_mag, steps)
    mean_values = np.zeros((steps, 3))
    stdd_values = np.zeros((steps, 3))
    expected_fields = np.zeros((steps, 3))
    all_curr_vals = np.zeros((steps, 3))
    directory = ''

    # create subdirectory to save measurements
    subDirBase = '{}_{}_field_meas'.format(theta, phi)
    folder, filePath = newMeasurementFolder(sub_dir_base=subDirBase)

    enableCurrents()
    # iterate through all possible steps
    for i in range(steps):
        # compute the currents theoretically needed to achieve an arbitrary magnetic field
        B_expected = tr.computeMagneticFieldVector(magnitudes[i], theta, phi)
        I_vector = tr.computeCoilCurrents(B_expected, windings, resistance)

        # set the computed currents on each channel
        for k in range(3):
            desCurrents[k] = I_vector[k]
        setCurrents(desCurrents, currDirectParam)
        # Let the field stabilize
        sleep(0.1)
        # collect measured and expected magnetic field (of the specified sensor in measurements)
        print('measurement nr. ', i+1)
        # collect measured magnetic field (of the specified sensor in measurements)
        mean_data, std_data = measure(filePath, N=15)
        mean_values[i] = mean_data
        stdd_values[i] = std_data
        expected_fields[i] = B_expected
        all_curr_vals[i] = I_vector

    # end of measurements
    disableCurrents()
    # plotting section
    saveDataPoints((all_curr_vals / 1000), mean_values,
                   stdd_values, expected_fields, directory)


def runCurrents(channels, t=0, direct=b'1'):
    """
    run arbitrary currents (less than maximum current) on each channel 


    Args:
        channels (np.array(list(int)), length 3): current values in [mA]. Make sure not to give more than 3!
        t (int, optional): time for which the magnetic field should be activated. 
            If zero, user can decide whether to change the current or deactivate it. Defaults to 0.
        direct (bytes, optional): current direct parameter (can usually be left alone). Defaults to b'1'.
    """
    global currDirectParam
    global desCurrents

    currDirectParam = direct
    # copy the computed current values (mA) into the desCurrents list (first 3 positions)
    # cast to int
    for i in range(len(channels)):
        if np.abs(channels[i]) > ECB_MAX_CURR:
            print("desired current exceeds limit")
            return
        desCurrents[i] = int(channels[i])

    # user specified time
    if t > 0:
        enableCurrents()
        setCurrents(desCurrents, currDirectParam)
        # prevent the connection with the ECB from timing out for long measurements.
        if t < 500:
            sleep(t)
        else:
            starttime = time()
            while time() - starttime < t:
                sleep(500)
                getCurrents()

        disableCurrents()

    # on until interrupted by user
    elif t == 0:
        enableCurrents()
        setCurrents(desCurrents, currDirectParam)
        # wait until user presses enter
        c1 = '0'
        while c1 != 'q':
            c1 = input(
                '[q] to disable currents\n[c]: get currents\n[r]: Set new currents\n[s]: get ECB status\n[t]: get coil temperature (useless)\n')
            if c1 == 'c':
                getCurrents()
            elif c1 == 'r':
                channels[0] = input('Channel 1 current: ')
                channels[1] = input('Channel 2 current: ')
                channels[2] = input('Channel 3 current: ')
                # handle new inputs
                for i in range(len(channels)):
                    if np.abs(channels[i]) > ECB_MAX_CURR:
                        print("desired current exceeds limit")
                        return
                    try:
                        desCurrents[i] = int(channels[i])
                    except:
                        print(
                            "non-integer value entered, setting channel {} to 0".format(i+1))
                        desCurrents[i] = 0

                setCurrents(desCurrents, currDirectParam)

            elif c1 == 's':
                print(getStatus())
            elif c1 == 't':
                getTemps()

        disableCurrents()
    else:
        return


def generateMagneticField(magnitude, theta, phi, t=0, direct=b'1'):
    """
    generate a magnetic field in an arbitrary direction and an arbitrary magnitude


    Args:
        magnitude (float): of the B-field, in [mT]
        theta (float): angle between desired field direction and z axis
        phi (float): azimuthal angle (measured from the x axis)
        t (int, optional): Time for which the magnetic field should be activated (if not 0). Defaults to 0.
        direct (bytes, optional): Current direct parameter. Defaults to b'1'.
    """
    global currDirectParam
    global desCurrents

    B_vector = tr.computeMagneticFieldVector(magnitude, theta, phi)
    I_vector = tr.computeCoilCurrents(B_vector, windings, resistance)

    currDirectParam = direct
    # copy the computed current values (mA) into the desCurrents list (first 3 positions)
    # cast to int
    for i in range(len(I_vector)):
        # make sure that the current is not too high
        if np.amax(I_vector) > ECB_MAX_CURR:
            print("desired current exceeds limit")
            return
        desCurrents[i] = int(I_vector[i])

    # user specified on time
    if t > 0:
        enableCurrents()
        setCurrents(desCurrents, currDirectParam)
        # prevent the connection with the ECB from timing out for long measurements.
        if t < 500:
            sleep(t)
        else:
            starttime = time()
            while time() - starttime < t:
                sleep(500)
                getCurrents()

        disableCurrents()
    # on until interrupted by user
    elif t == 0:
        enableCurrents()

        setCurrents(desCurrents, currDirectParam)

        # wait until user presses enter
        c1 = '0'
        while c1 != 'q':
            c1 = input(
                '[q] to disable currents\n[c]: get currents\n[s]: get ECB status\n[t]: get coil temperature\n')
            if c1 == 'c':
                getCurrents()
            elif c1 == 'r':
                inp1 = input('New B-Field magnitude: ')
                inp2 = input('Polar angle theta: ')
                inp3 = input('Azimuthal angle phi: ')
                try:
                    magnitude = int(inp1)
                except:
                    print("non-integer value entered, magnitude = 20")
                    magnitude = 20
                try:
                    theta = int(inp2)
                except:
                    print("non-integer value entered, theta = 0")
                    theta = 0
                try:
                    phi = int(inp3)
                except:
                    print("non-integer value entered, phi = 0")
                    phi = 0

                B_vector = tr.computeMagneticFieldVector(magnitude, theta, phi)
                I_vector = tr.computeCoilCurrents(
                    B_vector, windings, resistance)

                for i in range(len(I_vector)):
                    if np.abs(I_vector[i]) > ECB_MAX_CURR:
                        print("desired current exceeds limit")
                        return

                    desCurrents[i] = int(I_vector[i])

                setCurrents(desCurrents, currDirectParam)

            elif c1 == 's':
                print(getStatus())
            elif c1 == 't':
                getTemps()

        disableCurrents()
    else:
        return


def switchConfigsAndMeasure(config1, config2, dt=0, rounds=10):
    """
    Switch quickly between two current configurations and keep track of the measured fields over time. The time in each state is dt.

    Args:
        config1 (np.array): currents on coil 1, 2 and 3 as an array with 3 entries.
        config2 (np.array): currents on coil 1, 2 and 3 as an array with 3 entries.
        dt (float): Time to remain in each state. Defaults to 0s.
        rounds (int): Number of times to switch. Defaults to 10
    """
    global currDirectParam
    global desCurrents

    subDirBase = 'dynamic_field_meas'
    folder, filePath = newMeasurementFolder(sub_dir_base=subDirBase)

    enableCurrents()

    while rounds > 0:
        for i in range(len(config1)):
            # make sure that the current is not too high
            if np.amax(config1) > ECB_MAX_CURR:
                print("desired current exceeds limit")
                return
            desCurrents[i] = int(config1[i])

        setCurrents(desCurrents, currDirectParam)
        #measure(folder, True, N=10)
        sleep(dt)

        for i in range(len(config2)):
            # make sure that the current is not too high
            if np.amax(config2) > ECB_MAX_CURR:
                print("desired current exceeds limit")
                return
            desCurrents[i] = int(config2[i])

        setCurrents(desCurrents, currDirectParam)
        #measure(folder, True, N=10)
        sleep(dt)

        rounds = rounds-1

    disableCurrents()
