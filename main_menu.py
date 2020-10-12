"""
filename: main_menu.py

This script is meant to be used as an interface with the ECB 820. The user can choose from various
methods to set currents on the different channels ('coils' in Pantec's terminology) and thus communicate
with the ECB. The standard interface is the command line, but another option is to integrate this into 
a GUI for the best user experience.

Author: Maxwell Guerne-Kieferndorf (QZabre)
        gmaxwell@student.ethz.ch

Date: 09.10.2020
"""
########## Standard library imports ##########
import numpy as np
import math
from time import time, sleep
# from datetime import datetime
# import sys
# from os import getcwd, path
# from pathlib import Path
# import csv

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


def MainMenu(initialized):
    """
    Main menu for ECB/Vector magnet operation an arbitrary magnitude   

    Args:
    - initialized: if the ECB is initialized, this will be 0
    """

    # is there a connection?
    if initialized == 0:
        c1 = '0'
        while c1 != 'x':
            print('----------- Main Menu -----------')
            print('[x] to exit\n')
            print('[1]: sweep multiple current values and make measurment with cube (specify coil configuration)')
            print('[2]: sweep theoretical magnetic field magnitudes, measure actual values (specify polar and azimuthal angles, magnitude range)')
            print('[3]: set currents manually on 3 channels (in mA)')
            print(
                '[4]: generate magnetic field (specify polar and azimuthal angles, magnitude)\n')
            print('[s]: get ECB status\n[r] roll a die\n')

            c1 = input()
            if c1 == '1':
                inp1 = input(
                    'configuration (z or x-y direction), acceptable inputs:\n z, xy0...6 (for different directions in xy plane), r for random test = ')
                inp2 = input('starting current in mA = ')
                inp3 = input('ending current in mA = ')
                inp4 = input('# of steps: ')
                inp5 = input('How many measurement runs? (only if config r was chosen): ')

                try:
                    config = inp1
                except:
                    print('expected valid input (z or xy), defaulting to z')
                    config = 'z'
                try:
                    start_val = int(inp2)
                except:
                    print('expected numerical value, defaulting to 0')
                    start_val = 0
                try:
                    end_val = int(inp3)
                except:
                    print('expected numerical value, defaulting to 1')
                    end_val = 1
                try:
                    steps = int(inp4)
                except:
                    print('expected numerical value, defaulting to 1')
                    steps = 1
                try:
                    randomRuns = int(inp5)
                except:
                    print('expected numerical value, defaulting to 0')
                    randomRuns = 0
                    
                if randomRuns < 1:
                    sweepCurrents(config, start_val, end_val, steps)
                elif (randomRuns >= 1 and config == 'r'):
                    while randomRuns > 0:
                        sweepCurrents(config, start_val, end_val, steps)
                        randomRuns = randomRuns-1
                else:
                    print('please enter a valid combination of inputs.')
                
            elif c1 == '2':
                inp1 = input('Angle to z axis in deg = ')
                inp2 = input('Angle to x axis in deg = ')
                inp3 = input('starting magnitude in mT = ')
                inp4 = input('ending magnitude in mT = ')
                inp5 = input('# of steps: ')

                try:
                    theta = int(inp1)
                except:
                    print('expected numerical value, defaulting to 0')
                    theta = 0
                try:
                    phi = int(inp2)
                except:
                    print('expected numerical value, defaulting to 0')
                    phi = 0
                try:
                    start_mag = int(inp3)
                except:
                    print('expected numerical value, defaulting to 0')
                    start_mag = 0
                try:
                    end_mag = int(inp4)
                except:
                    print('expected numerical value, defaulting to 1')
                    end_mag = 1
                try:
                    steps = int(inp5)
                except:
                    print('expected numerical value, defaulting to 1')
                    steps = 1

                rampVectorField(theta, phi, start_mag, end_mag, steps)

            elif c1 == '3':
                inp1 = input('Channel 1: ')
                inp2 = input('Channel 2: ')
                inp3 = input('Channel 3: ')
                inp4 = input('timer (leave empty -> manual termination) = ')
                try:
                    coil1 = int(inp1)
                except:
                    print('expected numerical value, defaulting to 0')
                    coil1 = 0
                try:
                    coil2 = int(inp2)
                except:
                    print('expected numerical value, defaulting to 0')
                    coil2 = 0
                try:
                    coil3 = int(inp3)
                except:
                    print('expected numerical value, defaulting to 0')
                    coil3 = 0

                if inp4 == '':
                    runCurrents(coil1, coil2, coil3, t=0, direct=b'1')
                else:
                    try:
                        timer = int(inp4)
                    except:
                        print('expected numerical value, defaulting to 0')
                        timer = 0
                        
                runCurrents(coil1, coil2, coil3, timer, direct=b'1')

            elif c1 == '4':
                inp1 = input('Magnitude in mT = ')
                inp2 = input('Angle to z axis in deg = ')
                inp3 = input('Angle to x axis in deg = ')
                inp4 = input('timer (leave empty -> manual termination) = ')
                try:
                    mag = int(inp1)
                except:
                    print('expected numerical value')
                    mag = 0
                try:
                    theta = int(inp2)
                except:
                    print('expected numerical value')
                    theta = 0
                try:
                    phi = int(inp3)
                except:
                    print('expected numerical value')
                    phi = 0
                if inp4 == '':
                    generateMagneticField(mag, theta, phi, t=0, direct=b'1')
                else:
                    try:
                        timer = int(inp4)
                    except:
                        print('expected numerical value')
                        timer = 0
                        
                generateMagneticField(mag, theta, phi, timer, direct=b'1')

            elif c1 == 's':
                print(getStatus())
            elif c1 == 'r':
                print(np.random.randint(1, 7))

    else:
        print('not connected')
        return


def sweepCurrents(config='z', start_val=0, end_val=1, steps=5):
    """
    sweeps all currents in the (1,1,1) or (1,0,-1) configuration, meaning we have either a z field or an x-y-plane field, measures magnetic field 
    and stores the measured values in various arrays

    Args:
    -config: 
    -start/end_val: current to start and end with, mA
    -steps: number of steps

    Args:
        ch1,ch2,ch3 (float, optional): current ratios (to the maximum number entered)
        config (str, optional): Choose from multiple possible 'configurations' (coil1, coil2, coil3) of currents. Possible values:
            - 'z': coil currents all the same
            - 'xy0', 'xy1',...: one coil on in positive direction, one in negative direction and one off
            - 'r': randomly generated configuration.
        Defaults to 'z'.
        start_val (int, optional): [description]. Defaults to 0.
        end_val (int, optional): [description]. Defaults to 1.
        steps (int, optional): [description]. Defaults to 5.
    """

    # initialization of all arrays
    all_curr_steps = np.linspace(start_val, end_val, steps)
    mean_values = np.zeros((steps, 3))
    stdd_values = np.zeros((steps, 3))
    expected_fields = np.zeros((steps, 3))
    all_curr_vals = np.zeros((steps, 3))
    directory = ''

    current_direction = np.ndarray(3)
    if config == 'z':
        current_direction[0] = 1
        current_direction[1] = 1
        current_direction[2] = 1
    elif config == 'xy0':
        current_direction[0] = 1
        current_direction[1] = 0
        current_direction[2] = -1
    elif config == 'xy1':
        current_direction[0] = 1
        current_direction[1] = -1
        current_direction[2] = 0
    elif config == 'xy2':
        current_direction[0] = -1
        current_direction[1] = 1
        current_direction[2] = 0
    elif config == 'xy3':
        current_direction[0] = -1
        current_direction[1] = 0
        current_direction[2] = 1
    elif config == 'xy4':
        current_direction[0] = 0
        current_direction[1] = -1
        current_direction[2] = 1
    elif config == 'xy5':
        current_direction[0] = 0
        current_direction[1] = 1
        current_direction[2] = -1
    elif config == 'r':
        arg = np.random.uniform(size=3)
        max_val = np.amax(arg)
        current_direction[0] = arg[0] / max_val
        current_direction[1] = arg[1] / max_val
        current_direction[2] = arg[2] / max_val
    else:
        print('invalid input!')
        return

    subDirBase = '({}{}{})_field_meas'.format(current_direction[0],current_direction[1],current_direction[2])
    folder = newMeasurementFolder(sub_dir_base=subDirBase)

    enableCurrents()
    # iterate through all possible steps
    for i in range(steps):
        # set the current on each channel
        for k in range(3):
            desCurrents[k] = current_direction[k] * all_curr_steps[i]
            
        all_curr_vals[i] = current_direction * all_curr_steps[i]
        
        B_expected = tr.computeMagField(
            current_direction * all_curr_steps[i], windings)
        setCurrents(desCurrents, currDirectParam)
        sleep(0.8)

        # collect measured and expected magnetic field (of the specified sensor in measurements)
        print('measurement nr. ', i+1)
        # see measurements.py for more details
        mean_data, std_data, directory = measure(sub_dirname=folder)
        mean_values[i] = mean_data
        stdd_values[i] = std_data
        expected_fields[i] = B_expected

        sleep(0.1)

    disableCurrents()

    # plotting section (generate plots still only works for 1 current value)
    saveDataPoints((all_curr_vals / 1000), mean_values,
                   stdd_values, expected_fields, directory)


def rampVectorField(theta, phi, start_mag, finish_mag, steps):
    """  Ramps magnetic field from start_magn to finish_magn in a specified number of steps and over a specified duration (sort of analogous to sweepCurrent)
    Measure the magnetic field values, save them to a file.


    Args:
        theta (float): Give the direction of the magnetic field. Angle between field vector and z axis
        phi (float): angle between field vector projection on xy plane and x axis
        start_mag (float): Magnitude in mT
        finish_mag (float): Magnitude in mT, finish_mag > start_mag must hold!!!
        steps (int): number of measurement steps
    """

    magnitudes = np.linspace(start_mag, finish_mag, steps)
    mean_values = np.zeros((steps, 3))
    stdd_values = np.zeros((steps, 3))
    expected_fields = np.zeros((steps, 3))
    all_curr_vals = np.zeros((steps, 3))
    directory = ''

    subDirBase = '{}_{}_field_meas'.format(theta,phi)
    folder = newMeasurementFolder(sub_dir_base=subDirBase)

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
        sleep(0.8)

        print('measurement nr. ', i+1)
        # collect measured magnetic field (of the specified sensor in measurements)
        mean_data, std_data, directory = measure(sub_dirname=folder)
        mean_values[i] = mean_data
        stdd_values[i] = std_data
        expected_fields[i] = B_expected
        all_curr_vals[i] = I_vector

        sleep(0.1)

    disableCurrents()
    # plotting section
    saveDataPoints((all_curr_vals / 1000), mean_values,
                   stdd_values, expected_fields, directory)


def runCurrents(channels, t=0, direct=b'1'):
    """
    run arbitrary currents (less than maximum current) on each channel 
    
    
    Args:
        coils (list(int), length 3): current values in [mA]. Make sure not to give more than 3!
        t (int, optional): time for which the magnetic field should be activated. 
            If zero, user can decide whether to change the current or deactivate it. Defaults to 0.
        direct (bytes, optional): current direct parameter (can usually be left alone). Defaults to b'1'.
    """
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

        sleep(t)

        disableCurrents()

    # on until interrupted by user
    elif t == 0:
        enableCurrents()
        # TODO: chk slew rate
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

                for i in range(len(channels)):
                    if np.abs(channels[i]) > ECB_MAX_CURR:
                        print("desired current exceeds limit")
                        return
                    # TODO: chk slew rate
                    desCurrents[i] = int(channels[i])
                    
            elif c1 == 's':
                print(getStatus())
            elif c1 == 't':
                getTemps()

        disableCurrents()
    else:
        return


#TODO: implement this function
def chkSlewRate():
    pass


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

    B_vector = tr.computeMagneticFieldVector(magnitude, theta, phi)
    I_vector = tr.computeCoilCurrents(B_vector, windings, resistance)

    # make sure that the current is not too high
    if np.amax(I_vector) > ECB_MAX_CURR:
        print("desired current exceeds limit")
        return

    currDirectParam = direct
    # copy the computed current values (mA) into the desCurrents list (first 3 positions)
    # cast to int
    for i in range(len(I_vector)):
        desCurrents[i] = int(I_vector[i])

    # user specified on time
    if t > 0:
        enableCurrents()
        setCurrents(desCurrents, currDirectParam)
        print('Current on coils 1, 2 and 3: [{}, {}, {}]'.format(
            I_vector[0], I_vector[1], I_vector[2]))
        sleep(t)

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
            elif c1 == 's':
                print(getStatus())
            elif c1 == 't':
                getTemps()

        disableCurrents()
    else:
        return


if __name__ == '__main__':
    ecbInit = openConnection()
    while ecbInit != 0:
        ecbInit = openConnection()

    MainMenu(ecbInit)
    closeConnection()
