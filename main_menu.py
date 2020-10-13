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
from utility_functions import *
from main_comm import *


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
            print(
                '[1]: sweep multiple current values and make measurment with cube (specify coil configuration)')
            print('[2]: sweep theoretical magnetic field magnitudes, measure actual values (specify polar and azimuthal angles, magnitude range)')
            print('[3]: set currents manually on 3 channels (in mA)')
            print(
                '[4]: generate magnetic field (specify polar and azimuthal angles, magnitude)')
            print(
                '[5]: Switch from one operating point to another multiple times while measuring\n')
            print('[s]: get ECB status\n[r] roll a die\n')

            c1 = input()
            if c1 == '1':
                inp1 = input(
                    'configuration (z or x-y direction), acceptable inputs:\n z, xy0...6 (for different directions in xy plane), r for random test = ')
                inp2 = input('starting current in mA = ')
                inp3 = input('ending current in mA = ')
                inp4 = input('# of steps: ')
                inp5 = input(
                    'How many measurement runs? (only if config r was chosen): ')

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
                    (config, start_val, end_val, steps)
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
                    runCurrents(
                        np.array([coil1, coil2, coil3]), t=0, direct=b'1')
                else:
                    try:
                        timer = int(inp4)
                    except:
                        print('expected numerical value, defaulting to 0')
                        timer = 0

                    runCurrents(
                        np.array([coil1, coil2, coil3]), timer, direct=b'1')

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
                
            elif c1 == '5':
                inp1 = input('configuration 1\nChannel 1: ')
                inp2 = input('Channel 2: ')
                inp3 = input('Channel 3: ')
                inp4 = input('configuration 2\nChannel 1: ')
                inp5 = input('Channel 2: ')
                inp6 = input('Channel 3: ')
                inp7 = input('time in each state (add 1.5s): ')
                inp8 = input('how many times to switch: ')
                try:
                    a1 = int(inp1)
                    b1 = int(inp2)
                    c1 = int(inp3)
                    config1 = np.array([a1,b1,c1])
                except:
                    print('expected numerical value, defaulting to (0,0,1)')
                    config1 = np.array([0,0,1])
                try:
                    a2 = int(inp4)
                    b2 = int(inp5)
                    c2 = int(inp6)
                    config2 = np.array([a2,b2,c2])
                except:
                    print('expected numerical value(s), defaulting to (0,1,0)')
                    config2 = np.array([0,1,0])
                try:
                    dt = float(inp7)
                except:
                    print('expected numerical value(s), defaulting to 0')
                    dt = 0
                try:
                    rounds = int(inp8)
                except:
                    print('expected numerical value(s), defaulting to 10')
                    rounds = 10

                switchConfigsAndMeasure(config1, config2, dt, rounds)

            elif c1 == 's':
                print(getStatus())
            elif c1 == 'r':
                print(np.random.randint(1, 7))

    else:
        print('not connected')
        return


if __name__ == '__main__':
    ecbInit = openConnection()
    while ecbInit != 0:
        ecbInit = openConnection()

    MainMenu(ecbInit)
    closeConnection()
