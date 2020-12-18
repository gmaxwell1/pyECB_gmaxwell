"""
filename: main_menu.py

This script is meant to be used as an interface with the ECB 820. The user can choose from various
methods to set currents on the different channels ('coils' in Pantec's terminology) and thus communicate
with the ECB. The standard interface for now is the command line, but the ultimate goal is to integrate this into
QS3.

Author: Maxwell Guerne-Kieferndorf (QZabre)
        gmaxwell@student.ethz.ch

Date: 09.10.2020
"""
########## Standard library imports ##########
import numpy as np
import math
from time import time, sleep
import csv
import os
from scipy import stats

########## local imports ##########
from utility_functions import *
from main_comm import *
from measurements import gotoPosition
import feedback as fb
from MetrolabTHM1176.thm1176 import MetrolabTHM1176Node


def MainMenu(initialized):
    """
    Main menu for ECB/Vector magnet operation. An arbitrary combination of currents can be set with this menu, thus
    any magnetic field may be generated as well.

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
                '[1]: sweep current values and record measurements with sensor (specify coil configuration)')
            print('[2]: sweep theoretical magnetic field vectors, measure actual components '
                  '\n\t(specify polar and azimuthal angles, magnitude range or rotational axis)')
            print('[3]: set currents manually on 3 channels (in mA)')
            print(
                '[4]: generate magnetic field (specify polar and azimuthal angles, magnitude)')

            print('[h] do a hysteresis test.\n')

            c1 = input()

            if c1 == '1':
                inp0 = input(
                    'Automatic or manual input? (m or f or nothing): ')
                datadir = input(
                    'Enter a valid directory name to save measurement data in: ')
                c1 = input('Automatic exit after finish? (x for yes): ')
                if datadir == '':
                    callCurrentSweep(inp0)
                else:
                    callCurrentSweep(inp0, datadir)

            # tentative implementation. Actual I-to-B actuation matrix needed. Many other features not added yet.
            elif c1 == '2':
                c1 = input('Automatic exit after finish? (x for yes): ')
                callSweepVectorField()

            elif c1 == '3':
                c1 = input('Automatic exit after finish? (x for yes): ')
                callRunCurrents()

            # tentative implementation. Actual I-to-B actuation matrix needed.
            elif c1 == '4':
                c1 = input('Automatic exit after finish? (x for yes): ')
                callGenerateVectorField()

            # elif c1 == '5':
            #     c1 = input('Automatic exit after finish? (x for yes): ')
            #     callFunctionGen()

            # elif c1 == '6':
            #     c1 = input('Automatic exit after finish? (x for yes): ')
            #     feedbackMode()

            elif c1 == 'h':
                c1 = input('Automatic exit after finish? (x for yes): ')
                callHysteresisSweep()

    else:
        print('not connected')
        return


def callCurrentSweep(mode='m', datadir='test_measurements'):
    """
    Setup function to call the utility function 'sweepCurrents', see 'utility_functions.py'. Manages necessary inputs.

    Args:
        mode (str, optional): Decides whether a file with configurations or a manual input will be read. The third option is
        using default configurations, such as the x, y or z direction. Defaults to 'm'.
        datadir (str, optional): subdirectory to save measurements in. default: 'test_measurements'
    """
    if mode == 'f':
        # must be a .csv file!
        inpFile = input('Enter a valid configuration file path: ')
        inp1 = input('starting current in mA = ')
        inp2 = input('ending current in mA = ')
        inp3 = input('# of steps: ')
        # the values for each measurement run should be the same for consistent results
        try:
            start_val = int(inp1)
        except:
            print('expected numerical value, defaulting to -4500')
            start_val = -4500
        try:
            end_val = int(inp2)
        except:
            print('expected numerical value, defaulting to 4500')
            end_val = 4500
        try:
            steps = int(inp3)
        except:
            print('expected numerical value, defaulting to 200')
            steps = 200

        node = MetrolabTHM1176Node(block_size=20, range='0.3 T', period=0.01, average=1)# as node:
        gotoPosition()
        inp = input('Do grid sweep function? (y/n) ')
        if inp == 'y':
            gridSweep(node, inpFile, datadir=datadir, current_val=start_val, demagnetize=True, today=False)
        else:
            with open(inpFile, 'r') as f:
                contents = csv.reader(f)
                next(contents)
                for row in contents:
                    config = np.array(
                        [float(row[0]), float(row[1]), float(row[2])])
                        
                    sweepCurrents(config_list=config, start_val=start_val, datadir=datadir,
                                    end_val=end_val, steps=steps, node=node, today=False)

    elif mode == 'm':
        inp1 = input('Configuration:\nChannel 1: ')
        inp2 = input('Channel 2: ')
        inp3 = input('Channel 3: ')
        inp4 = input('starting current in mA = ')
        inp5 = input('ending current in mA = ')
        inp6 = input('# of steps: ')

        try:
            a1 = float(inp1)
            b1 = float(inp2)
            c1 = float(inp3)
            config = np.array([a1, b1, c1])
        except:
            print(
                'expected numbers (float or int), defaulting to (0,0,1)')
            config = np.array([0, 0, 1])
        try:
            start_val = int(inp4)
        except:
            print('expected numerical value, defaulting to 0')
            start_val = 0
        try:
            end_val = int(inp5)
        except:
            print('expected numerical value, defaulting to 1')
            end_val = 1
        try:
            steps = int(inp6)
        except:
            print('expected numerical value, defaulting to 1')
            steps = 1

        with MetrolabTHM1176Node(block_size=20, range='0.3 T', period=0.01, average=5) as node:
            gotoPosition()

            sweepCurrents(config_list=config, start_val=start_val, datadir=datadir,
                          end_val=end_val, steps=steps, node=node, today=False)

    else:
        print('Using preset current configuration.')
        config = input('Which config (z, xy or r): ')
        inp1 = input('How many measurement runs?: ')
        inp2 = input('starting current in mA = ')
        inp3 = input('ending current in mA = ')
        inp4 = input('# of steps: ')

        try:
            randomRuns = int(inp1)
        except:
            print('expected numerical value, defaulting to 0')
            randomRuns = 0
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

        while randomRuns > 0:
            sweepCurrents(config=config, start_val=start_val,
                          end_val=end_val, steps=steps, datadir=datadir,
                          node=MetrolabTHM1176Node(block_size=20, range='0.3 T', period=0.01, average=5))
            randomRuns = randomRuns-1


def callHysteresisSweep():
    """
    Setup function to call the utility function 'sweepCurrents', see 'utility_functions.py'. Manages necessary inputs.

    Args:
        mode (str, optional): Decides whether a file with configurations or a manual input will be read. The third option is
        using default configurations, such as the x, y or z direction. Defaults to 'm'.
        datadir (str, optional): subdirectory to save measurements in. default: 'test_measurements'
    """
    inp1 = input('Configuration:\nChannel 1: ')
    inp2 = input('Channel 2: ')
    inp3 = input('Channel 3: ')
    datadir = input('Which directory should the output file be saved in? ')
    inp5 = input('maximum current in mA = ')
    inp6 = input('# of steps: ')

    try:
        a1 = float(inp1)
        b1 = float(inp2)
        c1 = float(inp3)
        config = np.array([a1, b1, c1])
    except:
        print(
            'expected numbers (float or int), defaulting to (0,0,1)')
        config = np.array([0, 0, 1])
    try:
        end_val = int(inp5)
    except:
        print('expected numerical value, defaulting to 1')
        end_val = 1
    try:
        steps = int(inp6)
    except:
        print('expected numerical value, defaulting to 1')
        steps = 1

    if datadir != '':
        sweepHysteresis(config_list=config, datadir=datadir,
                        end_val=end_val, steps=steps, today=False)
    else:
        sweepHysteresis(config_list=config, end_val=end_val,
                        steps=steps, today=False)


def callSweepVectorField():
    """
    Setup function to call the utility function 'rampVectorField', see 'utility_functions.py'. Manages necessary inputs.
    """
    inp1 = input('Angle to z axis in deg = ')
    inp2 = input('Angle to x axis in deg = ')
    inp3 = input('starting magnitude in mT = ')
    rot = input(
        'Rotate constant magnitude field around a specified axis? (otherwise ramp field magnitude in a constant direction): ')

    if rot == 'y':
        axisang1 = input('Axis polar angle: ')
        axisang2 = input('Axis azimuthal angle: ')
        try:
            axis = (int(axisang1), int(axisang2))
        except:
            print('expected numerical value, defaulting to (0,0) or z-axis')
            axis = (0, 0)

    else:
        inp4 = input('ending magnitude in mT = ')
        try:
            end_mag = int(inp4)
        except:
            print('expected numerical value, defaulting to 1')
            end_mag = 1

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
        steps = int(inp5)
    except:
        print('expected numerical value, defaulting to 1')
        steps = 1

    with MetrolabTHM1176Node(block_size=20, range='0.3 T', period=0.01, average=5) as node:
        # node = MetrolabTHM1176Node(block_size=20, sense_range_upper="0.3 T", period=0.001)
        gotoPosition()

        if rot == 'y':
            rampVectorField(node, theta, phi, start_mag,
                            steps=steps, rotate=axis)
        else:
            rampVectorField(node, theta, phi, start_mag, end_mag, steps=steps)


def callRunCurrents():
    """
    Setup function to call the utility function 'runCurrents', see 'utility_functions.py'. Manages necessary inputs.
    """
    inp0 = input('timed mode? (y/n) ')
    configs = []
    timers = []
    char = ''
    while char != 'x':
        inp1 = input('configuration 1\nChannel 1: ')
        inp2 = input('Channel 2: ')
        inp3 = input('Channel 3: ')
        inp4 = input('timer duration: ')
        try:
            a1 = float(inp1)
            b1 = float(inp2)
            c1 = float(inp3)
            configs.append(np.array([a1, b1, c1]))
            timers.append(float(inp4))
        except:
            print('expected numerical value, defaulting to (0,0,1)')
            configs.append(np.array([0, 0, 1]))
            timers.append(0)

        if inp0 == 'y':
            char = input('another config (enter x to end)')
        else:
            char = 'x'
            
    inp5 = input('demagnetize afterwards? (y/n) ')
    
    if inp0 == '':
        subdir = input('Which subdirectory should measurements be saved to? ')
        # with MetrolabTHM1176Node(block_size=20, range='0.3 T', period=0.01, average=5) as node:
        #     # node = MetrolabTHM1176Node(block_size=20, sense_range_upper="0.3 T", period=0.001)
        #     char = input('Calibrate Metrolab sensor? (y/n): ')
        #     if char == 'y':
        #         calibration(node, calibrate=True)

        runCurrents(configs, direct=b'1', subdir=subdir, demagnetize=(inp5=='y'))

    runCurrents(configs, timers, direct=b'1', demagnetize=(inp5=='y'))


def callGenerateVectorField():
    """
    Setup function to call the utility function 'generateMagneticField', see 'utility_functions.py'. Manages necessary inputs.
    """
    inp1 = input('magnitude: ')
    inp2 = input('polar angle (theta): ')
    inp3 = input('azimuthal angle (phi): ')
    subdir = input('Which subdirectory should measurements be saved to? ')
    try:
        magnitude = float(inp1)
    except:
        print('expected numerical value, defaulting to 0')
        magnitude = 0
    try:
        theta = float(inp2)
    except:
        print('expected numerical value, defaulting to 0')
        theta = 0
    try:
        phi = float(inp3)
    except:
        print('expected numerical value, defaulting to 0')
        phi = 0

    # with MetrolabTHM1176Node(block_size=20, range='0.3 T', period=0.01, average=5) as node:
    #         # node = MetrolabTHM1176Node(block_size=20, sense_range_upper="0.3 T", period=0.001)
    #         char = input('Calibrate Metrolab sensor? (y/n): ')
    #         if char == 'y':
    #             calibration(node, calibrate=True)

    generateMagneticField(magnitude, theta, phi, subdir=subdir)


def callFunctionGen():
    """
    Setup function to call the utility function 'switchConfigsAndMeasure', see 'utility_functions.py'. Manages necessary inputs.
    """

    function = input('Function to generate (sqr or sin): ')

    configs = []
    char = ''
    while char != 'x':
        inp1 = input('configuration 1\nChannel 1: ')
        inp2 = input('Channel 2: ')
        inp3 = input('Channel 3: ')
        try:
            a1 = float(inp1)
            b1 = float(inp2)
            c1 = float(inp3)
            configs.append(np.array([a1, b1, c1]))
        except:
            print('expected numerical value, defaulting to (0,0,1)')
            configs.append(np.array([0, 0, 1]))

        if function == 'sqr':
            char = input('another config (enter x to end)')
        else:
            char = 'x'

    inp7 = input('current level (in mA): ')
    try:
        amplitude = int(inp7)
    except:
        print('expected numerical value(s), defaulting to 0')
        amplitude = 1000

    if function == 'sqr':
        inp8 = input('how many times to repeat: ')
    else:
        inp8 = input('Frequency: ')
    try:
        rounds = int(inp8)
    except:
        print('expected numerical value(s), defaulting to 10')
        rounds = 10

    if function == 'sin':
        inp1 = input('Finesse value (divisions per second): ')
        try:
            finesse = int(inp1)
        except:
            print('expected numerical value(s), defaulting to 10')
            finesse = 10
    else:
        finesse = 0

    inp2 = input('Duration? ')
    try:
        duration = float(inp2)
    except:
        print('expected numerical value(s), defaulting to 10*pi')
        duration = 10*np.pi

    measure = input('Measure (y for yes): ')
    if measure == 'y':
        measure = True
    else:
        measure = False
    if function == 'sin':
        measDur = 1.1*duration
    else:
        measDur = rounds * duration

    # functionGenerator(config1, config2, ampl=amplitude, function='sqr', duration=30, frequency=rounds, meas=True, measDur=2.05*rounds*30)
    functionGenerator(configs, ampl=amplitude, function=function, frequency=rounds,
                      finesse=finesse, duration=duration, meas=measure, measDur=measDur)


def feedbackMode():
    """
    Setup function to call the utility function 'callableFeedback', see 'feedback.py'. Manages necessary inputs.
    """
    import pandas as pd

    print('Enter the magnetic Field vector info:')
    source = input('With an input file? (y/n) ')
    if source != 'y':
        coordinates = input('coordinate system: ')
        B_0 = input('Component 1 = ')
        B_1 = input('Component 2 = ')
        B_2 = input('Component 3 = ')
        B_info_arr = [[coordinates, np.array(
            [float(B_0), float(B_1), float(B_2)])]]

    else:
        inpFile = input('Enter a valid configuration file path: ')
        B_info_arr = []
        with open(inpFile, 'r') as f:
            contents = csv.reader(f)
            next(contents)
            for row in contents:
                B_info_arr.append(
                    [row[0], np.array([float(row[1]), float(row[2]), float(row[3])])])

    BVectors = []
    currConfigs = []
    for k in range(len(B_info_arr)):
        B_info = B_info_arr[k]
        BVector, dBdI, cur = fb.callableFeedback(
            B_info, maxCorrection=20, threshold=1, calibrate=True, ECB_active=True)
        BVectors.append(BVector)
        currConfigs.append(cur)

    BVectors = np.array(BVectors)
    currConfigs = np.array(currConfigs)

    subdir = input('Which directory should the output file be saved in? ')
    filename = input('Enter a valid filename(no special chars): ')

    if subdir == '':
        subdir = r'data_sets\linearization_matrices'
    if filename == '' or filename[0] == ' ':
        filename = 'dataset1'

    now = datetime.now().strftime('%y_%m_%d_%H%M%S')
    filename = f'{now}-{filename}.csv'

    filepath = os.path.join(subdir, filename)

    df = pd.DataFrame({'expected Bx [mT]': BVectors[:, 0],
                       'expected By [mT]': BVectors[:, 1],
                       'expected Bz [mT]': BVectors[:, 2],
                       'channel 1 [A]': currConfigs[:, 0],
                       'channel 2 [A]': currConfigs[:, 1],
                       'channel 3 [A]': currConfigs[:, 2]})
    df.to_csv(filepath, index=False, header=True)


if __name__ == '__main__':
    ecbInit = openConnection()
    while ecbInit != 0:
        ecbInit = openConnection()

    MainMenu(ecbInit)
    closeConnection()
