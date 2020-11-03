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

########## local imports ##########
from utility_functions import *
from main_comm import *
from measurements import calibration
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
            # print(
                # '[4]: generate magnetic field (specify polar and azimuthal angles, magnitude)')
            print(
                '[5]: Generate a time varying field (sin and sqr waves are possible currently)')
            print('[6]: get ECB status\n')
            print('[r] roll a die\n')

            c1 = input()

            if c1 == '1':
                inp0 = input('Automatic or manual input? (m or nothing): ')
                datadir = input('Enter a valid directory name to save measurement data in: ')
                if datadir == '':
                    callCurrentSweep(inp0)
                else:
                    callCurrentSweep(inp0, datadir)

            # tentative implementation. Actual I-to-B actuation matrix needed. Many other features not added yet.
            elif c1 == '2':
                callSweepVectorField()
                
            elif c1 == '3':
                callRunCurrents()
               
            # tentative implementation. Actual I-to-B actuation matrix needed.
            # elif c1 == '4':
            #     callGenerateVectorField()

            elif c1 == '5':
                callFunctionGen()
                
            elif c1 == '6':
                print(getStatus())
                
            elif c1 == 'r':
                # just for fun :D
                c1 = np.random.randint(1, 7)
                print('doing option ', c1, ': ')
                c1 = str(c1)

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

        c1 = input('Automatic exit after finish? (x for yes): ')

        with MetrolabTHM1176Node(block_size=20, sense_range_upper="0.3 T", period=0.01, average=10) as node:
            char = input('Calibrate Metrolab sensor? (y/n): ')
            if char == 'y':
                calibration(node, calibrate=True)

            with open(inpFile, 'r') as f:
                contents = csv.reader(f)
                next(contents)
                for row in contents:
                    config = np.array(
                        [float(row[0]), float(row[1]), float(row[2])])
                    sweepCurrents(config_list=config, start_val=start_val, datadir=datadir,
                                  end_val=end_val, steps=steps, node=node)

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

        c1 = input('Automatic exit after finish? (x for yes): ')

        with MetrolabTHM1176Node(block_size=20, sense_range_upper="0.3 T", period=0.01, average=10) as node:
            char = input('Calibrate Metrolab sensor? (y/n): ')
            if char == 'y':
                calibration(node, calibrate=True)

            sweepCurrents(config_list=config, start_val=start_val, datadir=datadir,
                          end_val=end_val, steps=steps, node=node)

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

        c1 = input('Automatic exit after finish? (x for yes): ')

        while randomRuns > 0:
            sweepCurrents(config=config, start_val=start_val,
                          end_val=end_val, steps=steps, datadir=datadir,
                          node=MetrolabTHM1176Node(block_size=20, range='0.3 T', period=0.01, average=5))
            randomRuns = randomRuns-1
            
            
def callSweepVectorField():
    """
    Setup function to call the utility function 'rampVectorField', see 'utility_functions.py'. Manages necessary inputs.
    """
    inp1 = input('Angle to z axis in deg = ')
    inp2 = input('Angle to x axis in deg = ')
    inp3 = input('starting magnitude in mT = ')
    rot = input('Rotate constant magnitude field around a specified axis? (otherwise ramp field magnitude in a constant direction): ')
    
    if rot == 'y':
        axisang1 = input('Axis polar angle: ')
        axisang2 = input('Axis azimuthal angle: ')
        try:
            axis = (int(axisang1), int(axisang2))
        except:
            print('expected numerical value, defaulting to (0,0) or z-axis')
            axis = (0,0)
        
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

    with MetrolabTHM1176Node(block_size=20, sense_range_upper="0.3 T", period=0.01) as node:
    # node = MetrolabTHM1176Node(block_size=20, sense_range_upper="0.3 T", period=0.001)
        char = input('Calibrate Metrolab sensor? (y/n): ')
        if char == 'y':
            calibration(node, calibrate=True)
            
        if rot == 'y':
            rampVectorField(node, theta, phi, start_mag, steps=steps, rotate=axis)
        else:
            rampVectorField(node, theta, phi, start_mag, end_mag, steps=steps)
            

def callRunCurrents():
    """
    Setup function to call the utility function 'runCurrents', see 'utility_functions.py'. Manages necessary inputs.
    """
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
        subdir = input('Which subdirectory should measurements be saved to? ')
        runCurrents(
            np.array([coil1, coil2, coil3]), t=0, direct=b'1', subdir=subdir)
    else:
        try:
            timer = int(inp4)
            c1 = input(
                'Automatic Termination after timer? (x for yes): ')
        except:
            print('expected numerical value, defaulting to 0')
            timer = 0
        runCurrents(
            np.array([coil1, coil2, coil3]), timer, direct=b'1')
        
        
def callGenerateVectorField():
    """
    Setup function to call the utility function 'generateMagneticField', see 'utility_functions.py'. Manages necessary inputs.
    """
    pass
        
        
def callFunctionGen():
    """
    Setup function to call the utility function 'switchConfigsAndMeasure', see 'utility_functions.py'. Manages necessary inputs.
    """
    
    function = input('Function to generate (sqr or sin): ')
    
    configs = []
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
            configs = configs[0]
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

    # functionGenerator(config1, config2, ampl=amplitude, function='sqr', duration=30, frequency=rounds, meas=True, measDur=2.05*rounds*30)
    functionGenerator(configs, ampl=amplitude, function=function, frequency=rounds, finesse=finesse, duration=duration, meas=measure, measDur=1.1*duration)
    

if __name__ == '__main__':
    ecbInit = openConnection()
    while ecbInit != 0:
        ecbInit = openConnection()

    MainMenu(ecbInit)
    closeConnection()

