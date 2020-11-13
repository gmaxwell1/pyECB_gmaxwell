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
from datetime import datetime
import threading
import matplotlib.pyplot as plt

########## local imports ##########
import transformations as tr
from main_comm import *
from main_comm import _setCurrents_
from measurements import *
# from modules.analysis_tools import generate_plots
from MetrolabTHM1176.thm1176 import MetrolabTHM1176Node
from modules.general_functions import save_time_resolved_measurement as strm, ensure_dir_exists, sensor_to_magnet_coordinates

##########  Current parameters ##########

desCurrents = [0, 0, 0, 0, 0, 0, 0, 0]  # in milliamps
currDirectParam = b'1'

##########  Vector magnet properties ##########

windings = 508  # windings per coil
resistance = 0.47  # resistance per coil

########## list for storing measured values ##########
returnDict = {'Bx': 0, 'By': 0, 'Bz': 0, 'time': 0}
# order of data: Bx list, By list, Bz list, time list
# threadLock = threading.Lock()
flags = [1]


class myMeasThread(threading.Thread):
    """
    Start a new thread for measuring magnetic field over time with the Metrolab sensor.
    Thread has a name name and multiple member variables.

    kwargs:
        name (str): thread name (default: 'MeasureThread')
        period (float): trigger period, should be in the interval (122e-6, 2.79)
                        (default: 0.1)
        averaging (int): number of measured values to average over.
                            (default: 1)            
        block_size (int): number of measured values to fetch at once.
                            (default: 1)
        duration (float): duration of measurement series
                            (default: 10)
        tempData (bool): Fetch temperature data or not?
                            (default: False)
                            
        self.returnList is a tuple containing the returned values from the measurement.
    """
    def __init__(self, threadID,**kwargs):
        threading.Thread.__init__(self)
        keys = kwargs.keys()
        
        self.threadID = threadID

        if 'name' in keys:
            self.name = kwargs['name']
        else:
            self.name = 'MeasureThread'
            
        if 'period' in keys:
            self.period = kwargs['period']
        else:
            self.period = 0.1
        
        if 'averaging' in keys:
            self.averaging = kwargs['averaging']
        else:
            self.averaging = 1
        
        if 'block_size' in keys:
            self.block_size = kwargs['block_size']
        else:
            self.block_size = 1
            
        if 'duration' in keys:
            self.duration = kwargs['duration']
        else:
            self.duration = 10
            
        if 'tempData' in keys:
            self.tempData = kwargs['tempData']
        else:
            self.tempData = False

    def run(self):
        global returnDict
        
        # threadLock.acquire()
        print("Starting " + self.name)
        
        try:
            returnDict = timeResolvedMeasurement(period=self.period, average=self.averaging,
                                                block_size=self.block_size, duration=self.duration,
                                                return_temp_data=self.tempData)
        except Exception as e:
            print('There was a problem!')
            print(e)
        # threadLock.release()
        print("Finished measuring. {} exiting.".format(self.name))
        
        
class inputThread(threading.Thread):
   def __init__(self, threadID):
        """
        start a new thread with an ID, name and member variable.

        Args:
            threadID (int): An identifier number assigned to the newly created thread.
        """
        threading.Thread.__init__(self)
        self.threadID = threadID

   def run(self):
        global flags
        # Get lock to synchronize threads
        c = input("Hit Enter to quit.\n")
        flags.insert(0, c)
        print('exiting...')


def sweepCurrents(node: MetrolabTHM1176Node, config_list=None, config='z', datadir='config_tests', start_val=0, end_val=1, steps=5, today=True):
    """
    sweeps all currents in the (1,1,1) or (1,0,-1) configuration, meaning we have either a z field or an x-y-plane field, measures magnetic field 
    and stores the measured values in various arrays

    Args:
        config_list (numpy array of size 3, optional): current configuration entered by user, Defaults to np.array([0,0,1])
        config (str, optional): Choose from multiple possible 'configurations' (coil1, coil2, coil3) of currents. Possible values:
            - 'r': randomly generated configuration.
            - 'z': (1,1,1)
            - 'y': (0,1,0)
        Defaults to 'z'.
        datadir (str): directory to save measurements in
        start_val (int, optional): start ramp at. Defaults to 0.
        end_val (int, optional): end ramp at. Defaults to 1.
        steps (int, optional): number of steps. Defaults to 5.
    """
    global currDirectParam
    global desCurrents

    # initialization of all arrays
    all_curr_steps = np.linspace(start_val, end_val, steps)
    mean_values = np.zeros((steps, 3))
    stdd_values = np.zeros((steps, 3))
    expected_fields = np.zeros((steps, 3))
    all_curr_vals = np.zeros((steps, 3))

    # some possible configurations
    current_direction = np.ndarray(3)
    if config_list is not None:
        max_val = np.amax(abs(config_list))
        current_direction[0] = config_list[0] / max_val
        current_direction[1] = config_list[1] / max_val
        current_direction[2] = config_list[2] / max_val
    elif config == 'z':
        current_direction[0] = 1
        current_direction[1] = 1
        current_direction[2] = 1
    elif config == 'y':
        current_direction[0] = -1
        current_direction[1] = 2
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
    fileprefix = '({}_{}_{})_field_meas'.format(int(10*current_direction[0]), 
        int(10*current_direction[1]), int(10*current_direction[2]))
    # folder, 
    if today:
        now = datetime.now().strftime('%y_%m_%d')
        filePath = f'data_sets\{datadir}_{now}'
    else:
        filePath = f'data_sets\{datadir}'

    
    enableCurrents()
    # with MetrolabTHM1176Node() as node:
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
        mean_data, std_data = measure(node, N=5, average=True)
        mean_values[i] = mean_data
        stdd_values[i] = std_data
        expected_fields[i] = B_expected

    demagnetizeCoils()
    # end of measurements
    disableCurrents()
    # saving data section (prepared for plotting)
    saveDataPoints((all_curr_vals / 1000), mean_values,
                   stdd_values, expected_fields, filePath, fileprefix)


def rampVectorField(node: MetrolabTHM1176Node, theta=90, phi=0, start_mag=20, finish_mag=30, rotate=None, steps=10):
    """
    Ramps magnetic field from start_magn to finish_magn in a specified number of steps and over a specified duration (sort of analogous to sweepCurrent)
    or rotates (i.e. 360°) the field vector (start_mag, theta, phi) around the axis given by rotate[0] and rotate[1] (angles in degrees) in that number of steps.
    Measure the magnetic field values, save them to a file.


    Args:
        theta (float): Give the direction of the magnetic field. Angle between field vector and z axis
        phi (float): angle between field vector projection on xy plane and x axis
        start_mag (float): Magnitude in mT
        finish_mag (float): Magnitude in mT, finish_mag > start_mag must hold!!!
        rotate (float, float): tuple of angles which define the axis to rotate around (if at all). rotate[0] is the polar angle (to z axis),
                               rotate[1] is the azimuthal ('xy plane') angle. In degrees (°)
        steps (int): number of measurement steps
    """
    global currDirectParam
    global desCurrents

    if rotate is None:
        magnitudes = np.linspace(start_mag, finish_mag, steps)
    else:
        gamma, psi = rotate
        angles = np.linspace(0, 359.9, steps) # rotation around axis
        magnitude = start_mag
        BStartVector = tr.computeMagneticFieldVector(magnitude, theta, phi)

    mean_values = np.zeros((steps, 3))
    stdd_values = np.zeros((steps, 3))
    expected_fields = np.zeros((steps, 3))
    all_curr_vals = np.zeros((steps, 3))

    # create filename under which to save measurements
    if rotate is None:
        fileprefix = '({}_{})_ramp_meas'.format(int(theta), int(phi))
    else:
        fileprefix = '({}_{})_rotate_meas'.format(int(theta), int(phi))
    
    now = datetime.now().strftime('%y_%m_%d')
    filePath = r'.\data_sets\test_measurements_{}'.format(now)
    
    enableCurrents()
    # with MetrolabTHM1176Node(range="0.3 T", period=0.01) as node:
        # iterate through all possible steps
    for i in range(steps):
        # tentative estimation of resulting B field
        if rotate is None:
            B_expected = tr.computeMagneticFieldVector(magnitudes[i], theta, phi)
        else:
            B_expected = tr.rotationMatrix(BStartVector, psi=psi, theta=gamma, alpha=angles[i])
            
        current_direction = tr.computeCoilCurrents(B_expected, windings, resistance)
        
        # set the current on each channel
        for k in range(3):
            desCurrents[k] = current_direction[k]
        all_curr_vals[i] = current_direction

        setCurrents(desCurrents, currDirectParam)
        # Let the field stabilize
        sleep(0.1)
        # collect measured and expected magnetic field (of the specified sensor in measurements)
        print('measurement nr. ', i+1)
        # see measurements.py for more details
        mean_data, std_data = measure(node, N=5, average=True)
        mean_values[i] = mean_data
        stdd_values[i] = std_data
        expected_fields[i] = B_expected

    demagnetizeCoils()
    # end of measurements
    disableCurrents()
    # saving data section (prepared for plotting)
    saveDataPoints((all_curr_vals / 1000), mean_values,
                   stdd_values, expected_fields, filePath, fileprefix)


def runCurrents(channels, t=0, direct=b'1', subdir='serious_measurements_for_LUT'):
    """
    run arbitrary currents (less than maximum current) on each channel
    when running without a timer, the current can be changed in a menu and the magnetic field can
    be measured with the metrolab sensor.


    Args:
        channels (np.array(list(int)), length 3): current values in [mA]. Make sure not to give more than 3!
        t (int, optional): time for which the magnetic field should be activated. 
            If zero, user can decide whether to change the current or deactivate it. Defaults to 0.
        direct (bytes, optional): current direct parameter (can usually be left alone). Defaults to b'1'.
        subdir (str): if in non-timed mode, the subdirectory in which to save potential measurements
    """
    global currDirectParam
    global desCurrents
        
    currDirectParam = direct
    # copy the computed current values (mA) into the desCurrents list (first 3 positions)
    # cast to int
    for i in range(len(channels)):
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
                '[q] to disable currents\n[c]: get currents\n[r]: Set new currents\n[s]: monitor magnetic field\n[f]: get magnetic field measurement\n')
            if c1 == 'c':
                getCurrents()
            elif c1 == 'r':
                channels[0] = input('Channel 1 current: ')
                channels[1] = input('Channel 2 current: ')
                channels[2] = input('Channel 3 current: ')
                # handle new inputs
                for i in range(len(channels)):
                    try:
                        desCurrents[i] = int(channels[i])
                    except:
                        print(
                            "non-integer value entered, setting channel {} to 0".format(i+1))
                        desCurrents[i] = 0

                setCurrents(desCurrents, currDirectParam)

            elif c1 == 's':
                with MetrolabTHM1176Node(period=0.05, range='0.3 T', average=20) as node:
                    test_thread = inputThread(1)
                    test_thread.start()
                    while flags[0]:
                        newBMeasurement = sensor_to_magnet_coordinates(np.array(node.measureFieldmT()))
                        B_magnitude = np.linalg.norm(newBMeasurement)
                        theta = np.degrees(np.arccos(newBMeasurement[2]/B_magnitude))
                        phi = np.degrees(np.arctan2(newBMeasurement[1], newBMeasurement[0]))
                        if flags[0]:
                            print(f'\rMeasured B field: ({newBMeasurement[0]:.2f}, {newBMeasurement[1]:.2f}, '
                                f'{newBMeasurement[2]:.2f}) / In polar coordinates: ({B_magnitude:.2f}, '
                                f'{theta:.2f}°, {phi:.2f}°)    ', sep='' ,end='', flush=True)
                            
                    flags.insert(0,1)

            elif c1 == 'f':
                try:
                    duration = int(input('Duration of measurement (default is 10s): '))
                except:
                    duration = 10
                params = {'block_size': 20, 'period': 1e-2, 'duration': duration, 'averaging': 5}

                faden = myMeasThread(1,**params)
                faden.start()
                
                faden.join()
                strm(returnDict, r'.\data_sets\{}'.format(subdir), now=True)
                
        demagnetizeCoils()
        disableCurrents()
    else:
        return


def generateMagneticField(magnitude, theta, phi, subdir='serious_measurements_for_LUT'):
    
    global currDirectParam
    global desCurrents

    B_vector = tr.computeMagneticFieldVector(magnitude, theta, phi)
    print(B_vector)
    I_vector = tr.computeCoilCurrents(B_vector)

    # copy the computed current values (mA) into the desCurrents list (first 3 positions)
    # cast to int
    for i in range(len(I_vector)):
        desCurrents[i] = int(I_vector[i])
        
    print(f'Currents on each channel: ({desCurrents[0]}, {desCurrents[1]}, {desCurrents[2]})')

    enableCurrents()
    setCurrents(desCurrents, currDirectParam)
    # wait until user presses q
    c1 = '0'
    while c1 != 'q':
        c1 = input('[q] to disable currents\n[c]: Set new currents\n[r]: Set new mag field\n'
                   '[s]: measure magnetic field and show measurement\n[f]: get magnetic field time-resolved measurement series\n')
        if c1 == 'c':
            channels = [0,0,0]
            channels[0] = input('Channel 1 current: ')
            channels[1] = input('Channel 2 current: ')
            channels[2] = input('Channel 3 current: ')
            # handle new inputs
            for i in range(len(channels)):
                try:
                    desCurrents[i] = int(channels[i])
                except:
                    print(f"non-integer value entered, setting channel {i+1} to 0")
                    desCurrents[i] = 0
            print(f'Currents on each channel: ({desCurrents[0]}, {desCurrents[1]}, {desCurrents[2]})')
            setCurrents(desCurrents, currDirectParam)

        elif c1 == 'r':
            inp1 = input('New magnitude: ')
            inp2 = input('New polar angle (theta): ')
            inp3 = input('New azimuthal angle (phi): ')
            # handle new inputs
            try:
                magnitude = float(inp1)
            except:
                print("non-integer value entered, setting magnitude to 0")
                magnitude = 0
            try:
                theta = float(inp2)
            except:
                print("non-integer value entered, setting theta to 0")
                theta = 0
            try:
                phi = float(inp3)
            except:
                print("non-integer value entered, setting phi to 0")
                phi = 0

            B_vector = tr.computeMagneticFieldVector(magnitude, theta, phi)
            I_vector = tr.computeCoilCurrents(B_vector, windings, resistance)
            # copy the computed current values (mA) into the desCurrents list (first 3 positions)
            # cast to int
            for i in range(len(I_vector)):
                desCurrents[i] = int(I_vector[i])
            print(f'Currents on each channel: ({desCurrents[0]}, {desCurrents[1]}, {desCurrents[2]})')
            setCurrents(desCurrents, currDirectParam)

        elif c1 == 's':
            with MetrolabTHM1176Node(period=0.05, range='0.3 T', average=20) as node:
                test_thread = inputThread(1)
                test_thread.start()
                while flags[0]:
                    newBMeasurement = sensor_to_magnet_coordinates(np.array(node.measureFieldmT()))
                    B_magnitude = np.linalg.norm(newBMeasurement)
                    theta = np.degrees(np.arccos(newBMeasurement[2]/B_magnitude))
                    phi = np.degrees(np.arctan2(newBMeasurement[1], newBMeasurement[0]))
                    if flags[0]:
                        print(f'\rMeasured B field: ({newBMeasurement[0]:.2f}, {newBMeasurement[1]:.2f}, '
                              f'{newBMeasurement[2]:.2f}) / In polar coordinates: ({B_magnitude:.2f}, '
                              f'{theta:.2f}°, {phi:.2f}°)    ', sep='' ,end='', flush=True)
                flags.insert(0,1)

        elif c1 == 'f':
            try:
                duration = int(input('Duration of measurement (default is 10s): '))
            except:
                duration = 10
            params = {'block_size': 20, 'period': 1e-2, 'duration': duration, 'averaging': 5}

            faden = myMeasThread(1,**params)
            faden.start()
            
            faden.join()
            strm(returnDict, r'.\data_sets\{}'.format(subdir), now=True)
    
    demagnetizeCoils()
    disableCurrents()


def functionGenerator(config_list, ampl=1000, function='sin', frequency=1, finesse=10, duration=10*np.pi, meas=False, measDur=0):
    """
    Switch quickly between two current configurations and keep track of the measured fields over time. The time in each state is dt.


    Args:
        config_list (list(s)): List of configurations to change between if function is 'sqr'. Otherwise only the first config
                                       will be used.
        ampl (int, optional): amplitude of the current to be applied in any configuration
        function (str, optional): Either periodically changing between constant current ('sqr') or a sinusoidal current. Defaults to 'sin'.
        frequency (int, optional): sin wave frequency if 'sin' is the function chosen, otherwise the number of times to repeat the cycle.
        finesse (int, optional): Only for 'sin'. How many segments each second is divided into. Defaults to 10.
        duration (int, optional): Duration for which currents are on in sinusoidal mode. Time to stay constant in each state when 'sqr'
                                  is chosen. Defaults to 10*pi.
        meas (bool, optional): Flag for starting a measurement. Defaults to False.
    """
    global currDirectParam
    global desCurrents
    
    if meas:
        params = {'block_size': 30, 'period': 1e-2, 'duration': measDur, 'averaging': 5}
        faden = myMeasThread(1,**params)

    enableCurrents()
    
    if meas:
        faden.start()
                
    if function == 'sin':
        steps = int(duration) * finesse + 1
        tspan = np.linspace(0, duration, steps)
        # dt = round(tspan[1] - tspan[0], 2)
        func1 = ampl * config_list[0][0] * np.sin(2*np.pi*frequency * tspan) # channel 1 values
        func2 = ampl * config_list[0][1] * np.sin(2*np.pi*frequency * tspan) # channel 2 values
        func3 = ampl * config_list[0][2] * np.sin(2*np.pi*frequency * tspan) # channel 3 values
        
        sleep(1/finesse - time() * finesse % 1 / finesse)
        for j in range(len(tspan)):
            desCurrents[0] = int(func1[j])
            desCurrents[1] = int(func2[j])
            desCurrents[2] = int(func3[j])

            _setCurrents_(desCurrents, currDirectParam)
            
            sleep(1/finesse - time() * finesse % 1 / finesse)

    elif function == 'sqr':
        num = len(config_list)
        steps = num * frequency
        # tspan = np.linspace(0, duration, steps)
        # dt = duration/steps
        funcs = [ampl * config for config in config_list]
        
        sleep(duration - time() % duration)
        for j in range(steps):
            
            desCurrents[0] = int(funcs[j % num][0])
            desCurrents[1] = int(funcs[j % num][1])
            desCurrents[2] = int(funcs[j % num][2])

            setCurrents(desCurrents, currDirectParam)
            
            sleep(duration - time() % duration)

    demagnetizeCoils()
    
    if meas:
        faden.join()
        
    disableCurrents()

    strm(returnDict, r'.\data_sets\time_measurements_12_11', now=True)

        
    
if __name__ == "__main__":
    params = {'block_size': 20, 'period': 1e-2, 'duration': 30, 'averaging': 5}
   
    faden = myMeasThread(1,**params)
    faden.start()
    sleep(1)
    openConnection()
    enableCurrents()
    functionGenerator([1,1,1], ampl=3000, function='sin', frequency=1, finesse=20, duration=4*np.pi, meas=False, measDur=1.2*4*np.pi)
    demagnetizeCoils()
    disableCurrents()

    faden.join()
    closeConnection()

    
    strm(returnDict, r'.\data_sets\time_measurements_03_11', now=True)
