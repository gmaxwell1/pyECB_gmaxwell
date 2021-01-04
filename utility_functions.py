"""
filename: utility_functions.py

This collection of functions has functions that can be used to manipulate the currents on the ECB channels. 
The user can choose from various methods to set currents on the different channels and thus generate a magnetic 
field with the vector magnet. For example, we can sweep through current values on the different channels and
simultaneously measure the actual magnetic field to compare the theory with the actual results.

Author: Maxwell Guerne-Kieferndorf (QZabre)
        gmaxwell@student.ethz.ch

Date: 15.12.2020
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
from modules.arduinoPythonInterface import ArduinoUno, saveTempData


##########  Current parameters ##########
desCurrents = [0, 0, 0, 0, 0, 0, 0, 0]  # in milliamps
currDirectParam = b'1'

##########  Vector magnet properties ##########

windings = 508  # windings per coil
resistance = 0.47  # resistance per coil

########## list for storing measured values ##########
returnDict = {'Bx': 0, 'By': 0, 'Bz': 0, 'time': 0, 'temp': 0}
# order of data: Bx list, By list, Bz list, time list
threadLock = threading.Lock()
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
    """

    def __init__(self, threadID, **kwargs):
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

    def run(self):
        global returnDict

        # threadLock.acquire()
        print("Starting " + self.name)

        try:
            returnDict = timeResolvedMeasurement(period=self.period, average=self.averaging,
                                                 block_size=self.block_size, duration=self.duration)
        except Exception as e:
            print('There was a problem!')
            print(e)
        # threadLock.release()
        print("Finished measuring. {} exiting.".format(self.name))


class inputThread(threading.Thread):
    def __init__(self, threadID):
        """
        Waits for the user to press enter.

        Args:
            threadID (int): An identifier number assigned to the newly created thread.
        """
        threading.Thread.__init__(self)
        self.threadID = threadID

    def run(self):
        # global variable flags is modified and can then be read/modified
        # by other threads. Only this thread will append a zero to flags.
        global flags
        c = input("Hit Enter to quit.\n")
        #make sure there are no concurrency issues
        threadLock.acquire()
        flags.insert(0, c)
        threadLock.release()
        
        print('exiting...')
        
        
class timerThread(threading.Thread):
    def __init__(self, threadID, duration):
        """
        start a new thread with an ID, name and member variable. Counts down time.

        Args:
            threadID (int): An identifier number assigned to the newly created thread.
            duration (float): timer duration
        """
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.startTime = time()
        self.countdownDuration = duration

    def run(self):
        timeDiff = time() - self.startTime
        while timeDiff < self.countdownDuration:
            print(f'\rtime remaining: {int((self.countdownDuration - timeDiff))//3600} hours, {int((self.countdownDuration - timeDiff))//60 % 60} '
                  f'minutes and {(self.countdownDuration - timeDiff)%60:.0f} seconds', end='', sep='', flush=False)
            
            timeDiff = time() - self.startTime
            sleep(0.096)


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
        # max_val = np.amax(abs(config_list))
        current_direction[0] = config_list[0]
        current_direction[1] = config_list[1]
        current_direction[2] = config_list[2]
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
        sleep(0.5)
        # collect measured and expected magnetic field (of the specified sensor in measurements)
        print('measurement nr. ', i+1)
        # see measurements.py for more details
        mean_data, std_data = measure(node, N=7, average=True)
        mean_values[i] = mean_data
        stdd_values[i] = std_data
        expected_fields[i] = B_expected

    demagnetizeCoils(current_direction)
    # end of measurements
    disableCurrents()
    # saving data section (prepared for plotting)
    saveDataPoints((all_curr_vals / 1000), mean_values,
                   stdd_values, expected_fields, filePath, fileprefix)
    

def gridSweep(node: MetrolabTHM1176Node, inpFile=r'config_files\configs_numvals2_length4.csv', datadir='config_tests',
              current_val=0, demagnetize=False, today=True):
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
    # all_curr_steps = np.linspace(start_val, end_val, steps)
    mean_values = []
    stdd_values = []
    expected_fields = []
    all_curr_vals = []
    
    enableCurrents()

    with open(inpFile, 'r') as f:
        contents = csv.reader(f)
        next(contents)
        for i, row in enumerate(contents):
            demagnetizeCoils(current_config = 5000*np.ones(3))
            
            config = np.array(
                [float(row[0]), float(row[1]), float(row[2])])
            
            for k in range(3):
                desCurrents[k] = config[k]*current_val
            all_curr_vals.append(config*current_val)
            
            #  tentative estimation of resulting B field
            B_expected = tr.computeMagField(config*current_val, windings)
            setCurrents(desCurrents, currDirectParam)
            # Let the field stabilize
            sleep(0.5)
            # collect measured and expected magnetic field (of the specified sensor in measurements)
            print('measurement nr. ', i+1)
            # see measurements.py for more details
            mean_data, std_data = measure(node, N=7, average=True)
            mean_values.append(mean_data)
            stdd_values.append(std_data)
            expected_fields.append(B_expected)
 
    # create subdirectory to save measurements
    fileprefix = 'field_meas'
    # folder,
    if today:
        now = datetime.now().strftime('%y_%m_%d')
        filePath = f'data_sets\{datadir}_{now}'
    else:
        filePath = f'data_sets\{datadir}'


    if demagnetize:
        demagnetizeCoils()
    # end of measurements
    disableCurrents()
    # saving data section (prepared for plotting)
    saveDataPoints((np.array(all_curr_vals) / 1000), np.array(mean_values),
                   np.array(stdd_values), np.array(expected_fields), filePath, fileprefix)


def sweepHysteresis(config_list=None, datadir='hysteresis_tests', end_val=1, steps=5, today=False):
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
    half_curr_steps = np.linspace(0, end_val, steps)
    mean_values = np.zeros((5*steps, 3))
    stdd_values = np.zeros((5*steps, 3))
    expected_fields = np.zeros((5*steps, 3))
    all_curr_vals = np.zeros((5*steps, 3))

    # some possible configurations
    current_direction = np.ndarray(3)
    if config_list is not None:
        max_val = np.amax(abs(config_list))
        current_direction[0] = config_list[0] / max_val
        current_direction[1] = config_list[1] / max_val
        current_direction[2] = config_list[2] / max_val
    else:
        print('invalid input!')
        return

    # create subdirectory to save measurements
    fileprefix = '({}_{}_{})_hyst_meas'.format(int(10*current_direction[0]),
                                               int(10*current_direction[1]), int(10*current_direction[2]))
    # folder,
    if today:
        now = datetime.now().strftime('%y_%m_%d')
        filePath = f'data_sets\{datadir}_{now}'
    else:
        filePath = f'data_sets\{datadir}'

    enableCurrents()

    demagnetizeCoils()

    # with MetrolabTHM1176Node(average = 1,range = '0.3T',block_size = 5, period = 0.02) as node:
    node = MetrolabTHM1176Node(
        average=1, range='0.3T', block_size=5, period=0.02)
    # iterate through all possible steps
    sign = [1, 1, -1, -1, 1]
    for j in range(5):
        for i in range(steps):
            if j == 1 or j == 3:
                k = steps - 1 - i
            else:
                k = i
            # set the current on each channel
            desCurrents[0:3] = sign[j]*current_direction*half_curr_steps[k]
            all_curr_vals[i + j*steps] = sign[j] * \
                current_direction*half_curr_steps[k]

            # tentative estimation of resulting B field
            B_expected = tr.computeMagField(
                all_curr_vals[i + j*steps], windings)

            setCurrents(desCurrents, currDirectParam)
            # Let the field stabilize
            sleep(0.4)
            # collect measured and expected magnetic field (of the specified sensor in measurements)
            print('measurement nr. ', i+1)
            # see measurements.py for more details
            mean_data, std_data = measure(node, N=5, average=True)
            mean_values[i + j*steps] = mean_data
            stdd_values[i + j*steps] = std_data
            expected_fields[i + j*steps] = B_expected

    demagnetizeCoils(current_direction)
    # end of measurements
    disableCurrents()
    # saving data section (prepared for plotting)
    saveDataPoints((all_curr_vals / 1000), mean_values,
                   stdd_values, expected_fields, filePath, fileprefix)


def rampVectorField(node: MetrolabTHM1176Node, theta=90, phi=0, start_mag=20, finish_mag=30, rotate=None, steps=10, demagnetize=False):
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
        angles = np.linspace(0, 359.9, steps)  # rotation around axis
        magnitude = start_mag
        BStartVector = tr.computeMagneticFieldVector(magnitude, theta, phi)

    rounds = 80
    mean_values = np.zeros((rounds*steps, 3))
    stdd_values = np.zeros((rounds*steps, 3))
    expected_fields = np.zeros((rounds*steps, 3))
    all_curr_vals = np.zeros((rounds*steps, 3))

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
    for r in range(rounds):
        for i in range(steps):
            demagnetizeCoils(current_config = 5000*np.ones(3))
            
            # tentative estimation of resulting B field
            if rotate is None:
                B_expected = tr.computeMagneticFieldVector(
                    magnitudes[i], theta, phi)
            else:
                B_expected = tr.rotationMatrix(
                    BStartVector, psi=psi, theta=gamma, alpha=angles[i])

            current_direction = tr.computeCoilCurrents(
                B_expected, windings, resistance)

            # set the current on each channel
            for k in range(3):
                desCurrents[k] = current_direction[k]
            all_curr_vals[i] = current_direction

            setCurrents(desCurrents, currDirectParam)
            # Let the field stabilize
            sleep(0.1)
            # collect measured and expected magnetic field (of the specified sensor in measurements)
            print('measurement nr. ', i+1, r)
            # see measurements.py for more details
            mean_data, std_data = measure(node, N=5, average=True)
            mean_values[i+r*steps] = mean_data
            stdd_values[i+r*steps] = std_data
            expected_fields[i+r*steps] = B_expected
            
            

    if demagnetize:
        demagnetizeCoils()
    # end of measurements
    disableCurrents()
    # saving data section (prepared for plotting)
    saveDataPoints((all_curr_vals / 1000), mean_values,
                   stdd_values, expected_fields, filePath, fileprefix)


def runCurrents(config_list, t=[], direct=b'1', subdir='serious_measurements_for_LUT', demagnetize=False):
    """
    run arbitrary currents (less than maximum current) on each channel
    when running without a timer, the current can be changed in a menu and the magnetic field can
    be measured with the metrolab sensor.


    Args:
        config_list (np.array(list(int)), length 3): current values in [mA]. Make sure not to give more than 3!
        t (int, optional): timer duration list. multiple timers -> different currents will be set for different 
                           amounts of time.
            If zero, user can decide whether to change the current or deactivate it. Defaults to 0.
        direct (bytes, optional): current direct parameter (can usually be left alone). Defaults to b'1'.
        subdir (str): if in non-timed mode, the subdirectory in which to save potential measurements
    """
    global currDirectParam
    global desCurrents

    currDirectParam = direct
    # copy the computed current values (mA) into the desCurrents list (first 3 positions)
    # cast to int

    # user specified time
    enableCurrents()

    # on until interrupted by user
    if len(t) == 0 or t[0] == 0:
        channels = config_list[0]
        for i in range(len(channels)):
            desCurrents[i] = int(channels[i])
            
        setCurrents(desCurrents, currDirectParam)
        # wait until user presses enter
        c1 = '0'
        while c1 != 'q':
            c1 = input(
                '[q] to disable currents\n[c]: get currents\n[r]: Set new currents\n[s]: monitor magnetic field (does not work)\n[f]: get magnetic field measurement (does not work)\n')
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
            #     with MetrolabTHM1176Node(period=0.05, range='0.3 T', average=20) as node:
                test_thread = inputThread(1)
                test_thread.start()
                sleep(0.1)
                while flags[0]:
            #             newBMeasurement = sensor_to_magnet_coordinates(
            #                 np.array(node.measureFieldmT()))
                    newBMeasurement = np.random.randn((3)) * 10
                    B_magnitude = np.linalg.norm(newBMeasurement)
                    theta = np.degrees(
                        np.arccos(newBMeasurement[2]/B_magnitude))
                    phi = np.degrees(np.arctan2(
                        newBMeasurement[1], newBMeasurement[0]))
                    if flags[0]:
                        print(f'\rMeasured B field: ({newBMeasurement[0]:.2f}, {newBMeasurement[1]:.2f}, '
                                f'{newBMeasurement[2]:.2f}) / In polar coordinates: ({B_magnitude:.2f}, '
                                f'{theta:.2f}°, {phi:.2f}°)    ', sep='', end='', flush=True)
                    sleep(0.5)

                threadLock.acquire()
                flags.insert(0, 1)
                threadLock.release()

            # elif c1 == 's':
            #     with MetrolabTHM1176Node(period=0.05, range='0.3 T', average=20) as node:
            #         test_thread = inputThread(1)
            #         test_thread.start()
            #         while flags[0]:
            #             newBMeasurement = sensor_to_magnet_coordinates(
            #                 np.array(node.measureFieldmT()))
            #             B_magnitude = np.linalg.norm(newBMeasurement)
            #             theta = np.degrees(
            #                 np.arccos(newBMeasurement[2]/B_magnitude))
            #             phi = np.degrees(np.arctan2(
            #                 newBMeasurement[1], newBMeasurement[0]))
            #             if flags[0]:
            #                 print(f'\rMeasured B field: ({newBMeasurement[0]:.2f}, {newBMeasurement[1]:.2f}, '
            #                       f'{newBMeasurement[2]:.2f}) / In polar coordinates: ({B_magnitude:.2f}, '
            #                       f'{theta:.2f}°, {phi:.2f}°)    ', sep='', end='', flush=True)

            #         flags.insert(0, 1)

            # elif c1 == 'f':
            #     try:
            #         duration = int(
            #             input('Duration of measurement (default is 10s): '))
            #     except:
            #         duration = 10
            #     params = {'block_size': 20, 'period': 1e-2,
            #               'duration': duration, 'averaging': 5}

            #     faden = myMeasThread(1, **params)
            #     faden.start()

            #     faden.join()
            #     strm(returnDict, r'.\data_sets\{}'.format(subdir), now=True)
    else:
        for index, timer in enumerate(t):
            channels = config_list[index]
            for i in range(len(channels)):
                desCurrents[i] = int(channels[i])
            
            setCurrents(desCurrents, currDirectParam)
            # prevent the connection with the ECB from timing out for long measurements.
            if timer < 500:
                countdown = timerThread(0,timer)
                countdown.start()
                sleep(timer)
                countdown.join()
            else:
                countdown = timerThread(0,timer)
                countdown.start()
                starttime = time()
                while time() - starttime < timer:
                    pause = min(500, timer - (time() - starttime))
                    sleep(pause)
                    getCurrents()
                countdown.join()
      
    if demagnetize:
        demagnetizeCoils()
        
    disableCurrents()
        

def generateMagneticField(vectors, t=[], subdir='serious_measurements_for_LUT', demagnetize=False):
    """
    A magnetic field is generated in an arbitrary direction which is specified by the user. The currents
    set on the different channels are computed with a linear model. See transformations.py.

    Args:
        vectors (list, optional): Nx3 vector of float values, containing N magnetic field vectors
                                  in spherical coordinates. Order: (B,theta,phi), theta is the polar
                                  angle and phi is the azimuthal angle.
        t (list, optional): timer durations for each field. May not be longer than 'vectors'. Defaults to [].
        subdir (str, optional): Subdirectory where any measurements will be saved. Defaults to 'serious_measurements_for_LUT'.
        demagnetize (bool, optional): If True, coils will be demagnetized after the magnetic field is deactivated. Defaults to False.
    """
    global currDirectParam
    global desCurrents

    enableCurrents()
    if len(t) == 0 or t[0] == 0:
        B_Field = vectors[0]
        B_Field_cartesian = tr.computeMagneticFieldVector(B_Field[0], B_Field[1], B_Field[2])
        channels = tr.computeCoilCurrents(B_Field_cartesian)
        for i in range(len(channels)):
            desCurrents[i] = int(channels[i])
            
        print(f'Currents on each channel: ({desCurrents[0]}, {desCurrents[1]}, {desCurrents[2]})')
        setCurrents(desCurrents, currDirectParam)
        # wait until user presses q
        c1 = '0'
        while c1 != 'q':
            c1 = input('[q] to disable currents\n[c]: Set new currents\n[r]: Set new mag field\n'
                       '[s]: monitor magnetic field (does not work)\n'
                       '[f]: get magnetic field time-resolved measurement series (does not work)\n')
            
            if c1 == 'c':
                channels = [0, 0, 0]
                channels[0] = input('Channel 1 current: ')
                channels[1] = input('Channel 2 current: ')
                channels[2] = input('Channel 3 current: ')
                # handle new inputs
                for i in range(len(channels)):
                    try:
                        desCurrents[i] = int(channels[i])
                    except:
                        print(
                            f"non-integer value entered, setting channel {i+1} to 0")
                        desCurrents[i] = 0
                print(
                    f'Currents on each channel: ({desCurrents[0]}, {desCurrents[1]}, {desCurrents[2]})')
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
                print(
                    f'Currents on each channel: ({desCurrents[0]}, {desCurrents[1]}, {desCurrents[2]})')
                setCurrents(desCurrents, currDirectParam)

            elif c1 == 's':
            #     with MetrolabTHM1176Node(period=0.05, range='0.3 T', average=20) as node:
                test_thread = inputThread(1)
                test_thread.start()
                sleep(0.1)
                while flags[0]:
            #             newBMeasurement = sensor_to_magnet_coordinates(
            #                 np.array(node.measureFieldmT()))
                    newBMeasurement = np.random.randn((3)) * 10
                    B_magnitude = np.linalg.norm(newBMeasurement)
                    theta = np.degrees(
                        np.arccos(newBMeasurement[2]/B_magnitude))
                    phi = np.degrees(np.arctan2(
                        newBMeasurement[1], newBMeasurement[0]))
                    if flags[0]:
                        print(f'\rMeasured B field: ({newBMeasurement[0]:.2f}, {newBMeasurement[1]:.2f}, '
                                f'{newBMeasurement[2]:.2f}) / In polar coordinates: ({B_magnitude:.2f}, '
                                f'{theta:.2f}°, {phi:.2f}°)    ', sep='', end='', flush=True)
                    sleep(0.5)

                threadLock.acquire()
                flags.insert(0, 1)
                threadLock.release()

            # elif c1 == 'f':
            #     try:
            #         duration = int(
            #             input('Duration of measurement (default is 10s): '))
            #     except:
            #         duration = 10
            #     params = {'block_size': 20, 'period': 1e-2,
            #               'duration': duration, 'averaging': 5}

            #     faden = myMeasThread(1, **params)
            #     faden.start()

            #     faden.join()
            #     strm(returnDict, r'.\data_sets\{}'.format(subdir), now=True)
    else:
        for index, timer in enumerate(t):
            B_Field = vectors[index]
            B_Field_cartesian = tr.computeMagneticFieldVector(B_Field[0], B_Field[1], B_Field[2])
            channels = tr.computeCoilCurrents(B_Field_cartesian)
            for i in range(len(channels)):
                desCurrents[i] = int(channels[i])
            
            setCurrents(desCurrents, currDirectParam)
            # prevent the connection with the ECB from timing out for long measurements.
            if timer < 500:
                countdown = timerThread(0,timer)
                countdown.start()
                sleep(timer)
                countdown.join()
            else:
                countdown = timerThread(0,timer)
                countdown.start()
                starttime = time()
                while time() - starttime < timer:
                    pause = min(500, timer - (time() - starttime))
                    sleep(pause)
                    getCurrents()
                countdown.join()
                
    if demagnetize:
        demagnetizeCoils()
        
    disableCurrents()

def callTempFieldMeasurement(demagnetize=True, today=False, demagnet_current=1000):
    """
    Set currents to constant value for a defined duration, afterwards to zero for the same duration. 
    Measure magnetic field and temperature simultaneously.
    """
    global currDirectParam
    global desCurrents
    
    datadir = input(
        'Enter a valid directory name to save measurement data in: ')
    
    if datadir == '':
        datadir='test_measurements'
    
    inp1 = input('constant currents in mA = ')
    try:
        start_val = int(inp1)
    except:
        print('expected numerical value, defaulting to 1000 mA')
        start_val = 1000
    
    inp2 = input('duration [s] of constant currents (half the total duration) = ')
    try:
        duration = int(inp2)
    except:
        print('expected numerical value, defaulting to 100 s')
        duration = 100
    
    # initialize Matrolab sensor and ensure that it is at the correct position
    node = MetrolabTHM1176Node(block_size=20, range='0.3 T', period=0.01, average=1)# as node:
    gotoPosition()
    
    # initialize temperature sensor and measurement routine and start measuring
    arduino = ArduinoUno('COM7')
    measure_temp = threading.Thread(target=arduino.getTemperatureMeasurements)
    measure_temp.start()
    t_start_measurement = time()
    
    # initialization of all arrays
    mean_values = []
    stdd_values = []
    all_curr_vals = []
    time_values = []
    
    # enable currents at ECB
    enableCurrents()

    # start with demagnetizing, wait a second afterwards
    print('demagnetize')
    demagnetizeCoils(current_config = demagnet_current*np.ones(3))
    sleep(1)
    
    # first half with constant, nonzero currents
    print('start first half')
    t_start = time()
    config = start_val * np.ones(3)
    
    # set currents before starting first half 
    for k in range(3):
        desCurrents[k] = config[k]
    setCurrents(desCurrents, currDirectParam)
    sleep(0.2)
    
    while time() - t_start < duration:
        
        # save currents
        all_curr_vals.append(config)
        
        # print progress
        ratio = (time() - t_start) / duration
        left = int(ratio * 30)
        right = 30-left
        print('\rprogress: [' + '#' * left + ' ' * right + ']',
              f' {ratio * 100:.0f}%', sep='', end='', flush=True)
        
        # see measurements.py for more details
        mean_data, std_data = measure(node, N=7, average=True)
        time_values.append(time() - t_start_measurement)
        mean_values.append(mean_data)
        stdd_values.append(std_data)
    
    # disable and enable currents, thereby trying to avoid communication timeout with ECB
    disableCurrents()
    enableCurrents()
        
    # try demagnetizing in between both halfs
    print('\ndemagnetize')
    try:
        demagnetizeCoils(current_config = demagnet_current*np.ones(3))
    except TimeoutError:
        # if there is a communication error with ECB, skip demagnetizing and continue measuring second half 
        print('skip demagnetizing')
        
    finally:
        # disable currents again 
        disableCurrents()
        
        # new config: zero currents
        config = np.zeros(3)
        
        # second half with zero currents
        print('\nstart second half')
        t_start = time()
        while time() - t_start < duration:
            
            # save current configuration
            all_curr_vals.append(config)
        
            # print progress
            ratio = (time() - t_start) / duration
            left = int(ratio * 30)
            right = 30-left
            print('\rprogress: [' + '#' * left + ' ' * right + ']',
                f' {ratio * 100:.0f}%', sep='', end='', flush=True)
            
            # see measurements.py for more details
            mean_data, std_data = measure(node, N=7, average=True)
            time_values.append(time() - t_start_measurement)
            mean_values.append(mean_data)
            stdd_values.append(std_data)
        
        # enable currents at ECB  for final demagnetization
        enableCurrents()
        # try demagnetizing 
        print('\ndemagnetize')
        try:
            demagnetizeCoils(current_config = demagnet_current*np.ones(3))
        except TimeoutError:
            # if a communication error with ECB occures, simply skip demagnetizing
            print('skip demagnetizing')
            
        finally:
            # in any case, collect data and save them!
            
            # end of measurements
            disableCurrents()
        
            # saving data section (prepared for plotting)
            directory = '.\\data_sets\\temperature_drift'
            data_filename_postfix = 'measured_field_at_const_currents'
            ensure_dir_exists(directory, verbose=False)
            
            # after field measurements are over, also stop temperature measurement
            arduino.stop = True
            measure_temp.join()
            temp_file_name = 'measured_temp_at_const_currents'
            saveTempData(arduino.data_stack, directory = directory, filename_suffix = temp_file_name)
            
            # save measured fields and corresponding currents
            I = np.array(all_curr_vals) / 1000
            mean_data = np.array(mean_values)
            std_data = np.array(stdd_values)
            time_data = np.array(time_values)
        
            df = pd.DataFrame({'time [s]': time_data,
                                'channel 1 [A]': I[:, 0],
                                'channel 2 [A]': I[:, 1],
                                'channel 3 [A]': I[:, 2],
                                'mean Bx [mT]': mean_data[:, 0],
                                'mean By [mT]': mean_data[:, 1],
                                'mean Bz [mT]': mean_data[:, 2],
                                'std Bx [mT]': std_data[:, 0],
                                'std By [mT]': std_data[:, 1],
                                'std Bz [mT]': std_data[:, 2]})
            print('success!')

            now = datetime.now().strftime('%y_%m_%d_%H-%M-%S')
            output_file_name = f'{now}_{data_filename_postfix}.csv'
            file_path = os.path.join(directory, output_file_name)
            df.to_csv(file_path, index=False, header=True)
    
    
    
# def functionGenerator(config_list, ampl=1000, function='sin', frequency=1, finesse=10, duration=10*np.pi, meas=False, measDur=0):
    # """
    # Switch quickly between two current configurations and keep track of the measured fields over time. The time in each state is dt.


    # Args:
    #     config_list (list(s)): List of configurations to change between if function is 'sqr'. Otherwise only the first config
    #                                    will be used.
    #     ampl (int, optional): amplitude of the current to be applied in any configuration
    #     function (str, optional): Either periodically changing between constant current ('sqr') or a sinusoidal current. Defaults to 'sin'.
    #     frequency (int, optional): sin wave frequency if 'sin' is the function chosen, otherwise the number of times to repeat the cycle.
    #     finesse (int, optional): Only for 'sin'. How many segments each second is divided into. Defaults to 10.
    #     duration (int, optional): Duration for which currents are on in sinusoidal mode. Time to stay constant in each state when 'sqr'
    #                               is chosen. Defaults to 10*pi.
    #     meas (bool, optional): Flag for starting a measurement. Defaults to False.
    # """
    # global currDirectParam
    # global desCurrents

    # if meas:
    #     params = {'block_size': 30, 'period': 1e-2,
    #               'duration': measDur, 'averaging': 5}
    #     faden = myMeasThread(1, **params)

    # enableCurrents()

    # if meas:
    #     faden.start()

    # if function == 'sin':
    #     steps = int(duration) * finesse + 1
    #     tspan = np.linspace(0, duration, steps)
    #     # dt = round(tspan[1] - tspan[0], 2)
    #     # channel 1 values
    #     func1 = ampl * config_list[0][0] * np.sin(2*np.pi*frequency * tspan)
    #     # channel 2 values
    #     func2 = ampl * config_list[0][1] * np.sin(2*np.pi*frequency * tspan)
    #     # channel 3 values
    #     func3 = ampl * config_list[0][2] * np.sin(2*np.pi*frequency * tspan)

    #     sleep(1/finesse - time() * finesse % 1 / finesse)
    #     for j in range(len(tspan)):
    #         desCurrents[0] = int(func1[j])
    #         desCurrents[1] = int(func2[j])
    #         desCurrents[2] = int(func3[j])

    #         _setCurrents_(desCurrents, currDirectParam)

    #         sleep(1/finesse - time() * finesse % 1 / finesse)

    # elif function == 'sqr':
    #     num = len(config_list)
    #     steps = num * frequency
    #     # tspan = np.linspace(0, duration, steps)
    #     # dt = duration/steps
    #     funcs = [ampl * config for config in config_list]

    #     sleep(duration - time() % duration)
    #     for j in range(steps):

    #         desCurrents[0] = int(funcs[j % num][0])
    #         desCurrents[1] = int(funcs[j % num][1])
    #         desCurrents[2] = int(funcs[j % num][2])

    #         setCurrents(desCurrents, currDirectParam)

    #         sleep(duration - time() % duration)

    # demagnetizeCoils()

    # if meas:
    #     faden.join()

    # disableCurrents()


if __name__ == "__main__":
    # params = {'block_size': 40, 'period': 0.01, 'duration': 40, 'averaging': 1}
  
    # faden = myMeasThread(1, **params)
    # faden.start()
    countdown = timerThread(0, 7200)
    countdown.start()
    countdown.join()

    # # magnitude = 100
    # # theta = 120
    # # phi = 280
    # # B_vector = tr.computeMagneticFieldVector(magnitude, theta, phi)
    # # print(B_vector)
    # # I_vector = tr.computeCoilCurrents(B_vector)
    # # ampl = np.amax(np.abs(I_vector))
    # # print(I_vector)
    # # returnDict = timeResolvedMeasurement(period=period, average=averaging,
    # #                         block_size=block_size, duration=duration)
    # sleep(10)

    # # openConnection()
    # # enableCurrents()
    # sleep(10)
    # # demagnetizeCoils()
    
    # faden.join()

    # # closeConnection()
    
    # # print(returnDict)

    # strm(returnDict, r'data_sets\noise_measurements', 'zero_field_close_withoutECB', now=True)
