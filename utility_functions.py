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
import threading
import matplotlib.pyplot as plt

########## local imports ##########
import transformations as tr
from main_comm import *
from measurements import *
# from modules.analysis_tools import generate_plots
from MetrolabTHM1176.thm1176 import MetrolabTHM1176Node

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
    def __init__(self, **kwargs):
        threading.Thread.__init__(self)
        keys = kwargs.keys()
        
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
            returnDict = timeResolvedMeasurement(period=self.period, averaging=self.averaging,
                                                block_size=self.block_size, duration=self.duration,
                                                return_temp_data=self.tempData)
        except Exception as e:
            print('There was a problem!')
            print(e)
        # threadLock.release()
        print("Finished measuring. {} exiting.".format(self.name))


def sweepCurrents(node: MetrolabTHM1176Node, config_list=None, config='z', start_val=0, end_val=1, steps=5):
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
    elif config == 'xy':
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
    filePath = '.\data_sets'
    
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
        mean_data, std_data = measure(node, filePath, N=5, average=True)
        mean_values[i] = mean_data
        stdd_values[i] = std_data
        expected_fields[i] = B_expected

    # end of measurements
    disableCurrents()
    # saving data section (prepared for plotting)
    saveDataPoints((all_curr_vals / 1000), mean_values,
                   stdd_values, expected_fields, filePath, fileprefix)


def rampVectorField(node: MetrolabTHM1176Node, theta=90, phi=0, start_mag=0, finish_mag=20, steps=10):
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

     # create subdirectory to save measurements
    fileprefix = '({}_{})_field_meas'.format(int(theta), int(phi))
    # folder, 
    filePath = '.\data_sets'
    
    enableCurrents()
    # with MetrolabTHM1176Node(range="0.3 T", period=0.01) as node:
        # iterate through all possible steps
    for i in range(steps):
        # tentative estimation of resulting B field
        B_expected = tr.computeMagneticFieldVector(magnitudes[i], theta, phi)
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
        mean_data, std_data = measure(node, filePath, N=5, average=True)
        mean_values[i] = mean_data
        stdd_values[i] = std_data
        expected_fields[i] = B_expected

    # end of measurements
    disableCurrents()
    # saving data section (prepared for plotting)
    saveDataPoints((all_curr_vals / 1000), mean_values,
                   stdd_values, expected_fields, filePath, fileprefix)


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
                '[q] to disable currents\n[c]: get currents\n[r]: Set new currents\n[s]: get ECB status\n[f]: get magnetic field measurement\n')
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
            #elif c1 == 'f':
                # params = {'block_size': 5, 'period': 1e-3, 'duration': 10, 'averaging': 1}

                # faden = threading.Thread(timeResolvedMeasurement, **params)
                # faden.start()
                # for key in returnDict.keys:
                #     print(key, 'measured values: ', returnDict[key])
                

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
    
    
    
if __name__ == "__main__":
    params = {'block_size': 20, 'period': 1e-2, 'duration': 10, 'averaging': 3}
    
    faden = myMeasThread(**params)
    faden.start()
    openConnection()
    sleep(1)
    generateMagneticField(60, 90, 0, 4)
    closeConnection()
    faden.join()
    
    item_name = ['Bx', 'By', 'Bz']
    labels = ['Bx', 'By', 'Bz', 'T']
    curve_type = ['F', 'F', 'F', 'T']
    to_show = [True, True, True, False]
    
    plotdata = [returnDict[key] for key in item_name]
    timeline = returnDict['time']
    
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    plt.draw()

    # Setup colors
    NTemp = curve_type.count('T')
    cmap1 = plt.get_cmap('autumn')
    colors1 = [cmap1(i) for i in np.linspace(0, 1, NTemp)]

    NField = curve_type.count('F')
    cmap2 = plt.get_cmap('winter')
    colors2 = [cmap2(i) for i in np.linspace(0, 1, NField)]

    colors = []
    count1 = 0
    count2 = 0
    for ct in curve_type:
        if ct == 'T':
            colors.append(colors1[count1])
            count1 += 1
        else:
            colors.append(colors2[count2])
            count2 += 1

    # Create the matplotlib lines for each curve
    lines = []
    for k, flag in enumerate(to_show):
        if flag:
            data_to_plot = plotdata[k]
            if curve_type[k] == 'F':
                ln, = ax1.plot(timeline, data_to_plot, label=labels[k], color=colors[k])
            else:
                ln, = ax2.plot(timeline, data_to_plot, label=labels[k], color=colors[k])
            lines.append(ln)

    ax1.legend(lines, labels, loc='best')
    plt.show()
