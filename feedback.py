"""
filename: feedback.py

This script is meant to be used control the magnetic field using feedback from measurement and a simple correction algorithm.

Author: Maxwell Guerne-Kieferndorf (QZabre)
        gmaxwell@student.ethz.ch

Date: 10.11.2020
"""

########## Standard library imports ##########
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from time import time, sleep
from datetime import datetime
import os
from pathlib import Path
import csv
import pandas as pd
import threading

from main_comm import *
from measurements import calibration
from MetrolabTHM1176.thm1176 import MetrolabTHM1176Node
import transformations as tr
import modules.general_functions as gen


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


def print_time(thread: inputThread, delay):
    while flags[0]:
        sleep(delay)
        now = datetime.now()
        if flags[0]:
            print(
                f"\rTime on thread {thread.threadID}: {now: %H:%M:%S}", end='', flush=True)

flags = [1]

def simpFeedback(currentStep=20, maxCorrection=30, thresholdPercentage=0.02, calibrate=False, ECB_active=False):
    """
    Feedback control of vector magnet to obtain any desired magnetic field.

    Args:
        currentStep (int, optional): The maximum step size for derivative (system) estimation. Defaults to 20.
        maxCorrection (int, optional): The maximum step size for current correction. Defaults to 30.
        thresholdPercentage (float, optional): Percentage of total magnitude the error needs to exceed to be considered. Defaults to 0.02.
        calibrate (bool, optional): Tells the user whether or not to calibrate the sensor before measuring. Defaults to False.

    Returns:
        derivatives (np.ndarray(3,3)): The average derivative matrix used to control the current values set.
        printCurrents (np.array(3)): The average of the currents that generate the wanted magnetic field with a small enough error.
    """
    global flags # for handling input
    # currents to set on ECB
    desCurrents = [0] * 8
    # we want to find the 'best' current configuration for a certain direction ultimately.
    # So we save the current values where the error is below the threshold, at the three
    # coil values will be averaged.
    goodCurrentValues = []
    # queue with the last measured magnetic field vector
    fieldVectors = []

    # # future option: enter the field in either coordinate system
    coordinates = input('cartesian or spherical coordinates? ')
    if coordinates == 'c':
        print('Enter the desired magnetic Field vector:')
        magneticField_x = input('B_x = ')
        magneticField_y = input('B_y = ')
        magneticField_z = input('B_z = ')
    elif coordinates == 's':
        print('Enter the desired magnetic Field magnitude and direction:')
        magneticField = input('B = ')
        theta = input('theta = ')
        phi = input('phi = ')
    else:
        print('Enter the desired magnetic Field vector (cartesian):')
        magneticField_x = input('B_x = ')
        magneticField_y = input('B_y = ')
        magneticField_z = input('B_z = ')

    try:
        # desired magnetic (flux density) field used as reference
        if coordinates == 's':
            magneticField = float(magneticField)
            desiredBField = tr.computeMagneticFieldVector(
                magneticField, float(theta), float(phi))
        else:
            desiredBField = np.array([float(magneticField_x), float(
                magneticField_y), float(magneticField_z)])
            magneticField = np.linalg.norm(desiredBField)
    except:
        desiredBField = np.array([20, 0, 0])
        print("There was an issue setting the magnetic field. Using the default value (20mT in x direction).")

    # connect with ECB
    if not ECB_active:
        openConnection()
    # initialize sensor to start measuring
    with MetrolabTHM1176Node(period=0.05, range='0.3 T', average=20) as node:
    # node = MetrolabTHM1176Node(period=0.1, range='0.3 T', average=5)
        calibration(node, calibrate=calibrate)

        configCurrents = tr.computeCoilCurrents(desiredBField)
        # compute linearization Matrix around the current operating point
        derivatives = localLinearization(
            node, configCurrents=configCurrents, currentStep=currentStep, numSteps=10)
        print(derivatives)

        enableCurrents()
        desCurrents[0:3] = configCurrents.tolist()
        setCurrents(desCurrents=desCurrents, direct=b'1')
        goodCurrentValues.append(desCurrents[0:3])
        
        # 'Press Enter to quit' prompt
        test_thread = inputThread(1)
        test_thread.start()
        # Active monitoring & controlling section
        while flags[0]:
            newBMeasurement = gen.sensor_to_magnet_coordinates(
                np.array(node.measureFieldmT()))

            B_magnitude = np.linalg.norm(newBMeasurement)
            theta = np.degrees(np.arccos(newBMeasurement[2]/B_magnitude))
            phi = np.degrees(np.arctan2(
                newBMeasurement[1], newBMeasurement[0]))

            if flags[0]:
                print(f'\rMeasured B field: ({newBMeasurement[0]:.2f}, {newBMeasurement[1]:.2f}, {newBMeasurement[2]:.2f})'
                        f' / In polar coordinates: ({B_magnitude:.2f}, {theta:.2f}°, {phi:.2f}°)    ', sep='', end='', flush=True)

            # currentStep = currentStep * (decayFactor ** i)
            errorsB = newBMeasurement - desiredBField
            fieldDiff = np.zeros(3)
            # if the error of some component is greater than 2%, it needs to be decreased.
            if abs(errorsB[0]) >= abs(thresholdPercentage * magneticField):
                # decrease the B_x component (proportional to the error)
                fieldDiff[0] = -errorsB[0]
            else:
                goodCurrentValues.append(desCurrents[0:3])

            if abs(errorsB[1]) >= abs(thresholdPercentage * magneticField):
                # decrease the B_y component (proportional to the error)
                fieldDiff[1] = -errorsB[1]
            else:
                goodCurrentValues.append(desCurrents[0:3])

            if abs(errorsB[2]) >= abs(thresholdPercentage * magneticField):
                # decrease the B_z component (proportional to the error)
                fieldDiff[2] = -errorsB[2]
            else:
                goodCurrentValues.append(desCurrents[0:3])

            if len(goodCurrentValues) >= 50:
                goodCurrentValues.pop(0)

            deltaCurrent = np.linalg.inv(derivatives).dot(fieldDiff)

            # limit step size
            for i in range(len(deltaCurrent)):
                if abs(deltaCurrent[i]) > maxCorrection:
                    deltaCurrent = np.sign(deltaCurrent) * maxCorrection

            # The set current values are adapted to correct the magnetic field
            desCurrents[0:3] = [int(desCurrents[0]+deltaCurrent[0]), int(desCurrents[1]+deltaCurrent[1]),
                                int(desCurrents[2]+deltaCurrent[2])]
            setCurrents(desCurrents=desCurrents)
            # wait for stabilization
            sleep(0.4)

    # current configurations with little to no error are averaged -> hopefully a good estimate
    # of the required currents for a given direction.
    printCurrents = np.mean(np.array(goodCurrentValues), 0)

    print(
        f'Average currents on channels 1, 2, 3: ({printCurrents[0]:.0f}, {printCurrents[1]:.0f}, {printCurrents[2]:.0f})')
    demagnetizeCoils()
    disableCurrents()
    if not ECB_active:
        closeConnection()
    
    return derivatives, np.array(goodCurrentValues), printCurrents


def localLinearization(node: MetrolabTHM1176Node, configCurrents=np.array([1000, 1000, 1000]), currentStep=20, numSteps=10):
    """
    Compute a local linearization of the relationship between current and B-field using measurements and a fixed step size.

    Args:
        node (MetrolabTHM1176Node): Instance of the Metrolab sensor.
        configCurrents (np.array(3)), optional): Current configuration around which to linearize. Defaults to np.array([1000,1000,1000]).
        currentStep (int, optional): Step size by which current is changed to compute slope. If the step
                                     size is too small, the change in magnetic field strength may be 
                                     inaccurate due to noise. Defaults to 20.
        numSteps (int, optional): The number of measurements over which the result is averaged. Defaults to 10.

    Returns:
        derivatives (3x3 np.ndarray): matrix with 'derivatives' of field components wrt current in each coil.
    """
    desCurrents = [0] * 8
    fieldVectors = []

    currentRefs = np.ndarray((3, 3))
    currentRefs[0, :] = np.array([1, 0, 0])
    currentRefs[1, :] = np.array([0, 1, 0])
    currentRefs[2, :] = np.array([0, 0, 1])

    differenceQuotients = []
    derivatives = np.zeros((3, 3))

    # set currents for the first time
    desCurrents[0:3] = configCurrents.tolist()
    setCurrents(desCurrents=desCurrents, direct=b'1')
    print(
        f'Currents on channels 1, 2, 3: ({configCurrents[0]}, {configCurrents[1]}, {configCurrents[2]})')

    newBMeasurement = gen.sensor_to_magnet_coordinates(
        np.array(node.measureFieldmT()))
    fieldVectors.insert(0, newBMeasurement)

    for i in range(3):
        for j in range(numSteps):
            desCurrents[0:3] = [int(desCurrents[0] + currentStep * currentRefs[i, 0]),
                                int(desCurrents[1] +
                                    currentStep * currentRefs[i, 1]),
                                int(desCurrents[2] + currentStep * currentRefs[i, 2])]
            setCurrents(desCurrents=desCurrents)

            newBMeasurement = gen.sensor_to_magnet_coordinates(
                np.array(node.measureFieldmT()))
            print(
                f'\rMeasured B field: ({newBMeasurement[0]:.2f}, {newBMeasurement[1]:.2f}, {newBMeasurement[2]:.2f})', sep='', end='', flush=True)
            fieldVectors.insert(0, newBMeasurement)
            differenceQuotients.insert(
                0, (newBMeasurement - fieldVectors[1]) / currentStep)
            sleep(0.2)
            fieldVectors.pop()

        desCurrents[0:3] = configCurrents.tolist()
        setCurrents(desCurrents=desCurrents, direct=b'1')

        for j in range(numSteps):
            a = -currentStep  # constant current decrease
            desCurrents[0:3] = [int(desCurrents[0] + a * currentRefs[i, 0]),
                                int(desCurrents[1] + a * currentRefs[i, 1]),
                                int(desCurrents[2] + a * currentRefs[i, 2])]
            setCurrents(desCurrents=desCurrents)

            newBMeasurement = gen.sensor_to_magnet_coordinates(
                np.array(node.measureFieldmT()))
            print(
                f'\rMeasured B field: ({newBMeasurement[0]:.2f}, {newBMeasurement[1]:.2f}, {newBMeasurement[2]:.2f})', sep='', end='', flush=True)
            fieldVectors.insert(0, newBMeasurement)
            differenceQuotients.insert(
                0, (newBMeasurement - fieldVectors[1]) / a)
            sleep(0.2)
            fieldVectors.pop()

        print('')
        # save the latest row of the 'derivatives' in the matrix
        derivatives[:, i] = np.mean(np.array(differenceQuotients), 0)
        differenceQuotients.clear()

    return derivatives


########## test stuff out ##########
if __name__ == '__main__':
    _ = simpFeedback(maxCorrection=20, calibrate=True)

    # simpFeedback()
