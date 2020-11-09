"""
filename: feedback.py

This script is meant to be used control the magnetic field using feedback from measurement and a simple PD algorithm.

Author: Maxwell Guerne-Kieferndorf (QZabre)
        gmaxwell@student.ethz.ch

Date: 02.11.2020
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
from MetrolabTHM1176.thm1176 import MetrolabTHM1176Node
import transformations as tr
import modules.general_functions as gen


lockObj = threading.Lock()

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
        # Get lock to synchronize threads
        c = input("Hit Enter to quit.\n")
        flags.insert(0, c)
        print('exiting...')


def print_time(thread: inputThread, delay):
    while flags[0]:
        sleep(delay)
        now = datetime.now()
        if flags[0]:
            print(f"\rTime on thread {thread.threadID}: {now: %H:%M:%S}", end='', flush=True)


def simpFeedback(decayFactor = 0.9, currentStep = 30):
    # manual feedback control
    # queue with the last (10) set current changes in it, make sure the length is always (10)
    currentsDiffs = []
    # currents to set on ECB
    desCurrents = [0] * 8
    # queue with the last (10) measured magnetic field vectors
    fieldVectors = []
    fieldVectorsDiffs = []
    # computed with matlab
    feedbackMatrix = np.ndarray((3,3))
    feedbackMatrix[0,:] = np.array([34.82,-8.18,25.74])
    feedbackMatrix[1,:] = np.array([-29.17,39.79,-9.73])
    feedbackMatrix[2,:] = np.array([21.76,-0.22,-7.68])
    feedbackMatrix_inv = np.linalg.inv(feedbackMatrix) * 1000 # *1000 to get mA
    
    # currentRefs = np.ndarray((19,3))
    # currentRefs[0,:] = np.array([10,0,0])
    # currentRefs[1,:] = np.array([0,10,0])
    # currentRefs[2,:] = np.array([0,0,10])
    # currentRefs[3,:] = np.array([-10,0,0])
    # currentRefs[4,:] = np.array([0,-10,0])
    # currentRefs[5,:] = np.array([0,0,-10])
    # currentRefs[6,:] = np.array([100,0,0])
    # currentRefs[7,:] = np.array([0,100,0])
    # currentRefs[8,:] = np.array([0,0,100])
    # currentRefs[9,:] = np.array([50,0,50])
    # currentRefs[10,:] = np.array([50,0,-50])
    # currentRefs[11,:] = np.array([50,50,0])
    # currentRefs[12,:] = np.array([50,-50,0])
    # currentRefs[13,:] = np.array([-50,50,0])
    # currentRefs[14,:] = np.array(-50,0,50])
    # currentRefs[15,:] = np.array([0,-50,50])
    # currentRefs[16,:] = np.array([0,50,50])
    # currentRefs[17,:] = np.array([0,50,-50])
    # currentRefs[18,:] = np.array([50,50,50])    
    
    # # future option: enter the field in either coordinate system
    # coordinates = input('cartesian or spherical coordinates? ')
    # if coordinates == 'c':
    print('Enter the desired magnetic Field vector:')
    magneticField_x = input('B_x = ')
    magneticField_y = input('B_y = ')
    magneticField_z = input('B_z = ')
    # elif coordinates == 's':
    #     print('Enter the desired magnetic Field magnitude and direction:')
    #     magneticField = input('B = ')
    #     theta = input('theta = ')
    #     phi = input('phi = ')
    # else:
    #     print('Enter the desired magnetic Field vector (cartesian):')
    #     magneticField_x = input('B_x = ')
    #     magneticField_y = input('B_y = ')
    #     magneticField_z = input('B_z = ')
    
    try:
        # desired magnetic (flux density) field used as reference
        desiredBField = np.array([float(magneticField_x),float(magneticField_y),float(magneticField_z)])
    except:
        desiredBField = np.array([20,0,0])
        print("There was an issue setting the magnetic field. Using the default value (20mT in x direction).")
        
    openConnection()
    enableCurrents()    
    # initialize sensor to start measuring
    with MetrolabTHM1176Node(period=0.05, range='0.3 T', average=20) as node:
    # node = MetrolabTHM1176Node(period=0.1, range='0.3 T', average=20)      
        deltaCurrent = np.zeros(3)
        errorsB = np.zeros(3)
        # fieldVectors.insert(0,np.zeros(3))

        # for i in range(19):
        #     # deltaCurrent = currentRefs[i,:]
        #     # desCurrents[0:3] = [int(desCurrents[0]+deltaCurrent[0]), int(desCurrents[1]+deltaCurrent[1]), 
        #     #                     int(desCurrents[2]+deltaCurrent[2])]
        #     # setCurrents(desCurrents=desCurrents)
        currentsDiffs.insert(0,deltaCurrent)
        #     newBMeasurement = gen.sensor_to_magnet_coordinates(np.array(node.measureFieldmT()))
        #     fieldVectorsDiffs.insert(0,newBMeasurement-fieldVectors[0])
        #     fieldVectors.insert(0,newBMeasurement)
        #     fieldVectors.pop()

        # set currents for the first time
        deltaCurrent = tr.computeCoilCurrents(desiredBField)
        desCurrents[0:3] = deltaCurrent.tolist()
        setCurrents(desCurrents=desCurrents, direct=b'1')
        currentsDiffs.insert(0, deltaCurrent)
        print(f'Currents on channels 1, 2, 3: ({desCurrents[0]}, {desCurrents[1]}, {desCurrents[2]})')
        
        # fieldVectorsDiffs.insert(0, newBMeasurement-fieldVectors[0])
        # fieldVectors.insert(0,newBMeasurement)
        # fieldVectors.pop()
        # Press Enter to quit prompt
        test_thread = inputThread(1)
        test_thread.start()
        i = 0
        while flags[0]:
            newBMeasurement = gen.sensor_to_magnet_coordinates(np.array(node.measureFieldmT()))
            
            if flags[0]:
                print(f'\rMeasured B field: ({newBMeasurement[0]:.2f}, {newBMeasurement[1]:.2f}, {newBMeasurement[2]:.2f})',sep='' ,end='', flush=True)

            # currentStep = currentStep * (decayFactor ** i)
            errorsB = desiredBField - newBMeasurement
            # if the error of some component is greater than 1%, it needs to be decreased.
            B_mag = np.linalg.norm(desiredBField)
            if errorsB[0] < -0.02 * B_mag:
                # decrease the B_x component (proportional to the error)
                deltaCurrent = np.array([-currentStep,currentStep,currentStep])
            elif errorsB[0] >= 0.02 * B_mag:
                # increase the B_x component (proportional to the error)
                deltaCurrent = np.array([currentStep,-currentStep,-currentStep])
            if errorsB[1] < -0.02 * B_mag:
                # decrease the B_y component (proportional to the error)
                deltaCurrent += np.array([currentStep,-currentStep,currentStep])
            elif errorsB[1] >= 0.02 * B_mag:
                # increase the B_y component (proportional to the error)
                deltaCurrent += np.array([-currentStep,currentStep,-currentStep])
            if errorsB[2] < -0.02 * B_mag:
                # decrease the B_z component (proportional to the error)
                deltaCurrent += np.array([-currentStep,0,-currentStep])
            elif errorsB[2] >= 0.02 * B_mag:
                # increase the B_z component (proportional to the error)
                deltaCurrent += np.array([currentStep,0,currentStep])

            # print(f'\rMeasured B field: ({deltaCurrent[0]}, {deltaCurrent[1]}, {deltaCurrent[2]})', end='', flush=True)
            currentsDiffs.insert(0, deltaCurrent)
            currentsDiffs.pop()

            desCurrents[0:3] = [int(desCurrents[0]+deltaCurrent[0]), int(desCurrents[1]+deltaCurrent[1]), 
                                int(desCurrents[2]+deltaCurrent[2])]
            setCurrents(desCurrents=desCurrents)
            # wait for stabilization
            sleep(0.4)
            
            # fieldVectorsDiffs.insert(0, newBMeasurement-fieldVectors[0])
            # fieldVectors.insert(0, newBMeasurement)
            
            # fieldVectorsDiffs.pop()
            # fieldVectors.pop()
            
            i += 1
            
    print(currentsDiffs)
    print(errorsB)
                    
    demagnetizeCoils()
    disableCurrents()
    closeConnection()
            

########## test stuff out ##########
if __name__ == '__main__':
    flags = [1]
    simpFeedback()
    
    # simpFeedback()