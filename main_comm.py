"""
filename: main_comm.py

This code is meant to be used to set a specific desired current on the channels
1, 2 and 3 (and up to 8, but this is not really needed in our case) of the ECB-820, 
either by expilcitly setting current values, or by specifying the desired magnetic field vector.

Author: Maxwell Guerne-Kieferndorf (QZabre)
        gmaxwell@student.ethz.ch

Date: 28.09.2020
"""
########## Standard library imports ##########
import numpy as np
import math
from time import time,sleep
from datetime import datetime
import sys
from os import getcwd, path
from pathlib import Path
import csv

########## local imports ##########
from ECB import * # this is ok, since all the function names are pretty much unique and known.
import transformations as tr


##########  Connection parameters ##########

ECB_ADDRESS = "192.168.237.47"
ECB_PORT = "7070"

##########  ECB state values (mostly unnecessary) ##########

# ECB_CONNECTED = False
# ECB_INIT = False
ECB_ERR = 0
# Max current on any coil in [mA] (default value)
ECB_MAX_CURR = 19800
# ECB_ACT_CURRENTS = np.zeros(8, dtype=int)

# ECB_MAX_TEMP = 50
##### 7: 0x07 = 00000111 , that is, 3 sensors enabled #####
# ECB_TEMP_MASK = 7
# ECB_CURRENTS_ENABLED = False
# overCurrent = 0

##########  Current parameters ##########

desCurrents = [0, 0, 0, 0, 0, 0, 0, 0]  # in milliamps
currDirectParam = b'1'

##########  Vector magnet properties ##########

windings = 508  # windings per coil
resistance = 0.45  # resistance per coil



##########  Error messages from ECB, see ECB_API documentation (no exceptions needed): ##########
#           ..\ECB_820_Pantec\Pantec_API\pantec_tftp\ecb_api_20150203\doxygen\html     

def _chk(msg):
    """Check for errors from the C library."""
    if msg:
        if msg == 0:
            print("Succesful")
        if msg == 1:
            print("ECB uninitialized")
        if msg == 2:
            print("No valid ECB socket")
        if msg == 3:
            print(" ECB current on communication error")
        if msg == 4:
            print(" ECB current off communication error")
        if msg == 5:
            print(" ECB revision communication error")
        if msg == 6:
            print(" ECB revision value error")
        if msg == 20:
            print("ECB set current communication error")
        if msg == 21:
            print("ECB get current communication error")
        if msg == 22:
            print("ECB set desired current, maximum current exceeded")
        if msg == 23:
            print("ECB set maximum current, value too high")
        if msg == 24:
            print("ECB set maximal current communication error")
        if msg == 26:
            print("ECB get drive temperature communication error")
        if msg == 27:
            print("ECB cannot take new current values, since not enabled")
        if msg == 30:
            print("ECB get coil status communication error")
        if msg == 31:
            print("ECB coil has event")
        if msg == 40:
            print("ECB set maximal temperature communication error")
        if msg == 50:
            print("ECB set temperature sensor mask communication error")
        if msg == 60:
            print("ECB factory reset communication error")
        if msg == 70:
            print("ECB set hearbeat timeout communication error")
        if msg == 71:
            print("ECB get hearbeat timeout communication error")
        print(
            "Unhandled error number: {}. \
            See DCx_User_and_SDK_Manual.pdf for details".format(msg))


##########  Set maximum current value on each ECB channel ##########

def setMaxCurrent(maxValue):

    if maxValue < 19800
        ECB_ERR = setMaxCurrent(maxValue)
    else:
        print("specified current was too high")
    
    if ECB_ERR != 0:
        _chk(ECB_ERR)
        return

    return



##########  Get current values from each ECB channel, print them to the console ##########
#
#           returns: a list of all the currents

def getCurrents():
    (ECB_ERR, result) = getActCurrents()
    
    if ECB_ERR != 0:
        _chk(ECB_ERR)
        return
    else:
        print("Channel: 1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8")
        print("Current [mA]: {0} \t {1} \t {2} \t {3} \t {4} \t {5} \t {6} \t {7}"
            .format(result[0], result[1], result[2], result[3], result[4],
                    result[5], result[6], result[7]))

    return result



##########  gets current values writes the values into a .csv file named according ##########
#           to today's date/time
#           Args:
#           -measurements: an integer specifying how many values are wanted
#           -frequency: how many times per second (approximately!!, not reliable) to get the 
#                       values

def pollCurrents(measurements, frequency=1):
    dt = 1/frequency
    k = 0

    # format the output file/folder
    now = datetime.now().strftime('%H%M%S')
    today = datetime.now().strftime('%Y-%m-%d')
    filename = '{}_{}_{}'.format(today, now, 'meas.csv')
    cwd = getcwd()
    workDir = path.join(cwd, today)
    # create new directory (for current day) if necessary
    Path(workDir).mkdir(parents=True,exist_ok=True)

    with open(path.join(workDir, filename), 'x', newline='') as newFile:
        meas_writer = csv.writer(newFile)
        meas_writer.writerow(['Measurement time', 'Ch1 [mA]', 'Ch2 [mA]', 'Ch3 [mA]',
                              'Ch4 [mA]', 'Ch5 [mA]', 'Ch6 [mA]', 'Ch7 [mA]', 'Ch8 [mA]'])
        line = []

        startTime = time()

        while k < measurements:
            elapsedTime = time() - startTime
            k = k + 1
            (ECB_ERR, result) = getActCurrents()

            if ECB_ERR != 0:
                _chk(ECB_ERR)
                return
            else:
                line = result

            line.insert(0, elapsedTime)
            meas_writer.writerow(line)

            loop_time = time() - startTime - elapsedTime
            # print(loop_time)
            # not definitive, just to see how significant of an effect the execution
            # time has on the sampling period, especially for high frequencies
            sleep(dt)
    
    ECB_ERR = 0



##########  test method to write a .csv file ##########

def textCSV(str):

    # format the output file/folder
    now = datetime.now().strftime('%H%M%S')
    today = datetime.now().strftime('%Y-%m-%d')
    filename = '{}_{}_{}'.format(today, now, 'meas.csv')
    cwd = getcwd()
    workDir = path.join(cwd, today)
    # create new directory (for current day) if necessary
    Path(workDir).mkdir(parents=True,exist_ok=True)

    with open(path.join(workDir, filename), 'x', newline='') as newFile:
        _writer = csv.writer(newFile)
        _writer.writerow(['Hello, ' + str + '!'])

        elapsedTime = 0
        startTime = time()

        while elapsedTime < 1:
            elapsedTime = time() - startTime

            #loop_time = time() - startTime - elapsedTime

            _writer.writerow([elapsedTime])
            # not definitive, just to see how significant of an effect the execution
            # time has on the sampling period, especially for high frequencies
            sleep(0.1)



##########  generate a magnetic field in an arbitrary direction and ##########
#           an arbitrary magnitude     
#           Args:
#           -magnitude: of the B-field, in [mT]
#           -theta: angle between desired field direction and z axis
#           -phi: azimuthal angle (measured from the x axis)
#           -t: time for which the magnetic field should be activated (if not 0)
#           -direct: current direct parameter (can usually be left alone)

def generateMagField(magnitude, theta, phi, t=0, direct=b'1'):
    B_vector = tr.computeMagneticFieldVector(magnitude, theta, phi)
    print(B_vector)
    I_vector = tr.computeCoilCurrents(B_vector, windings, resistance)
    # make sure that the current is not too high
    if np.amax(I_vector) > ECB_MAX_CURR:
        ECB_ERR = 22
        _chk(ECB_ERR)
        ECB_ERR = 0
        return

    currDirectParam = direct
    # copy the computed current values (mA) into the desCurrents list (first 3 positions)
    # cast to int
    for i in range(len(I_vector)):
        desCurrents[i] = int(I_vector[i])
    print(desCurrents)
    # user specified on time
    if t > 0:
        ECB_ERR = enableECBCurrents()
        if ECB_ERR != 0:
            _chk(ECB_ERR)

        ECB_ERR = setDesCurrents(desCurrents, currDirectParam)
        if ECB_ERR != 0:
            _chk(ECB_ERR)
            return

        sleep(t)

        ECB_ERR = disableECBCurrents()
        if ECB_ERR != 0:
            _chk(ECB_ERR)
    # on until interrupted by user
    elif t == 0:
        ECB_ERR = enableECBCurrents()
        if ECB_ERR != 0:
            _chk(ECB_ERR)

        ECB_ERR = setDesCurrents(desCurrents, currDirectParam)
        if ECB_ERR != 0:
            _chk(ECB_ERR)
            return

        # wait until user presses enter
        input('Press Enter to disable.')

        ECB_ERR = disableECBCurrents()
        if ECB_ERR != 0:
            _chk(ECB_ERR)
    else:
        return



########## operate the ECB in the desired mode (test stuff out) ##########

if __name__ == '__main__':
    initECBapi(ECB_ADDRESS, ECB_PORT)
    #enableECBCurrents()
    generateMagField(magnitude=30,theta=0,phi=0)
    #setDesCurrents([2232,2232,2232,0,0,0,0,0], currDirectParam)

    #getCurrents()

    #pollCurrents(100,10)
    #sleep(10)
    #disableECBCurrents()
    exitECBapi()
