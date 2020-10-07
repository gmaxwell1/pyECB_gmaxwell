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
from time import time, sleep
import sys
from os import getcwd, path
from pathlib import Path
import csv

########## local imports ##########
from ECB import * # this is ok, since all the function names are pretty much unique and known.


__all__ = [
        'openConnection',
        'closeConnection',
        'enableCurrents',
        'disableCurrents',
        'setMaxCurrent',
        'setCurrents',
        'getCurrents',
        'getTemps',
        'getStatus',
        'ECB_MAX_CURR'
        ]

##########  Connection parameters ##########

ECB_ADDRESS = "192.168.237.47"
ECB_PORT = "7070"

##########  ECB state values (mostly unnecessary) ##########

# ECB_CONNECTED = False
ECB_ERR = 0
# Max current on any coil in [mA] (default value)
ECB_MAX_CURR = 19800
# ECB_ACT_CURRENTS = np.zeros(8, dtype=int)

# ECB_MAX_TEMP = 50
##### 7: 0x07 = 00000111 , that is, 3 sensors enabled #####
# ECB_TEMP_MASK = 7
ECB_CURRENTS_ENABLED = False
# overCurrent = 0


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

def openConnection(IPAddress=ECB_ADDRESS, port=ECB_PORT):
    ECB_ERR = initECBapi(IPAddress, port)

    if ECB_ERR != 0:
        _chk(ECB_ERR)
        return ECB_ERR

    return ECB_ERR


##########  Set maximum current value on each ECB channel ##########

def closeConnection():
    exitECBapi()


##########  Enable current controllers ##########

def enableCurrents():
    ECB_ERR = enableECBCurrents()

    if ECB_ERR != 0:
        _chk(ECB_ERR)
        return ECB_ERR

    ECB_CURRENTS_ENABLED = True
    return ECB_CURRENTS_ENABLED


##########  Disable current controllers ##########

def disableCurrents():
    ECB_ERR = disableECBCurrents()

    if ECB_ERR != 0:
        _chk(ECB_ERR)
        return ECB_ERR

    ECB_CURRENTS_ENABLED = False
    return ECB_CURRENTS_ENABLED


##########  Set maximum current value on each ECB channel ##########

def setMaxCurrent(maxValue=19000):

    if maxValue < 19800:
        ECB_ERR = setMaxCurrent(maxValue)
        if ECB_ERR != 0:
            _chk(ECB_ERR)
            return ECB_ERR
        ECB_MAX_CURR = maxValue
    else:
        print("specified current was too high")


##########  Get current values from each ECB channel, print them to the console ##########
#
#           returns: error code if an error occurs

def setCurrents(desCurrents=[0,0,0,0,0,0,0,0], direct=b'0'):
    ECB_ERR = setDesCurrents(desCurrents, direct)
    
    if ECB_ERR != 0:
        _chk(ECB_ERR)
        return ECB_ERR


##########  Get current values from each ECB channel, print them to the console ##########
#
#           returns: a list of all the currents or error code

def getCurrents():

    (ECB_ERR, result) = getActCurrents()

    if ECB_ERR != 0:
        _chk(ECB_ERR)
        return ECB_ERR
    else:
        print("Channel: \t 1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8")
        print("Current [mA]: {0} \t {1} \t {2} \t {3} \t {4} \t {5} \t {6} \t {7}"
            .format(result[0], result[1], result[2], result[3], result[4],
                    result[5], result[6], result[7]))

    return result


##########  Get temperature values from each sensor, print them to the console ##########
#
#           returns: a tuple with all of the values of interest if no error occurs
#                    otherwise an error code is returned

def getTemps():
    (ECB_ERR, result, hall_list, currents_list, coil_status) = getCoilValues()
    
    if ECB_ERR != 0:
        _chk(ECB_ERR)
        return ECB_ERR
    else:
        print("Channel: 1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8")
        print("Temperature [°C]: {0} \t {1} \t {2} \t {3} \t {4} \t {5} \t {6} \t {7}"
            .format(result[0], result[1], result[2], result[3], result[4],
                    result[5], result[6], result[7]))

    return (result, hall_list, currents_list, coil_status)


##########  Get temperature values from each sensor, print them to the console ##########
#
#           returns: a tuple with all of the values of interest if there is no error
#                    otherwise an error code is returned

def getStatus():
    (ECB_ERR, status) = getECBStatus()
    
    if ECB_ERR != 0:
        _chk(ECB_ERR)
        return ECB_ERR

    return status


########## operate the ECB in the desired mode (test stuff out) ##########
if __name__ == '__main__':
    initECBapi(ECB_ADDRESS, ECB_PORT)
    #enableECBCurrents()
    #generateMagField(magnitude=60,theta=0,phi=0)
    #setDesCurrents([2232,2232,2232,0,0,0,0,0], currDirectParam)

    #getCurrents()

    #pollCurrents(100,10)
    #sleep(10)
    #disableECBCurrents()
    exitECBapi()
