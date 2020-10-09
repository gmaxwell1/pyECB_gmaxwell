"""
filename: main_comm.py

This code is meant to bundle the communication with the ECB-820, by simplifying the most basic necessary functions.
There may be some more functions that are necessary in the future

Author: Maxwell Guerne-Kieferndorf (QZabre)
        gmaxwell@student.ethz.ch

Based on code by Moritz Reinders, Denis von Arx

Date: 09.10.2020
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
# this is ok, since all the function names are pretty much unique and known.
from ECB import *


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


##########  Error messages from ECB, Error messages from ECB, see ECB_API documentation (no exceptions needed): ##########
#           ..\ECB_820_Pantec\Pantec_API\pantec_tftp\ecb_api_20150203\doxygen\html

def _chk(msg):
    """
    Check for errors from the C library, see ECB_API documentation (no exceptions needed):
           ..\ECB_820_Pantec\Pantec_API\pantec_tftp\ecb_api_20150203\doxygen\html    
    """
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


def openConnection(IPAddress=ECB_ADDRESS, port=ECB_PORT):
    """
    Open a connection with the ECB on port ECB_PORT.

    Args:
    -IPAddress: IP address of the ECB. Should always be 192.168.237.47, unless the configuration was changed.
    -port: port on which TCP connection is opened. Should always be 7070, unless the configuration was changed.

    Returns: ECB error code
    """
    ECB_ERR = initECBapi(IPAddress, port)

    if ECB_ERR != 0:
        _chk(ECB_ERR)
        return ECB_ERR

    return ECB_ERR


def closeConnection():
    """
    Close the connection with the ECB.
    """
    exitECBapi()


def enableCurrents():
    """
    Enable current controllers.

    Returns: error code iff an error occurs, otherwise True (whether ECB currents are enabled)
    """
    ECB_ERR = enableECBCurrents()

    if ECB_ERR != 0:
        _chk(ECB_ERR)
        return ECB_ERR

    ECB_CURRENTS_ENABLED = True
    return ECB_CURRENTS_ENABLED


def disableCurrents():
    """
    Disable current controllers.

    Returns: error code iff an error occurs, otherwise False (whether ECB currents are enabled)
    """
    ECB_ERR = disableECBCurrents()

    if ECB_ERR != 0:
        _chk(ECB_ERR)
        return ECB_ERR

    ECB_CURRENTS_ENABLED = False
    return ECB_CURRENTS_ENABLED


def setMaxCurrent(maxValue=19000):
    """
    Set maximum current values for each ECB channel, as long as they are under the threshold specified in the API source code.
    Args:
    -maxValue

    Returns: error code iff an error occurs
    """
    if maxValue < 19800:
        ECB_ERR = setMaxCurrent(maxValue)
        if ECB_ERR != 0:
            _chk(ECB_ERR)
            return ECB_ERR
        ECB_MAX_CURR = maxValue
    else:
        print("specified current was too high")


def setCurrents(desCurrents=[0, 0, 0, 0, 0, 0, 0, 0], direct=b'0'):
    """
    Set current values for each ECB channel

    Args:
    -desCurrents: list of length 8 containing int values, where the '0th' value is the desired current on channel 1 (units of mA),
    the '1st' is the desired current on channel 2 and so on.
    -direct: if 1, the existing ECB buffer will be cleared and desCurrents will be directly applied.
    If 0, the desired currents will be appended to the buffer.

    Returns: error code iff an error occurs
    """
    ECB_ERR = setDesCurrents(desCurrents, direct)

    if ECB_ERR != 0:
        _chk(ECB_ERR)
        return ECB_ERR


def getCurrents():
    """
    Get current values from each ECB channel, print them to the console

    Returns: a list of all the currents or error code
    """
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


def getTemps():
    """
    Get temperature values from each sensor, print them to the console

    returns: a tuple with all of the values of interest if no error occurs otherwise an error code is returned
    """
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


def getStatus():
    """
    Get ECB status (Number)

    returns: status, or error code iff there is an error
    """
    (ECB_ERR, status) = getECBStatus()

    if ECB_ERR != 0:
        _chk(ECB_ERR)
        return ECB_ERR

    return status


########## operate the ECB in the desired mode (test stuff out) ##########
if __name__ == '__main__':
    print(initECBapi(ECB_ADDRESS, ECB_PORT))
    # enableECBCurrents()
    # generateMagField(magnitude=60,theta=0,phi=0)
    #setDesCurrents([2232,2232,2232,0,0,0,0,0], currDirectParam)

    # getCurrents()
    print(getStatus())
    # pollCurrents(100,10)
    sleep(10)
    # disableECBCurrents()
    exitECBapi()
