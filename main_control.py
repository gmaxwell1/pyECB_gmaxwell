"""
filename: main_control.py

This code is meant to be used to set a specific desired current on the channels
1, 2 and 3 (and up to 8, but this is not really needed in our case) of the ECB-820, 
either by expilcitly setting current values, or by specifying the desired magnetic field vector.

Author: Maxwell Guerne-Kieferndorf (QZabre)
        gmaxwell@student.ethz.ch

Date: 28.09.2020
"""
# Standard library imports
import numpy as np
import sys
import math
import time
# local imports
import ECB

# error messages from ECB
def _chk(msg):
    """Check for errors from the C library."""
    if msg:
        if msg == 0:
            print("Succesful")
        if msg == 1:
            raise RuntimeError("ECB uninitialized")
        if msg == 2:
            raise RuntimeError("No valid ECB socket")
        if msg == 3:
            raise RuntimeError(" ECB current on communication error")
        if msg == 4:
            raise RuntimeError(" ECB current off communication error")
        if msg == 5:
            raise RuntimeError(" ECB revision communication error")
        if msg == 6:
            raise RuntimeError(" ECB revision value error")
        if msg == 20:
            raise RuntimeError("ECB set current communication error")
        if msg == 21:
            raise RuntimeError("ECB get current communication error")
        if msg == 22:
            raise RuntimeError("ECB set desired current, maximum current exceeded")
        if msg == 23:
            raise RuntimeError("ECB set maximum current, value too high")
        if msg == 24:
            raise RuntimeError("ECB set maximal current communication error")
        if msg == 26:
            raise RuntimeError("ECB get drive temperature communication error")
        if msg == 27:
            raise RuntimeError("ECB cannot take new current values, since not enabled")
        if msg == 30:
            raise RuntimeError("ECB get coil status communication error")
        if msg == 31:
            raise RuntimeError("ECB coil has event")
        if msg == 40:
            raise RuntimeError("ECB set maximal temperature communication error")
        if msg == 50:
            raise RuntimeError("ECB set temperature sensor mask communication error")
        if msg == 60:
            raise RuntimeError("ECB factory reset communication error")
        if msg == 70:
            raise RuntimeError("ECB set hearbeat timeout communication error")
        if msg == 71:
            raise RuntimeError("ECB get hearbeat timeout communication error")
        raise RuntimeError(
            "Unhandled error number: {}. \
            See DCx_User_and_SDK_Manual.pdf for details".format(msg))

"""
Class that creates an ECB object and allows communication with the Pantec ECB-820 through
this object.
"""
class ECB820:
    # connection parameters
    ECB_ADDRESS = "192.168.237.47"
    ECB_PORT = "7070"

    # ECB state values
    ecb_state = {}
    ecb_state["ECB_CONNECTED"] = False
    ecb_state["ECB_INIT"] = False
    ecb_state["ECB_ERR"] = 0
    ecb_state["ECB_MAX_CURR"] = 19800  # Max current on any coil in [mA]
    ecb_state["ECB_ACT_CURRENTS"] = np.zeros(8, dtype=int)
    ecb_state["ECB_HEARTBEAT_TO"] = 0
    ecb_state["ECB_MAX_TEMP"] = 50
    ecb_state["ECB_TEMP_MASK"] = 7  # 7: 0x07 = 00000111 , that is, 3 sensors enabled
    ecb_state["ECB_CURRENTS_ENABLED"] = False



    # Current parameters
    desCurrents = np.zeros(8, dtype=int)  # in milliamps
    currDirectParam = b'0'
    # overCurrent = 0

    # Vector magnet properties
    # windings = 508  # windings per coil
    # resistance = 0.45  # resistance per coil

    # Store active magnetic field
    # magnetic_field_unit = np.zeros((3, 1))  # in mT. Unit vector
    # magnetic_field = np.zeros((3, 1))  # in mT. Unit vector scaled with magnitude

    # elapsed_time = 0


    # Initialization of ECB with predefined values. If anything goes wrong, the initialization will be canceled.
    def __init__(self):
        print("Initializing ECB")
        self.ecb_state["ECB_ERR"] = ECB.initECBapi(self.ECB_ADDRESS, self.ECB_PORT)
        if self.ecb_state["ECB_ERR"] != 0:
            print("Initialization failed")
            _chk(self.ecb_state["ECB_ERR"])
            return

        self.ecb_state["ECB_CONNECTED"] = True

        self.ecb_state["ECB_ERR"] = ECB.setHeartbeatTimeout(self.ecb_state["ECB_HEARTBEAT_TO"])
        if self.ecb_state["ECB_ERR"] != 0:
            print("Error in setting the heartbeat timeout")
            ECB.exitECBapi()
            self.ecb_state["ECB_CONNECTED"] = False
            _chk(self.ecb_state["ECB_ERR"])
            return

        self.ecb_state["ECB_ERR"] = ECB.setMaxTemp(self.ecb_state["ECB_MAX_TEMP"])
        if self.ecb_state["ECB_ERR"] != 0:
            print("Error in setting the maximum temperature")
            ECB.exitECBapi()
            self.ecb_state["ECB_CONNECTED"] = False
            _chk(self.ecb_state["ECB_ERR"])
            return

        self.ecb_state["ECB_ERR"] = ECB.setTempMask(self.ecb_state["ECB_TEMP_MASK"])
        if self.ecb_state["ECB_ERR"] != 0:
            print("Error in setting the temperature sensor mask")
            ECB.exitECBapi()
            self.ecb_state["ECB_CONNECTED"] = False
            _chk(self.ecb_state["ECB_ERR"])
            return

        self.ecb_state["ECB_INIT"] = True

    def checkECBState(self):
        for key in self.ecb_state:
            print(key, ": ", self.ecb_state[key])
        print("If any parameters have strange values, reset them or reset the ECB.")
        result = ECB.getCoilStatus()
        print(result)

    # by default, no currDirect
    def setCurrDirect(self, param):
        self.currDirectParam = param

    # set the current values on any number (between 1 and 8) of channels to the values specified
    # in channels (a list), or leave 0 if nothing is specified.
    def setCurrents(self, channels):

        for i in range(len(channels)):
            self.desCurrents[i] = channels[i]

        if np.amax(self.desCurrents) > self.ecb_state["ECB_MAX_CURR"]:
            print("Error: specified current is too high!")
            return
        print(self.desCurrents.tolist())
        self.ecb_state["ECB_ERR"] = ECB.setDesCurrents(self.desCurrents.tolist(), self.currDirectParam)
        if self.ecb_state["ECB_ERR"] != 0:
            print("Error in setting the desired currents")
            _chk(self.ecb_state["ECB_ERR"])
            return

        self.ecb_state["ECB_ERR"] = ECB.enableECBCurrents()
        if self.ecb_state["ECB_ERR"] != 0:
            print("Error enabling currents")
            _chk(self.ecb_state["ECB_ERR"])
            return

        self.ecb_state["ECB_CURRENTS_ENABLED"] = True
        return

    def disableCurrents(self):
        self.ecb_state["ECB_ERR"] = ECB.disableECBCurrents()
        if self.ecb_state["ECB_ERR"] != 0:
            print("Error disabling currents")
            _chk(self.ecb_state["ECB_ERR"])
        
        self.ecb_state["ECB_CURRENTS_ENABLED"] = False
        return


    # get current values on first 3 channels from ECB
    def getCurrents(self):
        result = ECB.getActCurrents()
        self.ecb_state["ECB_ACT_CURRENTS"] = result[1]
        print("Channel: 1 \t 2 \t 3 ")
        print("\t {0} \t {1} \t {2} ".format(self.ecb_state["ECB_ACT_CURRENTS"][0], self.ecb_state["ECB_ACT_CURRENTS"][1], self.ecb_state["ECB_ACT_CURRENTS"][2]))


    # close connection with ECB
    def turnOffECB(self):
        self.desCurrents = np.zeros(8, dtype=int)
        self.setCurrDirect(b'0')
        self.ecb_state["ECB_ERR"] = ECB.setDesCurrents(self.desCurrents.tolist(), self.currDirectParam)
        if self.ecb_state["ECB_ERR"] != 0:
            print("Error in setting the desired currents")
            _chk(self.ecb_state["ECB_ERR"])
        ECB.exitECBapi()
        self.ecb_state["ECB_CONNECTED"] = False
        return

"""def toFlatList(inList):
    result = []
    if not isinstance(inList, list):
        raise RuntimeError("Error: invalid input for function toFlatList")
    for i in inList:
        if isinstance(i,inList):
            toFlatList(i)
        else:
            result.append(i)
    return result
"""

def main():
    ECB1 = ECB820()
    ECB1.checkECBState()
    ECB1.setCurrDirect(b'1')
    ECB1.setCurrents([3000,0,1000])
    period = 15
    polls = 5
    dt = period/polls
    while polls > 0:
        polls = polls - 1
        ECB1.getCurrents()
        time.sleep(dt)

    ECB1.disableCurrents()
    ECB1.checkECBState()
    time.sleep(5)

    ECB1.turnOffECB()


if __name__ == '__main__':
    main()

