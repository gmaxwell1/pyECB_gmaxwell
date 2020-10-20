#! /usr/bin/env python

# Tesla - A ROS-based framework for performing magnetic manipulation
#
# Copyright 2018 Multi Scale Robotics Lab
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import usbtmc
import numpy as np
import logging
from time import sleep

class MetrolabTHM1176Node(object):

    def __init__(self):

        logging.basicConfig(filename='metrolab.log', level=logging.DEBUG)

        self.sensor = usbtmc.Instrument(0x1bfa, 0x0498)
        ret = self.sensor.ask("*IDN?")
        # if ret != "Metrolab Technology SA,THM1176-MF":
        #         raise(IOError, "Error opening Metrolab device")

        self.average_count = 10
        self.unit = "MT"
        self.sense_range_upper = "0.1 T"

        # Write settings to device
        #self.sensor.write(":format:data default")
        #self.sensor.write(":average:count " + str(self.average_count))
        self.sensor.write(":unit " + self.unit)
        self.sensor.write(":sense:range:upper 0.1T")
        self.sensor.write(":sense:flux:range:auto off")
       
        logging.debug('range upper %s', self.sensor.ask(":sense:range:upper?"))

        self.posx = 0.0
        self.posy = 0.0
        self.posz = 0.0
        self.frequency = 5.0

        logging.info('... End init')

    def calibrate(self):
        self.sensor.write(":calibration:initiate")
        self.sensor.write(":calibration:state on")
    
    def setAveragingCount(self, num_meas):
        avg_max = int(self.sensor.ask(":average:count? maximum"))
        
        if num_meas <= avg_max and num_meas > 0:
            self.sensor.write(":average:count " + str(num_meas))
            return True
        else:
            print("MetrolabTHM1176:setAveringCount: value has to be between 1 and " + str(avg_max))
            return False
    
    def getNumMeasurements(self):
        ret = self.sensor.ask(':CALC:AVER:COUNT?')
        return int(ret)
            
    def measureFieldmT(self):
        Bx = float(self.sensor.ask(':measure:scalar:flux:x? 0.01T,5').strip('MT'))
        By = float(self.sensor.ask(":measure:y? 0.01T,5").strip('MT'))
        Bz = float(self.sensor.ask(":measure:z? 0.01T,5").strip('MT'))
        
        return [Bx, By, Bz]
        
    def measureFieldArraymT(self, num_meas):
        #self.setAveragingCount(num_meas)
        ret = self.sensor.ask(":measure:array:x? " + str(num_meas) + ", 0.01T,5")
        Bx_str = ret.split(",")
        Bx = []
        for val in Bx_str:
            Bx.append(float(val.strip('MT')))
        
        ret = self.sensor.ask(":measure:array:y? " + str(num_meas) + ", 0.01T,5")
        By_str = ret.split(",")
        By = []
        for val in By_str:
            By.append(float(val.strip('MT')))
            
        ret = self.sensor.ask(":measure:array:z? " + str(num_meas) + ", 0.01T,5")
        Bz_str = ret.split(",")
        Bz = []
        for val in Bz_str:
            Bz.append(float(val.strip('MT')))
        
        if (len(Bx) != num_meas or len(By) != num_meas or len(Bz) != num_meas):
            raise ValueError("length of Bx, By, Bz do not match num_meas")
        
        return [Bx, By, Bz]
        
    def getAvailableUnits(self):
        unit_str = self.sensor.ask(":UNIT:ALL?")
        return unit_str.split(',')
    
    def getUnit(self):
        units_str = self.sensor.ask(":UNIT?")
        return units_str
    
    def getAvailableSenseRangeUpper(self):
        upper_str = self.sensor.ask(':SENS:ALL?')
        return upper_str.split(',')
        
    def getSenseRangeUpper(self):
        upper_str = self.sensor.ask(":SENSE:RANGE:UPPER?")
        return upper_str
    
    def setAutoRangeEnabled(self, on):
        str = 'ON' if on else 'OFF'
        self.sensor.write(':SENS:AUTO ' + str)
    
    def isAutoRangeEnabled(self):
        ret = self.sensor.ask(':SENS:AUTO?')
        return ret == 'ON'

if __name__ == '__main__':

    node = MetrolabTHM1176Node()
<<<<<<< HEAD

    print(node.getSenseRangeUpper())
    print(node.getAvailableSenseRangeUpper())

=======
    #node.calibrate()
>>>>>>> cf8ca2f050158fd16bdeccc81c1b64dff013f3b1
    #Continuously measure and publish mag field until shutdown
    k = 0
    while k < 10:
        field=node.measureFieldmT()
        print(node.getUnit())
        k = k+1
        
        print("Measured Field: [%f,%f,%f]\n"%(field[0],field[1],field[2]))
        
        sleep(1)