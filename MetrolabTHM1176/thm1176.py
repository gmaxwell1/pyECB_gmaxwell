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
from time import sleep, time

class MetrolabTHM1176Node(object):
    """
    Class representing the metrolab THM1176-MF magnetometer. Can be used to read out data and adapt the sensor's settings.
    """
    ranges = ['0.1T', '0.3T', '1T', '3T']
    trigger_period_bounds = (122e-6, 2.79)
    base_fetch_cmd = ':FETC:ARR:'
    axes = ['X', 'Y', 'Z']
    # field_axes = ['Bx', 'By', 'Bz']
    # fetch_kinds = ['Bx', 'By', 'Bz', 'Timestamp',
    #                'Temperature']  # Order matters, this is linked to the fetch command that is sent to retrived data
    # 
    defaults = {'block_size': 10, 'period': 0.5, 'range': '0.1T', 'average': 1, 'unit' : 'MT', 'n_digits' : 5}
    # id_fields = ['manufacturer', 'model', 'serial', 'version']
    
    def __init__(self, *args, **kwargs):
        
        logging.basicConfig(filename='metrolab.log', level=logging.DEBUG)
        
        self.sensor = usbtmc.Instrument(0x1bfa, 0x0498)
        
        self.max_transfer_size = 49216  # (4096 samples * 3 axes * 4B/sample + 64B for time&temp&...
        # Show sensor name in Terminal
        ret = self.sensor.ask("*IDN?")
        print(ret)

        self.average_count = defaults['average']
        self.unit = defaults['unit']
        self.range = defaults['range'] # can be 0.1, 0.3, 1 or 3
        self.n_digits = defaults['n_digits']
        self.period = defaults['period']
        self.block_size = defaults['block_size']

        # Write settings to device
        #self.sensor.write(":format:data default")
        self.sensor.write(":unit " + self.unit)
        self.sensor.write(":sense:range:upper " + self.range)

        logging.debug('range upper %s', self.sensor.ask(":sense:range:upper?"))

        # self.period = 5

        logging.info('... End init')
        
        
    def setup(self, **kwargs):
        '''
        :param kwargs:
        :return:
        '''
        keys = kwargs.keys()

        if 'block_size' in keys:
            self.block_size = kwargs['block_size']

        if 'period' in keys:
            if self.trigger_period_bounds[0] <= kwargs['period'] <= self.trigger_period_bounds[1]:
                self.period = kwargs['period']
            else:
                print('Invalid trigger period value.')
                print('Setting to default...')
                self.period = self.defaults['period']

        if 'range' in keys:
            if kwargs['range'] in self.ranges:
                self.range = kwargs['range']

        if 'average' in keys:
            self.average_count = kwargs['average']

        if 'n_digits' in keys:
            if kwargs['n_digits'] < 5:
                self.n_digits = kwargs['n_digits']
            
        if 'unit' in keys:
            self.unit = kwargs['unit']
            
        
        self.sensor.write(":UNIT " + self.unit)
        self.sensor.write(":SENS " + self.range)
        self.sensor.write(":SENS:AUTO OFF")
        self.setAveragingCount()
        self.set_periodic_trigger()

        cmd = ''
        for axis in self.axes:
            cmd += self.base_fetch_cmd + axis + '? {},{};'.format(self.block_size, self.n_digits)
        cmd += ':FETC:TIM?;:FETCH:TEMP?;*STB?'
        self.fetch_cmd = cmd
        

    def calibrate(self):
        self.sensor.write(":CAL:INIT")
        self.sensor.write(":CAL:STAT ON")
        
    
    def setAveragingCount(self):
        avg_max = int(self.sensor.ask(":AVER:COUN? MAX"))
        
        if average_count <= avg_max and average_count > 0:
            self.sensor.write(":AVER:COUN {}".format(str(average_count)))
            return True
        else:
            print("MetrolabTHM1176:setAveragingCount: value has to be between 1 and " + str(avg_max))
            return False

    
    def getNumMeasurements(self):
        ret = self.sensor.ask(':CALC:AVER:COUNT?')
        return int(ret)
    
            
    def measureFieldmT(self):
        Bx = float(self.sensor.ask(':measure:scalar:flux:x? 0.01T,' + self.n_digits).strip('MT'))
        By = float(self.sensor.ask(':measure:y? 0.01T,' + self.n_digits).strip('MT'))
        Bz = float(self.sensor.ask(':measure:z? 0.01T,' + self.n_digits).strip('MT'))
        
        return [Bx, By, Bz]
    
    def measureFieldArraymT(self, num_meas=10):
        """
        Make a certain number of measurements of each field direction.

        Args:
            num_meas (int, optional): How many times to measure each component. Defaults to 10.
            
        Raises:
            ValueError: if not all measurements were made correctly.

        Returns:
        """
        ret = self.sensor.ask(":READ:array:x? " + str(num_meas) + ", 0.01T,5")
        Bx_str = ret.split(",")
        Bx = []
        for val in Bx_str:
            Bx.append(float(val.strip('MT')))
        
        ret = self.sensor.ask(":READ:array:y? " + str(num_meas) + ", 0.01T,5")
        By_str = ret.split(",")
        By = []
        for val in By_str:
            By.append(float(val.strip('MT')))
            
        ret = self.sensor.ask(":READ:array:z? " + str(num_meas) + ", 0.01T,5")
        Bz_str = ret.split(",")
        Bz = []
        for val in Bz_str:
            Bz.append(float(val.strip('MT')))
        
        if (len(Bx) != num_meas or len(By) != num_meas or len(Bz) != num_meas):
            raise ValueError("length of Bx, By, Bz do not match num_meas")
        
        return [Bx, By, Bz]
    
    # added by gmaxwell
    def set_periodic_trigger(self):
        """
        Set the probe to run in periodic trigger mode with a given period, continuously
        - param period

        Returns:
        """
        self.sensor.write(':TRIG:SOUR TIM')
        self.sensor.write(':TRIG:TIM {:f}S'.format(self.period))
        self.sensor.write(':TRIG:COUN {}'.format(self.block_size))
        self.sensor.write(':INIT:CONT ON')
        
        
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
    
    
    # context manager to ba able to use a with...as... statement    
    def __enter__(self):
        if not self.sensor.connected:
            self.sensor.open()
            
        return self
    
    
    def __exit__(self, type, value, traceback):
        if self.sensor.connected:
            self.sensor.close()
            return not self.sensor.connected
        else:
            return isinstance(value, TypeError)
        
        
if __name__ == '__main__':
    num_meas = 10
    with MetrolabTHM1176Node(sense_range_upper="0.1 T") as node:
        # node.sensor.write(":TRIG:TIM " + arg)
        # node.sensor.write(":TRIG:SOUR TIM")
        
        ret = node.sensor.ask(":READ:array:x? " + str(num_meas) + ", 0.01T,5")
        Bx_str = ret.split(",")
        Bx = []
        for val in Bx_str:
            Bx.append(float(val.strip('MT')))

        print(Bx)