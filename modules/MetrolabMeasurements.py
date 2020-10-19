""" 
filename: MetrolabMeasurements.py

The following functions are used for measurements of the magnetic field using the Metrolab THM1176 sensor. 

Date: 19.10.2020
"""

#%%
import numpy as np
import os 
from time import time, sleep
from datetime import datetime
import pandas as pd

try:
    from modules.serial_reader import MeasurementError, ensure_dir_exists
    from MetrolabTHM1176.thm1176 import MetrolabTHM1176Node
except ModuleNotFoundError:
    import sys
    sys.path.insert(1, os.path.join(sys.path[0], '..'))
    from modules.serial_reader import MeasurementError, ensure_dir_exists
    from MetrolabTHM1176.thm1176 import MetrolabTHM1176Node

#%%
def readoutMetrolabSensor(node: MetrolabTHM1176Node, measure_runs=1, fname_postfix='data_sets',
                    directory = './data_sets', verbose=False, save_data=True):
    """
    Read measurement outcomes of Metrolab THM1176 sensor and return the estimated B-field [mT], also save data if desired.

    Args:
    - cube (serial.Serial): represents the magnetic field sensor
    - measure_runs (int): Number of samples per fuction call per sensor
    - save_data (bool): if True, the results are stored in a csv-file
    - directory (string): valid path of the folder where data should be stored. The default name of the data file is
    'yy_mm_dd_hh-mm-ss_{fname_postfix}.csv'
    - fname_postfix (string): postfix of data file (csv).

    Return:
    - meas_time (ndarray of length=measure_runs): Contains the time of each measurement relative 
    to start of measurements 
    - meas_data (ndarray of shape=(measure_runs, 3)): Contains the measured field components

    Exceptions: 
    Raise a MeasurementError if something went wrong during measurements. 
    Possible reasons are that a sensor was skipped or that an incomplete message was received from the sensor. 

    """
    # initialize ndarray to store the measurement data and time
    meas_data = np.zeros((measure_runs, 3))
    meas_time = np.zeros(measure_runs)

    # get current time before starting
    t_start = time()

    # read in data from sensor for all measurement runs and sensors 
    for run in range(measure_runs):

        # read the current output of sensor and save measured magnetic field
        meas_data[run, :] = node.measureFieldmT()

        # save the current time
        current_time = time() - t_start
        meas_time[run] = current_time

    # save data if desired
    if save_data:

        # Measurement Time (s), Sensor Number, X-Axis (mT), Y-Axis (mT), Z-Axis (mT)
        df = pd.DataFrame({ 'Measurement Time (s)': meas_time, 
                            'Bx [mT]': meas_data[:, 0],
                            'By [mT]': meas_data[:, 1],
                            'Bz [mT]': meas_data[:, 2]})

        if directory is None:
            directory = os.getcwd()
        ensure_dir_exists(directory)

        now = datetime.now().strftime("%y_%m_%d_%H-%M-%S") 
        output_file_name = "{}_{}.csv".format(now, fname_postfix)
        data_filepath = os.path.join(directory, output_file_name)

        try:
            df.to_csv(data_filepath, index=False, header=True)
        except FileNotFoundError:
            data_filepath = os.path.join(os.getcwd(), output_file_name)
            df.to_csv(data_filepath, index=False, header=True)
    
    return meas_time, meas_data


#%%
if __name__ == '__main__':

    sensor = MetrolabTHM1176Node()
    # sensor.calibrate()


    meas_times, field = readoutMetrolabSensor(sensor, measure_runs=2, directory='./test_data')

    print(field)




# %%
