""" 
filename: MetrolabMeasurements.py

The following functions are used for measurements of the magnetic field using the Metrolab THM1176 sensor. 

Author: Nicholas Meinhardt (QZabre)
        nmeinhar@student.ethz.ch

Date: 19.10.2020
"""

#%%
# standard library imports 
import numpy as np
import os 
from time import time, sleep
from datetime import datetime
import pandas as pd

# local imports 
try:
    from modules.serial_reader import MeasurementError, ensure_dir_exists
except ModuleNotFoundError:
    import sys
    sys.path.insert(1, os.path.join(sys.path[0], '..'))
    from modules.serial_reader import MeasurementError, ensure_dir_exists
finally:
    from MetrolabTHM1176.thm1176 import MetrolabTHM1176Node
    from modules.conexcc_control import save_in_dir
    from modules.general_functions import sensor_to_magnet_coordinates


#%%
def readoutMetrolabSensor(node: MetrolabTHM1176Node, measure_runs=1, fname_postfix='data_sets',
                          directory = './data_sets', verbose=False, save_data=True):
    """
    Read measurement outcomes of Metrolab THM1176 sensor and return the estimated B-field [mT] in magnet coordinates, 
    also save data if desired.

    Magnet coordinates: coil 2 is mounted along -y axis.

    Args:
    - node (MetrolabTHM1176Node): represents the Metrolab THM 1176 sensor
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
    if measure_runs == 1:
        # initialize ndarray to store the measurement data and time
        meas_data = np.zeros(3)
        meas_time = 0

        # get current time before starting
        t_start = time()

        # read the current output of sensor and save measured magnetic field
        meas_data = np.array(node.measureFieldmT())

        # save the current time
        meas_time = time() - t_start

    elif measure_runs > 1:
        # get current time before starting
        # t_start = time()

        # read the current output of sensor and save measured magnetic field
        meas_data = np.array(node.measureFieldArraymT(measure_runs)).swapaxes(0, 1)

        # assume that each measurement took the same time, such that the times of measurements 
        # are equally spaced between initial and final time and offset such that t_start = 0
        meas_time = np.linspace(0, time() - t_start, num=len(meas_data))
        
    else:
        raise ValueError("measure_runs must be >= 1!")

    # due to the setup, transform sensor coordinates to magnet coordinates
    meas_data = sensor_to_magnet_coordinates(meas_data)

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

def get_mean_dataset_MetrolabSensor(node: MetrolabTHM1176Node, sampling_size, verbose=False, max_num_retrials=5):
    """
    Estimate field vectors with Metrolab sensor sampling_size-times and return the mean and std 
    as 1d-ndarrays of length 3. 

    Args: 
    - node (MetrolabTHM1176Node): represents the Metrolab THM 1176 sensor
    - sampling_size (int): sampling size to estimate mean magnetic field vector and std, 
    i.e. number of times the sensor is read out in series before averaging 
    - verbose (bool): switching on/off print-statements for displaying progress

    Return: 
    - mean_data, std_data (ndarrays of of shape (number sensors, 3)): Mean magnetic field and 
    its standard deviation as a vector. 

    Raises MeasurementError if no valid measurement data could be aquired after max_num_retrials tries. 
    """
    # perform measurement and collect the raw data 
    for _ in range(max_num_retrials):
        try:
            meas_data = np.array(node.measureFieldArraymT(sampling_size)).swapaxes(0, 1)
        except:
            pass
        else:
            break
    
    try:
        # due to the setup, transform sensor coordinates to magnet coordinates
        meas_data = sensor_to_magnet_coordinates(meas_data)
        
    # if it was not possible to obtain valid measurement results after max_num_retrials, raise MeasurementError, too
    except UnboundLocalError:
        raise MeasurementError

    # estimate the mean and std from raw data for each sensor
    mean_data = np.mean(meas_data, axis=0)
    std_data = np.std(meas_data, axis=0)
    
    # # save raw data if desired
    # if save_mean_data:
    #     # if no directory is provided, just take current working directory
    #     if directory is None:
    #         directory = os.getcwd()
    #     save_in_dir(meas_data, directory, 'data', now=True)

    # return mean and std fields either for only this sensor or for all sensors
    return mean_data, std_data


#%%
if __name__ == '__main__':

    with MetrolabTHM1176Node() as sensor:
        # first try with single measurements using measureFieldmT method of sensor
        meas_times, field = readoutMetrolabSensor(sensor, measure_runs=10, directory='./test_data')

        durations = meas_times[1:] - meas_times[:-1]
        print('average duration (single) = {:.5f} s +- {:.5f} s'.format(np.mean(durations), np.std(durations)))
        print(field)

        # then try with array of measurements using measureFieldArraymT method of sensor
        # meas_times, field = readoutMetrolabSensor(sensor, measure_runs=10, directory='./test_data',
        #                                             single_measurements=False)
        # durations = meas_times[1:] - meas_times[:-1]
        # print('average duration (arrays) = {:.5f} s +- {:.5f} s'.format(np.mean(durations), np.std(durations)))
        # print(field)



# %%
