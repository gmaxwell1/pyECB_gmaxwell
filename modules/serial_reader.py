""" 
filename: serial_reader.py

The following functions are used for communicating with the Hall sensor cube using a serial data port (USB)
and reading/writing the received data. 

Author: Jona Buehler 2020

Documentation and Updates by Nicholas Meinhardt, Maxwell Guerne-Kieferndorf (QZabre)
        					 nmeinhar@student.ethz.ch, gmaxwell@student.ethz.ch

Date: 09.10.2020
"""

########## Standard library imports ##########
import numpy as np
import serial
import os
from time import time
from datetime import datetime


def open_port(numport='4', baudrate=256000, timeout=2):
    """
    Open the COM-numport port and return instance of serial.Serial class.

    Args:
    - numport is the number of COM port that should be opened
    - baudrate
    - timeout [s] before raising exception
    """
    try:
        # open serial port; baud rate = 256000
        port = serial.Serial('COM{}'.format(numport), 256000, timeout=2)
        return port
    except:
        print('Failed to establish connection.')
        return


def ensure_dir_exists(directory, access_rights=0o755, purpose_text='', verbose=False):
    """
    Ensure that the directory exists and create the respective folders if it doesn't.

    Prints message that either a folder has been created or that it already exists. 
    Note that only read and write options can be set under Windows, rest ignored.

    Args:
    - directory is a path which should exist after calling this function
    - access_rights: set rights for reading and writing.
    - purpose_text (str): add more information what the dictionary was created for

    Return:
    - 0 if directory needed to be created
    - 1 if directory already exists
    - -1 if there was an exception
    """
    try:
        os.mkdir(directory, access_rights)
        os.chmod(directory, access_rights)
        if verbose:
            print('Created directory {}: {}'.format(
                purpose_text, os.path.split(directory)[1]))
        return 0
    except FileExistsError:
        if verbose:
            print('Folder already exists, no new folder created.')
        return 1
    except Exception as e:
        # if verbose:
        # this is important and should always be printed - gmaxwell, 8.10.2020
        print('Failed to create new directory due to {} error'.format(e))
        return -1


def get_new_data_set(interactive=False, numport='4', measure_runs=int(1), fname_postfix='measured', runtimelimit_per_run=5,
                                                dirname='data_sets', sub_dirname=None, cube: serial.Serial = None, verbose=False, no_enter=False,
                                                on_stage=False, specific_sensor=None, omit_64=False):
    """
    Read data from cube with 64 Hall sensors and either return (if specific sensor chosen) 
    or write the estimated B-field to an output csv-file (specific_sensor=None,  default).

    Note: The 64th sensor might be broken, hence it can be ignored using the omit_64 flag!

    Args:
    - interactive (bool): Flag to allow parameter adjustment over terminal window. 
      In this mode, it is not necessary to provide `cube`, `numport` and `measure_runs` parameters
    - numport (string): The COM port number of the sensor/serial device
    - measure_runs: INTEGER! Number of samples per fuction call per sensor
    - runtimelimit_per_run (float): in seconds, waiting time for sensor before error message 
      per measurment run

    - dirname: name of folder where data is stored in general. Can be left alone unless you need a new main data folder.
    - sub_dirname: name of folder where measured data files are stored. Make sure to always set this specifically!
    - fname_postfix (string): postfix of data file (csv).

    - cube: instance of serial.Serial class, representing the magnetic field sensor
    - no_enter (bool): if True, measurement starts automatically, else the user is asked to press enter to start.
      This flag only matters if specific_sensor == None!
    - on_stage (bool): flag used to set the action upon occurence an error when reading a measurement outcome 
      from the sensor. If False, continue measuring and write a "'Read Failure', 0,0,0,0"-line to file. 
      If True, the measurement is stopped entirely and the output file is deleted. 
    - specific_sensor (int in [1,64]): ID of a specific Hall sensor of the whole cube. 
      If this argument is provided, only the output of this sensor will be considered. 
      All remaining sensors are neglected

    - verbose: switching on/off print-statements for displaying progress

    Return (several possibilities):
    - if specific_sensor=None and 
     - on_stage = False: work_dir+output_file_name
     - on_stage = True: 0, str(work_dir), str(output_file_name),
    - specific_sensor != None: return either a 1-d ndarray containing the 3 measured magnetic field components
    or return 1 if measurement fails
    """
    # a = serial.Serial()
    # a.readline()

    # if specific sensor is selected
    if specific_sensor is not None:
        meas_time = 0
        start_time = time()

        if cube is None:
            cube = open_port(numport=numport)

        # print('Test output:', cube.readline()) # first line is typically incomplete, hence trash
        # Actually not needed, but still nice to check that everything is working:
        # the first line can certainly be split. If the first entry is not a number,
        # the comparison raises exception. if length less than 4 nothing is returned anyway

        while meas_time < (measure_runs*runtimelimit_per_run):
            try:
                meas_time = time() - start_time
                data = cube.readline().split(b',')
                if int(data[0].decode('ascii')) == int(specific_sensor) and len(data) == 4:
                    B_field = np.array([float(data[1].decode('ascii')),
                                        float(data[2].decode('ascii')),
                                        float(data[3].decode('ascii'))])
                    return B_field  # note that only a single measurement run is performed here!
            except:
                continue
        if verbose:
            print("Unable to read sensor ", specific_sensor,
                  " during ", measure_runs*runtimelimit_per_run, "s")
        return 1

    # consider all sensors now
    else:
        meas_time = 0

        # set the last sensor ID, depending on whether sensor 64 is omitted or not
        if omit_64:
            number_sensors = 63
        else:
            number_sensors = 64

        # initialize cube interactively, automatically or just take the passed cube argument if provided.
        if interactive:
            print(
                'Make sure to have the calibration cube positioned correctly (axis alignment)')
            while cube == None:
                numport = str(
                    input("Please enter the calibration cube port number X (COM): "))  # 4
                if numport == 'quit':
                    exit()
                cube = open_port(numport=numport)

                try:
                    measure_runs = int(
                        input('Please enter number of measurement runs: '))
                except:
                    print('Invalid Input: number of runs must be integer')
                fname_postfix = input(
                    'Please enter a name for your output file: ')
        elif cube is None:
            cube = open_port(numport=numport)

        # create a subfolder called dirname, if it already exists print a message -gmaxwell, 8.10.2020
        cwd = os.getcwd()
        work_dir = os.path.join(cwd, dirname)
        access_rights = 0o755
        ensure_dir_exists(work_dir, purpose_text='to save data',
                          access_rights=access_rights, verbose=verbose)
        # if specified, create subfolder to store data from this measurement run -gmaxwell, 8.10.2020
        if sub_dirname is not None:
            work_dir = os.path.join(work_dir, sub_dirname)
            ensure_dir_exists(work_dir, purpose_text='for specific measurement run',
                              access_rights=access_rights, verbose=verbose)

        os.chmod(work_dir, access_rights)
        now = datetime.now().strftime("%y_%m_%d_%H-%M-%S")  # Jona used "%d_%m_%y_%H-%M-%S"
        # name of the csv-file will be in the format: dd_mm_yy_hh-mm-ss_fname_postfix,
        # where fname_postfix is the input the user gave in the line above
        output_file_name = "{}_{}".format(now, fname_postfix + '.csv')

        # initialize array for for storing measurement values of all 64 sensors of cube
        measurements = []
        for i in range(64):  # this is the same as [[]]*64 !
            measurements.append([])

        with open(os.path.join(work_dir, output_file_name), 'w') as f:
            #os.chmod(work_dir + output_file_name+'\\', access_rights)
            # header line
            f.write(
                'Measurement Time (s), Sensor Number, X-Axis (mT), Y-Axis (mT), Z-Axis (mT)\r\n')

            if verbose:
                # first line is typically incomplete, hence trash
                print('Test output:', cube.readline())
            # Actually not needed, but still nice to check that everything is working:
            # the first line can certainly be split. If the first entry is not a number,
            # the comparison raises exception. if length less than 4 nothing is returned anyway

            if not no_enter:
                lets_go = input("Press Enter to start measurement!")

            start_time = time()
            sens_num_to_read = 1
            finished_runs = 0
            failure = False

            while meas_time < (measure_runs*runtimelimit_per_run):
                try:
                    meas_time = time() - start_time
                    data = cube.readline().split(b',')
                    #print(data[0].decode('ascii'), sens_num_to_read)
                    if int(data[0].decode('ascii')) == int(sens_num_to_read):
                        if len(data) == 4:
                            # = sens_num_to_read, else the upper if condition is not met
                            s_number = int(data[0])
                            field = (meas_time, float(data[1]), float(
                                data[2]), float(data[3]))
                            measurements[s_number-1].append(field)
                            f.write('{},{},{},{},{}\n'.format(meas_time, int(data[0]),
                                                              float(data[1]), float(data[2]), float(data[3])))
                            sens_num_to_read += 1
                except ValueError as e:
                    # this is a normal error during the first readout of a run,
                    # when only an incomplete message is received from the sensor.
                    # The following messages should be fine again
                    if 'invalid literal for int() with base 10:' == e.args[0][:39]:
                        pass
                    else:
                        print(e)
                        print('current output from sensor: {}'.format(
                            [data[i].decode('ascii') for i in range(len(data))]))
                        print("Have to restart measurement at this position.")
                        f.write('failure')
                        failure = True
                        return 1, str(work_dir), str(output_file_name)
                except Exception as e:
                    if not on_stage:
                        print('Line failed.')
                        print(e)
                        # this line appears several times during the measurement and we are not sure why.
                        # just ignore it, if it appears only a couple of times (5-10 times max per one run of 120s)
                        f.write('{},{},{},{},{}'.format(
                            'Read Failure', 0, 0, 0, 0) + '\r\n')
                    else:
                        print(e)
                        print('current raw output from sensor: {}'.format(
                            [data[i] for i in range(len(data))]))
                        print('current output from sensor: {}'.format(
                            [data[i].decode('ascii') for i in range(len(data))]))
                        print("Have to restart measurement at this position.")
                        f.write('failure')
                        failure = True
                        return 1, str(work_dir), str(output_file_name)

                # break while-loop if measure_runs runs have been performed
                # The 64th sensor seems to be broken, thus leave it out!
                if int(sens_num_to_read) > int(number_sensors):
                    # Note that not every sensor is evaluated measure_runs-times!
                    # Only the 63th sensor is certainly measured that often
                    finished_runs += 1
                    if finished_runs == measure_runs:
                        break
                    else:
                        sens_num_to_read = 1

        # print('Results of sensor 64:')
        # [print(measurements[63][i]) for i in range(len(measurements[63])) ]

        # delete the measurement file if an error occured
        if failure:
            os.remove(os.path.join(work_dir, output_file_name))

        if verbose:
            print('Number of finished runs: {}\n'.format(finished_runs))
        if not on_stage:
            return str(work_dir+output_file_name)
        else:
            return 0, str(work_dir), str(output_file_name)


# %%
if __name__ == "__main__":
    _ = get_new_data_set(interactive=True)

    # ------------------------Testing area---------------------------------------------------
    # %%
    # test the whole measurement setup

    # establish permanent connection to calibration cube
    port = 'COM4'
    # open serial port; baud rate = 256000
    t1 = time()
    with serial.Serial(port, 256000, timeout=2) as cube:
        results = get_new_data_set(
            numport=4, measure_runs=2000, cube=cube, on_stage=False, no_enter=True)
    duration = time() - t1
    print('duration: {:.4f} s'.format(duration))
    print(results)

    # %%
    # test duration and compare sensor 64 with others, use rod magnet to see a difference
    port = 'COM4'
    t1 = time()
    for _ in range(5):
        # open serial port; baud rate = 256000
        with serial.Serial(port, 256000, timeout=2) as cube:
            # results = get_new_data_set(numport=4, specific_sensor=1, measure_runs=int(2), cube=cube)
            # print('sensor 4: {}'.format(results))

            results = get_new_data_set(
                numport=4, specific_sensor=64, measure_runs=int(2), cube=cube)
            print('sensor 64: {}'.format(results))
    duration = time() - t1

    print(duration/100)
