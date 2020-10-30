""" 
filename: time_related_plots.py

This file contains functions that are used for plotting data extracted from the Metrolab THM1176-MF Sensor. 
Mainly for plots of the magnetic field vs. time or Fourier analysis.

Author: Maxwell Guerne-Kieferndorf (Qzabre)
        gmaxwell@student.ethz.ch
        
Date: 27.10.2020
"""

import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, MaxNLocator
from datetime import datetime


def extract_time_dependence(filepath, sensorIsMetrolab=True, omit_64=False, ):
    """
    Extract and return time and field data from the provided file.

    Args: 
    - filepath (string): valid path of the data file
    - sensorIsMetrolab (bool): if True, the data originate from Metrolab THM1176 sensor, 
    else from Calibration Cube
    - omit_64 (bool): flag to omit sensor 64 if True (only reasonable if sensorIsMetrolab=False)

    Return:
    - times (1d-ndarray of length measure_runs): containing the time estimates of measurements.
    If sensorIsMetrolab=False, times is an ndarray of shape (number_sensors, measure_runs)
    - B_fields (ndarray of shape (measure_runs, 3)): contains measured x,y,z-components of magnetic field. 
    If , B_fields is an ndarray of shape (number_sensors, measure_runs, 3)
    """
    # import the measurement data from csv file 
    dataD = pd.read_csv(filepath)
    if sys.version_info[0] == 3:
        data = dataD.to_numpy()
    else:
        data = dataD.values

    if sensorIsMetrolab:
        # estimate number of measurement rounds
        measure_runs = len(data)

        # collect results for each sensor and save mean, std and abs(std/mean)     
        times = data[:,0]
        B_fields = data[:, 1:4]

    else:
        # adapt number of sensors depending on omit_64 flag
        if omit_64:
            number_sensors = 63
        else:
            number_sensors = 64
        
        # estimate number of measurement rounds
        measure_runs = len(data) // number_sensors
    
        # initialize arrays for measurement outcomes
        times = np.zeros((number_sensors, measure_runs))
        B_fields = np.zeros((number_sensors, measure_runs, 3))

        # collect results for each sensor and save mean, std and abs(std/mean)
        for i in range(number_sensors):
            # collect data for sensor i
            sensor_i_data = np.zeros((measure_runs, 5))
            for k in range(measure_runs):
                if (data[i+k*number_sensors,:].dtype == 'float64'):
                    sensor_i_data[k,:] = data[i+k*number_sensors,:]
                else:
                    print("could not convert data properly! wrong data type: ", data[i+k*number_sensors,:].dtype)
                    sensor_i_data[k,:] = 0
            
            times[i,:] = sensor_i_data[:,0]
            B_fields[i,:,:] = sensor_i_data[:, 2:5]
        
    return  times, B_fields


def generateAndSavePlot(filepath=r'.\data_sets\time_measurements_23_10\20_10_23_15-31-57_time_resolved.csv', plot_components='xyz',
                        show_image=True, save_image=False, save_dir=None, show_dev_from_mean=False, image_name_postfix='B_Field_Time'):
    """
    Generate a plot of field components on the y-axis vs. time on the x-axis.

    Args: 
    - times (1d ndarray of floats): contains the time values [s]
    - fields (ndarray of shape (len(time), 3)): contains the Bx, By and Bz components of the measured magnetic field
    - plot_component (string): contains letters as flags, where valid flags are 'x','y', 'z', 'm' and 'p'
    The number of letters may vary, but at least one valid letter should be contained. For each flag,
    the according quantity is added to the plot. Each quantity is plotted at most once.  
        - 'x', 'y', 'z': single component of magnetic field, p.e. B_z
        - 'm': magnitude of magnetic field
        - 'p': in-plane magnitude |B_xy| of magnetic field in xy-plane
    - show_image (bool): If True, plt.show() is executed at the end
    - show_dev_from_mean (bool): If False, the values of the desired components are plotted.
    If True, the deviation from their mean values is plotted. 

    Return:
    - fig (plt.Figure): figure instance of the plot 
    - axs (plt.Axes): axes instance of the plot
    """
    # extract data and convert to ndarray
    raw_data = pd.read_csv(filepath).to_numpy()
    
    times = raw_data[:,0]
    fields = raw_data[:,1:4]
    
    # generate Figure and Axis instances
    fig, ax = plt.subplots()
    
    # plot the desired contents
    if show_dev_from_mean:
        ax.set_ylabel('deviation from mean, $\Delta B$=$B$ - $\overline{B}$ [mT]')
        ax.hlines
        if 'x' in plot_components:
            ax.plot(times, fields[:,0] - np.mean(fields[:,0]), label = '$\Delta$ $B_x$')
        if 'y' in plot_components:
            ax.plot(times, fields[:,1] - np.mean(fields[:,1]), label = '$\Delta$ $B_y$')
        if 'z' in plot_components:
            ax.plot(times, fields[:,2] - np.mean(fields[:,2]), label = '$\Delta$ $B_z$')
        if 'm' in plot_components:
            magn_B_fields = np.linalg.norm(fields, axis=1)
            ax.plot(times, magn_B_fields - np.mean(magn_B_fields), label = '$\Delta$ $|B|$')
        if 'p' in plot_components:
            inplane_B_fields = np.linalg.norm(fields[:,:2], axis=1)
            ax.plot(times, inplane_B_fields - np.mean(inplane_B_fields), label = '$\Delta$ $|B_{xy}|$')
    else:
        ax.set_ylabel('magnetic field, $B$ [mT]')
        if 'x' in plot_components:
            ax.plot(times, fields[:,0], label = '$B_x$')
        if 'y' in plot_components:
            ax.plot(times, fields[:,1], label = '$B_y$')
        if 'z' in plot_components:
            ax.plot(times, fields[:,2], label = '$B_z$')
        if 'm' in plot_components:
            magn_B_fields = np.linalg.norm(fields, axis=1)
            ax.plot(times, magn_B_fields, label = '$|B|$')
        if 'p' in plot_components:
            inplane_B_fields = np.linalg.norm(fields[:,:2], axis=1)
            ax.plot(times, inplane_B_fields, label = '|$B_{xy}$|')

    # label axes
    ax.set_xlabel('time, $t$ [s]')

    # show legend
    ax.legend()
    
    # save image
    if save_image:
        # set the directory name and current datetime if not passed as argument
        if save_dir is None:
            save_dir = os.getcwd()
        now = datetime.now().strftime('%y_%m_%d_%H-%M-%S')

        output_file_name = '{}_{}.png'.format(now, image_name_postfix)
        file_path = os.path.join(save_dir, output_file_name)
        fig.savefig(file_path, dpi=300)

    if show_image:
        plt.show()
    
    return fig, ax


def spectralAnalysis(filepath=r'.\data_sets\time_measurements_23_10\20_10_23_15-31-57_time_resolved.csv', plot_components='xyz',
                      height_per_plot=2, save_image=False, save_dir=None, image_name_postfix='Power_Spectral_Density'):
    """
    Generate a plot of field components on the y-axis vs. time on the x-axis.

    Args: 
    - times (1d ndarray of floats): contains the time values [s]
    - fields (ndarray of shape (len(time), 3)): contains the Bx, By and Bz components of the measured magnetic field
    - plot_component (string): contains letters as flags, where valid flags are 'x','y', 'z', 'm' and 'p'
    The number of letters may vary, but at least one valid letter should be contained. For each flag,
    the according quantity is added to the plot. Each quantity is plotted at most once.  
        - 'x', 'y', 'z': single component of magnetic field, p.e. B_z
        - 'm': magnitude of magnetic field
    - show_image (bool): If True, plt.show() is executed at the end

    Return:
    - fig (plt.Figure): figure instance of the plot 
    - axs (plt.Axes): axes instance of the plot
    """
     # extract data and convert to ndarray
    raw_data = pd.read_csv(filepath).to_numpy()
    
    times = raw_data[:,0]
    dt = times[1]
    fs = 1 / dt # sampling frequency for FFT
    fields = raw_data[:,1:4]

    # if an empty string is passed as plot_components, set it to 'x'
    if len(plot_components) == 0:
        plot_components = 'x'

    # create a simple plot with as many axes as letters in plots_yaxis
    number_plots = len(plot_components)
    fig, axs = plt.subplots(number_plots, 2, sharex='col')
    fig.set_size_inches(10, number_plots * height_per_plot)
    
    ylabels = [[],[]]
    # exp averaging parameter
    for flag in plot_components:
        # magnetic field in x-direction
        if flag == 'x':
            ylabels[0].append('$B_x$ [mT]')
            coeff = np.fft.fft(fields[:,0], n=1024)
            ylabels[1].append('x direction PSD')
        # magnetic field in y-direction
        elif flag == 'y':
            ylabels[0].append('$B_y$ [mT]')
            ylabels[1].append('y direction PSD')
        # magnetic field in z-direction
        elif flag == 'z':
            ylabels[0].append('$B_z$ [mT]')
            ylabels[1].append('z direction PSD')
        else:
            raise ValueError(
                '{} is not a valid flag, it should be in [\'x\', \'y\', \'z\']!'.format(flag))
    
    # label axes
    # xticks = np.arange(times[0], times[-1], 25)
    axs[-1,0].set_xlabel('Time [s]')
    # axs[-1,0].set_xticks(xticks)
    # AutoMinorLocator()
        
    for i in range(len(axs)):
        axs[i,0].set_ylabel(ylabels[0][i])
        axs[i,0].plot(times, fields[:,i])
        axs[i,1].psd(fields[:,i], NFFT=512, Fs=fs, pad_to=len(times), c='g')
        axs[i,1].set_ylabel(ylabels[1][i])
        
    axs[0,0].set_title(r'Measured magnetic field component', fontsize=16)
    axs[0,1].set_title(r'Power Spectral density (units: $dB_{10}/Hz$)', fontsize=16)
    
    plt.show()
    
    # save image
    if save_image:
        # set the directory name and current datetime if not passed as argument
        if save_dir is None:
            save_dir = os.getcwd()
        now = datetime.now().strftime('%y_%m_%d_%H-%M-%S')

        output_file_name = '{}_{}.png'.format(now, image_name_postfix)
        file_path = os.path.join(save_dir, output_file_name)
        fig.savefig(file_path, dpi=300)
    
    return fig, axs, 


if __name__ == "__main__":
    # generateAndSavePlot(filepath=r'.\data_sets\time_measurements_27_10\20_10_27_11-10-20_time_resolved.csv', save_image=True, save_dir=r'.\data_sets\time_measurements_27_10', image_name_postfix='stability_test')
    
    _, axs = spectralAnalysis(filepath=r'.\data_sets\time_measurements_23_10\20_10_23_15-29-43_time_resolved.csv', plot_components='xyz',
                      save_image=True, save_dir=r'.\data_sets\time_measurements_23_10')
    
    plt.psd