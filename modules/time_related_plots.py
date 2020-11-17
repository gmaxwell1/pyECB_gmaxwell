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


def generateAndSavePlot(filepath=r'.\data_sets\time_measurements_23_10\20_10_23_15-31-57_time_resolved.csv', plot_components='xyz', separate=False,
                        show_image=True, save_image=False, output_file_name='B_vs_t', save_dir=None, show_dev_from_mean=False, statistics=False):
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
    - output_file_name(str): what to call the output file
    - savedir (str): where to save the output file
    - show_dev_from_mean (bool): If False, the values of the desired components are plotted.
    If True, the deviation from their mean values is plotted. 
    - statistics (bool): show mean and standard deviation of each component, angles and magnitude.

    Return:
    - fig (plt.Figure): figure instance of the plot 
    - axs (plt.Axes): axes instance of the plot
    - times: array with timestamps
    - fields: field data 
    - plot_components: string with plot flags
    """
    # extract data and convert to ndarray
    raw_data = pd.read_csv(filepath).to_numpy()
    
    times = raw_data[:,0]
    fields = raw_data[:,1:4]
    
    if not separate:
        # generate Figure and Axis instances
        fig, ax = plt.subplots()
        ax = np.array([ax])
    else:
        fig, ax = plt.subplots(len(plot_components))
        fig.set_size_inches(6, len(plot_components) * 3)
        
    
    # plot the desired contents
    if show_dev_from_mean:
        # ax[0].hlines
        if 'x' in plot_components:
            ax[0].plot(times, fields[:,0] - np.mean(fields[:,0]), label='$\Delta$ $B_x$', color='C0')
        if 'y' in plot_components:
            ax[1].plot(times, fields[:,1] - np.mean(fields[:,1]), label='$\Delta$ $B_y$', color='C1')
        if 'z' in plot_components:
            ax[2].plot(times, fields[:,2] - np.mean(fields[:,2]), label='$\Delta$ $B_z$', color='C2')
        # if 'm' in plot_components:
        #     magn_B_fields = np.linalg.norm(fields, axis=1)
        #     ax.plot(times, magn_B_fields - np.mean(magn_B_fields), label = '$\Delta$ $|B|$')
        # if 'p' in plot_components:
        #     inplane_B_fields = np.linalg.norm(fields[:,:2], axis=1)
        #     ax.plot(times, inplane_B_fields - np.mean(inplane_B_fields), label = '$\Delta$ $|B_{xy}|$')
    else:
        ax[0].set_ylabel('magnetic flux density, $B$ [mT]')
        if 'x' in plot_components:
            if separate:
                ax[0].plot(times, fields[:,0], label = '$B_x$', color='C0')
                ax[0].set_ylabel('magnetic field, $B_x$ [mT]')                
            else:
                ax[0].plot(times, fields[:,0], label = '$B_x$', color='C0')
        if 'y' in plot_components:
            if separate:
                ax[1].plot(times, fields[:,1], label = '$B_y$', color='C1')
                ax[1].set_ylabel('magnetic field, $B_y$ [mT]')
            else:
                ax[0].plot(times, fields[:,1], label = '$B_y$', color='C1')
        if 'z' in plot_components:
            if separate:
                ax[2].plot(times, fields[:,2], label = '$B_z$', color='C2')
                ax[2].set_ylabel('magnetic field, $B_z$ [mT]')
            else:
                ax[0].plot(times, fields[:,2], label = '$B_z$', color='C2')
        # if 'm' in plot_components:
        #     magn_B_fields = np.linalg.norm(fields, axis=1)
        #     ax.plot(times, magn_B_fields, label = '$|B|$')
        # if 'p' in plot_components:
        #     inplane_B_fields = np.linalg.norm(fields[:,:2], axis=1)
        #     ax.plot(times, inplane_B_fields, label = '|$B_{xy}$|')

    # label axes
    ax[-1].set_xlabel('time, $t$ [s]')

    # show legend
    ax[0].legend()
    
    if statistics:
        mag_x = round(np.mean(fields[:,0]),2)
        mag_y = round(np.mean(fields[:,1]),2)
        mag_z = round(np.mean(fields[:,2]),2)    

        std_x = round(np.std(fields[:,0]),2)
        std_y = round(np.std(fields[:,1]),2)
        std_z = round(np.std(fields[:,2]),2)
        
        mag = round(np.sqrt(mag_x ** 2 + mag_y ** 2 + mag_z ** 2),2)
        theta = round(np.degrees(np.arccos(mag_z/mag)),2)
        phi = round(np.degrees(np.arctan2(mag_y, mag_x)),2)
            
        ax[0].set_title('$B_{{x,avg}}$ = {0} $\pm$ {1} $mT$\t$|B|$ = {2} $mT$\n$B_{{y,avg}}$ = {3} $\pm$ {4} $mT$\t$\\theta$ = {5}°\n$B_{{z,avg}}$ = {6} $\pm$ {7} $mT$\t$\\phi$ = {8}°'
                        .format(mag_x, std_x, mag, mag_y, std_y, theta, mag_z, std_z, phi), fontsize=16)

        plt.tight_layout()
    
    # save image
    if save_image:
        # set the directory name and current datetime if not passed as argument
        if save_dir is None:
            save_dir = os.getcwd()
        # now = datetime.now().strftime('%y_%m_%d_%H-%M-%S')

        img_path = os.path.join(save_dir, output_file_name)
        fig.savefig(img_path, dpi=300)

    if show_image:
        plt.show()
    
    return fig, ax, times, fields, plot_components

def add_insets_time_plots(axs, x_vals, plot_data, zoom_component, begin_idx=0, end_idx=5, inset_x = 0.15, inset_y = 0.15,
                          inset_ylim_factor = 0.25, manual_inset_ylim=None, color=None):
    """
    Add insets into angular plots, such that a line is resolved more precisely.

    Args:
    - axs (ndarray of plt.Axes objectes): represent the subplots, returned by generate_functions
    - x_vals (ndarray of length = #measurements): contain the x-values of all plots, 
    returned by generate_functions
    - plot_data, plot_std_data, plot_exp_data (nbdarrays): plot data that are returned 
    by the generate_plots function, they are of shape (#plots, #measurements, 1)
    - flags_yaxis (str): contains letters as flags, which indicate which quantity is plotted 
    on which axis. Hence, the len(flags_yaxis) = len(axs) = #plots. For this function, 
    only the two following flags are relevant: 
    't': angle theta with respect to z-axis, 
    'f': in-plane angle phi with respect to x-axis
    - inset_ylim_factor (float): this factor can be used to set the ylims of the inset axes. 
    The inset plot will show all data that are in [min - inset_ylim_factor*yrange, max + inset_ylim_factor*yrange],
    where min (max) is the minimum (maximum) value on the chosen side and yrange = max-min.
    - which_side (str): either 'left' or 'right', decides for which side the insets are created
    - manual_inset_ylim (None, tuple or list): the limits of y-axis of the inset axis 
    are set to the provided tuple(s). If a list is passed, it must have the same length as axs and each 
    entry contains the limits for the corresponding subplot. Default is None, in this case the limits 
    are estimated automatically using the min/max y-values
    and the inset_ylim_factor for all subplots. 
    - mask_selected_plots (ndarray of bools or None): Mask of the same length as axs, only for the 
    subplots where mask_selected_plots is True insets are added. This allows to avoid insets
    for certain subplots. Default is None, in which case insets are added for all subplots. 
    """
    if zoom_component == 'x':
        color='C0'
    elif zoom_component == 'y':
        color='C1'
    elif zoom_component == 'z':
        color='C2'

    data_in_area = plot_data[begin_idx:end_idx+1]
    
    x_area = x_vals[begin_idx:end_idx+1]
    
    # compute mean and rms values
    mean = np.mean(data_in_area)
    stdd = np.std(data_in_area)

    # inset axis
    axins = axs.inset_axes([inset_x, inset_y, 0.4, 0.3])
    # axins.errorbar(x_vals, plot_data[i], yerr=plot_std_data[i],
    #                             linestyle='', marker='.', capsize = 2
    axins.plot(x_area, data_in_area, color)
    axins.set_ylabel('$B_{}$'.format(zoom_component))

    # define which sub-regions of original plot should be shown by setting the limits of inset axis
    x_range = x_area[-1] - x_area[0]
    axins.set_xlim(x_area[0]-0.1*x_range, x_area[-1] + 0.1*x_range)        

    # unless manual limits for inset-y-axis are provided, choose them automatically
    y_range_inset = np.max(data_in_area) - np.min(data_in_area)
    axins.set_ylim( np.min(data_in_area) - inset_ylim_factor*y_range_inset, 
                    np.max(data_in_area) + inset_ylim_factor*y_range_inset)
    
    
    # add boxes to indicate which region of original plot is shown in the inset
    axs.indicate_inset_zoom(axins)
    
    return axs, round(mean,3), round(stdd,3)


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
    
    if number_plots == 1:
        axs = np.array([axs])
    #     axs[1] = np.array([axs[1]])
        
    plot_data = []
    ylabels = [[],[]]
    # exp averaging parameter
    for flag in plot_components:
        # magnetic field in x-direction
        if flag == 'x':
            ylabels[0].append('$B_x$ [mT]')
            # coeff = np.fft.fft(fields[:,0], n=1024)
            ylabels[1].append('x direction PSD')
            plot_data.append(fields[:,0])
        # magnetic field in y-direction
        elif flag == 'y':
            ylabels[0].append('$B_y$ [mT]')
            ylabels[1].append('y direction PSD')
            plot_data.append(fields[:,1])
        # magnetic field in z-direction
        elif flag == 'z':
            ylabels[0].append('$B_z$ [mT]')
            ylabels[1].append('z direction PSD')
            plot_data.append(fields[:,2])
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
        axs[i,0].plot(times, plot_data[i])
        a = axs[i,1].psd(plot_data[i], NFFT=512, Fs=fs, pad_to=len(times), c='g')
        axs[i,1].set_ylabel(ylabels[1][i])
        
    # axs[0,0].set_title(r'Measured magnetic field component', fontsize=16)
    axs[0,1].set_title(r'Power Spectral density (units: $dB_{10}/Hz$)', fontsize=16)
    
    # plt.show()
    # save image
    if save_image:
        # set the directory name and current datetime if not passed as argument
        if save_dir is None:
            save_dir = os.getcwd()
        now = datetime.now().strftime('%y_%m_%d_%H-%M-%S')

        output_file_name = '{}_{}.png'.format(now, image_name_postfix)
        file_path = os.path.join(save_dir, output_file_name)
        fig.savefig(file_path, dpi=300)
    
    return fig, axs, times, plot_data


if __name__ == "__main__":
    
    data_directory = r'data_sets\time_measurements_12_11'
    # # files = [ fi for fi in os.listdir(data_directory) if fi.endswith(".csv") ]
    # # for item in files:
    filepath = os.path.join(data_directory, '20_11_12_15-11-40_time_resolved.csv')
    # img_name = filepath.strip(data_directory).strip('_time_resolved.csv').strip('\\') + 'sinusoidal_3A'
    generateAndSavePlot(filepath=filepath, show_image=True, plot_components='xyz', save_image=False, save_dir=data_directory,
                        separate=False, statistics=True)
    
    # fig, ax, times, plot_data = spectralAnalysis(r'data_sets\time_measurements_03_11\20_11_03_15-31-01_time_resolved.csv', 'z', 2.5)
    
    # _, mean, std = add_insets_time_plots(ax[0], times, fields[:,2], 'z', begin_idx=1000, end_idx=1900, inset_x = 0.4, inset_y = 0.3,
    #                       inset_ylim_factor = 0.1, manual_inset_ylim=None, color=None)
    # ampl = np.amax(np.abs(plot_data))
    # print(np.abs(plot_data))
    # ampl_list = []
    # for item in np.abs(plot_data[0]):
    #     if item >= 0.95 * ampl:
    #         ampl_list.append(item)
    # ampl = round(np.mean(np.array(ampl_list)),2)
    # std_ampl = round(np.std(np.array(ampl_list)),2)
    # std = round(np.std(plot_data),2)
    
    # freq = 1/(times[1]-times[0])
    
    # ax[0,0].set_title('Sine frequency: 0.11  $Hz$, measured at 100 $Hz$\nAmplitude: {} $\pm$ {} $mT$ \nRMS: {} $mT_{{rms}}$'.format(ampl, std_ampl, std), fontsize=16)
    # ax[0,0].set_title('Sine frequency: 1 $Hz$, measured at {} $Hz$\nAmplitude: {} $\pm$ {} $mT$'.format(freq, ampl, std_ampl, std), fontsize=16)
    
    # mag_x = round(np.mean(fields[:,0]),2)
    # mag_y = round(np.mean(fields[:,1]),2)
    # mag_z = round(np.mean(fields[:,2]),2)    

    # std_x = round(np.std(fields[:,0]),2)
    # std_y = round(np.std(fields[:,1]),2)
    # std_z = round(np.std(fields[:,2]),2)
    
    # mag = round(np.sqrt(mag_x ** 2 + mag_y ** 2 + mag_z ** 2),2)
    # theta = round(np.degrees(np.arccos(mag_z/mag)),2)
    # phi = round(np.degrees(np.arctan2(mag_y, mag_x)),2)
        
    # ax[0].set_title('$B_{{x,avg}}$ = {0} $\pm$ {1} $mT$\t$|B|$ = {2} $mT$\n$B_{{y,avg}}$ = {3} $\pm$ {4} $mT$\t$\\theta$ = {5}°\n$B_{{z,avg}}$ = {6} $\pm$ {7} $mT$\t$\\phi$ = {8}°'
    #                 .format(mag_x, std_x, mag, mag_y, std_y, theta, mag_z, std_z, phi), fontsize=16)

    # plt.tight_layout()
    # # print(fields[1:6,0])
    
    # plt.show()