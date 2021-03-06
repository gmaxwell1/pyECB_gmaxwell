"""
filename: time_related_plots.py

This file contains functions that are used for plotting data extracted from various sensors
(Metrolab, ADT temp sensors).
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


# low pass filter with cutoff at 1 Hz to cut out noise from 100 Hz sampling rate measurements
lowPass100Hz = np.array([-0.0000,  0.0012,  0.0013,  0.0020,  0.0030,  0.0041,  0.0055,  0.0072,  0.0091,  0.0113,
                         0.0137,  0.0164,  0.0192,  0.0222, 0.0253,  0.0284,  0.0315,  0.0344,  0.0372,  0.0398,
                         0.0420,  0.0439,  0.0453,  0.0463,  0.0468,  0.0468,  0.0463,  0.0453, 0.0439,  0.0420,  0.0398,
                         0.0372,  0.0344,  0.0315,  0.0284,  0.0253,  0.0222,  0.0192,  0.0164,  0.0137,  0.0113,  0.0091,
                         0.0072,  0.0055,  0.0041,  0.0030,  0.0020,  0.0013,  0.0012, -0.0000])
# low pass filter with cutoff at 1 Hz to cut out noise from 100 Hz sampling rate measurements
lowPass100Hz_1 = np.array([0.0002, 0.0005, 0.0010, 0.0017, 0.0027, 0.0040, 0.0058, 0.0080, 0.0108, 0.0140, 0.0176, 0.0217,
                           0.0260, 0.0305, 0.0351, 0.0395, 0.0436, 0.0473, 0.0503, 0.0526, 0.0540, 0.0544, 0.0540, 0.0526,
                           0.0503, 0.0473, 0.0436, 0.0395, 0.0351, 0.0305, 0.0260, 0.0217, 0.0176, 0.0140, 0.0108, 0.0080,
                           0.0058, 0.0040, 0.0027, 0.0017, 0.0010, 0.0005, 0.0002])

# low pass filter with cutoff at 1 Hz to cut out noise from 20 Hz sampling rate measurements
lowPass20Hz = np.array([-0.0184, 0.0246, 0.1335, 0.2671,
                       0.3285, 0.2671, 0.1335, 0.0246, -0.0184])


def extract_time_dependence(filepath, sensorIsMetrolab=True, omit_64=False):
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
        times = data[:, 0]
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
                if (data[i+k*number_sensors, :].dtype == 'float64'):
                    sensor_i_data[k, :] = data[i+k*number_sensors, :]
                else:
                    print("could not convert data properly! wrong data type: ",
                          data[i+k*number_sensors, :].dtype)
                    sensor_i_data[k, :] = 0

            times[i, :] = sensor_i_data[:, 0]
            B_fields[i, :, :] = sensor_i_data[:, 2:5]

    return times, B_fields


def generateAndSaveTempPlot(times, temps, plot_components='123', separate=False, show_image=True, save_image=False,
                            output_file_name='Temp_vs_t', save_dir=None, show_dev_from_mean=False, statistics=False):
    """
    Generate a plot of temperature measurements on the y-axis vs. time on the x-axis.

    Args:
    - times: np.array with N entries (timeline)
    - temps: np.ndarray with 3xN entries for each sensor
    - plot_component (string): contains letters as flags, where valid flags are '1','2','3'
    The number of letters may vary, but at least one valid letter should be contained. For each flag,
    the according quantity is added to the plot. Each quantity is plotted at most once.
        - '1', '2', '3': temperature sensor 1,2 or 3 data
    - separate (bool): plot each sensor's data separately
    - show_image (bool): If True, plt.show() is executed at the end
    - output_file_name(str): what to call the output file
    - save_dir (str): where to save the output file
    - show_dev_from_mean (bool): If False, the values of the desired components are plotted.
    If True, the deviation from their mean values is plotted.
    - statistics (bool): show mean and standard deviation of each component, angles and magnitude.

    Return:
    - fig (plt.Figure): figure instance of the plot
    - axs (plt.Axes): axes instance of the plot
    - times: array with timestamps
    - temps: temperature data from each sensor
    - plot_components: string with plot flags
    """
    
    numplots = 1
    if not separate:
        # generate Figure and Axis instances
        fig, ax = plt.subplots(numplots, sharex=True)
        if numplots == 1:
            ax = np.array([ax])
        fig.set_size_inches(8, 5)
    else:
        numplots = len(plot_components)
        fig, ax = plt.subplots(numplots, sharex=True)
        fig.set_size_inches(8, numplots * 3)

    # plot the desired contents
    i = 0
    for c in plot_components:
        if show_dev_from_mean:
            if c == '1':
                if separate:
                    ax[0].plot(times, temps[:, 0] - np.mean(temps[:, 0]),
                                 label='$\Delta$ $T_1$', color='C0')
                    ax[0].set_ylabel('$\Delta$ $T_1$ [°C]')
                else:
                    ax[0].plot(times, temps[:, 0] - np.mean(temps[:, 0]),
                                 label='$\Delta$ $T_1$', color='C0')
                    ax[0].set_ylabel('Temperature Deviation $\Delta$ $T$ [°C]')
            elif c == '2':
                if separate:
                    ax[i].plot(times, temps[:, 1] - np.mean(temps[:, 1]),
                                 label='$\Delta$ $T_2$', color='C1')
                    ax[i].set_ylabel('$\Delta$ $T_2$ [°C]')
                else:
                    ax[0].plot(times, temps[:, 1] - np.mean(temps[:, 1]),
                                 label='$\Delta$ $T_2$$', color='C1')
                    ax[0].set_ylabel('Temperature Deviation $\Delta$ $T$ [°C]')
            elif c == '3':
                if separate:
                    ax[i].plot(times, temps[:, 2] - np.mean(temps[:, 2]),
                                 label='$\Delta$ $T_3$', color='C2')
                    ax[i].set_ylabel('$\Delta$ $T_3$ [°C]')
                else:
                    ax[0].plot(times, temps[:, 2] - np.mean(temps[:, 2]),
                                 label='$\Delta$ $T_3$', color='C2')
                    ax[0].set_ylabel('Temperature Deviation $\Delta$ $T$ [°C]')
        else:
            if c == '1':
                if separate:
                    ax[0].plot(times, temps[:, 0], label='$T_1$', color='C0')
                    ax[0].set_ylabel('$T_1$ [°C]')
                else:
                    ax[0].plot(times, temps[:, 0], label='$T_1$', color='C0')
                    ax[0].set_ylabel('Temperature Deviation $T$ [°C]')
            elif c == '2':
                if separate:
                    ax[i].plot(times, temps[:, 1], label='$T_2$', color='C1')
                    ax[i].set_ylabel('$T_2$ [°C]')
                else:
                    ax[0].plot(times, temps[:, 1], label='$T_2$', color='C1')
                    ax[0].set_ylabel('Temperature Deviation $T$ [°C]')
            elif c == '3':
                if separate:
                    ax[i].plot(times, temps[:, 2], label='$T_3$', color='C2')
                    ax[i].set_ylabel('$T_3$ [°C]')
                else:
                    ax[0].plot(times, temps[:, 2], label='$T_3$', color='C2')
                    ax[0].set_ylabel('Temperature Deviation $T$ [°C]')
                    
        i = i+1

    # label axes
    ax[-1].set_xlabel('time, $t$ [s]')
    # show legend
    if not separate:
        ax[0].legend()

    # save image
    if save_image:
        # set the directory name and current datetime if not passed as argument
        if save_dir is None:
            save_dir = os.getcwd()
        now = datetime.now().strftime('%y_%m_%d_%H-%M-%S')
        output_file_name = now + output_file_name
        img_path = os.path.join(save_dir, output_file_name)
        fig.savefig(img_path, dpi=300)

    if show_image:
        plt.show()

    return fig, ax, times, temps, plot_components


def generateAndSavePlot(filepath=r'.\data_sets\time_measurements_23_10\20_10_23_15-31-57_time_resolved.csv', plot_components='xyz', separate=False,
                        show_image=True, save_image=False, save_dir=None, output_file_name='B_vs_t', show_dev_from_mean=False, statistics=False):
    """
    Generate a plot of field components on the y-axis vs. time on the x-axis.

    Args: 
    - filepath: path of csv file to be read/plotted
    - plot_component (string): contains letters as flags, where valid flags are 'x','y', 'z', 'm' and 'p'
    The number of letters may vary, but at least one valid letter should be contained. For each flag,
    the according quantity is added to the plot. Each quantity is plotted at most once.  
        - 'x', 'y', 'z': single component of magnetic field, p.e. B_z
        - 'm': magnitude of magnetic field
        - 'p': in-plane magnitude |B_xy| of magnetic field in xy-plane
    - separate (bool): plot each field component separately
    - save_image (bool): if true, save image in save_dir, a directory path
    - show_image (bool): If True, plt.show() is executed at the end
    - output_file_name(str): what to call the output file
    - show_dev_from_mean (bool): If False, the values of the desired components are plotted.
    If True, the deviation from their mean values is plotted. 
    - statistics (bool): show mean and standard deviation of each separate component, angles and magnitude.

    Return:
    - fig (plt.Figure): figure instance of the plot 
    - axs (plt.Axes): axes instance of the plot
    - times: array with timestamps
    - fields: field data 
    - plot_components: string with plot flags
    """
    # extract data and convert to ndarray
    raw_data = pd.read_csv(filepath).to_numpy()

    times = raw_data[:, 0]
    fields = raw_data[:, 1:4]
    if len(raw_data[0]) == 5:
        temp = raw_data[:, 4]
    numplots = 1

    if not separate:
        # generate Figure and Axis instances
        if 't' in plot_components:
            numplots = 2
        fig, ax = plt.subplots(numplots, sharex=True)
        if numplots == 1:
            ax = np.array([ax])
        fig.set_size_inches(8, 5)
    else:
        numplots = len(plot_components)
        fig, ax = plt.subplots(numplots, sharex=True)
        fig.set_size_inches(8, numplots * 3)

    # plot the desired contents
    i = 0
    for c in plot_components:
        if show_dev_from_mean:
            if c == 'x':
                if separate:
                    ax[0].plot(times, fields[:, 0] - np.mean(fields[:, 0]),
                                 label='$\Delta$ $B_x$', color='C0')
                    ax[0].set_ylabel('$\Delta$ $B_x$ [mT]')
                else:
                    ax[0].plot(times, fields[:, 0] - np.mean(fields[:, 0]),
                                 label='$\Delta B_x$', color='C0')
                    ax[0].set_ylabel('Magnetic Field Deviation $\Delta$ $B$ [mT]')
            elif c == 'y':
                if separate:
                    ax[i].plot(times, fields[:, 1] - np.mean(fields[:, 1]),
                                 label='$\Delta B_y$', color='C1')
                    ax[i].set_ylabel('$\Delta$ $B_y$ [mT]')
                else:
                    ax[0].plot(times, fields[:, 1] - np.mean(fields[:, 1]),
                                 label='$\Delta B_y$', color='C1')
                    ax[0].set_ylabel('Magnetic Field Deviation $\Delta$ $B$ [mT]')
            elif c == 'z':
                if separate:
                    ax[i].plot(times, fields[:, 2] - np.mean(fields[:,2]),
                                 label='$\Delta B_z$', color='C2')
                    ax[i].set_ylabel('$\Delta$ $B_z$ [mT]')
                else:
                    ax[0].plot(times, fields[:, 2] - np.mean(fields[:,2]),
                                 label='$\Delta B_z$', color='C2')
                    ax[0].set_ylabel('Magnetic Field Deviation $\Delta$ $B$ [mT]')
            elif c == 't':
                if separate:
                    ax[numplots-1].plot(times, temp - np.mean(temp),
                                          label='$\Delta T$', color='C3')
                    ax[numplots-1].set_ylabel('Temp. Deviation $\Delta$ $T$ [no unit]')
                else:
                    ax[numplots-1].plot(times, temp - np.mean(temp),
                                          label='$\Delta T$', color='C3')
                    ax[numplots-1].set_ylabel('Temp. Deviation $\Delta$ $T$ [no unit]')
                    ax[numplots-1].legend()
        else:
            if c == 'x':
                if separate:
                    ax[0].plot(times, fields[:, 0], label='$B_x$', color='C0')
                    ax[0].set_ylabel('magnetic field, $B_x$ [mT]')
                else:
                    ax[0].plot(times, fields[:, 0], label='$B_x$', color='C0')
                    ax[0].set_ylabel('magnetic flux density, $B$ [mT]')
            elif c == 'y':
                if separate:
                    ax[i].plot(times, fields[:, 1], label='$B_y$', color='C1')
                    ax[i].set_ylabel('magnetic field, $B_y$ [mT]')
                else:
                    ax[0].plot(times, fields[:, 1], label='$B_y$', color='C1')
                    ax[0].set_ylabel('magnetic flux density, $B$ [mT]')
            elif c == 'z':
                if separate:
                    ax[i].plot(times, fields[:, 2], label='$B_z$', color='C2')
                    ax[i].set_ylabel('magnetic field, $B_z$ [mT]')
                else:
                    ax[0].plot(times, fields[:, 2], label='$B_z$', color='C2')
                    ax[0].set_ylabel('magnetic flux density, $B$ [mT]')
            elif c == 't':
                if separate:
                    ax[numplots-1].plot(times, fields[:, 1], label='$T$', color='C3')
                    ax[numplots-1].set_ylabel('Absolute temperature, $T$ [no unit]')
                else:
                    ax[numplots-1].plot(times, temp, label='$T$', color='C3')
                    ax[numplots-1].set_ylabel('Absolute temperature, $T$ [no unit]')
                    ax[numplots-1].legend()
                    
        i = i+1

    # label axes
    ax[-1].set_xlabel('time, $t$ [s]')

    # show legend
    if not separate:
        ax[0].legend()

    if statistics:
        mag_x = round(np.mean(fields[:, 0]), 2)
        mag_y = round(np.mean(fields[:, 1]), 2)
        mag_z = round(np.mean(fields[:, 2]), 2)

        std_x = round(np.std(fields[:, 0]), 2)
        std_y = round(np.std(fields[:, 1]), 2)
        std_z = round(np.std(fields[:, 2]), 2)

        mag = round(np.sqrt(mag_x ** 2 + mag_y ** 2 + mag_z ** 2), 2)
        theta = round(np.degrees(np.arccos(mag_z/mag)), 2)
        phi = round(np.degrees(np.arctan2(mag_y, mag_x)), 2)
        
        # delta_temp = np.amax(temp) - np.amin(temp)

        ax[0].set_title(f'$B_{{x,avg}}$ = {mag_x} $\pm$ {std_x} $mT$\t$|B|$ = {mag} $mT$\n$B_{{y,avg}}$ = {mag_y} '
                          f'$\pm$ {std_y} $mT$\t$\\theta$ = {theta}°\n$B_{{z,avg}}$ = {mag_z} $\pm$ {std_z} $mT$'
                          f'\t$\\phi$ = {phi}°', fontsize=16)

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



def add_insets_time_plots(axs, x_vals, plot_data, zoom_component, begin_idx=0, end_idx=5, inset_x=0.15, inset_y=0.15,
                          inset_ylim_factor=0.25, manual_inset_ylim=None, color=None):
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
        color = 'C0'
    elif zoom_component == 'y':
        color = 'C1'
    elif zoom_component == 'z':
        color = 'C2'

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
    axins.set_ylim(np.min(data_in_area) - inset_ylim_factor*y_range_inset,
                   np.max(data_in_area) + inset_ylim_factor*y_range_inset)

    # add boxes to indicate which region of original plot is shown in the inset
    axs.indicate_inset_zoom(axins)

    return axs, round(mean, 3), round(stdd, 3)


def spectralAnalysis(filepath=r'.\data_sets\time_measurements_23_10\20_10_23_15-31-57_time_resolved.csv', plot_components='xyz',
                     height_per_plot=2, save_image=False, save_dir=None, image_name_postfix='Power_Spectral_Density'):
    """
    Generate a plot of the Power spectral density of field measurements on the y-axis vs. frequency on the x-axis.
    Includes corresponding time domain plots as well.

    Args: 
    - filepath: path of file containing data to be plotted (must be csv file)
    - plot_component (string): contains letters as flags, where valid flags are 'x','y', 'z', 'm'
    The number of letters may vary, but at least one valid letter should be contained. For each flag,
    the according quantity is added to the plot. Each quantity is plotted at most once.  
        - 'x', 'y', 'z': single component of magnetic field, p.e. B_z
        - 'm': magnitude of magnetic field
    - height_per_plot: the height of each plot in inches
    - save_image (bool): If True, image is saved in save_dir, a directory path
    - image_name_postfix: this name is added to the end of the file.
    
    Return:
    - fig (plt.Figure): figure instance of the plot 
    - axs (plt.Axes): axes instance of the plot
    - times, plot_data: time stamps, data to be plotted
    """
    # extract data and convert to ndarray
    raw_data = pd.read_csv(filepath).to_numpy()

    times = raw_data[:, 0]
    dt = times[1]
    fs = 1 / dt  # sampling frequency for FFT
    fields = raw_data[:, 1:4]

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
    ylabels = [[], []]
    # exp averaging parameter
    for flag in plot_components:
        # magnetic field in x-direction
        if flag == 'x':
            ylabels[0].append('$B_x$ [mT]')
            # coeff = np.fft.fft(fields[:,0], n=1024)
            ylabels[1].append('x direction PSD')
            plot_data.append(fields[:, 0])
        # magnetic field in y-direction
        elif flag == 'y':
            ylabels[0].append('$B_y$ [mT]')
            ylabels[1].append('y direction PSD')
            plot_data.append(fields[:, 1])
        # magnetic field in z-direction
        elif flag == 'z':
            ylabels[0].append('$B_z$ [mT]')
            ylabels[1].append('z direction PSD')
            plot_data.append(fields[:, 2])
        else:
            raise ValueError(
                '{} is not a valid flag, it should be in [\'x\', \'y\', \'z\']!'.format(flag))

    # label axes
    # xticks = np.arange(times[0], times[-1], 25)
    axs[-1, 0].set_xlabel('Time [s]')
    # axs[-1,0].set_xticks(xticks)
    # AutoMinorLocator()

    for i in range(len(axs)):
        axs[i, 0].set_ylabel(ylabels[0][i])
        axs[i, 0].plot(times, plot_data[i])
        a = axs[i, 1].psd(plot_data[i], NFFT=512, Fs=fs,
                          pad_to=len(times), c='g')
        axs[i, 1].set_ylabel(ylabels[1][i])

    # axs[0,0].set_title(r'Measured magnetic field component', fontsize=16)
    axs[0, 1].set_title(
        r'Power Spectral density (units: $dB_{10}/Hz$)', fontsize=16)

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
    data_directory = r'C:\Users\Magnebotix\Desktop\Qzabre_Vector_Magnet\2_Misc_Code\Temperature Sensors\ADT7410_temperature_measurements\Measurement_over_time'
    # files = [ fi for fi in os.listdir(data_directory) if fi.endswith(".csv") ]
    # for item in files:
    filename = '20_12_15_13-03-56_Siglent_3A_1coil_2base_3pole.csv'

    filepath = os.path.join(data_directory, filename)
    raw_data = pd.read_csv(filepath).to_numpy()
    # in case columns are swapped for some reason
    # times = raw_data[:, 3]
    # times = (times * 128) / 1000
    # raw_data[:, 0] = (raw_data[:, 0] * 1000) / 128
    times = raw_data[:, 0]
    
    dt = times[5] - times[4]
    N = len(times)
    times_new = np.arange(0, N * dt, dt)
    # times_new = times
    temps = raw_data[:, 1:4]
    # temps = np.array([raw_data[:, 0],raw_data[:, 1],raw_data[:, 2]])
    # temps = np.swapaxes(temps, 0, 1)

    # # digitally filter data measured before plotting
    # LPF = lowPass100Hz_1
    # filtered_data = []
    # # global lowPass100Hz
    # for n in range(len(temps)):
    #     # remove clearly wrong spikes:
    #     if n > 10 and abs(temps[n, 0] - filtered_data[n-1][0]) > 30:
    #         temps[n, 0] = filtered_data[n-1][0]
    #     if n > 10 and abs(temps[n, 1] - filtered_data[n-1][1]) > 30:
    #         temps[n, 1] = filtered_data[n-1][1]
    #     if n > 10 and abs(temps[n, 2] - filtered_data[n-1][2]) > 30:
    #         temps[n, 2] = filtered_data[n-1][2]
            
    #     avg_data_point = np.array([0, 0, 0])
    #     # if 0 < n <= len(LPF) - 1:
    #     #     for m in range(n):
    #     #         avg_data_point = avg_data_point + temps[m, :] * LPF[n-m]
    #     if n > len(LPF) - 1:
    #         index_start = n - len(LPF) + 1
    #         for m in range(len(LPF)):
    #             k = index_start + m
    #             avg_data_point = avg_data_point + temps[k, :] * LPF[n-k]
    #     else:
    #         avg_data_point = temps[n, :]
        
    #     filtered_data.append(avg_data_point)
        
    # fieldFiltered = np.array(filtered_data)

        
    fig1, ax1, times1, fields1, plot_components1 = generateAndSaveTempPlot(times_new, temps, show_image=False, plot_components='123', save_image=False,
                                                                           save_dir=data_directory, output_file_name='Temp_vs_t', separate=True, show_dev_from_mean=False,
                                                                           statistics=False)

    # fig, ax, times, fields = spectralAnalysis(filepath, 'xyz', 2.5)

    # _, mean, std = add_insets_time_plots(ax1[0], times1, fields1[:,0]-np.mean(fields1[:,0]), 'x', begin_idx=470, end_idx=1270, inset_x = 0.5, inset_y = 0.6,
    #                       inset_ylim_factor = 0.1, manual_inset_ylim=None, color=None)
    # # ampl = np.amax(np.abs(plot_data))
    # print(np.abs(plot_data))
    # ampl_list = []
    # for item in np.abs(plot_data[0]):
    #     if item >= 0.95 * ampl:
    #         ampl_list.append(item)
    # ampl = round(np.mean(np.array(ampl_list)),2)
    # std_ampl = round(np.std(np.array(ampl_list)),2)
    # std = round(np.std(plot_data),2)

    # freq = 1/(times[1]-times[0])

    # ax1[0].plot(times1[40:],fieldFiltered[40:,0],'C3',label='moving avg.')
    # ax1[1].plot(times1[40:],fieldFiltered[40:,1],'C3',label='moving avg.')
    # ax1[2].plot(times1[40:],fieldFiltered[40:,2],'C3',label='moving avg.')
    
    # ax1[0].legend()
    # ax1[1].legend()
    # ax1[2].legend()
    # p2p_x_rise = np.amax(fields1[:, 0]) - np.amin(fields1[:, 0])
    # p2p_y_rise = np.amax(fields1[:, 1]) - np.amin(fields1[:, 1])
    # p2p_z_rise = np.amax(fields1[:, 2]) - np.amin(fields1[:, 2])
    # field_dev1 = fields1 - np.mean(fields1,axis=0)
    
    # p2p_x_prerise = np.amax(fields1[:,0]) - np.amin(fields1[:,0])
    # p2p_y_prerise = np.amax(fields1[:,1]) - np.amin(fields1[:,1])
    # p2p_z_prerise = np.amax(fields1[:,2]) - np.amin(fields1[:,2])

    # p2p_x_postrise = np.amax(fields1[745:1820, 0]) - np.amin(fields1[745:1820, 0])
    # p2p_y_postrise = np.amax(fields1[745:1820, 1]) - np.amin(fields1[745:1820, 1])
    # p2p_z_postrise = np.amax(fields1[745:1820, 2]) - np.amin(fields1[745:1820, 2])
    
    # p2p_x = np.mean([p2p_x_prerise,p2p_x_postrise])
    # p2p_y = np.mean([p2p_y_prerise,p2p_y_postrise])
    # p2p_z = np.mean([p2p_z_prerise,p2p_z_postrise])
    
    # mag_x = round(np.mean(fields1[:, 0]),2)
    # mag_y = round(np.mean(fields1[:, 1]),2)
    # mag_z = round(np.mean(fields1[:, 2]),2)

    # std_x = np.std(fields1[:,0])
    # std_y = np.std(fields1[:,1])
    # std_z = np.std(fields1[:,2])
    
    # std_x_0 = np.std(fields1[745:1820,0])
    # std_y_0 = np.std(fields1[745:1820,1])
    # std_z_0 = np.std(fields1[745:1820,2])
    
    # print(std_x,', ',std_y,', ',std_z)

    # mag = round(np.sqrt(mag_x ** 2 + mag_y ** 2 + mag_z ** 2),2)
    # theta = round(np.degrees(np.arccos(mag_z/mag)),2)
    # phi = round(np.degrees(np.arctan2(mag_y, mag_x)),2)

    # ax1[0].set_title(f'\n$\Delta B_{{x,pp,max}}$ = {p2p_x_prerise:.2f} $mT$\t$\Delta B_{{x1,RMS}}$ = {std_x:.2f} $mT$'
    #                  f'\n$\Delta B_{{y,pp,max}}$ = {p2p_y_prerise:.2f} $mT$\t$\Delta B_{{y1,RMS}}$ = {std_y:.2f} $mT$'
    #                  f'\n$\Delta B_{{z,pp,max}}$ = {p2p_z_prerise:.2f} $mT$\t$\Delta B_{{z1,RMS}}$ = {std_z:.2f} $mT$')

    # ax1[0].set_xticks([500,1500,2500,3500,4500,5500,6500,7500,8500],minor=True)
    # ax1[0].set_xticks([1000,2000,3000,4000,5000,6000,7000,8000,9000])
    ax1[0].grid()
    ax1[1].grid()
    ax1[2].grid()
    ax1[0].set_title('Temperature measured on coil1 ($T_1$), base ($T_2$) and pole ($T_3$)')#\nwith ECB current 5A in all coils, (0.8,8,8)A and (1,1,8)A')
    
    plt.tight_layout()

    plt.show()
