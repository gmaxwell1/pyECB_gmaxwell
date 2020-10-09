""" 
filename: analysis_tools.py

This file contains functions that are used for plotting data extracted from the Hall Sensor Cube. Mainly for 2D plots of the magnetic field.

Author: Nicholas Meinhardt (Qzabre)
        nmeinhar@student.ethz.ch
        
Date: 09.10.2020
"""

########## Standard library imports ##########
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from datetime import datetime


# %%
def estimate_std_theta(mean_values, std_values):
    """
    Estimate the standard deviation of the estimated angle theta (wrt. z-axis), 
    assuming that x,y,z components of magnetic field are independent variables.

    Args: mean_data_specific_sensor, std_data_specific_sensor are ndarrays of shape (#measurements, 3).

    Returns: ndarray of length #measurements containing standard deviations of theta
    """
    # estimate magnitudes
    mag = np.linalg.norm(mean_values, axis=1)

    # estimate partial derivatives
    deriv_arccos = 1 / np.sqrt(1 - mean_values[:, 2]**2 / mag**2)

    pderiv_x = deriv_arccos * mean_values[:, 2]*mean_values[:, 0] / mag**3
    pderiv_y = deriv_arccos * mean_values[:, 2]*mean_values[:, 1] / mag**3
    pderiv_z = deriv_arccos * (mean_values[:, 2]**2 / mag**3 - 1/mag)

    var_theta = (pderiv_x*std_values[:, 0])**2 + (pderiv_y *
                                                  std_values[:, 1])**2 + (pderiv_z*std_values[:, 2])**2
    return np.sqrt(var_theta)

def estimate_std_phi(mean_values, std_values):
    """
    Estimate the standard deviation of the estimated in-plane angle phi (wrt. x-axis), 
    assuming that x,y,z components of magnetic field are independent variables.

    Args: mean_data_specific_sensor, std_data_specific_sensor are ndarrays of shape (#measurements, 3).

    Returns: ndarray of length #measurements containing standard deviations of phi
    """
    # estimate partial derivatives
    deriv_arctan = np.cos(np.arctan2(mean_values[:,1], mean_values[:,0]))**2

    pderiv_x = deriv_arctan * (-mean_values[:,1]/ mean_values[:,0]**2)
    pderiv_y = deriv_arctan / mean_values[:,0]

    var_phi = (pderiv_x*std_values[:,0])**2 + (pderiv_y*std_values[:,1])**2 
    return np.sqrt(var_phi)

def estimate_std_magnitude(mean_values, std_values):
    """
    Estimate the standard deviation of the estimated magnitude |B|, 
    assuming that x,y,z components of magnetic field are independent variables.

    Args: mean_values, std_values are ndarrays of shape (#measurements, 3).

    Returns: ndarray of length #measurements containing standard deviations of |B|
    """
    # estimate magnitudes
    mag = np.linalg.norm(mean_values, axis=1)

    return np.sqrt(np.einsum('ij,ij->i', mean_values**2, std_values**2)) / mag


def estimate_std_inplane(mean_values, std_values):
    """
    Estimate the standard deviation of the estimated in-plane magnitude |B_xy|, 
    assuming that x,y,z components of magnetic field are independent variables.

    Args: mean_values, std_values are ndarrays of shape (#measurements, 3).

    Returns: ndarray of length #measurements containing standard deviations of |B_xy|
    """
    # estimate magnitudes
    inplane_mag = np.linalg.norm(mean_values[:, 0:2], axis=1)

    return np.sqrt(np.einsum('ij,ij->i', mean_values[:, 0:2]**2, std_values[:, 0:2]**2)) / inplane_mag


def extract_raw_data_from_file(filepath):
    """
    Extract current, mean and std values of measured magnetic field and expected magnetic field from file.

    Also checks whether the current is provided as single float value or as vector
    """
    # extract data and convert to ndarray
    raw_data = pd.read_csv(filepath).to_numpy()

    # current can be single-valued or a vector, differentiate these cases
    dimension_I = len(raw_data[0]) - 9
    if dimension_I == 1:
        I = raw_data[:, 0]
        mean_data_specific_sensor = raw_data[:, 1:4]
        std_data_specific_sensor = raw_data[:, 4:7]
        expected_fields = raw_data[:, 7:10]
    else:
        I = raw_data[:, 0:3]
        mean_data_specific_sensor = raw_data[:, 3:6]
        std_data_specific_sensor = raw_data[:, 6:9]
        expected_fields = raw_data[:, 9:12]

    return I, mean_data_specific_sensor, std_data_specific_sensor, expected_fields

def plot_I_vs_B(I, mean_values, std_values, expected_values, directory, save_image = True, 
                        show_labels=True, ylim=None, xlim=None, image_name_postfix='I_vs_B',
                        remove_first_half=True):
    """
    Generate plot of the currents in all three coils vs the magnetude of the magnetic field.


    Args:   
    - I, mean_values, std_values, expected_values are ndarrays of shape (#measurements, 3), 
    containing applied current, experimentally estimated/expected mean values and standard deviations 
    for x,y,z-directions.
    - directory: valid path of directory where the image should be saved
    - save_image (bool): flag to save or not save the image
    - show_labels (bool): flag to switch on/off labels of plots
    - ylim: tuple containing (ymin, ymax). If None, the default limits are taken
    - xlim: tuple containing (xmin, xmax). If None, the default limits are taken
    - image_name_postfix: The image is saved as '%y_%m_%d_%H-%M-%S_'+ image_name_postfix +'.png'
    - remove_first_half (bool): ignore first half of all input data arrays, 
    which would correspond to a magnetic field pointing towards the negative direction
    compared to the desired direction
    """
    if remove_first_half:
        # if measurements range from negative to positive direction of B-field
        # filter out the first half of all measurements and only keep rest for uniqueness
        number_data = I.shape[0]
        I = I[number_data//2 + 1:,:]
        mean_values = mean_values[number_data//2 + 1:,:]
        std_values = std_values[number_data//2 + 1:,:]
        expected_values = expected_values[number_data//2 + 1:,:]

    # calculate magnitudes of measured and expected magnetic field and their standard deviations
    mean_magnitudes = np.linalg.norm(mean_values, axis=1)
    mean_std_magnitudes = estimate_std_magnitude(mean_values, std_values)
    expected_magnitudes = np.linalg.norm(expected_values, axis=1)

    # initialize colorbar to have measurement and simulation data in a similar color for each coil
    cmap = cm.get_cmap('tab20')

    # generate plots
    fig, ax = plt.subplots()
    for i in range(3):
        ax.errorbar(mean_magnitudes, I[:,i], xerr=mean_std_magnitudes, label = '$I_{}$ measured'.format(i+1),
                        linestyle='', marker='.', capsize = 2, color = cmap(i*0.1))
        ax.plot(expected_magnitudes, I[:,i], linestyle='--', label = '$I_{}$ simulation'.format(i+1),
                        color = cmap(i*0.1+0.05))
    
    # set labels of axes
    ax.set_ylabel('$I$ [A]')
    ax.set_xlabel('$|B|$ [mT]')

    # show legend if desired
    if show_labels:
        ax.legend()

    # set axis limits
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # save image
    if save_image:
        # set the directory name and current datetime if not passed as argument
        if directory is None:
            directory = os.getcwd()
        now = datetime.now().strftime('%y_%m_%d_%H-%M-%S')

        output_file_name = '{}_{}.png'.format(now, image_name_postfix)
        file_path = os.path.join(directory, output_file_name)
        fig.savefig(file_path, dpi=300)

    fig.show()



def generate_plots(I, mean_values, std_values, expected_values, flag_xaxis = 'I1', flags_yaxis = 'zma',
                        plot_delta_sim = False, directory= None, image_name_postfix = 'B_vs_I', 
                        height_per_plot = 2, save_image = True, distance=3.0, xlim = None, 
                        ylim_field_abs = None, ylim_field_z = None, show_labels=True, remove_first_half=True):
    """
    Generate plots of B vs I containing errorbars with mean values and standard deviations. 

    Note that one can decide on which quantities should be plotted on the x and y axes, 
    as well as how many plots should be generated, by adapting the flags_yaxis and flag_xaxis parameters. 

    Args: 
    - I, mean_values, std_values, expected_values are ndarrays of shape (#measurements, 3), 
    containing applied current, experimentally estimated/expected mean values and standard deviations 
    for x,y,z-directions.
    - flag_xaxis (string): Switch between the following quantities displayed on the x-axis:
        - 'I1', 'I2', 'I3' for currents in coils 1,2,3 
        - 'P' for power
        - 'B' for magnitude of magnetic field
    All axes share the same x-axis, and only the lowest plot has a label and ticks (currently)! 
    If invalid (other than the listed letters) flags are passed, the default 'I1' is assumed.
    - flag_yaxis: string containing the letters as flags, where valid letters are mentioned below.
    The number of letters may vary, but at least one valid letter should be contained. For each flag,
    a plot is generated with the according quantity plotted on the y-axis. The plots are generated
    in the same order as the flag-letters are passed. 
        - 'z': Bz-component of magnetic field
        - 'm': magnitude of magnetic field
        - 'p': in-plane magnitude |B_xy| of magnetic field in xy-plane
        - 't': angle theta with respect to z-axis
        - 'f': in-plane angle phi with respect to x-axis
    - plot_delta_sim: If False, both the measured and expected values are plotted. 
    If True, the difference between simulation and measurement is plotted.
    - directory: valid path of directory where the image should be saved
    - image_name_postfix: The image is saved as '%y_%m_%d_%H-%M-%S_'+ image_name_postfix +'.png'
    - height_per_plot: height of each plot in inches. Default for several plots is 2.
    The usual height of a single plot is 4.
    - save_image (bool): flag to save or not save the image
    - distance (float): distance between sensor and tip [mm], which is added as plot label 
    - xlim: tuple containing (xmin, xmax) of plots. If None, the default limits are taken
    - ylim_field_abs, ylim_field_z: tuple containing (ymin, ymax) for the plots of magnetic field,
    where '_abs' indicates magnitudes (only positive) and '_z' indicates the field along z (positive and negative)
    - show_labels (bool): flag to switch on/off labels of plots
    - remove_first_half (bool): ignore first half of all input data arrays, 
    which would correspond to a magnetic field pointing towards the negative direction
    compared to the desired direction

    Return: 
    - x_vals: ndarray of length = #measurements, containing the x-values of all plots
    - plot_mean_data, plot_std_data, plot_expected_data
    which are ndarrays of shape (#plots, #measurements, 1).
    The number of plots (#plots) is set by the length of the flags_yaxis parameter. 
    For each plot, the y-data are contained in plot_mean_data, and the errorbars in plot_std_data. 
    The values that are expected based on previous simulations are returned as plot_expected_data for each plot.
    Note that if plot_delta_sim==False, the actual measurement data are returned. 
    Else, if plot_delta_sim==True, plot_mean_data is the difference between expected and measured data, 
    i.e. plot_mean_data = plot_expected_data - (measurement data).   
    """
    if remove_first_half:
        # if measurements range from negative to positive direction of B-field
        # filter out the first half of all measurements and only keep rest for uniqueness
        number_data = I.shape[0]
        I = I[number_data//2 + 1:,:]
        mean_values = mean_values[number_data//2 + 1:,:]
        std_values = std_values[number_data//2 + 1:,:]
        expected_values = expected_values[number_data//2 + 1:,:]

    # if an empty string is passed as flags_yaxis, set it to 'm'
    if len(flags_yaxis) == 0:
        flags_yaxis = 'm'

    # create a simple plot with as many axes as letters in plots_yaxis
    number_plots = len(flags_yaxis)
    fig, axs = plt.subplots(number_plots, 1, sharex=True)
    fig.set_size_inches(6, number_plots * height_per_plot)

    # if number_plots=1, axs is returned as AxesSubplot class instead of an ndarray containing
    # instances of this class. Since the following requires a ndarray, ensure to have an ndarray!
    if number_plots == 1:
        axs = np.array([axs])

    # calculate magnitudes of measured and expected magnetic field
    mean_magnitudes = np.linalg.norm(mean_values, axis=1)
    expected_magnitudes = np.linalg.norm(expected_values, axis=1)

    # collect plot data.
    # Note: errorbars display std, estimate errors for magnitudes (and angle) using propagation of uncertainty,
    # assuming that the measured fields in x,y,z direction are independent variables
    plot_mean_data = []
    plot_std_data = []
    plot_expected_data = []
    ylabels = []
    for flag in flags_yaxis:
        # magnetic field in z-direction
        if flag == 'z':
            plot_mean_data.append(mean_values[:, 2])
            plot_std_data.append(std_values[:, 2])
            plot_expected_data.append(expected_values[:, 2])
            ylabels.append('$B_z$ [mT]')
        # magnitude of magnetic field
        elif flag == 'm':
            plot_mean_data.append(mean_magnitudes)
            plot_std_data.append(estimate_std_magnitude(mean_values, std_values))
            plot_expected_data.append(expected_magnitudes)
            ylabels.append('$|B|$ [mT]')
        # angle theta (wrt to z-axis)
        elif flag == 't':
            plot_mean_data.append(np.degrees(
                np.arccos(mean_values[:, 2]/mean_magnitudes)))
            plot_std_data.append(np.degrees(
                estimate_std_theta(mean_values, std_values)))
            plot_expected_data.append(np.degrees(
                np.arccos(expected_values[:, 2]/expected_magnitudes)))
            ylabels.append('$\\theta$ [°]')
        # angle phi (wrt to x-axis)
        elif flag == 'f':
            plot_mean_data.append(np.degrees(np.arctan2(mean_values[:,1], mean_values[:,0])))
            plot_std_data.append(np.degrees(estimate_std_phi(mean_values, std_values)))
            plot_expected_data.append(np.degrees(np.arctan2(expected_values[:,1], expected_values[:,0])))
            ylabels.append('$\\phi$ [°]')
        # in-plane component (along xy) of magnetic field
        elif flag == 'p':
            plot_mean_data.append(
                np.sqrt(mean_values[:, 0]**2 + mean_values[:, 1]**2))
            plot_std_data.append(estimate_std_inplane(mean_values, std_values))
            plot_expected_data.append(
                np.sqrt(expected_values[:, 0]**2 + expected_values[:, 1]**2))
            ylabels.append('$|B_{xy}|$ [mT]')
        # account for invalid flags:
        else:
            raise ValueError(
                '{} is not a valid flag, it should be in [\'z\', \'m\', \'a\', \'p\']!'.format(flag))

    # plot current ('I'), power ('P') or magnitude of magnetic field ('B') on xaxis, depending on flag_xaxis
    if flag_xaxis == 'P':
        R = 0.47    # resistance [Ohm] of each coil
        x_vals = R * np.sum(I**2, axis=1)  # sum over all three coils
        axs[-1].set_xlabel('$P$ [W]')
    elif flag_xaxis == 'B':
        x_vals = mean_magnitudes
        axs[-1].set_xlabel('$|B|$ [mT]')
    elif flag_xaxis == 'I2':
        x_vals = I[:, 1]
        axs[-1].set_xlabel('$I_2$ [A]')
    elif flag_xaxis == 'I3':
        x_vals = I[:, 2]
        axs[-1].set_xlabel('$I_3$ [A]')
    else:
        x_vals = I[:, 0]
        axs[-1].set_xlabel('$I_1$ [A]')

    # actual plotting
    for i in range(len(axs)):
        axs[i].set_ylabel(ylabels[i])

        # plot either both measured and expected values or only the differences between both.
        if plot_delta_sim:
            axs[i].errorbar(x_vals, plot_expected_data[i]-plot_mean_data[i], yerr=plot_std_data[i],
                            linestyle='', marker='.', capsize=2,
                            label='$\\Delta$ = simulation-measurements')
        else:
            axs[i].errorbar(x_vals, plot_mean_data[i], yerr=plot_std_data[i],
                                    linestyle='', marker='.', capsize = 2, 
                                    label = 'measured @ {:.1f} mm'.format(distance))
            axs[i].plot(x_vals, plot_expected_data[i], linestyle='--', 
                                    label = 'simulation @ 3 mm')
            if show_labels:
                axs[i].legend()

    # add a Delta at the front of each label if differences should be plotted
    if plot_delta_sim:
        ylabels = ['$\\Delta$ {}'.format(label) for label in ylabels]

    # settings that are different for plots of angle and plots of field
    
    
    for i in range(len(flags_yaxis)):
        if flags_yaxis[i] == 't':
            # if the angle is plotted on one of the axes, set the limits to (-5, 185)
            if not plot_delta_sim:
                axs[i].set_ylim((-5,185))
                axs[i].yaxis.set_ticks(np.arange(0, 210, 30))
        elif flags_yaxis[i] == 'f':
            # if the angle is plotted on one of the axes, set the limits to (-5, 365)
            if not plot_delta_sim:
                axs[i].set_ylim((-190,190))
                axs[i].yaxis.set_ticks(np.arange(-180, 240, 60))
        elif flags_yaxis[i] == 'z':
            axs[i].set_ylim(ylim_field_z)
        else:
            axs[i].set_ylim(ylim_field_abs)
    
    # set limits for x-axis 
    axs[-1].set_xlim(xlim)

    # further adjustments
    plt.tight_layout()

    # save image
    if save_image:
        # set the directory name and current datetime if not passed as argument
        if directory is None:
            directory = os.getcwd()
        now = datetime.now().strftime('%y_%m_%d_%H-%M-%S')

        output_file_name = '{}_{}.png'.format(now, image_name_postfix)
        file_path = os.path.join(directory, output_file_name)
        fig.savefig(file_path, dpi=300)

    fig.show()

    return x_vals, plot_mean_data, plot_std_data, plot_expected_data


def sigmoid(x, k, a):
    """
    Sigmoid function with growth rate k and maximum value a, rescaled such that center point is at origin.
    """
    return a * (1 - np.exp(-k*x)) / (1 + np.exp(-k*x))


def abs_sigmoid(x, k, a):
    """
    Return absolute value of sigmoid function with growth rate k and maximum value a, 
    rescaled such that center point is at origin.
    """
    return np.abs(sigmoid(x, k, a))


def coth(x):
    """Return cosh(x)/sinh(x)"""
    return np.cosh(x) / np.sinh(x)


def brillouin_fct(x, J, a):
    """
    Implement the Brillouin function, which is used to describe paramagnet. 
    """
    s = 1/(2*J)
    return a * ((1+s) * coth((1+s)*x) - s * coth(s*x))


def abs_brillouin_fct(x, J, a):
    """
    Implement the absolute value of Brillouin function, which is used to describe paramagnet. 
    """
    return np.abs(brillouin_fct(x, J, a))


def lin_and_const(x, x_kink, a):
    """
    Implement a function that raises linearly as a*x until |x|= x_kink and remains constant afterwards.

    Note that x_kink >= 0 is required.
    """
    try:
        len(x)
    except TypeError:
        return a*x if np.abs(x) <= x_kink else a*x_kink*np.sign(x)
    else:
        return np.array([a*x[i] if np.abs(x[i]) <= x_kink else a*x_kink*np.sign(x[i]) for i in range(len(x))])


def abs_lin_and_const(x, x_kink, a):
    """
    Implement a function that raises linearly as |a*x| until |x|= x_kink and remains constant afterwards.

    Note that x_kink >= 0 is required and that the returned value is >= 0.
    """
    return np.abs(lin_and_const(x, x_kink, a))

