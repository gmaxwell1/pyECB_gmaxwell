# Author Nicholas Meinhardt 2020

#%%
import numpy as np
import pandas as pd 
import os
import matplotlib.pyplot as plt
from datetime import datetime


#%%
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
    deriv_arccos = 1/ np.sqrt(1- mean_values[:,2]**2 / mag**2)

    pderiv_x = deriv_arccos * mean_values[:,2]*mean_values[:,0] / mag**3
    pderiv_y = deriv_arccos * mean_values[:,2]*mean_values[:,1] / mag**3
    pderiv_z = deriv_arccos * (mean_values[:,2]**2 / mag**3 - 1/mag)

    var_theta = (pderiv_x*std_values[:,0])**2 + (pderiv_y*std_values[:,1])**2 + (pderiv_z*std_values[:,2])**2
    return np.sqrt(var_theta)

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
    inplane_mag = np.linalg.norm(mean_values[:,0:2], axis=1)

    return np.sqrt(np.einsum('ij,ij->i', mean_values[:,0:2]**2, std_values[:,0:2]**2)) / inplane_mag

def extract_raw_data_from_file(filepath):
    """
    Extract current, mean and std values of measured magnetic field and expected magnetic field from file.
    """
    raw_data = pd.read_csv(filepath).to_numpy()
    I = raw_data[:,0]
    mean_data_specific_sensor = raw_data[:,1:4]
    std_data_specific_sensor = raw_data[:,4:7]
    expected_fields = raw_data[:,7:10]

    return I, mean_data_specific_sensor, std_data_specific_sensor, expected_fields

def generate_plots(I, mean_values, std_values, expected_values, flag_xaxis = 'I', flags_yaxis = 'zma',
                        plot_delta_sim = False, directory= None, data_filename_postfix = 'B_vs_I', 
                        height_per_plot = 2, save_image = True):
    """
    Generate plots of B vs I containing errorbars with mean values and standard deviations. 
    
    Note that one can decide on which quantities should be plotted on the x and y axes, 
    as well as how many plots should be generated, by adapting the flags_yaxis and flag_xaxis parameters. 

    Args: 
    - I, mean_values, std_values, expected_values are ndarrays of shape (#measurements, 3), 
    containing applied current, experimentally estimated/expected mean values and standard deviations 
    for x,y,z-directions.
    - flag_xaxis (string): Switch between the following quantities displayed on the x-axis:
        - 'I' for current, 
        - 'P' for power, 
        - 'B' for magnitude of magnetic field.
    All axes share the same x-axis, and only the lowest plot has a label and ticks (currently)! 
    If invalid (other than the listed letters) flags are passed, the default 'I' is assumed.
    - flag_yaxis: string containing the letters as flags, where valid letters are mentioned below.
    The number of letters may vary, but at least one valid letter should be contained. For each flag,
    a plot is generated with the according quantity plotted on the y-axis. The plots are generated
    in the same order as the flag-letters are passed. 
        - 'z': Bz-component of magnetic field
        - 'm': magnitude of magnetic field
        - 'p': in-plane magnitude |B_xy| of magnetic field in xy-plane
        - 'a': angle theta with respect to z-axis
    - plot_delta_sim: If False, both the measured and expected values are plotted. 
    If True, the difference between simulation and measurement is plotted.
    - directory: valid path of directory where the image should be saved
    - data_filename_postfix: The image is saved as '%y_%m_%d_%H-%M-%S_'+ data_filename_postfix +'.png'
    - height_per_plot: height of each plot in inches. Default for several plots is 2.
    The usual height of a single plot is 4.
    - save_image: flag to save or not save the image

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
    # if an empty string is passed as flags_yaxis, set it to 'm'
    if len(flags_yaxis) == 0:
        flags_yaxis = 'm'

    # create a simple plot with as many axes as letters in plots_yaxis
    number_plots = len(flags_yaxis)
    fig, axs = plt.subplots(number_plots, 1, sharex=True)
    fig.set_size_inches(6, number_plots * height_per_plot)

    # if number_plots=1, axs is returned as AxesSubplot class instead of an ndarray containing 
    # instances of this class. Since the following requires a ndarray, ensure to have an ndarray!
    if number_plots ==1:
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
            plot_mean_data.append(mean_values[:,2])
            plot_std_data.append(std_values[:,2])
            plot_expected_data.append(expected_values[:,2])
            ylabels.append('$B_z$ [mT]')
        # magnitude of magnetic field
        elif flag == 'm':
            plot_mean_data.append(mean_magnitudes)
            plot_std_data.append(estimate_std_magnitude(mean_values, std_values))
            plot_expected_data.append(expected_magnitudes)
            ylabels.append('$|B|$ [mT]')
        # angle theta (wrt to z-axis)
        elif flag == 'a':
            plot_mean_data.append(np.degrees(np.arccos(mean_values[:,2]/mean_magnitudes)))
            plot_std_data.append(np.degrees(estimate_std_theta(mean_values, std_values)))
            plot_expected_data.append( np.degrees(np.arccos(expected_values[:,2]/expected_magnitudes)))
            ylabels.append('$\\theta$ [°]')
        # in-plane component (along xy) of magnetic field
        elif flag == 'p':
            plot_mean_data.append(np.sqrt(mean_values[:,0]**2 + mean_values[:,1]**2))
            plot_std_data.append(estimate_std_inplane(mean_values, std_values))
            plot_expected_data.append(np.sqrt(expected_values[:,0]**2 + expected_values[:,1]**2))
            ylabels.append('$|B_{xy}|$ [mT]')
        # account for invalid flags:
        else:
            raise ValueError('{} is not a valid flag, it should be in [\'z\', \'m\', \'a\', \'p\']!'.format(flag))
    

    # plot current ('I'), power ('P') or magnitude of magnetic field ('B') on xaxis, depending on flag_xaxis
    if flag_xaxis == 'P':
        R = 0.47    # resistance [Ohm] of each coil]
        x_vals = R * I**2
        x_vals *= 3 # multiply by the number of coils to get overall power
        axs[-1].set_xlabel('$P$ [W]')
    elif flag_xaxis == 'B':
        x_vals = mean_magnitudes
        axs[-1].set_xlabel('$|B|$ [mT]')
    else:
        x_vals = I
        axs[-1].set_xlabel('$I$ [A]')

    # actual plotting 
    for i in range(len(axs)):
        axs[i].set_ylabel(ylabels[i])

        # plot either both measured and expected values or only the differences between both. 
        if plot_delta_sim:
            axs[i].errorbar(x_vals, plot_expected_data[i]-plot_mean_data[i], yerr=plot_std_data[i],
                                linestyle='', marker='.', capsize = 2, label = '$\\Delta$ = simulation-measurements')
        else:
            axs[i].errorbar(x_vals, plot_mean_data[i], yerr=plot_std_data[i],
                                    linestyle='', marker='.', capsize = 2, label = 'measured @ ~ 3.5 mm')
            axs[i].plot(x_vals, plot_expected_data[i], linestyle='--', marker='.', label = 'simulation @ 3 mm')
            axs[i].legend()

    # add a Delta at the front of each label if differences should be plotted
    if plot_delta_sim:
        ylabels = ['$\\Delta$ {}'.format(label) for label in ylabels]

    # if the angle is plotted on one of the axes, set the limits to (-5, 185)
    if not plot_delta_sim:
        for i in range(len(flags_yaxis)):
            if flags_yaxis[i] == 'a':
                axs[i].set_ylim((-5,185))
                axs[i].yaxis.set_ticks(np.arange(0, 210, 30))
    
    # further adjustments
    plt.tight_layout()

    # save image
    if save_image :
        # set the directory name and current datetime if not passed as argument
        if directory is None:
            directory = os.getcwd()
        now = datetime.now().strftime('%y_%m_%d_%H-%M-%S')

        output_file_name = '{}_{}.png'.format(now, data_filename_postfix) 
        file_path = os.path.join(directory, output_file_name)
        fig.savefig(file_path, dpi = 300)
    
    fig.show()

    return x_vals, plot_mean_data, plot_std_data, plot_expected_data
