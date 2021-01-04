""" 
filename: analysis_tools.py

This file contains functions that are used for plotting data extracted from the Hall Sensor Cube. 
Mainly for 2D plots of the magnetic field.

Author: Nicholas Meinhardt (Qzabre)
        nmeinhar@student.ethz.ch
        
Date: 09.10.2020
"""

########## Standard library imports ##########
import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, MaxNLocator
from datetime import datetime
from itertools import product
from scipy.optimize import curve_fit
from scipy.spatial.transform import Rotation as R

from modules.general_functions import angle_wrt_z, inplane_angle_wrt_x
from modules.interpolation_tools import find_start_of_saturation
from modules.general_functions import estimate_RMS_error


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
    return np.degrees(np.sqrt(var_theta))

def estimate_std_phi(mean_values, std_values):
    """
    Estimate the standard deviation of the estimated in-plane angle phi (wrt. x-axis) in degrees, 
    assuming that x,y,z components of magnetic field are independent variables.

    Args: mean_data_specific_sensor, std_data_specific_sensor are ndarrays of shape (#measurements, 3).

    Returns: ndarray of length #measurements containing standard deviations of phi in degrees
    """
    # estimate partial derivatives
    deriv_arctan = np.cos(np.arctan2(mean_values[:,1], mean_values[:,0]))**2

    pderiv_x = deriv_arctan * (-mean_values[:,1]/ mean_values[:,0]**2)
    pderiv_y = deriv_arctan / mean_values[:,0]

    var_phi = (pderiv_x*std_values[:,0])**2 + (pderiv_y*std_values[:,1])**2 
    return np.degrees(np.sqrt(var_phi))

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

def extract_raw_data_from_2d_scan(filepath):
    """
    Extract positions and measured magnetic fields from file that was created after 2d-scan.

    Returns: positions, B_field (both are arrays of shape (#points, 3))
    """
    # extract data and convert to ndarray
    raw_data = pd.read_csv(filepath).to_numpy()

    # current can be single-valued or a vector, differentiate these cases
    dimension_I = len(raw_data[0]) - 9
    
    positions = raw_data[:, 0:3]
    B_field = raw_data[:, 3:6]


    return positions, B_field

def plot_2d_scan(positions, B_field, origin=None, Cont=False, Scat_Mag=False, Show=True,
             title_on=True, cmin=None, cmax=None, plot_component='m', center_position=None, levels=None):
    """
    Generate 3d plot of 2d-sweep data in the sensor coordinate system.

    Note that sensor and stage coordinate systems differ, the provided positions are expected to be in 
    sensor coordinates already.

    Args:
    - positions_corrected (ndarray) is of shape ((grid_number+1)**2, 3) and contains the stage positions 
    in sensor coordinates
    - B_field  (ndarray): is of shape ((grid_number+1)**2, 3) containing the measured fields
    which is the final distance between gird points
    - origin (None or 1d-array of length 3): 
    If None, the origin is set automatically such that the x and y data range between 0 and a positive number.
    If an array is passed as offset position, note that a position (p.e. the chosen center position) has to be 
    transformed to sensor coordinates before passing it as argument. 
    - cmin, cmax (float): min and max of colormap for displaying magnetic field strength. 
    If both any of them is None, the scale is chosen automatically based on the data
    - Scat_Mag (bool): flag to switch on/off scatter plot (only points) of magnetic field
    - Cont (bool): flag to switch on/off plotting contours for z=const, i.e. flat surfaces with interpolation 
    between data points and the surfaces are stacked on top of each other.
    - Show (bool): flag to switch on/off plt.show() command 
    - plot_component (string): flag to choose how data should be plotted, possible values are:
    'm' (or any invalid value): euclidian norm/magnitude of magnetic field vector. 
    'x', 'y', 'z': plot the according vector component.
    'xy': euclidian norm of in-plane component of field (i.e. x,y components)
    'theta': angle wrt z-axis in degrees
    'phi': inplane angle wrt to x-axis in degrees
    - center_position (None or 1d array of length 2): If an array [x,y] is provided, a red dot is plotted
    at this position. The passed positions are in stage coordinates and will be transformed to 
    sensor coordinates. 
    - levels (int or 1d array of increasing numbers): level parameter passed to contourf, 
    sets the number of color levels

    Return: fig, ax: Figure and Axis objects of the created plot. 
    """
    # if origin is None, set it automatically. Note that x,y positions in stage coordinates (which are positive)
    # are transformed to sensor coordinates, resulting in negative coordinates.
    if origin is None:
        origin = -1 * np.array([np.min(positions[:,0]), np.min(positions[:,1]), np.min(positions[:,2])])

    # store positions in arrays x and y, which can be used to generate plots later
    points_per_side = int(np.sqrt(len(positions)))
    x = np.zeros((points_per_side, points_per_side))
    y = np.zeros((points_per_side, points_per_side))
    mag = np.zeros((points_per_side, points_per_side))

    # for i in range(points_per_side):
    #     for j in range(points_per_side):
    #         x[j,i] = positions[i+points_per_side*j, 0] + origin[0]
    #         y[j,i] = positions[i+points_per_side*j, 1] + origin[1]
    x = positions[:,0].reshape((points_per_side, points_per_side)) + origin[0]
    y = positions[:,1].reshape((points_per_side, points_per_side)) + origin[1]
    
    # choose which data should be plotted. 
    # abbriviate notation using itertools.product
    if plot_component == 'x':
        label = 'field component $B_x$ [mT]'
        for i,j in product(range(points_per_side), range(points_per_side)):
            mag[j,i] = B_field[i+points_per_side*j, 0]
    elif plot_component == 'y':
        label = 'field component $B_y$ [mT]'
        for i,j in product(range(points_per_side), range(points_per_side)):
            mag[j,i] = B_field[i+points_per_side*j, 1]
    elif plot_component == 'z':
        label = 'field component $B_z$ [mT]'
        for i,j in product(range(points_per_side), range(points_per_side)):
            mag[j,i] = B_field[i+points_per_side*j, 2]
    elif plot_component == 'xy':
        label = 'in-plane magnitude $|B_{xy}|$ [mT]'
        for i,j in product(range(points_per_side), range(points_per_side)):
            mag[j,i] = np.sqrt(B_field[i+points_per_side*j, 0]**2 + B_field[i+points_per_side*j, 1]**2)
    elif plot_component == 'theta':
        label = 'angle wrt. z-axis $\\theta$ [°]'
        for i,j in product(range(points_per_side), range(points_per_side)):
            mag[j,i] = np.degrees(angle_wrt_z(B_field[i+points_per_side*j, :]))
    elif plot_component == 'phi':
        label = 'in-plane angle wrt. x-axis $\phi$ [°]'
        for i,j in product(range(points_per_side), range(points_per_side)):
            mag[j,i] = np.degrees(inplane_angle_wrt_x(B_field[i+points_per_side*j, :]))
    else:
        label = 'total magnitude $|B|$ [mT]'
        for i,j in product(range(points_per_side), range(points_per_side)):
            mag[j,i] = np.linalg.norm(B_field[i+points_per_side*j,:])
            
    # due to the structure of the 2d-grid run, the direction of every second sweep along x is reversed.
    # for the scatter plot, this does not matter, while it makes a difference for contour plot. 
    # account for this by inverting order
    for j in range(points_per_side):
        if j % 2 == 1:
            x[j,:] = np.flip(x[j,:])
            y[j,:] = np.flip(y[j,:])
            mag[j,:] = np.flip(mag[j,:])

    # set limits for colormap if they are not provided yet
    if cmin == None or cmax == None:
        cmin = np.amin(mag) 
        cmax = np.amax(mag) 

    # create figure 
    fig, ax = plt.subplots()

    # countour plot
    if Cont:
        cf = ax.contourf(x, y, mag, vmin=cmin, vmax=cmax, levels=levels)
        fig.colorbar(cf, ax=ax, boundaries=(cmin, cmax), label=label)
    # scatter plot
    if Scat_Mag:  
        mag = mag.flatten()
        cf = ax.scatter(x, y, mag, c=mag, vmin=cmin, vmax=cmax)
        fig.colorbar(cf, ax=ax, label=label)

    # if provided, plot the center position after converting it to sensor coordintates
    if center_position is not None:
        center_position = [-center_position[1], -center_position[0]]
        ax.plot(center_position[0]+ origin[0], center_position[1]+ origin[1], 
                                marker='+', color='r', label='determined center')  

    # final axis settings
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")

    return fig, ax

def plot_I_vs_B(I, mean_values, std_values, expected_values, directory, save_image = True, 
                        show_labels=True, ylim=None, xlim=None, image_name_postfix='I_vs_B',
                        remove_half=0):
    """
    Generate plot of the currents in all three coils vs the magnetude of the magnetic field.

    Args:   
    - I, mean_values, std_values, expected_values are ndarrays of shape (#measurements, 3), 
    containing applied current, experimentally estimated/expected mean values and standard deviations 
    for x,y,z-directions.
    - directory: valid path of directory where the image should be saved
    - save_image (bool): flag to save or not save the image
    - show_labels (bool): flag to switch on/off labels of plots
    - ylim: tuple of floats containing (ymin, ymax). If None, the default limits are taken
    - xlim: tuple of floats containing (xmin, xmax). If None, the default limits are taken
    - image_name_postfix: The image is saved as '%y_%m_%d_%H-%M-%S_'+ image_name_postfix +'.png'
    - remove_half (int): if 1, ignore first half of all input data arrays, 
    which would correspond to a magnetic field pointing towards the negative direction
    compared to the desired direction. If 2, ignore second half. Default is 0, 
    no data are ignored.
    """
    number_data = I.shape[0]
    if remove_half == 1:
        # if measurements range from negative to positive direction of B-field
        # filter out the first half of all measurements and only keep rest for uniqueness
        I = I[number_data//2 + 1:,:]
        mean_values = mean_values[number_data//2 + 1:,:]
        std_values = std_values[number_data//2 + 1:,:]
        expected_values = expected_values[number_data//2 + 1:,:]
    elif remove_half == 2:
        # if measurements range from negative to positive direction of B-field
        # filter out the second half of all measurements and only keep rest for uniqueness
        I = I[:number_data//2 + 1,:]
        mean_values = mean_values[:number_data//2 + 1,:]
        std_values = std_values[:number_data//2 + 1,:]
        expected_values = expected_values[:number_data//2 + 1,:]

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

    plt.show()

def get_phi(values, cut_phi_at_0=False):
    """
    Return the in-plane angle phi with respect to x-axis in degrees.
    """
    angles = np.degrees(np.arctan2(values[...,1], values[...,0]))

    # if cut should be at 0 degrees, add 360 degrees to all negative values
    if cut_phi_at_0:
        mask = angles < 0
        angles[mask] = 360 + angles[mask]

    return angles

def get_theta(values):
    """
    Return the angle theta between z axis and the provided values in degrees.
    """
    mag = np.linalg.norm(values, axis=-1)
    angles = np.degrees(np.arccos(values[...,2]/mag))

    return angles

def generate_plots(I, mean_values, std_values, expected_values, flag_xaxis = 'I1', flags_yaxis = 'zmt',
                        plot_delta_sim = False, directory = None, image_name_postfix = 'B_vs_I', 
                        height_per_plot = 2, save_image = True, distance = 3.0, xlim = None, 
                        ylim_field_abs = None, ylim_field_single = None, ylim_theta = None, ylim_phi = None,
                        show_labels = True, remove_half = 0, ygrid = False, cut_phi_at_0 = False,
                        show_image = True, show_expected_data=True):
    """
    Generate plots of B vs I containing errorbars with mean values and standard deviations. 

    Note that one can decide on which quantities should be plotted on the x and y axes, 
    as well as how many plots should be generated, by adapting the flags_yaxis and flag_xaxis parameters. 

    Args: 
    - I, mean_values, std_values, expected_values (ndarrays) are of shape (#measurements, 3), 
    containing applied current, experimentally estimated/expected mean values and standard deviations 
    for x,y,z-directions.
    - flag_xaxis (string): Switch between the following quantities displayed on the x-axis:
        - 'I1', 'I2', 'I3' for currents in coils 1,2,3 
        - 'P' for power
        - 'B' for magnitude of magnetic field
    All axes share the same x-axis, and only the lowest plot has a label and ticks (currently)! 
    If invalid (other than the listed letters) flags are passed, the default 'I1' is assumed.
    - flag_yaxis (string): contains letters as flags, where valid letters are mentioned below.
    The number of letters may vary, but at least one valid letter should be contained. For each flag,
    a plot is generated with the according quantity plotted on the y-axis. The plots are generated
    in the same order as the flag-letters are passed. 
        - 'x', 'y', 'z': single component of magnetic field, p.e. B_z
        - 'm': magnitude of magnetic field
        - 'p': in-plane magnitude |B_xy| of magnetic field in xy-plane
        - 't': angle theta with respect to z-axis
        - 'f': in-plane angle phi with respect to x-axis
    - plot_delta_sim: If False, both the measured and expected values are plotted. 
    If True, the difference between simulation and measurement is plotted.
    - directory (string): valid path of directory where the image should be saved
    - image_name_postfix (string): The image is saved as '%y_%m_%d_%H-%M-%S_'+ image_name_postfix +'.png'
    - height_per_plot (float): height of each plot [inches]. Default for several plots is 2.
    The usual height of a single plot is 4.
    - save_image (bool): flag to save or not save the image
    - distance (float): distance between sensor and tip [mm], which is added as plot label 
    - xlim (float): tuple containing (xmin, xmax) of plots. If None, the default limits are taken
    - ylim_field_abs, ylim_field_z, ylim_theta, ylim_phi: tuple of floats containing (ymin, ymax) 
    for the plots of magnetic field or angle, where '_abs' indicates magnitudes (only positive),
    '_z' indicates the field along z (positive and negative) and '_phi' and '_theta' the respective angles
    - show_labels (bool): flag to switch on/off labels of plots
    - remove_half (int): if 1, ignore first half of all input data arrays, 
    which would correspond to a magnetic field pointing towards the negative direction
    compared to the desired direction. If 2, ignore second half. Default is 0, 
    no data are ignored.
    - ygrid (bool): flag for switching on/off grid lines for y-axis
    - cut_phi_at_0 (bool): if True, the values of phi are in [0,360], else in [-180, +180]
    - show_image (bool): if True, plt.show() is added at the end
    - show_expected_data (bool): if True, the expected data are plotted as dashed orange line. 

    Return: fig, axs, x_vals, plot_mean_data, plot_std_data, plot_expected_data

    - fig, axs: plt.Figure object and array of plt.Axes objects that represent the figure and subplots, respectively
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
    if remove_half==1:
        # if measurements range from negative to positive direction of B-field
        # filter out the first half of all measurements and only keep rest for uniqueness
        number_data = I.shape[0]
        I = I[number_data//2 + 1:,:]
        mean_values = mean_values[number_data//2 + 1:,:]
        std_values = std_values[number_data//2 + 1:,:]
        expected_values = expected_values[number_data//2 + 1:,:]
    elif remove_half == 2:
        # if measurements range from negative to positive direction of B-field
        # filter out the second half of all measurements and only keep rest for uniqueness
        number_data = I.shape[0]
        I = I[:number_data//2 + 1,:]
        mean_values = mean_values[:number_data//2 + 1,:]
        std_values = std_values[:number_data//2 + 1,:]
        expected_values = expected_values[:number_data//2 + 1,:]

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
        # magnetic field in x-direction
        if flag == 'x':
            plot_mean_data.append(mean_values[:, 0])
            plot_std_data.append(std_values[:, 0])
            plot_expected_data.append(expected_values[:, 0])
            ylabels.append('$B_x$ [mT]')
        # magnetic field in y-direction
        elif flag == 'y':
            plot_mean_data.append(mean_values[:, 1])
            plot_std_data.append(std_values[:, 1])
            plot_expected_data.append(expected_values[:, 1])
            ylabels.append('$B_y$ [mT]')
        # magnetic field in z-direction
        elif flag == 'z':
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
            plot_mean_data.append(get_theta(mean_values))
            plot_std_data.append(estimate_std_theta(mean_values, std_values))
            plot_expected_data.append(get_theta(expected_values))
            ylabels.append('$\\theta$ [°]')
        # angle phi (wrt to x-axis)
        elif flag == 'f':
            plot_mean_data.append(get_phi(mean_values, cut_phi_at_0=cut_phi_at_0))
            plot_std_data.append(estimate_std_phi(mean_values, std_values))
            plot_expected_data.append(get_phi(expected_values, cut_phi_at_0=cut_phi_at_0))
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
                '{} is not a valid flag, it should be in [\'x\', \'y\', \'z\', \'m\', \'a\', \'p\']!'.format(flag))

    # plot current ('I'), power ('P') or magnitude of magnetic field ('B') on xaxis, depending on flag_xaxis
    if flag_xaxis == 'P':
        R = 0.47    # resistance [Ohm] of each coil
        x_vals = R * np.sum(I**2, axis=1)  # sum over all three coils
        axs[-1].set_xlabel('overall power $P$ [W]')
    elif flag_xaxis == 'B':
        x_vals = mean_magnitudes
        axs[-1].set_xlabel('total magnitude $|B|$ [mT]')
    elif flag_xaxis == 'I2':
        x_vals = I[:, 1]
        axs[-1].set_xlabel('current in coil 2, $I_2$ [A]')
    elif flag_xaxis == 'I3':
        x_vals = I[:, 2]
        axs[-1].set_xlabel('current in coil 3, $I_3$ [A]')
    else:
        x_vals = I[:, 0]
        axs[-1].set_xlabel('current in coil 1, $I_1$ [A]')

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
            if show_expected_data:
                axs[i].plot(x_vals, plot_expected_data[i], linestyle='--', 
                                    label = 'linear model @ 3mm')
    if show_labels:
        axs[0].legend()

    # add a Delta at the front of each label if differences should be plotted
    if plot_delta_sim:
        ylabels = ['$\\Delta$ {}'.format(label) for label in ylabels]

    # settings that are different for plots of angle and plots of field
    for i in range(len(flags_yaxis)):
        if flags_yaxis[i] == 't':
            # if the angle is plotted on one of the axes, set the limits to (-5, 185)
            if not plot_delta_sim:
                if ylim_theta is None:
                    data_min = np.min(plot_mean_data[i])
                    data_max = np.max(plot_mean_data[i])
                    ylim_lower = -0.05 if data_min < 10 else 0.9 * data_min
                    ylim_upper =  1.1 * data_max
                    axs[i].set_ylim(ylim_lower, ylim_upper)
                else:
                    axs[i].set_ylim((-5,185))
                    axs[i].yaxis.set_ticks(np.arange(0, 210, 30))
            max_angle_diff = np.max(plot_mean_data[i]) - np.min(plot_mean_data[i]) 
            # if max_angle_diff > 60:
            #     axs[i].yaxis.set_major_locator(MultipleLocator(30))
            # elif max_angle_diff > 30:

            axs[i].yaxis.set_major_locator(MaxNLocator(nbins='auto', steps=[1,3,6,10]))
            axs[i].yaxis.set_minor_locator(AutoMinorLocator(3))
            
        elif flags_yaxis[i] == 'f':
            # if the angle is plotted on one of the axes, set the limits to (-5, 365)
            if not plot_delta_sim:
                if ylim_phi is None:
                    data_min = np.min(plot_mean_data[i])
                    data_max = np.max(plot_mean_data[i])
                    ylim_lower = 1.1 * data_min if data_min < 0 else 0.9 * data_min
                    ylim_upper = 0.9 * data_max if data_max < 0 else 1.1 * data_max
                    axs[i].set_ylim(ylim_lower, ylim_upper)
                else:
                    axs[i].set_ylim((-190,190))
                    axs[i].yaxis.set_ticks(np.arange(-180, 240, 60))
            # axs[i].yaxis.set_major_locator(MultipleLocator(90))
            axs[i].yaxis.set_major_locator(MaxNLocator(nbins='auto', steps=[1,3,6,10]))
            axs[i].yaxis.set_minor_locator(AutoMinorLocator(3))

        elif flags_yaxis[i] in ['x', 'y', 'z']:
            # either match limits and data or take the passed argument
            if (ylim_field_single is None) and (not plot_delta_sim):
                data_min = np.min(plot_mean_data[i])
                data_max = np.max(plot_mean_data[i])
                ylim_lower = 1.1 * data_min if data_min < 0 else 0.9 * data_min
                ylim_upper = 0.9 * data_max if data_max < 0 else 1.1 * data_max
                axs[i].set_ylim(ylim_lower, ylim_upper)
            else:
                axs[i].set_ylim(ylim_field_single)
            axs[i].yaxis.set_minor_locator(AutoMinorLocator())
        else:
            # either match limits and data or take the passed argument
            if (ylim_field_abs is None) and (not plot_delta_sim):
                data_min = np.min(plot_mean_data[i])
                data_max = np.max(plot_mean_data[i])
                ylim_lower = -0.05 if data_min < 1 else 0.9 * data_min
                ylim_upper = 1.05* data_max
                axs[i].set_ylim(ylim_lower, ylim_upper)
            else:
                axs[i].set_ylim(ylim_field_abs)
            axs[i].yaxis.set_minor_locator(AutoMinorLocator())


        # add grid lines if desired
        if ygrid:
            axs[i].yaxis.grid()
            

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

    if show_image:
        plt.show()

    return fig, axs, x_vals, np.array(plot_mean_data), np.array(plot_std_data), np.array(plot_expected_data)

def add_annotations(axs, x_vals, plot_data, flags_yaxis, mask_selected_plots = None, 
                fontweight = 'bold', fraction_arrow_length = 0.18, fraction_alignment = 0.05):
    """ 
    Add annotations for the maximum and minimum values as arrows pointing at these values.

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
    - mask_selected_plots (ndarray of bools or None): Mask of the same length as axs, only for the 
    subplots where mask_selected_plots is True annotations are added. This allows to avoid annotations
    for certain subplots. Default is None, in which case annotations are added for all subplots. 
    - fontweight (str): valid argument for fontweight of the label containing the actual value of a 
    maximum (minimum), is passed to plt.annotate()
    - fraction_arrow_length (float): ratio of arrow length and the plotted yrange in data coordinates
    - fraction_alignment (float): if a maximum (minimum) value is close to the left or right side of the 
    plot, the label is aligned to the left/right, such that it does not overlap with the axes. 
    Here, fraction_alignment determines what 'close' means: if a current is within a distance of
    fraction_alignment*xrange from the axes, the alignment is set to right or left. 
    In all other cases, the alignment is set to 'center'.
    """
    if mask_selected_plots is None:
        mask_selected_plots = np.ones(len(axs), dtype=bool)

    for i in range(len(axs)):
        # for plots of total or in-plane magnetudes, only add one arrow where the maximum value is reached
        if mask_selected_plots[i] and flags_yaxis[i] in ['p','m']:
            # find the maximum value first
            i_max = np.argmax(plot_data[i])
            max_B_val = plot_data[i, i_max]
            corresponding_current = x_vals[i_max]
            
            # on left side of plot align towards right, on the other side towards left and keep the center in the middle
            xlim = axs[i].get_xlim()
            x_range = xlim[1] - xlim[0]
            if corresponding_current < xlim[0] + fraction_alignment*x_range:
                text_alignment = 'left'
            elif corresponding_current > xlim[1] - fraction_alignment*x_range:
                text_alignment = 'right'
            else:
                text_alignment = 'center'

            # find a good length for arrow
            ylim = axs[i].get_ylim()
            arrow_length = abs(fraction_arrow_length * (ylim[1] - ylim[0]))

            # draw the arrow and text
            axs[i].annotate('{:.1f} mT'.format(max_B_val), xy=(corresponding_current, max_B_val), 
                    xytext=(corresponding_current, max_B_val - arrow_length), 
                    horizontalalignment=text_alignment, fontweight = fontweight,
                    arrowprops=dict(facecolor='black', shrink=0.05, width= 2, headwidth=8))
        
        # for plots of single magnetic field components, add two arrows 
        elif mask_selected_plots[i] and flags_yaxis[i] in ['x','y','z']:
            # find the max and min value first
            i_min = np.argmin(plot_data[i])
            i_max = np.argmax(plot_data[i])
            min_B_val = plot_data[i, i_min]
            max_B_val = plot_data[i, i_max]
            corresponding_current_min = x_vals[i_min]
            corresponding_current_max = x_vals[i_max]
            
            # on left side of plot align towards right, on the other side towards left and keep the center in the middle
            xlim = axs[i].get_xlim()
            x_range = xlim[1] - xlim[0]
            if corresponding_current_min < xlim[0] + fraction_alignment*x_range:
                text_alignment_min = 'left'
            elif corresponding_current_min > xlim[1] - fraction_alignment*x_range:
                text_alignment_min = 'right'
            else:
                text_alignment_min = 'center'

            if corresponding_current_max < xlim[0] + fraction_alignment*x_range:
                text_alignment_max = 'left'
            elif corresponding_current_max > xlim[1] - fraction_alignment*x_range:
                text_alignment_max = 'right'
            else:
                text_alignment_max = 'center'

            # find a good length for arrow
            ylim = axs[i].get_ylim()
            arrow_length = abs(fraction_arrow_length * (ylim[1] - ylim[0]))

            # draw the arrow and text for minimum value 
            axs[i].annotate('{:.1f} mT'.format(min_B_val), xy=(corresponding_current_min, min_B_val), 
                    xytext=(corresponding_current_min, min_B_val + arrow_length), 
                    horizontalalignment=text_alignment_min, fontweight = fontweight,
                    verticalalignment = 'bottom',
                    arrowprops=dict(facecolor='black', shrink=0.05, width= 2, headwidth=8))

            # draw the arrow and text for maximum value 
            axs[i].annotate('{:.1f} mT'.format(max_B_val), xy=(corresponding_current_max, max_B_val), 
                    xytext=(corresponding_current_max, max_B_val - arrow_length), 
                    horizontalalignment=text_alignment_max, fontweight = fontweight,
                    verticalalignment = 'top',
                    arrowprops=dict(facecolor='black', shrink=0.05, width= 2, headwidth=8))


def add_insets_for_angular_plots(axs, x_vals, plot_data, plot_std_data, plot_exp_data, flags_yaxis, 
                        inset_ylim_factor = 0.25, which_side = 'right', manual_inset_ylim=None,
                        mask_selected_plots = None):
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
    if mask_selected_plots is None:
        mask_selected_plots = np.ones(len(axs), dtype=bool)

    for i in range(len(axs)):
        # for plots of total or in-plane magnetudes, only add one arrow where the maximum value is reached
        if flags_yaxis[i] in ['t','f'] and mask_selected_plots[i]:
            # Preselect the data on the desired side. Note that due to different configurations, 
            # the xvals can be sorted in ascending or decreasing order, such that a comparison 
            # of first and last value is required to pick the desired side. 
            len_data = len(plot_data[i])
            if which_side == 'left':
                inset_x = 0.1
                if x_vals[0] < x_vals[-1]:
                    # sorted in ascending order
                    data_on_side = plot_data[i,:len_data//2]
                else:
                    # sorted in decreasing order
                    data_on_side = plot_data[i,len_data//2:]
            else:
                inset_x = 0.65
                if x_vals[0] < x_vals[-1]:
                    # sorted in ascending order
                    data_on_side = plot_data[i,len_data//2:]
                else:
                    # sorted in decreasing order
                    data_on_side = plot_data[i,:len_data//2]

            # estimate mean value of angle on the chosen side. 
            mean_value_on_side = np.mean(data_on_side)

            # depending on the mean value, set the coordinates of inset axis in upper or lower half of the plot
            if mean_value_on_side > np.mean(axs[i].get_ylim()):
                inset_y = 0.15
            else:
                inset_y = 0.55

            # inset axis
            axins = axs[i].inset_axes([inset_x, inset_y, 0.33, 0.4])
            axins.errorbar(x_vals, plot_data[i], yerr=plot_std_data[i],
                                        linestyle='', marker='.', capsize = 2)
            axins.plot(x_vals, plot_exp_data[i], linestyle='--')

            # define which sub-regions of original plot should be shown by setting the limits of inset axis
            xlim = axs[i].get_xlim()
            x_range = xlim[1] - xlim[0]
            if which_side == 'left':
                axins.set_xlim(np.min(x_vals), xlim[0]+ 0.5*x_range)
            else:
                axins.set_xlim(xlim[0]+ 0.5*x_range, np.max(x_vals))

            # unless manual limits for inset-y-axis are provided, choose them automatically
            if manual_inset_ylim is None:
                y_range_inset = np.max(data_on_side) - np.min(data_on_side)
                axins.set_ylim( np.min(data_on_side) - inset_ylim_factor*y_range_inset, 
                                np.max(data_on_side) + inset_ylim_factor*y_range_inset)
            else:
                if isinstance(manual_inset_ylim, list):
                    if manual_inset_ylim[i] is None:
                        y_range_inset = np.max(data_on_side) - np.min(data_on_side)
                        axins.set_ylim( np.min(data_on_side) - inset_ylim_factor*y_range_inset, 
                                        np.max(data_on_side) + inset_ylim_factor*y_range_inset)
                    else: 
                        axins.set_ylim(manual_inset_ylim[i])
                elif isinstance(manual_inset_ylim, tuple):
                    axins.set_ylim(manual_inset_ylim[i])
                else:
                    raise ValueError('manual_inset_ylim must be None, list or tuple')

            # add boxes to indicate which region of original plot is shown in the inset
            axs[i].indicate_inset_zoom(axins)

def add_power_axis_on_top(axs, flag_xaxis, I, color_twinsy = 'firebrick'):
    """
    Add overall power as twiny axis on top in a different color and with nice tickmarks.
    
    Note that this is only reasonable if flag_xaxis is set to 'I1', 'I2' or 'I3'.

    Args:
    - axs (ndarray of plt.Axes objectes): represent the subplots, returned by generate_functions
    - flag_xaxis (str): indicates which coil the current plotted on x-axis refers to
    - I (ndarray of shape (#measurements, 3)): Contains the currents in all three coils, which is used to 
    estimate the ratio between the currents in all three coils. Returned by extract_raw_data_from_file.
    - color_twinsy (str): valid color name from matplotlib for the twin axis, including the ticks and labels
    """
    # determine the coil of which the current is plotted
    if flag_xaxis[0] == 'I':
        coil_number = int(flag_xaxis[1])

    # add twin axes for each plot
    ax_twinsy = np.array([axs[i].twiny() for i in range(len(axs))])

    # add x-axis on top to display total power via setting the limits with negative and positive powers
    # and update to only positive labels later (advantage is that tick positions are chosen automatically,
    # which usually is nicer than translating the ticks of current to power values):
    for i in range(len(axs)):
        # get the limits of current axis and transform to power, however with the same sign as the currents
        Ilimits = axs[i].get_xlim()
        Pmin, Pmax = [np.sign(value)* estimate_power(I[0,:], value, coil_number=coil_number) for value in Ilimits]
        
        # set limits of twin axis accordingly
        ax_twinsy[i].set_xlim(Pmin, Pmax)

        # set color of ticks
        ax_twinsy[i].tick_params(axis='x', colors=color_twinsy)
        ax_twinsy[i].spines['top'].set_color(color_twinsy)
        
        if i == 0:
            # update the tick labels, such that absolute values are displayed. Also keep the automatically chosen 
            # format (int or float) of ticks
            power_ticks = ax_twinsy[i].get_xticks()
            if '.' in str(ax_twinsy[i].get_xticklabels()[0]):
                ticks_type = float
            else:
                ticks_type = int
            corrected_power_ticks = [ticks_type(abs(tick)) for tick in power_ticks]

            # set tick labels and axis label
            ax_twinsy[i].set_xticklabels(corrected_power_ticks, color=color_twinsy)
            ax_twinsy[i].set_xlabel('overall power, $P$ [W]', color=color_twinsy)

        # for subsequent plots, hide tick labels of twinx axes but keep the ticks by setting tick labels to empty strings
        else:
            labels = [item.get_text() for item in ax_twinsy[i].get_xticklabels()]
            empty_string_labels = ['']*len(labels)
            ax_twinsy[i].set_xticklabels(empty_string_labels)

def estimate_power(ratios, value, coil_number=1, R = 0.47 ):
    """
    Return the total power [W] of all three coils when the current in the passed coil (1, 2 or 3) is value [A].

    Args:
    - ratios (nonzero 1d-ndarray of length 3): ratios of the currents in the three coils. This can be the first element of 
    the current array of shape (number measurements, 3). Note that the entry corresponding to the provided coil_number
    must not be 0, because this coil would have zero current for the whole measurement series. 
    - value (float): current value [A] in coil with number coil_number
    - coil_number (int): number of the coil for which value is provided, can be 1, 2 or 3.
    - R (float): resistance [Ohm] of a coil 

    Return: power (float)
    """
    # check reasonability of inputs
    if np.all(ratios == 0):
        raise ValueError('ratios must be non-zero!')
    if coil_number not in [1,2,3]:
        raise ValueError('coil_number must be in \{1,2,3\}, not {}!'.format(coil_number))
    if ratios[coil_number-1] == 0:
        raise ValueError('The provided ratios are invalid, since coil {} has zero current.'.format(coil_number))
    
    # normalize ratios, such that the desired coil has a factor 1, and multiply by value
    currents = value * ratios/ratios[coil_number-1]

    return R * np.sum(currents**2)  # sum over all three coils

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


def find_closest_measured_field(pos, field, a, b):
    """
    Return the measured magnetic field at the closest measurement position next to [a,b]
    """
    distances = np.sqrt((pos[:,0]-a)**2 + (pos[:,1]-b)**2)
    index_closest = np.argmin(distances)
    
    return field[index_closest]

def distance_to_center(pos, field, a, b):
    """
    Estimate |B_i - B_ab|_2 for magnetic field vectors B_i at position i and B_ab at position 
    [a,b] for all i in range(len(i)).
    """
    # take the measured field at the position that is closest to [a,b]
    closest_field = find_closest_measured_field(pos, field, a, b)
    return  np.sqrt((field[:,0] - closest_field[0])**2 + (field[:,1] - closest_field[1])**2 + (field[:,2]- closest_field[2])**2)

def weighted_distances(pos, field, a, b):
    """
    Estimate |B_i - B_ab|_2 for magnetic field vectors B_i at position i and B_ab at position 
    [a,b] for all i in range(len(i)).
    """
    # take the measured field at the position that is closest to [a,b]
    distances = np.sqrt((pos[:,0]-a)**2 + (pos[:,1]-b)**2)
    inplane_fields = np.sqrt((field[:,0])**2 + (field[:,1])**2)
    return distances*inplane_fields**2

def find_center_of_mass(pos, field):
    """
    Estimate the resulting weights for each position in pos when assuming that this position is the center, 
    using the distance_to_center function to estimate the weights. 
    Return the position with smalles weight. 
    """
    objectives = np.zeros(len(pos))
    for i in range(len(pos)):
        objectives[i] = np.sum(weighted_distances(pos, field, pos[i,0], pos[i,1]))

    index_min = np.argmin(objectives)
    return pos[index_min, :2]

def find_center(pos, field):
    """
    Estimate the resulting weights for each position in pos when assuming that this position is the center, 
    using the distance_to_center function to estimate the weights. 
    Return the position with smalles weight. 
    """
    objectives = np.zeros(len(pos))
    for i in range(len(pos)):
        objectives[i] = np.sum(distance_to_center(pos, field, pos[i,0], pos[i,1]))

    index_min = np.argmin(objectives)
    return pos[index_min, :2]

def weights_of_rays(pos, field, phi, a, b, length=20, num=20):
    """
    For each point in pos, generate three rays with this point at center, where the rays are at 120 degrees 
    angle to each other and the first one has an angle phi wrt to x-axis.
    For each point on these rays, find the closest position in pos at which a measurement was performed 
    and save the norm of the in-plane field component of the magnetic field measured at this point. 
    Return an ndarray of shape (3, num), where each entry contains the in-plane component of 
    the closest measured point.
    """
    # generate rays at 120 degrees between each other, starting from a, b
    rays, radius_rays = generate_lines(phi, a, b, length=length, num=num)
    
    field_distances = np.ones((3,num))
    for i in range(3):
        for j in range(num):
            # take the measured field at the position that is closest to a point on one of the rays
            closest_field = find_closest_measured_field(pos, field, rays[i, j, 0], rays[i, j, 1])
            field_distances[i,j] = np.linalg.norm(closest_field[:2])
    return field_distances

def find_center_of_rays(pos, field, phi=30, length=20, num=20):
    """
    For each point in pos, estimate the weight using weights_of_rays of this point as center. 
    The weights are returned as ndarray of shape (3, num), average over the squares of the second axis 
    as a weighted sum, where the latter weights are set originate from a Gaussian function.
    Return the point with the smalles weight. 

    Note that the result is heavily dependent on the gaussian function used (i.e. the prefactor and 
    the factor inside the exponential).
    """
    objectives = np.zeros(len(pos))
    for i in range(len(pos)):
        # increase weights of points that are closer to center
        gaussian = lambda x: np.exp(-x/(num/10))
        weights_along_ray = gaussian(np.linspace(0,10,num=20))
        objectives_per_ray = np.average(weights_of_rays(pos, field, phi, pos[i,0], pos[i,1])**2, 
                                            axis=1, weights=weights_along_ray)
        objectives[i] = np.sum(objectives_per_ray)

    index_min = np.argmin(objectives)

    return pos[index_min, :2]

def generate_lines(phi, center_x, center_y, length=20, num=20):
    """
    Generate 3 ndarrays of shape (length, 2) that correspond to three lines in xy-plane starting 
    center position [center_x, center_y] and moving outwards, all at an angle of 120 degrees.
    One of the lines is at angle phi wrt to x-axis.

    Args: 
    - phi (float) is angle wrt to x-axis in degrees
    """
    radius = np.linspace(0, length, num=num)

    angles = np.array([0, 120, -120]) + phi
    angles = np.radians(angles)

    rays = np.zeros((3, num, 2))
    for i in range(3):
        rays[i,:,0] = np.cos(angles[i])*radius[:] + center_x
        rays[i,:,1] = np.sin(angles[i])*radius[:] + center_y

    return rays, radius

def distance_to_lines(pos, field, phi, a, b, length=20, num=20):
    # generate rays at 120 degrees between each other, starting from a, b
    rays, radius_rays = generate_lines(phi, a, b, length=length, num=num)
    
    field_distances = np.ones(len(pos))
    for i in range(len(pos)):
        # find distance from pos[i] to (a,b)
        radius = np.linalg.norm(pos[i,:2]- np.array([a, b]))

        # find three points on ray at the same radius from center 
        i_same_radius = np.argmin(np.abs(radius-radius_rays))
        
        # find the closest point of the three
        i_closest_point = np.argmin(np.linalg.norm(rays[:,i_same_radius] - pos[i,:2], axis=1))

        # take the measured field at the position that is closest to the closest point on one of the rays
        closest_field = find_closest_measured_field(pos, field, 
                            rays[i_closest_point, i_same_radius, 0], rays[i_closest_point, i_same_radius, 1])
        
        field_distances[i] = np.sqrt((field[i,0] - closest_field[0])**2 + (field[i,1] - closest_field[1])**2 )

    return field_distances

def find_center_of_rays(pos, field, phi=30, length=20, num=20):
    """
    For each point in pos, estimate the weight using weights_of_rays of this point as center. 
    The weights are returned as ndarray of shape (3, num), average over the squares of the second axis 
    as a weighted sum, where the latter weights are set originate from a Gaussian function.
    Return the point with the smalles weight. 

    Note that the result is heavily dependent on the gaussian function used (i.e. the prefactor and 
    the factor inside the exponential).
    """
    objectives = np.zeros(len(pos))
    for i in range(len(pos)):
        # increase weights of points that are closer to center
        gaussian = lambda x: np.exp(-x/(num/10))
        weights_along_ray = gaussian(np.linspace(0,10,num=20))
        objectives_per_ray = np.average(weights_of_rays(pos, field, phi, pos[i,0], pos[i,1])**2, 
                                            axis=1, weights=weights_along_ray)
        objectives[i] = np.sum(objectives_per_ray)

    index_min = np.argmin(objectives)

    return pos[index_min, :2]

def find_center_of_rays_all(pos, field, phi=30, length=20, num=20):
    """
    For each point in pos, estimate the weight using weights_of_rays of this point as center. 
    The weights are returned as ndarray of shape (3, num), average over the squares of the second axis 
    as a weighted sum, where the latter weights are set originate from a Gaussian function.
    Return the point with the smalles weight. 

    Note that the result is heavily dependent on the gaussian function used (i.e. the prefactor and 
    the factor inside the exponential).
    """
    objectives = np.zeros(len(pos))
    for i in range(len(pos)):
        # increase weights of points that are closer to center
    
        objectives[i] = np.sum(distance_to_lines(pos, field, phi, pos[i,0], pos[i,1])**2)

    index_min = np.argmin(objectives)

    return pos[index_min, :2]


def plot_rotation_plane(I, mean_values, std_values, expected_values, distance = 1.5, show_labels = True,
                        rotated_expection=None):
    """
    Generate a plot of the two magnetic field components within the plane of rotation.

    Args: 
    - I, mean_values, std_values, expected_values (ndarrays) are of shape (#measurements, 3), 
    containing applied current, experimentally estimated/expected mean values and standard deviations 
    for x,y,z-directions.
    - distance (float): distance between sensor and tip [mm], which is added as plot label 
    - show_labels (bool): flag to switch on/off labels of plots

    Return: fig, axs
    """
    # create figure
    fig, ax = plt.subplots()
    fig.set_size_inches(4, 4)
    
    # find the rotation axis
    if np.all(np.isclose(expected_values[:,0], expected_values[0,0])):
        rot_axis = 0 # 'x'
    elif np.all(np.isclose(expected_values[:,1], expected_values[0,1])):
        rot_axis = 1 # 'y'
    elif np.all(np.isclose(expected_values[:,2], expected_values[0,2])):
        rot_axis = 2 # 'z'
    else:
        raise NotImplementedError('Rotations about general axes are not implemented yet.')

    print('rotation about {} axis'.format(rot_axis))

    # set which components of magnetic field should be plotted -> the ones which are rotated
    plot_components = np.delete([0,1,2], rot_axis)

    # plot components orthogonal to rotation axis
    ax.errorbar(mean_values[:, plot_components[0]], mean_values[:, plot_components[1]], 
                            xerr = std_values[:, plot_components[0]], yerr = std_values[:, plot_components[1]],
                            linestyle='', marker='.', capsize = 2, 
                            label = 'measured @ {:.1f} mm'.format(distance))
    ax.plot(expected_values[:, plot_components[0]], expected_values[:, plot_components[1]], linestyle='--', 
                                    label = 'linear model @ 3mm')
    if rotated_expection is not None:
        ax.plot(rotated_expection[:, plot_components[0]], rotated_expection[:, plot_components[1]], linestyle='--', 
                                    label = 'rotated linear model')

    # set axis labels
    components = ['$B_x$ [mT]', '$B_y$ [mT]', '$B_z$ [mT]']
    ax.set_xlabel(components[plot_components[0]])
    ax.set_ylabel(components[plot_components[1]])

    # set aspect ratio to one, such that a circle actually looks round 
    ax.set_aspect('equal')

    # if desired, switch on legend
    if show_labels:
        ax.legend()

    plt.tight_layout()

    return fig, ax

def swap_components(a, index_last):
    """
    Permute the indices of the last axis of a, such that index_last is the last one. 
    The final array is sorted such that the order of the second axis becomes 'yzx', 'zxy' or 'xyz',
    so the permutation has positive sign.

    Args:
    - a (ndarray with length 3 along last dimension): array for which the indices of the last axis 
    should be permuted
    - index_last (int): current index of the array, which should become the last index = 2
    """
    # if the order is already correct, return a as it is
    if index_last == 2:
        return a
    else:
        b = np.zeros_like(a)
        for i in range(3):
            b[..., i] = a[..., (i-2+index_last)%3]

        return b

def plot_vs_rotation_angle(I, mean_values, std_values, expected_values, height_per_plot = 2,
                        flags_yaxis = 'a', distance = 1.5, show_labels = True,
                        plot_delta_sim=False, rotated_expection=None, allow_negative_at_beginning=False):
    """
    Generate plots of the desired rotation angle on x-axis vs measured rotation angle and the 
    third field component on the y-axis. 

    Args: 
    - I, mean_values, std_values, expected_values (ndarrays) are of shape (#measurements, 3), 
    containing applied current, experimentally estimated/expected mean values and standard deviations 
    for x,y,z-directions.
    - flags_yaxis (string): contains letters as flags, where valid letters are mentioned below.
    The number of letters may vary, but at least one valid letter should be contained. For each flag,
    a plot is generated with the according quantity plotted on the y-axis. The plots are generated
    in the same order as the flag-letters are passed. 
        - 'x', 'y', 'z': single component of magnetic field, p.e. B_z
        - 'o': off-plane component of magnetic field, p.e. B_z if the rotation is about z-axis
        - 'm': total magnitude of magnetic field
        - 'p': in-plane magnitude of magnetic field within plane of rotation 
        - 'r': in-plane rotation angle within plane of rotation
        - 'a': angle wrt plane of rotation
    - height_per_plot (float): height of each plot [inches]. Default for several plots is 2.
    The usual height of a single plot is 4.
    - distance (float): distance between sensor and tip [mm], which is added as plot label 
    - show_labels (bool): flag to switch on/off labels of plots
    - allow_negative_at_beginning (bool): it can happen that a few angles would be negative in the beginning, 
    which are set to values slightly below 360 degrees. Correct for this by resetting those to negative values
    if True is passed, else it is ignored.

    Return: fig, axs
    """
    # if an empty string is passed as flags_yaxis, set it to 'a'
    if len(flags_yaxis) == 0:
        flags_yaxis = 'a'

    # create a simple plot with as many axes as letters in plots_yaxis
    number_plots = len(flags_yaxis)
    fig, axs = plt.subplots(number_plots, 1, sharex=True)
    fig.set_size_inches(6, number_plots * height_per_plot)
    
    # find the rotation axis
    if np.all(np.isclose(expected_values[:,0], expected_values[0,0])):
        rot_axis = 0 # 'x'
    elif np.all(np.isclose(expected_values[:,1], expected_values[0,1])):
        rot_axis = 1 # 'y'
    elif np.all(np.isclose(expected_values[:,2], expected_values[0,2])):
        rot_axis = 2 # 'z'
    else:
        raise NotImplementedError('Rotations about general axes are not implemented yet.')

    # estimate total magnitude and magnitude within the plane of ration
    mean_magnitudes = np.linalg.norm(mean_values, axis=1)
    expected_magnitudes = np.linalg.norm(expected_values, axis=1)
    if rotated_expection is not None:
        rotated_magnitudes = np.linalg.norm(rotated_expection, axis=1)

    # if number_plots=1, axs is returned as AxesSubplot class instead of an ndarray containing
    # instances of this class. Since the following requires a ndarray, ensure to have an ndarray!
    if number_plots == 1:
        axs = np.array([axs])

    rot_plane_mags = np.linalg.norm(np.delete(mean_values, rot_axis, axis=1), axis=1)
    rot_plane_mags_expected = np.linalg.norm(np.delete(expected_values, rot_axis, axis=1), axis=1)
    if rotated_expection is not None:
        rot_plane_mags_rotated = np.linalg.norm(np.delete(rotated_expection, rot_axis, axis=1), axis=1)

    # collect plot data.
    # Note: errorbars display std, estimate errors for magnitudes (and angle) using propagation of uncertainty,
    # assuming that the measured fields in x,y,z direction are independent variables
    plot_mean_data = []
    plot_std_data = []
    plot_expected_data = []
    plot_rotated_data = []
    ylabels = []
    for flag in flags_yaxis:
        # magnetic field in x-direction
        if flag == 'x':
            plot_mean_data.append(mean_values[:, 0])
            plot_std_data.append(std_values[:, 0])
            plot_expected_data.append(expected_values[:, 0])
            if rotated_expection is not None:
                plot_rotated_data.append(rotated_expection[:, 0])
            ylabels.append('$B_x$ [mT]')
        # magnetic field in y-direction
        elif flag == 'y':
            plot_mean_data.append(mean_values[:, 1])
            plot_std_data.append(std_values[:, 1])
            plot_expected_data.append(expected_values[:, 1])
            if rotated_expection is not None:
                plot_rotated_data.append(rotated_expection[:, 1])
            ylabels.append('$B_y$ [mT]')
        # magnetic field in z-direction
        elif flag == 'z':
            plot_mean_data.append(mean_values[:, 2])
            plot_std_data.append(std_values[:, 2])
            plot_expected_data.append(expected_values[:, 2])
            if rotated_expection is not None:
                plot_rotated_data.append(rotated_expection[:, 2])
            ylabels.append('$B_z$ [mT]')
        # field component parallel to rotation axis
        elif flag == 'o':
            plot_mean_data.append(mean_values[:, rot_axis])
            plot_std_data.append(std_values[:, rot_axis])
            plot_expected_data.append(expected_values[:, rot_axis])
            if rotated_expection is not None:
                plot_rotated_data.append(rotated_expection[:, rot_axis])
            components = ['$B_x$ [mT]', '$B_y$ [mT]', '$B_z$ [mT]']
            ylabels.append(components[rot_axis])
        # magnitude of magnetic field
        elif flag == 'm':
            plot_mean_data.append(mean_magnitudes)
            plot_std_data.append(estimate_std_magnitude(mean_values, std_values))
            plot_expected_data.append(expected_magnitudes)
            if rotated_expection is not None:
                plot_rotated_data.append(rotated_magnitudes)
            ylabels.append('$|B|$ [mT]')
        # magnetude of magnetic field within plane of ration 
        elif flag == 'p':
            plot_mean_data.append(rot_plane_mags)
            plot_std_data.append(estimate_std_inplane(swap_components(mean_values,rot_axis), 
                                                    swap_components(std_values,rot_axis)))
            plot_expected_data.append(rot_plane_mags_expected)
            if rotated_expection is not None:
                plot_rotated_data.append(rot_plane_mags_rotated)
            if rot_axis == 0:
                ylabels.append('$|B_{yz}|$ [mT]')
            elif rot_axis == 1:
                ylabels.append('$|B_{xz}|$ [mT]')
            elif rot_axis == 2:
                ylabels.append('$|B_{xy}|$ [mT]')
        # angle within plane of rotation
        elif flag == 'r':
            swapped_mean_values = swap_components(mean_values, rot_axis)
            swapped_std_values = swap_components(std_values, rot_axis)
            swapped_exp_values = swap_components(expected_values, rot_axis)
            measured_phi = get_phi(swapped_mean_values, cut_phi_at_0=True)
            measured_error_phi = estimate_std_phi(swapped_mean_values, swapped_std_values)
            expected_phi = get_phi(swapped_exp_values, cut_phi_at_0=True)
            # allow negative values in the beginning if provided:
            if allow_negative_at_beginning:
                for i in range(len(mean_values) // 2):
                    if measured_phi[i] > 300:
                        measured_phi[i] -= 360
                    if expected_phi[i] > 300:
                        expected_phi[i] -= 360
            plot_mean_data.append(measured_phi)
            plot_std_data.append(measured_error_phi)
            plot_expected_data.append(expected_phi)
            if rotated_expection is not None:
                swapped_rotated_values = swap_components(rotated_expection, rot_axis)
                rotated_phi = get_phi(swapped_rotated_values, cut_phi_at_0=True)
                # allow negative values in the beginning if provided:
                if allow_negative_at_beginning:
                    for i in range(len(mean_values) // 2):
                        if rotated_phi[i] > 300:
                            rotated_phi[i] -= 360
                plot_rotated_data.append(rotated_phi)
            ylabels.append('rotation angle, $\\Phi\'$ [°]')
        # angle wrt plane of rotation
        elif flag == 'a':
            swapped_mean_values = swap_components(mean_values, rot_axis)
            swapped_std_values = swap_components(std_values, rot_axis)
            plot_mean_data.append(90 - np.degrees(
                np.arccos(mean_values[:, rot_axis]/mean_magnitudes)))
            plot_std_data.append(estimate_std_theta(swapped_mean_values, swapped_std_values))
            plot_expected_data.append(90 - np.degrees(
                np.arccos(expected_values[:, rot_axis]/expected_magnitudes)))
            if rotated_expection is not None:
                plot_rotated_data.append(90 - np.degrees(
                    np.arccos(rotated_expection[:, rot_axis]/rotated_magnitudes)))
            ylabels.append('offset angle, $\\Theta$ [°]')
        # account for invalid flags:
        else:
            raise ValueError(
                '{} is not a valid flag, it should be in [\'x\', \'y\', \'z\', \'m\', \'a\', \'p\']!'.format(flag))

    # estimate the desired/theoretical rotation angle wrt to the first espected vector. 
    x_vals = get_phi(swap_components(expected_values, rot_axis), cut_phi_at_0=True)
    x_vals -= get_phi(swap_components(expected_values[0,:].reshape((1,3)), rot_axis), cut_phi_at_0=True)
    x_vals = x_vals % 360
   
    # it is possible that x_vals is not sorted, which results in a horizontal line for the linear model in the plots,
    # which is not desired. Thus, sort x_vals and shuffle the plot data accordingly. 
    # Note that this does not alter the plotted values at all, it only corrects for the horizontal line 
    # as an artifact from an unsorted x-axis
    i_sort = np.argsort(x_vals)
    x_vals = x_vals[i_sort]
    plot_mean_data = np.array(plot_mean_data)[:, i_sort]
    plot_std_data = np.array(plot_std_data)[:, i_sort]
    plot_expected_data = np.array(plot_expected_data)[:, i_sort]
    plot_rotated_data = np.array(plot_rotated_data)[:, i_sort]
    
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
                                    label = 'linear model @ 3mm')
            if rotated_expection is not None:
                axs[i].plot(x_vals, plot_rotated_data[i], linestyle='--', 
                                        label = 'rotated linear model')
            
    if show_labels:
        axs[0].legend()

    # set axis labels
    axs[-1].set_xlabel('desired rotation angle $\Phi$ [°]')
    axs[-1].yaxis.set_major_locator(MaxNLocator(nbins='auto', steps=[1,3,6,10]))


    # if desired, switch on legend
    if show_labels:
        axs[0].legend()

    plt.tight_layout()

    return fig, axs

def fit_rotation_matrix(mean_values, expected_values, convention = 'xzx'):
    """
    Fits a rotation that yields the minimum sum of square distances between mean_values and rotated
    expected_values. This method uses the curve_fit method of scipy.optimize, which makes use of the 
    least_squares method.

    Args: 
    - mean_values, expected_values (ndarrays of shape (#measurements, 3)): contain the measured and expected data
    - convention (str): the convention used to describe a general 3d-rotation by Euler angles. The resulting 
    angles depend on the convention, while the calculated rotation remains the same.

    Returns:
    - p, pcov (list of floats): the estimated Euler angles that achieve the best fitting rotation.
    - rotated_expections (ndarray of shape (#measurements, 3)): the estimated values when applying the
    final rotation to the expected_values 
    """
    p, pcov = curve_fit(lambda x, alpha, beta, gamma: apply_rotation(x, alpha, beta, gamma, convention=convention).flatten(), 
                        expected_values, mean_values.flatten(), p0=[5,90,0])
    return p, pcov, apply_rotation(expected_values, *p)

def apply_rotation(x, alpha, beta, gamma, convention = 'xzx'):
    """ 
    Applies a rotation with Euler angles alpha, beta and gamma to an input vector

    Args: 
    - x (ndarray of length 3 or shape (N,3)): Input vector(s) which should be rotated
    - alpha, beta, gamma (flaot): Euler angles of the rotation
    - convention (str): the convention used to describe a general 3d-rotation by Euler angles.

    Returns the rotated input vectors
    """
    r = R.from_euler(convention, [alpha, beta, gamma], degrees=True)
    return r.apply(x)

def rotation_on_basis_vectors(alpha, beta, gamma, convention = 'xzx', verbose=True):
    """ 
    Estimate the effect of a rotation with Euler angles alpha, beta, gamma on the coordinate axes.
    Note that considering the coordinate axes actually involves the inverse of the rotation.

    Args:
    - alpha, beta, gamma (flaot): Euler angles of the rotation
    - convention (str): the convention used to describe a general 3d-rotation by Euler angles.
    - verbose (bool): Flag to switch on/off additional printing of effects

    Returns:
    - rotated (ndarray of shape (3,3)): Contains coordinate axes of the rotated system in the original coordinates.
    The first axis covers the three axes, the second the respective coordinates.
    - delta_phi (ndarray of length 3): differences between the polar angles of the original and 
    transformed axes (expressed in original coordinates)
    - delta_phi (ndarray of length 3): differences between the azimuthal angles of the original and 
    transformed axes (expressed in original coordinates)
    """
    r = R.from_euler(convention, [alpha, beta, gamma], degrees=True).inv()
    basis_vectors = np.array([[1,0,0], [0,1,0], [0,0,1]])
    rotated = np.array([r.apply(basis_vectors[i]) for i in range(3)])

    delta_phi = get_phi(basis_vectors) - get_phi(rotated)
    delta_theta = get_theta(basis_vectors) - get_theta(rotated)

    if verbose:
        print('effect on axes (note that this involves inverse of rotation):')
        axes = ['x', 'y', 'z']
        for i in range(3):
            print('{}-axis -> {} (Dphi = {:.2f}°, Dtheta = {:.2f}°)'.format(axes[i], np.round(rotated[i],2),
                                                                        delta_phi[i], delta_theta[i]))
    
    return rotated, delta_phi, delta_theta


def get_remanent_B(H_falling, H_rising, B_falling, B_rising):
    # find zero values of H
    i_fall = np.argmin(np.abs(H_falling))
    i_rise = np.argmin(np.abs(H_rising))
    # return y-axis offsets
    return B_falling[i_fall], B_rising[i_rise]

def linearly_fit_root(a, b):
    return -a[1]*(b[0]-a[0])/(b[1]-a[1]) + a[0]

def get_coercivity(H_falling, H_rising, B_falling, B_rising):
    # falling branch
    if np.min(np.abs(B_falling)) == 0:
        # check for zero B values, although it is rather unlikely that an experiments produces exactly zero field 
        i_fall = np.argmin(np.abs(B_falling))
        H_coer_fall = H_falling[i_fall]
    else:
        # take two entries with lowest absolut values, draw line through them and estimate shift of origin
        i1, i2 = np.argsort(np.abs(B_falling))[:2]
        H_coer_fall = linearly_fit_root([H_falling[i1], B_falling[i1]], [H_falling[i2], B_falling[i2]])

    # rising branch
    if np.min(np.abs(B_rising)) == 0:
        # check for zero B values, although it is rather unlikely that an experiments produces exactly zero field 
        i_rise = np.argmin(np.abs(B_rising))
        H_coer_rise = H_rising[i_rise]
    else:
        # take two entries with lowest absolut values, draw line through them and estimate shift of origin
        i1, i2 = np.argsort(np.abs(B_falling))[:2]
        H_coer_rise = linearly_fit_root([H_rising[i1], B_rising[i1]], [H_rising[i2], B_rising[i2]])

    return H_coer_fall, H_coer_rise

def get_direction_vector(vecs):
    """
    Estimate the overall normalized direction vector within the xy-plane for a set of vectors, 
    pointing from (about) along the curve defined by the vectors and projected onto the xy-plane.

    Args:
    - vecs (ndarray of shape (N,3) or (N,2)): Input data vectors (for example of measured 
    magnetic field vectors). The third dimension is ignored entirely. 

    Returns: direction (1d-ndarray of length 2): normalized direction vector
    """
    # first check that x-data are not too close to 0, take 5 mT as boundary
    if not np.all(np.abs(vecs[:,0]) < 5 ):
        # fit data projected onto xy-plane with linear function plus an offset
        lin_fct = lambda x, a, b: a*x +b
        p, _ = curve_fit(lin_fct, vecs[:,0], vecs[:,1], p0 = [1,0])
        
        # take the x-value with largest magnitude to correctly set positive direction along provided data
        x = vecs[np.argmax(np.abs(vecs[:,0])), 0]
        direction = np.array([x, lin_fct(x,*p)-lin_fct(0,*p)])

    # else do the same thing as before, but exchange x and y when fitting:
    else:
        lin_fct = lambda x, a, b: a*x +b
        p, _ = curve_fit(lin_fct, vecs[:,1], vecs[:,0], p0 = [1,0])
        
        # take the y-value with largest magnitude to correctly set positive direction along provided data
        y = vecs[np.argmax(np.abs(vecs[:,1])), 1]
        direction = np.array([lin_fct(y,*p)-lin_fct(0,*p), y])

    # normalize direction before returning
    return direction / np.linalg.norm(direction)

def get_relative_inplane_angles(fields, verbose=False):
    """
    Estimate the two angles between the directions of first + second entry along first dimension 
    and between directions of first and third entry along first dimension, each projected onto the xy-plane. 
    This function is supposed to be used with data that originate from ramping the current in 
    the three coils individually, while the remaining coils are switched off.  

    Args:
    - fields (ndarray of shape (3, N, 3)): estimated fields with the three field directions on last axis
    and the data originating from ramping the coils individually along first dimension.
    - verbose (bool): If True, print the resulting three directions

    Return angles (1d-ndarray of length 2), containing the angles between virgin hysteresis curves of
    first and second coil, first and third coil, and second and third coil. 
    """
    # get normalized direction vectors first
    directions = np.zeros((3,2))
    for i in range(3):
        directions[i] = get_direction_vector(fields[i])
    
    if verbose:
        [print('direction for coil {}: {}'.format(i+1, directions[i])) for i in range(3)]

    # estimate angles using dot product
    angles = np.arccos(np.array([np.dot(directions[0], directions[1]),
                                np.dot(directions[0], directions[2]),
                                np.dot(directions[1], directions[2]) ]))

    # convert to degrees before returning
    return np.degrees(angles)


def evaluate_performance(measured, fitted):
    """
    Evaluate the performance of the fit by estimating different measured of 
    deviation between fitted and measured parameters. 
    Estimated parameters are RMS error, and the angular accuracy in terms of 
    RMS angular error as well as mean, std, min, max and median angular error.

    Args:
    - measured, fitted (ndarrays of shape (N,3)): Estimated and fitted values. 
    """
    # estimate RMS errors 
    RMSE = estimate_RMS_error(measured.flatten(), fitted.flatten())

    # evaluate angular accuracy
    dot = np.array([np.dot(measured[i], fitted[i]) for i in range(len(measured))])
    norms_measured = np.linalg.norm(measured, axis=1)
    norms_fits = np.linalg.norm(fitted, axis=1)
    alphas = np.degrees(np.arccos(dot / (norms_measured * norms_fits)))

    # print all measures
    print(f'RMS error fit: {RMSE:.2f} mT')
    print('RMS angular error: {:.2f}°'.format(estimate_RMS_error(alphas, np.zeros_like(alphas))))
    print('mean angular error: {:.2f}°, std: {:.2f}°'.format(np.mean(alphas), np.std(alphas)))
    print('min / max angular error: {:.2f}° / {:.2f}°'.format(np.min(alphas), np.max(alphas)))
    print('median angular error: {:.2f}°'.format(np.median(alphas)))

def collectAndExtract(directory, B_min, remove_saturation = True,
                        verbose=False, fraction_cutoff = 0.02,
                        affine_fct = lambda x, a, b: a*x + b):
    """
    Extract data from a measurement series, where all measurement runs are stored in 
    various files in a directory. All data below the threshold B_min are ignored.
    If desired, data above saturation are ignored, too.

    Args:
    - directory (str): valid path of directory that contains the data
    - B_min (float): threshold value that sets the minimum considered magnetic field magnitude 
    - remove_saturation (bool): if True, data above saturation are removed as well. 
    For estimating the boundaries of the linear regime and where saturation starts the function
    find_start_of_saturation is used, which takes the remaining optional arguments fraction_cutoff
    and affine_fct. This algorithm is not perfect and requires some fine-tuning, but it seems 
    to work qualitatively well, so the estimated boundaries are not far from what one would 
    pick as start of saturation. 
    - fraction_cutoff (float): paramter passed to find_start_of_saturation 
    - affine_fct (function): paramter passed to find_start_of_saturation 

    Return:
    - currents, B_measured, B_expected (ndarrays of shape (N, 3)): All currents, measured and expected
    fields extracted from the directory that are above the threshold magnitude (and below saturation
    if desired). 
    """
    # collect all csv-files in this directory
    filenames = []
    [filenames.append(file) for file in os.listdir(directory) if file.endswith(".csv")]
    
    # if in list, remove the linear_fits file from previous fits
    try:
        filenames.remove('linear_fits.csv')
    except ValueError:
        pass
    print(f'files considered: {len(filenames)}')

    # loop through all csv files in a dictionary and fit the data
    for i in range(len(filenames)):
        if verbose:
            print(filenames[i])

        # read in raw measurment data
        data_filepath = os.path.join(directory, filenames[i])
        I, mean_data, std_data, expected_fields = extract_raw_data_from_file(data_filepath)

        # -> this could be used to exclude data above saturation, but could be left out
        # estimate minimum and maximum indices of region within which the linear relation between current and field holds  
        # even though find_start_of_saturation offers the possibility to specify the considered component, 
        # keep the default stting, which detects the field component that has the greatest absolute field values. 
        # This should work fine for situations, where one component is dominating. 
        i_min, i_max = find_start_of_saturation(I, mean_data, std_data, fraction_cutoff=fraction_cutoff,
                                fitting_fct = affine_fct)

        # estimate field magnitudes
        magnitudes = np.linalg.norm(mean_data, axis=1)

        # set up mask to only keep data with magnitudes larger than B_min
        mask_keep = magnitudes >= B_min

        # optionally: also remove potentially saturated part
        if remove_saturation:
            mask_keep[:i_min+1] = False
            mask_keep[i_max:] = False

        # collect all relevant data
        if i == 0:
            B_measured = mean_data[mask_keep]
            B_expected = expected_fields[mask_keep]
            currents = I[mask_keep]
        else:
            B_measured = np.append(B_measured, mean_data[mask_keep], axis=0)
            B_expected = np.append(B_expected, expected_fields[mask_keep], axis=0)
            currents = np.append(currents, I[mask_keep], axis=0)

    print(f'final shape of considered array: {B_measured.shape}')
    return currents, B_measured, B_expected