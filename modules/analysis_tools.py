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

from modules.general_functions import angle_wrt_z, inplane_angle_wrt_x


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
    Return the in-plane angle phi with respect to x-axis.
    """
    angles =  np.degrees(np.arctan2(values[:,1], values[:,0]))

    # if cut should be at 0 degrees, add 360 degrees to all negative values
    if cut_phi_at_0:
        mask = angles < 0
        angles[mask] = 360 + angles[mask]

    return angles

def generate_plots(I, mean_values, std_values, expected_values, flag_xaxis = 'I1', flags_yaxis = 'zmt',
                        plot_delta_sim = False, directory = None, image_name_postfix = 'B_vs_I', 
                        height_per_plot = 2, save_image = True, distance = 3.0, xlim = None, 
                        ylim_field_abs = None, ylim_field_single = None, ylim_theta = None, ylim_phi = None,
                        show_labels = True, remove_half = 0, ygrid = False, cut_phi_at_0 = False,
                        show_image = True):
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
        # magnetic field in z-direction
        if flag == 'x':
            plot_mean_data.append(mean_values[:, 0])
            plot_std_data.append(std_values[:, 0])
            plot_expected_data.append(expected_values[:, 0])
            ylabels.append('$B_x$ [mT]')
        elif flag == 'y':
            plot_mean_data.append(mean_values[:, 1])
            plot_std_data.append(std_values[:, 1])
            plot_expected_data.append(expected_values[:, 1])
            ylabels.append('$B_y$ [mT]')
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
            plot_mean_data.append(np.degrees(
                np.arccos(mean_values[:, 2]/mean_magnitudes)))
            plot_std_data.append(np.degrees(
                estimate_std_theta(mean_values, std_values)))
            plot_expected_data.append(np.degrees(
                np.arccos(expected_values[:, 2]/expected_magnitudes)))
            ylabels.append('$\\theta$ [°]')
        # angle phi (wrt to x-axis)
        elif flag == 'f':
            plot_mean_data.append(get_phi(mean_values, cut_phi_at_0=cut_phi_at_0))
            plot_std_data.append(np.degrees(estimate_std_phi(mean_values, std_values)))
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
            axs[i].plot(x_vals, plot_expected_data[i], linestyle='--', 
                                    label = 'actuation matrix method @ 3mm')
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

# %%
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