""" 
filename: analysis_tools.py

This file contains functions that are used for fitting and interpolating the measurement data.

Author: Nicholas Meinhardt (Qzabre)
        nmeinhar@student.ethz.ch
        
Date: 09.11.2020
"""

#%%
# standard library imports
import numpy as np
from numpy.linalg import norm
import os
from scipy.optimize import curve_fit
from scipy.spatial import ConvexHull
from math import isnan
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# local imports
from modules.analysis_tools import *

#%%
def find_start_of_saturation(I, mean_values, std_values, component = None, verbose = False,
                            fitting_fct = lambda x, a, b: a*x+b, fraction_cutoff = 0.02):
    """
    Find the indices that represent the start of saturation for both current directions. 

    Args:
    - I, mean_values, std_values, expected_values (ndarrays) are of shape (#measurements, 3), 
    containing applied current, experimentally estimated/expected mean values and standard deviations 
    for x,y,z-directions.
    - component (None or int): field-component that should be considered for the fit, i.e. 0 for the x-component,
    1 for y-component and 2 for z-component. 
    - fitting_fct (function): function that should be used for fitting. Default is an affine transformation of the 
    form y = ax + b. 
    - fraction_cutoff (float): fraction of the maximum absolute value of measured fields along the provided component
    that should be taken as a maximum acceptable cutoff for differing between linear and non-linear regime. 
    Measurement values that deviate from the linear fit by more than fraction_cutoff * (maximum absolute value) are
    classified as "out of the linear regime", while fields with smaller deviations are classified as 
    "still inside the linear regime". The larger the fraction, the more data included in the 'linear regime', 
    resulting in larger saturation values
    - verbose (bool): flag to switch on/off printing additional information on the estimated fitting parameters

    Return: i_min, i_max (int): indices where the linear regime starts and ends
    """

    # if component is not passed, take the component (along 2nd axis) that contains the maximum absolute value
    if component is None:
        component = np.argmax(np.max(np.abs(mean_values), axis=0))

    # for not having to deal with different currents in different coils, consider equally distributed current data in 
    # a normalized range of [-1, 1]
    N = len(I)
    x = np.linspace(-1,1, N)

    # make initial guess for fitting dependent on number of arguments the fitting function takes
    if fitting_fct.__code__.co_argcount == 2:
        p0 = [10]
    elif fitting_fct.__code__.co_argcount < 2:
        raise ValueError
    else:
        p0 = np.append([10], np.zeros(fitting_fct.__code__.co_argcount-2))

    # starting from middle plus minus 10 data points, fit measurement data and estimate the sum of the square errors on
    # the selected data. Repeat this scheme while adding more data to both sides, such that equally many data points are
    # considerd on left and right of middle. Later, the fit with smallest sum of square errors will be chosen.
    params = []
    errors = []
    distances =[]
    for aditional_data in range(10, N//2):
        # set the currently considered range of data points
        i_min, i_max = N//2 - aditional_data, N//2 + aditional_data+1
        # fit the data within this range
        p, pcov = curve_fit(fitting_fct, x[i_min:i_max], mean_values[i_min:i_max, component], 
                                p0 = p0, sigma = std_values[i_min:i_max, component])
        # collect fitting parameters and compute sum of square distances (divided by number of considered data points)
        params.append(p)
        distances.append(np.sum((mean_values[i_min:i_max,component]-fitting_fct(x[i_min:i_max],*p))**2))
        # to speed up next guess, set the initial guess for next fitting to the current parameters
        p0 = p

    # transform lists to ndarrays
    params = np.array(params)
    distances = np.array(distances)

    # find the fit with smallest sum of squares errors (objective of least-squares optimization)
    i_chosen = np.argmin(distances / (2*np.arange(10, N//2)+1))
    
    # estimate deviation of prediction and measurements
    deviations = np.abs((mean_values[:, component] - fitting_fct(x[:],*params[i_chosen])))

    # find the values on right and left for which the deviations reach 1 mT
    cutoff = fraction_cutoff * np.max(np.abs(mean_values[:,component]))
    i_min = np.argmax(deviations < cutoff)
    i_max = len(deviations) - np.argmax(np.flip(deviations) < cutoff) - 1

    if verbose:
        print('chosen index: {}, data points considered: {}'.format(i_chosen, 2*(10+i_chosen)+1))
        print('fit parameters chosen: {}'.format(params[i_chosen]))
        print('start saturation at: {:.2f} A and {:.2f} A'.format(x[i_min], x[i_max]))
        print('in between there are {} data points'.format(i_max-i_min+1))
                            
    return i_min, i_max

def apply_fitting_to_files_in_directory(directory, verbose = False, show_plots = False, fraction_cutoff = 0.02,
                                    affine_fct = lambda x, a, b: a*x + b):
    """
    Apply the linear fitting scheme of estimate_linear_fit_params_and_boundaries
    to all data files in the provided directory and collect the relevant values, such as 
    field directions and fitted slopes. 

    Args:
    - directory (string): directory that contains the data files in csv-format
    - verbose (bool): If True, the fitted slopes and offsets are printed for each data file
    - show_plots (bool): If True, plots of the field components versus one of the coils are generated,
    including the fitted affine functions. 
    - fraction_cutoff (float): Passed to find_start_of_saturation, consult documentation 
    of this function for details 
    - affine_fct (function): function used for fitting the data, passed to find_start_of_saturation

    Return:
    - list_filenames (1d-ndarray of length num_files): Filenames of all considered data files
    - list_slopes, list_slopes_std (ndarrays of shape (num_files, 3, 3)):
    contain the fitted slopes as well as their respective standard deviations. 
    - list_offsets, list_offsets_std (ndarrays of shape (num_files, 3)):
    contain the fitted offsets as well as their respective standard deviations. 
    - list_B_directions (ndarray of shape (num_files, 3)): contains the respective 
    field directions as vectors. 
    - list_phi, list_phi_std, list_theta, list_theta_std (1d-ndarrays of length num_files): contain 
    the respective field directions as polar and azimuthal angles, including their standard deviations.
    - list_B_max_1, list_B_max_2 (1d-ndarrays of length num_files): contain magnetic fields 
    at the left (1) and right (2) boundary of the linear regime

    Note: 
    - The directions provided as vector or as angles correspond to the slopes 
    for the first valid coil (i.e. with not-nan entries). Thus, the sign of the direction depends 
    on the considered coil. 
    - The standard deviations are estimated when fitting the data by an affine function.
    """
    # collect all csv-files in this directory
    filenames = []
    [filenames.append(file) for file in os.listdir(directory) if file.endswith(".csv")]

    # if in list, remove the linear_fits file from previous fits
    try:
        filenames.remove('linear_fits.csv')
    except ValueError:
        pass

    # initialize lists to store the relevant data
    list_filenames = []
    list_B_directions = []
    list_phi = []
    list_phi_std = []
    list_theta = []
    list_theta_std = []
    list_slopes = []
    list_slopes_std = []
    list_offsets = []
    list_offsets_std = []
    list_B_max_1 = []
    list_B_max_2 = []

    # loop through all csv files in a dictionary and fit the data
    for data_filename in filenames:
        if verbose:
            print(data_filename)

        # fit data
        results = estimate_linear_fit_params_and_boundaries(directory, data_filename, 
                                                        verbose = verbose, fraction_cutoff = fraction_cutoff,
                                                        affine_fct = affine_fct)
        B_direction, phi, phi_std, theta, theta_std, slopes, slopes_std, offsets, offsets_std, B_max_1, B_max_2 = results

        # collect fit parameters
        list_filenames.append(data_filename)
        list_B_directions.append(B_direction)
        list_phi.append(phi)
        list_phi_std.append(phi_std)
        list_theta.append(theta)
        list_theta_std.append(theta_std)
        list_slopes.append(slopes)
        list_slopes_std.append(slopes_std)
        list_offsets.append(offsets)
        list_offsets_std.append(offsets_std)
        list_B_max_1.append(B_max_1)
        list_B_max_2.append(B_max_2)

        # plot the results
        data_filepath = os.path.join(directory, data_filename)
        I, mean_data, std_data, expected_fields = extract_raw_data_from_file(data_filepath)
        i_min, i_max = find_start_of_saturation(I, mean_data, std_data, fraction_cutoff=fraction_cutoff,
                                fitting_fct = affine_fct)
        mask = np.array([not isnan(slopes[0, coil]) for coil in range(3)])
        valid_coil = np.arange(3)[mask][0]  
        x = I[:, valid_coil]

        if show_plots:
            fig, axs = plt.subplots(3)
            fig.set_size_inches(6, 8)
            for component in range(3):
                axs[component].plot(x, mean_data[:,component], linestyle=None, marker= '.', label = 'measured')
                axs[component].plot(x, x* slopes[component,valid_coil] + offsets[component], linestyle='--', 
                                label='fitted: ({:.2f} mT/A)*I + {:.2f} mT'.format(slopes[component,valid_coil], offsets[component]))
                ymin, ymax = axs[component].get_ylim()
                axs[component].vlines([x[i_min], x[i_max]], ymin, ymax, alpha=0.5, linewidth=1, 
                            label='linear regime in [{:.2f} A, {:.2f} A]'.format(I[i_min, 0], I[i_max, 0]))
                axs[component].hlines([B_max_1[component], B_max_2[component]], x[i_min], x[i_max], alpha=0.5, linewidth=1, 
                            label='linear regime in [{:.2f} mT, {:.2f} mT]'.format(B_max_1[component], B_max_2[component]))

                axs[component].set_ylabel('$B_{}$ [mT]'.format(['x', 'y', 'z'][component]))
                axs[component].set_xlabel('$I_1$ [A]')
                axs[component].set_ylim( 1.1*np.min(mean_data[:,component]), 1.1*np.max(mean_data[:,component]))
                axs[component].legend()
            fig.suptitle(data_filename)
            plt.tight_layout()
            plt.show()

    # convert lists to arrays
    list_B_directions = np.array(list_B_directions)
    list_slopes = np.array(list_slopes)
    list_slopes_std = np.array(list_slopes_std)
    list_offsets = np.array(list_offsets)
    list_offsets_std = np.array(list_offsets_std)
    list_phi = np.array(list_phi).flatten()
    list_phi_std = np.array(list_phi_std).flatten()
    list_theta = np.array(list_theta).flatten()
    list_theta_std = np.array(list_theta_std).flatten()
    list_B_max_1 = np.array(list_B_max_1)
    list_B_max_2 = np.array(list_B_max_2)

    return [list_filenames, list_slopes, list_slopes_std, list_offsets, list_offsets_std, 
            list_B_directions, list_phi, list_phi_std, list_theta, list_theta_std,
            list_B_max_1, list_B_max_2]


def estimate_linear_fit_params_and_boundaries(directory, filename, verbose = False, fraction_cutoff = 0.02,
                                    affine_fct = lambda x, a, b: a*x + b):
    """
    Extract data from provided file in directory and fit the data with linear function (and offset, so an affine transformation).
    Estimate the linear slopes for all components x,y,z and wrt. all coils 1,2,3, resulting in overall 9 parameters.
    Additionally, standard deviation of these slopes are estimated. When no slopes can be fitted, because the 
    current in one coil remains constant, the according entries are set to NaN. From the linear fit (and after ignoring
    the offsets) the direction of magnetic field is estimated. The accoring boundaries, at which
    the relation between field and current leaves the linear regime, are estimated and returned as well.

    Args: 
    - directory, filename (str): Valid directory and file that contains measurement data.
    - verbose (bool): if True, the fitted slopes and offsets are printed for all components (x,y,z) 
    and coils (1,2,3).
    - fraction_cutoff (float): Passed to find_start_of_saturation, consult documentation 
    of this function for details 
    - affine_fct (function): function used for fitting the data, passed to find_start_of_saturation

    Return:
    - slopes, slopes_std (ndarray of shape (3,3)): Contain the fitted slopes of 
    affine transformation, or the according standard deviation, respectively. 
    - offsets, offsets_std (ndarray of length 3): Contain the fitted offsets of 
    affine transformation, or the according standard deviation, respectively. Note that they are the same
    for all coils.
    - B-direction: normalized vector pointing towards the same direction 
    - B_max_1, B_max_2 (ndarrays of length 3): magnetic field at the left (1) and right (2) boundary 
    of the linear regime
    """
    # read in raw measurment data
    data_filepath = os.path.join(directory, filename)
    I, mean_data, std_data, expected_fields = extract_raw_data_from_file(data_filepath)

    # estimate minimum and maximum indices of region within which the linear relation between current and field holds  
    # even though find_start_of_saturation offers the possibility to specify the considered component, 
    # keep the default stting, which detects the field component that has the greatest absolute field values. 
    # This should work fine for situations, where one component is dominating. 
    i_min, i_max = find_start_of_saturation(I, mean_data, std_data, fraction_cutoff=fraction_cutoff,
                            fitting_fct = affine_fct)

    # initialize array to collect the linear coeffiecients (and offsets) for all three components (and coils)
    slopes = np.zeros((3,3))
    slopes_std = np.zeros((3,3))
    offsets = np.zeros((3,3))
    offsets_std = np.zeros((3,3))

    # loop over the three currents in the coils 
    for coil in range(3):

        # if the current in one coil remains constant, it is useless to fit a linear slope with this current as x-value
        # thus, set the slopes and offsets to nan, to distinguish it from actual fits
        if np.all(I[i_min:i_max, coil] == I[i_min, coil]):
            slopes[:, coil] =  np.ones(3) * np.nan
            slopes_std[:, coil] = np.ones(3) * np.nan
            offsets[:, coil] = np.ones(3) * np.nan
            offsets_std[:, coil] = np.ones(3) * np.nan

        else: 
            # loop over the three components and estimate slope and offset each
            for component in range(3):
                # fit the data within i_min and i_max
                p, pcov = curve_fit(affine_fct, I[i_min:i_max, coil], mean_data[i_min:i_max, component], 
                                            p0=[10,0], sigma=std_data[i_min:i_max, component])
                perr = np.sqrt(np.diag(pcov))

                slopes[component, coil] = p[0]
                slopes_std[component, coil] = perr[0]
                offsets[component, coil] = p[1]
                offsets_std[component, coil] = perr[1]

    # estimate field direction from fitted slopes (i.e. not the nan-values):
    mask = np.array([not isnan(slopes[0, coil]) for coil in range(3)])
    valid_coil = np.arange(3)[mask][0]
    B_direction = slopes[:, valid_coil] 

    # note that offsets in one component are the same for all coils, thus reduce to 1d-arrays
    offsets = offsets[:, valid_coil]
    offsets_std = offsets_std[:, valid_coil]

    # estimate direction in spherical coordinates and the maximum magnetic field for which the linear fit works
    phi = get_phi(B_direction)
    phi_std = estimate_std_phi(B_direction.reshape(1,3), slopes_std[:, valid_coil].reshape(1,3))
    theta = get_theta(B_direction)
    theta_std = estimate_std_theta(B_direction.reshape(1,3), slopes_std[:, valid_coil].reshape(1,3))

    # if desired, print all relevant outputs to terminal
    if verbose:
        # print slopes for all components and all coils
        for component in range(3):
            axis = ['x', 'y', 'z'][component]
            print('along {}: a = {} +- {}'.format(axis, np.round(slopes[component], 3), np.round(slopes_std[component], 3)))
        # print offsets for all components but only a single coil
        for component in range(3):
            axis = ['x', 'y', 'z'][component]
            print('along {}: b = {} +- {}'.format(axis, np.round(offsets[component],3), 
                                                        np.round(offsets_std[component],3)))
        # print field direction both as vector and as angular values
        print('direction (vector): {}'.format(np.round(B_direction, 3)))
        print('direction (angular): theta = {:.2f}°, phi = {:.2f}°'.format(theta, phi))

    # estimate the maximum magnetic field (strength) for within which the linear regime holds
    B_max_1 = B_direction * I[i_min, valid_coil]
    B_max_2 = B_direction * I[i_max, valid_coil]

    return B_direction, phi, phi_std, theta, theta_std, slopes, slopes_std, offsets, offsets_std, B_max_1, B_max_2

def delaunay_triangulation_spherical_surface(phi, theta):
    """
    Delaunay triangulation on the surface of the sphere of the points defined by theta and phi. 
    This is simply the 3D convex hull of the points, since all vectors have the same length.
    
    Args:
    - phi, theta (1d ndarrays of same length n_data): polar and azimuthal angles of the data points. 
    
    Return:
    - vertices_triangulation (ndarray of shape (nfacet, 3, 2)): Represents the vertices of the 
    Delaunay triangulation on the sphere, where nfacet is the number of simplices in the 
    simplicial complex built from the provided set of points. The second dimension contains
    the three vertices of a simplex, and the last dimension covers the corresponding azimuthal angle
    theta as first entry and the planar angle phi as second entry. 
    - indices_triangulation (ndarray of shape (nfacet, 3)): Integer array containing the indeces 
    of the data points that form the vertices of the simplices.
    - points (ndarray of shape (n_data, 3)): All points on the sphere that correspond to the provided
    angles in Carthesian coordinates. 
    """
    # generate vectors on unit sphere from angles theta and phi
    x = np.sin(np.radians(theta)) * np.cos(np.radians(phi))
    y = np.sin(np.radians(theta)) * np.sin(np.radians(phi))
    z = np.cos(np.radians(theta))
    points = np.swapaxes([x, y, z], 0,1)
    
    # generate convex hull from all points
    hull = ConvexHull(points)
    
    # loop through simplices on surface of convex hull, which are supposed to be trianlges 
    # and collect the corresponding angles theta and phi
    vertices_triangulation = [] 
    indices_triangulation = []
    for simplex in hull.simplices:
        vertices_triangulation.append([theta[simplex], phi[simplex]])
        vertices_triangulation.append([theta[simplex], phi[simplex]])
        indices_triangulation.append(simplex)
    vertices_triangulation = np.swapaxes(vertices_triangulation, 1, 2)

    return vertices_triangulation, np.array(indices_triangulation), points

def add_reverse_directions(phis, phis_std, thetas, thetas_std, slopes, slopes_std, offsets, offsets_std):
    """
    Add the reverse directions to the angle, slope and offset arrays and make sure they do not contain these
    reverse directions already. 

    Args:
    - phis, phis_std, thetas, thetas_std (1d-ndarrays of length num_files): contain 
    the respective field directions as polar and azimuthal angles, including their standard deviations.
    - slopes, slopes_std, offsets, offsets_std (ndarrays of shape (num_files, 3, 3)):
    contain the fitted slopes and offsets as well as their respective standard deviations. 

    Return:
    - phis, phis_std, thetas, thetas_std (1d-ndarrays of length 2*num_files): contain 
    the respective field directions as polar and azimuthal angles, including their standard deviations,
    which are now complemented by the reverse directions: 
    phi -> 180° + phi (note that cut is at 180° by default)
    theta -> 180° - theta
    - slopes, slopes_std, offsets, offsets_std (ndarrays of shape (2*num_files, 3, 3)):
    contain the fitted slopes and offsets as well as their respective standard deviations,
    which are now complemented by the (negative) slopes and (same) offsets for the reverse directions. 
    """
    # check whether the reversed direction of first entry in thetas is already contained in thetas
    # in this way, it is prevented to add the reverse directions twice, while assuming that the 
    # measured values are never exactly the opposite of one of the other values.
    if 180 - thetas[0] in thetas:
        print('Most likely, the arrays have been complemeted by reverse directions already')
        return phis, phis_std, thetas, thetas_std, slopes, slopes_std, offsets, offsets_std

    # estimate number of measurements
    N = len(phis)

    # initialize complemented arrays by repeating the original ones
    # note that angles are 1d arrays, while slopes and offsets are of shape (N,3,3)
    phis_filled = np.tile(phis, 2)
    phis_std_filled = np.tile(phis_std, 2)
    thetas_filled = np.tile(thetas, 2)
    thetas_std_filled = np.tile(thetas_std, 2)

    slopes_filled = np.tile(slopes, [2,1,1])
    slopes_std_filled = np.tile(slopes_std, [2,1,1])
    offsets_filled = np.tile(offsets, [2,1,1])
    offsets_std_filled = np.tile(offsets_std, [2,1,1])
    
    # update the second half of filled up arrays by the reverse directions (stds remain the same)
    for i in range(N):
        # simply invert the angles. To stay consistent with cut of phi at +-180°, 
        # conditional setting required (even though it would not be necessary)
        phis_filled[N+i] = 180 + phis[i] if phis[i] <= 0 else -180 + phis[i]
        thetas_filled[N+i] = 180 - thetas[i]

        # while offsets and standard deviations of both slopes and offsets remain the same, 
        # the slopes get an overall negative sign
        slopes_filled[N+i] = - slopes[i] 

    return phis_filled, phis_std_filled, thetas_filled, thetas_std_filled, slopes_filled, slopes_std_filled, offsets_filled, offsets_std_filled


def add_triangles_to_3dplot(ax, pts, inidces_simplices, spherical = True):
    """
    Plot the triangles of Delaunay triangulation by connecting the three vertices of each triangle
    with each user. These connections can either be parts of great circles on the sphere (spherical=True)
    or straight lines (spherical=False), which are convex combinations of the points they connect. 
    Note that in the latter case, the lines correspond to the edges of the convex hull of all points.
  
    Args:
    - ax (Axes3D): Axis of plot with 3d projection
    - pts (ndarray of shape (N, 3)): all vertices, which are points on the unit sphere
    - indices_triangulation (ndarray of shape (nfacet, 3)): Integer array containing the indeces 
    of the data points that form the vertices of the simplices.
    - spherical (bool): Flag to switch between straight lines and great circles as edges of the triangles.
    """
    fractions = np.linspace(0, 1, 20, endpoint=False)

    for s in inidces_simplices:
        # Cycle back to the first coordinate to plot the whole triangle
        s = np.append(s, s[0])  
        
        if spherical:
            # initialize array to store all positions along the path
            lines = []
            for i in range(len(s)-1):
                direction = pts[s[i+1]] - pts[s[i]]
                [lines.append(pts[s[i]] + r*direction) for r in fractions]
            # append the first point to close the loop 
            lines.append(lines[0])
            
            # normalize all points along path, such that they point towards the sphere
            lines = lines/ np.linalg.norm(lines, axis=1).reshape(-1,1)

            ax.plot(lines[:, 0], lines[:, 1], lines[:, 2], 'r-', alpha=0.4)

        else:
            ax.plot(pts[s, 0], pts[s, 1], pts[s, 2], 'r-', alpha=0.4)


def save_linear_slopes_to_file(filenames, slopes, slopes_std, offsets, offsets_std, B_directions, 
                                phis, phis_std, thetas, thetas_std, B_max_1, B_max_2, save_directory=None,
                                output_file_name = 'linear_fits.csv'):
    """
    Save all data that were collected from a series of measurements to a single csv file. 

    Args:
    - filenames (ndarray of length N): file names, where estimated slopes, offsets,... originate from
    - slopes, slopes_std (ndarrays of shape (N,3,3)): Fitted slopes and their standard deviations
    - offsets, offsets_std (ndarrays of shape (N,3)): Fitted offsets and their  standard deviations
    - B_directions (ndarray of shape (N,3)): vectorial direction of B-field
    - phis, phis_std, thetas, thetas_std (ndarray of length N): contain the angular direction of 
    magnetic field and the respective standard deviations
    - B_max_1, B_max_2 (ndarray of length N): magnetic fields at the left (1) and right (2) boundary 
    of the linear regime, estimated for the first valid coil (i.e. with nonzero current)
    - output_file_name (str): Name of csv file where the result are stored
    """
    # save data to file 
    data_dict = {'Filename': filenames, 
                'phi [°]': phis,
                'std phi [°]': phis_std,
                'theta [°]': thetas,
                'std theta [°]': thetas_std}

    # add all slopes
    component_labels = ['x', 'y', 'z']
    for component in range(3):
        for coil in range(3):
            label ='dB{}/dI{} [mT/A]'.format(component_labels[component], coil+1)
            data_dict[label] = slopes[:, component, coil]
            data_dict['std {} [mT]'.format(label)] = slopes_std[:, component, coil]

    # add all offsets
    for component in range(3):  
        label = 'offset Delta B{} [mT]'.format(component_labels[component])
        data_dict[label] = offsets[:, component]
        data_dict['std {} [mT]'.format(label)] = offsets_std[:, component]

    # add the direction of magnetic field as vector
    for component in range(3):  
        label = 'direction B{} [mT]'.format(component_labels[component])
        data_dict[label] = B_directions[:, component]
    
    # add the magnetic field at left boundary of linear regime
    for component in range(3):  
        label = 'max B{} left boundary [mT]'.format(component_labels[component])
        data_dict[label] = B_max_1[:, component]

    # add the magnetic field at right boundary of linear regime
    for component in range(3):  
        label = 'max B{} right boundary [mT]'.format(component_labels[component])
        data_dict[label] = B_max_2[:, component]

    # save data to csv file
    df = pd.DataFrame(data_dict)

    if save_directory is None:
        save_directory = os.getcwd()
    data_filepath = os.path.join(save_directory, output_file_name)
    df.to_csv(data_filepath, index=False, header=True)


def extract_fit_parameters_from_file(filepath):
    """
    Extract all fit paramters from a csv-file, which originates from an already 
    analyzed set of measurements.

    Args: filepath (str): valid path of the csv-file

    Return: 
    - filenames (ndarray of length N): file names, where estimated slopes, offsets,... originate from
    - slopes, slopes_std, offsets, offsets_std (ndarrays of shape (N,3,3)): Fitted slopes and offsets, 
    and the respective standard deviations
    - B_directions (ndarray of shape (N,3)): vectorial direction of B-field
    - phis, phis_std, thetas, thetas_std (ndarray of length N): contain the angular direction of 
    magnetic field and the respective standard deviations
    - B_max_1, B_max_2 (ndarray of length N): magnetic fields at the left (1) and right (2) boundary 
    of the linear regime, estimated for the first valid coil (i.e. with nonzero current)
    - output_file_name (str): Name of csv file where the result are stored
    """
    # extract data and convert to ndarray
    raw_data = pd.read_csv(filepath).to_numpy()

    # assign values to correct variables
    filenames = raw_data[:, 0]

    # conversion to np.float needed, because np.degrees(phis) would fail else
    phis = raw_data[:, 1].astype(float)
    phis_std = raw_data[:, 2].astype(float)
    thetas = raw_data[:, 3].astype(float)
    thetas_std = raw_data[:, 4].astype(float)

    slopes = raw_data[:, 5:23:2].reshape(-1, 3,3)
    slopes_std = raw_data[:, 6:23:2].reshape(-1, 3,3)

    offsets = raw_data[:, 23:29:2]
    offsets_std = raw_data[:, 24:29:2]
    
    B_directions = raw_data[:, 29:32]
    B_max_1 = raw_data[:, 32:35]
    B_max_2 = raw_data[:, 35:38]

    return filenames, slopes, slopes_std, offsets, offsets_std, B_directions, phis, phis_std, thetas, thetas_std, B_max_1, B_max_2

def add_triangles_to_3dplot(ax, pts, inidces_simplices, spherical = True, colored_triangles = False):
    """
    Plot the triangles of Delaunay triangulation by connecting the three vertices of each triangle
    with each user. These connections can either be parts of great circles on the sphere (spherical=True)
    or straight lines (spherical=False), which are convex combinations of the points they connect. 
    Note that in the latter case, the lines correspond to the edges of the convex hull of all points.
  
    Args:
    - ax (Axes3D): Axis of plot with 3d projection
    - pts (ndarray of shape (N, 3)): all vertices, which are points on the unit sphere
    - indices_triangulation (ndarray of shape (nfacet, 3)): Integer array containing the indeces 
    of the data points that form the vertices of the simplices.
    - spherical (bool): Flag to switch between straight lines and great circles as edges of the triangles.
    - colored_triangles (bool): if True, colored triangles are drawn onto the sphere. This only works if 
    spherical=True, too. Default is False, where no surfaces are drawn.
    """
    fractions = np.linspace(0, 1, 10, endpoint=False)
    np.random.seed(100)

    for s in inidces_simplices:
        # Cycle back to the first coordinate to plot the whole triangle
        s = np.append(s, s[0])  
        
        if spherical:
            # initialize array to store all positions along the path
            lines = []
            for i in range(len(s)-1):
                direction = pts[s[i+1]] - pts[s[i]]
                [lines.append(pts[s[i]] + r*direction) for r in fractions]
            # append the first point to close the loop 
            lines.append(lines[0])
            
            # normalize all points along path, such that they point towards the sphere
            lines = lines/ np.linalg.norm(lines, axis=1).reshape(-1,1)
    
            ax.plot(lines[:, 0], lines[:, 1], lines[:, 2], color='gray', 
                                        linestyle='-', alpha=0.4, lw=0.5)
            
            if colored_triangles:
                random_color = colors.rgb2hex(np.random.rand(3))
                polygon = Poly3DCollection([lines], alpha=1.0)
                polygon.set_color(random_color)
                ax.add_collection3d(polygon)

        else:
            ax.plot(pts[s, 0], pts[s, 1], pts[s, 2], 'r-', alpha=0.4)