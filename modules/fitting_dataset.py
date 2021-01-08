""" 
filename: estimate_derivatives.py

This script is meant to read in measurement data from a set of measurements, 
p.e. config_tests_20_11_10, extract all data above a threshold magnitude of 10 mT 
and below saturation and fit the remaining data.

Author: Nicholas Meinhardt (QZabre)
        nmeinhar@student.ethz.ch
        
Date: 16.12.2020
"""


#%%
# standard library imports
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.optimize import curve_fit, least_squares
from scipy.special import binom

# local imports
try:
    from modules.analysis_tools import *
except ModuleNotFoundError:
    import sys
    sys.path.insert(1, os.path.join(sys.path[0], '..'))
finally:
    from modules.analysis_tools import *
    from modules.interpolation_tools import *
    from modules.general_functions import estimate_RMS_error


#%%
def fit_data(xdata, ydata, degree):
    """ 
    Fit the provided data with a tri-linear or tri-quadratic polynomial and 
    return the fitting paramters and fitted y-values. Fitting is performed 
    component-wise and the estimated fitting parameters as well as the fitted 
    results are combined afterwards.

    Note: so far only trilinear and tri-quadratic fitting is implemented. 

    Args: 
    - currents (ndarray of shape (N, 3)): input data containing the currents 
    for all three coils.
    - ydata (ndarray of shape (N, 3)): output data containing all three components 
    of the magnetic field.
    - degree (int): degree of polynomial used for fitting.

    Return:
    - A (ndarray of shape (3,m)): fitting paramters for the three field components.
    The number m depends on the degree of the polynomial: 
     - m=3 for degree 1, 
     - m=9 for degree 2, where the first 3 entries correspond to the linear terms 
     and the rest to the quadratic terms.
    """
    if degree == 1:
        p0 = np.array([10,10,10])

        # define actuation matrix, first dimension for field components, second for coils
        A = np.zeros((3,3))

        # estimate linear actuation matrix component-wise along (first) dimension for field component
        for component in range(3):
            results = least_squares(lambda a: residuals(xdata, ydata[:,component], a), p0)
            A[component] = results.x

    elif degree == 2:
        p0 = np.zeros(9)
        p0[:3] = 10

        # define actuation matrix, first dimension for field components, second for coils
        A = np.zeros((3,9))

        # estimate quadratic actuation matrix component-wise along (first) dimension for field component
        for component in range(3):
            results = least_squares(lambda a: residuals(xdata, ydata[:,component], a), p0)
            A[component] = results.x

    elif degree == 3:
        p0 = np.zeros(19)
        p0[:3] = 10

        # define actuation matrix, first dimension for field components, second for coils
        A = np.zeros((3,19))

        # estimate quadratic actuation matrix component-wise along (first) dimension for field component
        for component in range(3):
            results = least_squares(lambda a: residuals(xdata, ydata[:,component], a), p0)
            A[component] = results.x

    else:
        raise NotImplementedError
    
    # estimate fitting results based on the fitting paramters
    fits = evaluate_fit(A, xdata)

    return A, fits

def residuals(xdata, ydata_1d, a):
    """
    Estimate the residuals, i.e. the differences between the actual ydata and 
    the results estimated from xdata by using parameters a. This function is intended 
    to be used as argument for scipy.optimize.least_squares as:
    lambda a: residuals(applied_currents, measured_fields, a)

    Note that so far only linear and quadratic fitting are implemented.
    Also note that ydata should be a 1d array, so only one component of the 
    magnetic field is considered here. 

    Args:
    - xdata (ndarray of shape (N, 3)): input data
    - ydata_1d (ndarray of length N): output data for a single component of magnetic field
    - a (1d ndarray): Fitting paramters, which are optimized by least_squares to minimize
    the sum of square residuals. Pass an array of length 3 to apply a linear function to 
    xdata, and an array of length 9 for a quadratic function

    Return:
    - res = B_computed - ydata_1d (ndarray of length N): Differences between computed 
    and actual results.

    Raise NotImplementedError if nonvalid paramteres a are passed.
    """   
    # linear fit
    if len(a) == 3:
        fit = xdata @ a

    # quadratic fit including cross terms
    elif len(a) == 9:
        # linear part
        fit_linear = xdata @ a[:3]

        # build upper-triangular matrix from remaining parameters
        A2 = np.zeros((3, 3))
        A2[np.triu_indices(3)] = a[3:]
        fit_quadratic = np.diag(xdata @ A2 @ xdata.T)

        fit = fit_linear + fit_quadratic

    # cubic fit including cross terms
    elif len(a) == 19:
        # linear part
        fit_linear = xdata @ a[:3]

        # build upper-triangular matrix from remaining parameters
        A2 = np.zeros((3, 3))
        A2[np.triu_indices(3)] = a[3:9]
        fit_quadratic = np.diag(xdata @ A2 @ xdata.T)

        # cubic part 
        # unique combinations of products, in total 10 for cubic 
        combis = np.array([[0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 1, 1], [0, 1, 2], 
                [0, 2, 2], [1, 1, 1], [1, 1, 2], [1, 2, 2], [2, 2, 2]])
        fit_cubic = np.zeros(len(xdata))
        for i in range(10):
            fit_cubic += a[i+9] * xdata[:,combis[i,0]] * xdata[:,combis[i,1]] * xdata[:,combis[i,2]]

        fit = fit_linear + fit_quadratic + fit_cubic

    else:
        raise NotImplementedError

    return fit - ydata_1d

def evaluate_fit(A, xdata):
    """
    Given the paramters A, estimate the outputs for the provided xdata.

    Args:
    - A (ndarray of shape (k,3)): fitting paramters for the three components individually. 
    The size k of the first dimension depends on the fitting degree. Currently implemented are
    k=3,9,19 only. 
    - xdata (ndarray of shape (N,3)): Input data for which the corresponding output data should be computed.

    Returns:
    - fits (ndarray of shape (N,3)): Estimated outputs based on the inputs and fitting paramters
    """
    # initialize array for fits
    fits = np.zeros_like(xdata)

    # linear fit
    if A.shape[1] == 3:
        for i in range(len(xdata)):
            fits[i] = A @ xdata[i]

    elif A.shape[1] == 9:
        # estimate expected fields based on fits 
        
        for i in range(len(xdata)):
            # linear part
            fits[i] = A[:,:3] @ xdata[i]

            # quadratic part 
            fits_quadratic = np.zeros(3)
            for component in range(3):
                A_tri = np.zeros((3, 3))
                A_tri[np.triu_indices(3)] = A[component, 3:]
                fits_quadratic[component] = xdata[i].reshape(1,3) @ A_tri @ xdata[i].reshape(3,1)

            fits[i] += fits_quadratic

    # cubic fit including cross terms
    elif A.shape[1] == 19:
        for i in range(len(xdata)):
            # linear part
            fits[i] = A[:,:3] @ xdata[i]

            # quadratic part 
            fits_quadratic = np.zeros(3)
            for component in range(3):
                A_tri = np.zeros((3, 3))
                A_tri[np.triu_indices(3)] = A[component, 3:9]
                fits_quadratic[component] = xdata[i].reshape(1,3) @ A_tri @ xdata[i].reshape(3,1)
            fits[i] += fits_quadratic

            # cubic part 
            combis = np.array([[0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 1, 1], [0, 1, 2], 
                [0, 2, 2], [1, 1, 1], [1, 1, 2], [1, 2, 2], [2, 2, 2]])
            fits_cubic = np.zeros(3)
            for component in range(3):
                for k in range(10):
                    fits_cubic[component] += A[component,k+9] * xdata[i, combis[k,0]] * \
                                                xdata[i, combis[k,1]] *   \
                                                xdata[i, combis[k,2]]
            fits[i] += fits_cubic

    return fits

def save_fitting_parameters(A, origin, degree, direction='B2I'):
    """
    Save fitting paramters A to a file in './fitting_paramters' as csv-file.

    Args:
    - A (ndarray of shape (3,k)): fitting paramters for the three field components or currents
    - origin (str): name of the folder that was used for fitting. The passed name will be the first
    part of the file name. 
    - degree (int): fitting degree, it will be used as second part of the file name
    - direction (str): direction of fitting, it will be used as last part of the file name.
    Reasonable values are 'B2I' and 'I2B', depending on which quantity is considered as output/input. 
    Accordingly, the unit of the fitting paramters change 
    """
    param_file_name = f'{origin}_degree{degree}_{direction}.csv'
    data_filepath = os.path.join('../fitting_parameters', param_file_name)

    df = pd.DataFrame({ 'A_x':  A[0, :], 
                        'A_y':  A[1, :], 
                        'A_z':  A[2, :]})
    

    df.to_csv(data_filepath, index=False, header=True)

def save_predicted_currents(currents, filepath):
    """
    Save predicted currents to a file in './fitting_paramters' as csv-file.

    Args:
    - currents (ndarray of shape (N,3)): estimated currents
    - filepath (str): path of created file
    """
    df = pd.DataFrame({ 'I_1 [A]':  currents[:, 0], 
                        'I_2 [A]':  currents[:, 1],
                        'I_3 [A]':  currents[:, 2]})
    df.to_csv(filepath, index=False, header=True)

def read_fitting_parameters(filepath):
    """
    Extract fitting paramters A from file.
    """
    # extract data and convert to ndarray
    A = pd.read_csv(filepath).to_numpy().T

    return A

def read_test_set(filepath):
    """
    Extract test set consisting of magnetic field vectors from provided file.
    """
    # extract data and convert to ndarray
    test_vectors = pd.read_csv(filepath).to_numpy()

    return test_vectors[:,:3]

