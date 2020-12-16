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


#%%
# estimate number of independent entries in partial derivative tensors, 
# since symmetry of partial derivatives yields symmetric tensors
k = 3
print(binom(3+k-1, k))


#  apply fitting to entire dataset ----------------------------------------------------
#%%
# read in data from measurement series and discard data below threshold or above saturation
directory = '../test_data/config_tests_20_11_10'
threshold = 10
currents, B_measured, B_expected = collectAndExtract(directory, threshold, remove_saturation=True)

#%%
# estimate actuation matrix by linearly fitting the data directly
A = np.zeros((3,3))
A_err = np.zeros((3,3))
for component in range(3):
    for coil in range(3):
        p, pcov = curve_fit(lambda x, a: a*x, currents[:, coil], B_measured[:, component], p0=[10])
        perr = np.sqrt(np.diag(pcov))
        A[component, coil] = p[0]
        A_err[component, coil] = perr[0]

print(A)

# estimate the predicted fields based on linear fitting
linear_fits = np.zeros_like(B_measured)

for i in range(len(currents)):
    linear_fits[i] = A @ currents[i]

# estimate and print performances
print('expected fields:')
evaluate_performance(B_measured, B_expected)

print('\nlinear fits:')
evaluate_performance(B_measured, linear_fits)




# %%
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
        B = xdata @ a

    # quadratic fit including cross terms
    elif len(a) == 9:
        # linear part
        B_linear = xdata @ a[:3]

        # build upper-triangular matrix from remaining parameters
        A2 = np.zeros((3, 3))
        A2[np.triu_indices(3)] = a[3:]
        B_quadratic = np.diag(xdata @ A2 @ xdata.T)

        B = B_linear + B_quadratic

    else:
        raise NotImplementedError

    return B - ydata_1d

# test fitting for a single component
x = np.ones((5,3))
for i in range(5):
    x[i] = i
y = x @ np.array([1,2,1])
# print(y)

p0 = 10 * np.zeros(9)
p0[:3] = 10

# p = least_squares(lambda a: residuals(x, y, a), p0)
results = least_squares(lambda a: residuals(currents, B_measured[:,0], a), p0)
print(results.x)
# print(x @ results.x)

    
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
        for i in range(3):
            results = least_squares(lambda a: residuals(xdata, ydata[:,i], a), p0)
            A[i] = results.x

        # estimate expected fields based on linear fits 
        fits = np.zeros_like(ydata)
        for i in range(len(currents)):
            fits[i] = A @ currents[i]


    elif degree == 2:
        p0 = np.zeros(9)
        p0[:3] = 10

        # define actuation matrix, first dimension for field components, second for coils
        A = np.zeros((3,9))

        # estimate quadratic actuation matrix component-wise along (first) dimension for field component
        for i in range(3):
            results = least_squares(lambda a: residuals(xdata, ydata[:,i], a), p0)
            A[i] = results.x

        # estimate expected fields based on linear fits 
        fits = np.zeros_like(ydata)
        for i in range(len(currents)):
            # linear part
            fits[i] = A[:,:3] @ currents[i]

            # quadratic part 
            B_quadratic = np.zeros(3)
            for component in range(3):
                A_tri = np.zeros((3, 3))
                A_tri[np.triu_indices(3)] = A[component, 3:]
                B_quadratic[component] = currents[i].reshape(1,3) @ A_tri @ currents[i].reshape(3,1)

            fits[i] += B_quadratic

    else:
        raise NotImplementedError

    return A, fits


# fit measurement data 
degree = 2
print(f'degree: {degree}')

A, B_fits = fit_data(currents, B_measured, degree)
print(A)

evaluate_performance(B_measured, B_fits)



#%%
# testing stuff

# I = np.array([[1,2,3], [1,1,4]])
I = np.array([1,2,3])
B = np.array([  [1,0,0],
                [0,1,1],
                [0,0,1]])

# res = I @ B @ I.T
res = I.reshape(1,3) @ B @ I.reshape(3,1) 

print(np.diag(res))

# %%
