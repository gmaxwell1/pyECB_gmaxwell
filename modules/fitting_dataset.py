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
# estimate number of independent entries in partial derivative tensors, 
# since symmetry of partial derivatives yields symmetric tensors
k = 3
print(binom(3+k-1, k))



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


#  apply fitting to entire dataset ----------------------------------------------------
#%%
# read in data, either from folder containing multiple measurement files or from a single measurement series

# # data from measurement series and discard data below threshold or above saturation
# directory = '../test_data/config_tests_20_12_19'
# threshold = 5
# currents, B_measured, B_expected = collectAndExtract(directory, threshold, remove_saturation=True)

# data from a single measurements series saved in one file
directory = '../predictions_on_test_set/measure_predictions'
filename = '20_12_27_20-03-07_config_tests_20_12_19_degree3_B2I_vectors_rng10_0-50mT_size1000_demagnet5A.csv'
data_filepath = os.path.join(directory, filename)

currents, B_measured, B_measured_std, B_expected = extract_raw_data_from_file(data_filepath)


#%%
# evaluate angular accuracy
dot = np.array([np.dot(B_measured[i], B_expected[i]) for i in range(len(B_measured))])
norms_measured = np.linalg.norm(B_measured, axis=1)
norms_fits = np.linalg.norm(B_expected, axis=1)
alphas = np.degrees(np.arccos(dot / (norms_measured * norms_fits)))

n, bins, patches = plt.hist(alphas, 50)
plt.xlabel('Angular Error, $\\alpha$ [°]')
plt.ylabel('Counts')
plt.text(15, 200, '$\\mu$ = {:.2f}°,\n'.format(np.mean(alphas))+ 
                '$\\sigma$ = {:.2f}°,\n'.format(np.std(alphas))+ 
                'min = {:.2f}°,\n'.format(np.min(alphas))+ 
                'max = {:.2f}°,\n'.format(np.max(alphas))+ 
                'median = {:.2f}°,\n'.format(np.median(alphas)))
plt.show()

#%%
# mask for selecting large angular deviations
threshold_angle = 15
mask = alphas >= threshold_angle
indices = np.nonzero(mask)[0]
# for i in indices:
#     print(f'i={i}: |B_exp| = {np.linalg.norm(B_expected[i]):.1f} mT ' +
#                 f'|B_meas| = {np.linalg.norm(B_measured[i]):.1f} mT ' +
#                 f'alpha = {alphas[i]:.1f}°')

n, bins, patches = plt.hist(np.linalg.norm(B_expected[indices], axis=1), 50)
plt.xlabel('expected B field, $B_{exp}$ [mT]')
plt.ylabel('Counts')
plt.show()

n, bins, patches = plt.hist(np.linalg.norm(B_measured[indices], axis=1), 50)
plt.xlabel('measured B field, $B_{meas}$ [mT]')
plt.ylabel('Counts')
plt.show()

n, bins, patches = plt.hist(get_theta(B_measured[indices]), 50)
plt.xlabel('theta, $\\theta$ [°]')
plt.ylabel('Counts')
plt.show()

n, bins, patches = plt.hist(get_phi(B_measured[indices]), 50)
plt.xlabel('theta, $\\phi$ [°]')
plt.ylabel('Counts')
plt.show()

#%%
plt.figure()
plt.plot(get_theta(B_measured), np.linalg.norm(B_measured, axis=1), marker='.', linestyle='')
plt.xlabel('$\\theta$ [°]')
plt.ylabel('measured B field, $B_{meas}$ [mT]')
plt.show()

plt.figure()
plt.plot(get_phi(B_measured), np.linalg.norm(B_measured, axis=1), marker='.', linestyle='')
plt.xlabel('$\\phi$ [°]')
plt.ylabel('measured B field, $B_{meas}$ [mT]')
plt.show()


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

print('actuation matrix from linear fits, individually estimated for each component and coil')
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



# --------------------------------------------------------------------------
# %%
# test fitting for a single component
x = np.ones((5,3))
for i in range(5):
    x[i] = i
y = x @ np.array([1,2,1])
# print(y)

p0 = np.zeros(9)
p0[:3] = 10

# p = least_squares(lambda a: residuals(x, y, a), p0)
results = least_squares(lambda a: residuals(currents, B_measured[:,0], a), p0)
print(results.x)
# print(x @ results.x)

    
# --------------------------------------------------------------------------
#%%
# fit measurement data with currents as inputs and magnetic fields as outputs
degree = 3
print(f'degree: {degree}')

A, B_fits = fit_data(currents, B_measured, degree)
print(A)

evaluate_performance(B_measured, B_fits)

# print(np.linalg.inv(A))

# --------------------------------------------------------------------------
#%%
# try opposite: fit currents as outputs for measured B_fields as inputs
degree = 3
print(f'degree: {degree}')

A, currents_fitted = fit_data(B_measured, currents, degree)
print(A)

# estimate RMS errors 
RMSE = estimate_RMS_error(currents.flatten(), currents_fitted.flatten())
print(f'RMS error fit: {RMSE:.4f} A')

mean = np.mean(np.abs(currents.flatten()- currents_fitted.flatten()))
std = np.std(np.abs(currents.flatten()- currents_fitted.flatten()))
median = np.median(np.abs(currents.flatten()- currents_fitted.flatten()))
print(f'mean +- std: {mean:.4f} A +- {std:.4f} A')
print(f'median: {median:.4f} A')

# save fitting paramters in file
save_fitting_parameters(A, os.path.split(directory)[-1], degree, direction='B2I')


# --------------------------------------------------------------------------
#%%
# make predictions for suitable currents for test set vectors using fits

# read in paramters from previous fits
# filename_fit_params = 'config_tests_20_12_19_degree3_B2I'
filename_fit_params = 'measure_predictions_degree3_B2I'
filepath_fit_params = f'../fitting_parameters/{filename_fit_params}.csv'
A = read_fitting_parameters(filepath_fit_params)
print(A)

# read in test set
filename_testset = 'walk_on_grid_max10_PointsPerDim7'
filepath_testset = f'../predictions_on_test_set/{filename_testset}.csv'
test_vectors = read_test_set(filepath_testset)
print(test_vectors.shape)

# apply fitting paramters/function to vectors of test set
predicted_currents = evaluate_fit(A, test_vectors)
filepath_currents = f'../predictions_on_test_set/predictions/{filename_fit_params}_{filename_testset}.csv'
save_predicted_currents(predicted_currents, filepath_currents)

#%%
# second option: test rotations
num = 100
rot_axis = 2
radius = 50
test_vectors = np.zeros((num, 3))
angles = np.linspace(0, 2*np.pi, num)
if rot_axis == 0:
    test_vectors[:,1] = radius*np.cos(angles)
    test_vectors[:,2] = radius*np.sin(angles)
elif rot_axis == 1:
    test_vectors[:,0] = radius*np.cos(angles)
    test_vectors[:,2] = radius*(-np.sin(angles))
elif rot_axis == 2:
    test_vectors[:,0] = radius*np.cos(angles)
    test_vectors[:,1] = radius*np.sin(angles)

# apply fitting paramters/function to vectors of test set
predicted_currents = evaluate_fit(A, test_vectors)
rot_axis_label = ['x','y','z'][rot_axis]
filepath_currents = f'../predictions_on_test_set/predictions/{filename_fit_params}_rotation_{rot_axis_label}_{radius}mT.csv'
save_predicted_currents(predicted_currents, filepath_currents)

#%%
# third option: test sweeps along main axes
num = 100
sweep_axis = 2
max_field = 50
test_vectors = np.zeros((num, 3))
test_vectors[:, sweep_axis] = np.linspace(-max_field, max_field, num)

# apply fitting paramters/function to vectors of test set
predicted_currents = evaluate_fit(A, test_vectors)
sweep_axis_label = ['x','y','z'][sweep_axis]
filepath_currents = f'../predictions_on_test_set/predictions/{filename_fit_params}_sweep_{sweep_axis_label}_axis_{max_field}mT.csv'
save_predicted_currents(predicted_currents, filepath_currents)


#%%
# plot predicted currents to check reasonability 
fig = plt.figure(figsize = plt.figaspect(1.))
ax = fig.add_subplot(111, projection='3d')

ax.scatter(predicted_currents[:, 0], predicted_currents[:, 1], predicted_currents[:, 2] )

# axis settings
ax.set_xlabel('$I_1$ [A]')
ax.set_ylabel('$I_2$ [A]')
ax.set_zlabel('$I_3$ [A]')

plt.tight_layout()
plt.show()


# ---------------------------------------------------------------------
#%%
# evaluate the performance of the fit on the training set

# read in measurement results
data_directory = '../predictions_on_test_set/measure_predictions'
data_filename = '21_01_05_03-45-27_measure_predictions_degree3_B2I_vectors_rng2222_0-50mT_size2000_demagnet5A.csv'
data_filepath = os.path.join(data_directory, data_filename)

I, B_measured, B_measured_std, expected = extract_raw_data_from_file(data_filepath)

# read in test set
filename_testset = 'vectors_rng2222_0-50mT_size2000'
filepath_testset = f'../predictions_on_test_set/{filename_testset}.csv'
test_vectors = read_test_set(filepath_testset)


# evaluate the performance
print('performance of cubic fit')
evaluate_performance(B_measured, test_vectors)

print('\nperformance of currently implemented model')
evaluate_performance(B_measured, expected)

# %%
# plot angular errors as histogram for all data and for the data >= 30 mT
dot = np.array([np.dot(B_measured[i], test_vectors[i]) for i in range(len(B_measured))])
norms_measured = np.linalg.norm(B_measured, axis=1)
norms_test = np.linalg.norm(test_vectors, axis=1)
alphas = np.degrees(np.arccos(dot / (norms_measured * norms_test)))

# generate histogram for all data
fig, ax = plt.subplots()
n, bins, patches = ax.hist(alphas, 50)
ax.set_xlabel('Angular Error, $\\alpha$ [°]')
ax.set_ylabel('Counts')
ax.text(0.7, 0.55, f'$\\mu$ = {np.mean(alphas):.2f}°,\n'+ 
                f'$\\sigma$ = {np.std(alphas):.2f}°,\n'+ 
                f'median = {np.median(alphas):.2f}°\n'+
                f'RMS $\\alpha$ = {estimate_RMS_error(alphas, np.zeros_like(alphas)):.2f}°\n'+
                f'min = {np.min(alphas):.2f}°,\n'+ 
                f'max = {np.max(alphas):.2f}°,\n', transform=ax.transAxes)

image_file_name = f'{os.path.splitext(data_filename)[0]}_performance_histogram.png'
image_path = os.path.join(data_directory, image_file_name)
fig.savefig(image_path, dpi=300)
plt.show()


# only consider data >= 30 mT
mask = norms_measured >= 30

fig, ax = plt.subplots()
n, bins, patches = ax.hist(alphas[mask], 50)
ax.set_xlabel('Angular Error for |B| $\\geq$ 30 mT, $\\alpha$ [°]')
ax.set_ylabel('Counts')
ax.text(0.7, 0.55, f'$\\mu$ = {np.mean(alphas[mask]):.2f}°,\n'+ 
                f'$\\sigma$ = {np.std(alphas[mask]):.2f}°,\n'+ 
                f'median = {np.median(alphas[mask]):.2f}°\n'+
                f'RMS $\\alpha$ = {estimate_RMS_error(alphas[mask], np.zeros_like(alphas[mask])):.2f}°\n'+
                f'min = {np.min(alphas[mask]):.2f}°,\n'+ 
                f'max = {np.max(alphas[mask]):.2f}°,\n', transform=ax.transAxes)

image_file_name = f'{os.path.splitext(data_filename)[0]}_performance_histogram_30-50mT.png'
image_path = os.path.join(data_directory, image_file_name)
fig.savefig(image_path, dpi=300)
plt.show()

#%%
print(np.linalg.norm([20,20,20]))
#%%
# plot angular error vs polar angle, azimuthal angle and field magnitude
fig, axs = plt.subplots(3, figsize = (6,8))

axs[0].plot(get_theta(B_measured), alphas, marker='.', linestyle='')
axs[0].set_xlabel('polar angle, $\\theta$ [°]')
axs[0].set_ylabel('Angular Error, $\\alpha$ [°]')

axs[1].plot(get_phi(B_measured), alphas, marker='.', linestyle='')
axs[1].set_xlabel('azimuthal angle, $\\phi$ [°]')
axs[1].set_ylabel('Angular Error, $\\alpha$ [°]')

axs[2].plot(np.linalg.norm(B_measured, axis=1), alphas, marker='.', linestyle='')
axs[2].set_xlabel('measured B field, $B_{meas}$ [mT]')
axs[2].set_ylabel('Angular Error, $\\alpha$ [°]')

[axs[i].set_ylim(bottom=0) for i in range(len(axs))]

plt.tight_layout()

image_file_name = f'{os.path.splitext(data_filename)[0]}_performance.png'
image_path = os.path.join(data_directory, image_file_name)
fig.savefig(image_path, dpi=300)
plt.show()



#  ---------------------------------------------------------------------------------------
# create plots for sweeps along main axes
#%%
# read in data 
data_directory = '../predictions_on_test_set/measure_predictions'
data_filename = '21_01_04_13-04-21_measure_predictions_degree3_B2I_sweep_z_axis_50mT_demagnet5A.csv'
data_filepath = os.path.join(data_directory, data_filename)

_, B_measured, B_measured_std, expected = extract_raw_data_from_file(data_filepath)

sweep_axis = 2
print('sweep along {} axis'.format(sweep_axis))

# define which components should be plotted - only plot combinations containing sweep_axis!
components = np.array([[sweep_axis, (sweep_axis+1) % 3], [sweep_axis, (sweep_axis+2) % 3]])

# generate a plot of both rotations
fig, axs = plt.subplots(ncols=2)
fig.set_size_inches(8, 3)

for i in range(2):
    # plot components orthogonal to rotation axis 
    # axs[i].plot(expected[:, components[i, 0]], expected[:, components[i, 1]], 
    #                         linestyle='', marker='.', label = 'linear model')

    axs[i].errorbar(B_measured[:, components[i, 0]], B_measured[:, components[i, 1]], 
                            xerr = B_measured_std[:, components[i, 0]], 
                            yerr = B_measured_std[:, components[i, 1]],
                            linestyle='', marker='.', capsize = 2, 
                            label = 'cubic model')

    # set axis labels
    labels_components = ['$B_x$ [mT]', '$B_y$ [mT]', '$B_z$ [mT]']
    axs[i].set_xlabel(labels_components[components[i, 0]])
    axs[i].set_ylabel(labels_components[components[i, 1]])

    # set aspect ratio to one, such that a circle actually looks round 
    # axs[i].legend()

# equalize the ranges
xmin = np.min([axs[i].get_xlim()[0] for i in range(len(axs))])
xmax = np.max([axs[i].get_xlim()[1] for i in range(len(axs))])
ymin = np.min([axs[i].get_ylim()[0] for i in range(len(axs))])
ymax = np.max([axs[i].get_ylim()[1] for i in range(len(axs))])
for i in range(len(axs)):
    axs[i].set_xlim(xmin, xmax)
    axs[i].set_ylim(ymin, ymax)

plt.tight_layout()

sweep_ax_label = ['x', 'y', 'z'][sweep_axis]
image_file_name = f'along_mainAxis_{sweep_ax_label}.png'
image_path = os.path.join(data_directory, image_file_name)
fig.savefig(image_path, dpi=300)

plt.show()

# %%
# evaluate performance along main axes
num = 100
max_field = 50
test_vectors = np.zeros((num, 3))
test_vectors[:, sweep_axis] = np.linspace(-max_field, max_field, num)

print('\nperformance of currently implemented model')
evaluate_performance(B_measured, test_vectors)


# estimate angular errors
dot = np.array([np.dot(B_measured[i], test_vectors[i]) for i in range(len(B_measured))])
norms_measured = np.linalg.norm(B_measured, axis=1)
norms_test = np.linalg.norm(test_vectors, axis=1)
alphas = np.degrees(np.arccos(dot / (norms_measured * norms_test)))

# plot angular errors versus desired field strength along main axis 
fig, ax = plt.subplots()
ax.plot(test_vectors[:, sweep_axis], alphas, linestyle='', marker='.')


labels_components = ['$B_x$ [mT]', '$B_y$ [mT]', '$B_z$ [mT]']
ax.set_xlabel(f'desired field {labels_components[sweep_axis]}')
ax.set_ylabel('Angular Error, $\\alpha$ [°]')
ax.set_ylim(bottom=0, top=17)

# add text containing statistics
plt.text(0.05, 0.6, '$\\langle \\alpha \\rangle$ = {:.2f}°,\n'.format(np.mean(alphas))+ 
                '$\\sigma (\\alpha)$ = {:.2f}°,\n'.format(np.std(alphas))+ 
                'min ($\\alpha$) = {:.2f}°,\n'.format(np.min(alphas))+ 
                'max ($\\alpha$) = {:.2f}°,\n'.format(np.max(alphas))+ 
                'median ($\\alpha$) = {:.2f}°,\n'.format(np.median(alphas)), transform=ax.transAxes)

plt.tight_layout()

sweep_ax_label = ['x', 'y', 'z'][sweep_axis]
image_file_name = f'along_mainAxis_{sweep_ax_label}_angular_error.png'
image_path = os.path.join(data_directory, image_file_name)
fig.savefig(image_path, dpi=300)

plt.show()




# %%

import pickle

filename = './fitting_parameters/new_model.sav'
test_vectors = np.array([[0,0,50]])

# load the model from disk
[loaded_model, loaded_poly] = pickle.load(open(filename, 'rb'))

# preprocess test vectors, st. they have correct shape for model
test_vectors_ = loaded_poly.fit_transform(test_vectors) 

# estimate prediction
predictions_new_sweep = loaded_model.predict(test_vectors_)
