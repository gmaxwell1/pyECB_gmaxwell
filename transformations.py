""" 
filename: transfomations.py

The following helper functions provide the calculation of magnetic fields and associated currents,
using a cubic model for the relation between the current and magnetic field (and vice versa).

Author: Moritz Reinders
        moritzr1@gmail.com

Edited by: Maxwell Guerne-Kieferndorf (QZabre)
           gmaxwell@student.ethz.ch

Date: 09.10.2020
latest update: 06.01.2021
"""

########## Standard library imports ##########
import math
import numpy as np
import pandas as pd
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model 
import pickle
import joblib



# Lookup table for current/magnetic field values
# LookupTable = {}
# filepathLUT = r'data_sets\linearization_matrices\20_11_16_180159-86LUT.csv'
# def buildLUT(path=filepathLUT):
#     global LookupTable
    
#     dataLUT = pd.read_csv(path).to_numpy()
#     directionsB = dataLUT[:,0:3]
#     currentConfigs = dataLUT[:,3:6]
    
#     for i in range(len(directionsB)):
#         unit = np.around(directionsB[i]/20, 3)
#         unitCurrents = np.around(currentConfigs[i]/20, 3)
        
#         LookupTable[f'direction {i+1}'] = (unit.tolist(), unitCurrents.tolist())
        
def read_fitting_parameters(filepath):
    """
    Extract fitting paramters A from file.
    """
    # extract data and convert to ndarray
    A = pd.read_csv(filepath).to_numpy().T

    return A


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
    

def computeMagneticFieldVector(magnitude, theta, phi):
    """
    Compute the cartesian coordinates of a magnetic field with an arbitrary direction and an arbitrary magnitude, given the spherical coordinates.

    Args:
    - magnitude: of the B-field, units: [mT]
    - theta: polar angle, between desired field direction and z axis
    - phi: azimuthal angle (angle measured counter clockwise from the x axis)

    Returns:
    Vector of 3 B field components (Bx,By,Bz), as a np.array, units: [mT]
    """

    x = math.sin(math.radians(theta)) * math.cos(math.radians(phi))
    y = math.sin(math.radians(theta)) * math.sin(math.radians(phi))
    z = math.cos(math.radians(theta))

    unitVector = np.array((x, y, z))
    unitVector = unitVector / np.linalg.norm(unitVector)

    return np.around(unitVector * magnitude, 3)


def computeCoilCurrents(B_fieldVector, windings=508, resistance=0.5):
    """
    Compute coil currents (in mA) required to generate the desired magnetic field vector. Actuation matrix derived from simulations so far

    Args:
    - B_fieldVector: np.array containing the B-field, in cartesian coordinates, magnitude units: [mT]
    x windings: number of windings per coil (508 for the QZabre electromagnet) not needed anymore
    x resistance: resistance per coil (ohms) not needed anymore

    Returns: 
    Vector of 3 current values, as a np.array, units: [mA]
    """
       # deprecated linear model
        # actMatrix = np.zeros((3, 3))
        # These values are still just estimates in the current config. But based on measurements. (06.11.) ~1.5mm above poles
        # dB_x/dI_{1,2,3}
        # actMatrix[0, 0] = 44
        # actMatrix[0, 1] = -14
        # actMatrix[0, 2] = -34
        # dB_y/dI_{1,2,3}
        # actMatrix[1, 0] = -35
        # actMatrix[1, 1] = 50
        # actMatrix[1, 2] = -13
        # dB_z/dI_{1,2,3}
        # actMatrix[2, 0] = 29
        # actMatrix[2, 1] = -5
        # actMatrix[2, 2] = -13
        # actMatrix_inverse = np.linalg.inv(actMatrix)
        # print('Inverse actuation matrix: \n', act_matrix_inverse)
        
    # read in paramters from previous fits, this is 'our' method for fitting the data points
    
    # filename_fit_params = 'measure_predictions_degree3_B2I'
    # filepath_fit_params = f'fitting_parameters\\{filename_fit_params}.csv'
    # A = read_fitting_parameters(filepath_fit_params)
    # B_fieldVector_reshape = B_fieldVector.reshape((1,3))
    # currVector = evaluate_fit(A, B_fieldVector_reshape)  # in amps
    filename = r'fitting_parameters\model_poly3_newSensorData_B2I.sav'

    # load the model from disk
    [loaded_model, loaded_poly] = pickle.load(open(filename, 'rb'))

    # preprocess test vectors, st. they have correct shape for model
    B_fieldVector_reshape = B_fieldVector.reshape((1,3))
    test_vectors_ = loaded_poly.fit_transform(B_fieldVector_reshape) 

    # estimate prediction
    currVector = loaded_model.predict(test_vectors_)
    
    currVector = currVector.reshape(3) * 1000  # in milliamps
    currVector = np.rint(currVector)  # round to nearest milliamp
    
    return currVector


def computeMagField(currVector, windings=508):
    """
    Compute magnetic field vector generated by the currents. Uses cubic regression model.

    Args:
    - currVector: Vector of 3 current values, as a np.array, units: [mA]
    x windings: number of windings per coil (508 for the QZabre electromagnet) not needed anymore

    Returns: 
    Vector of 3 B field components (Bx,By,Bz), as a np.array, units: [mT]
    """
    
    # actMatrix = np.zeros((3, 3))   
    # actMatrix = actMatrix * windings
    # dB_x/dI_{1,2,3}
    # actMatrix[0, 0] = 44
    # actMatrix[0, 1] = -14
    # actMatrix[0, 2] = -34
    # dB_y/dI_{1,2,3}
    # actMatrix[1, 0] = -35
    # actMatrix[1, 1] = 50
    # actMatrix[1, 2] = -13
    # dB_z/dI_{1,2,3}
    # actMatrix[2, 0] = 29
    # actMatrix[2, 1] = -5
    # actMatrix[2, 2] = -13 # z-field is positive when I3 is negative and increases when I3 decreases

    filename = r'fitting_parameters\model_poly3_newSensorData_I2B.sav'

    # load the model from disk
    [loaded_model, loaded_poly] = pickle.load(open(filename, 'rb'))

    # preprocess test vectors, st. they have correct shape for model
    currVector = currVector / 1000  # in amps

    currVector_reshape = currVector.reshape((1,3))
    test_vectors_ = loaded_poly.fit_transform(currVector_reshape) 

    # estimate prediction
    B_fieldVector = loaded_model.predict(test_vectors_)
    
    # currVector = np.rint(currVector)  # round to nearest milliamp
    # B_fieldVector = actMatrix.dot(currVector)  # in millitesla

    return B_fieldVector.reshape(3)


def rotationMatrix(inVector=np.array([1,0,0]), theta=90, psi=0, alpha=10.1):
    """
    Rotate a vector around any axis (defined by angles psi and theta), rotation around axis given by alpha

    Args:
    - inVector: arbitrary cartesian R^3 vector
    - psi: azimuthal angle, in degrees
    - theta: polar angle (measured downward from z axis), in degrees
    - alpha: counterclockwise rotation amount around the specified axis, in degrees

    Returns: 
    Vector of 3 B field components (Bx,By,Bz), as a np.array, units: [mT]
    """

    a_x = math.sin(math.radians(theta)) * math.cos(math.radians(psi))
    a_y = math.sin(math.radians(theta)) * math.sin(math.radians(psi))
    a_z = math.cos(math.radians(theta))

    axis = np.array([a_x, a_y, a_z])
    #a_unit = a_unit / np.linalg.norm(a_unit)

    #print('Axis of rotation (Cartesian) : \n ', a_unit)

    rot_matrix = np.zeros((3, 3))
    rot_matrix[0, 0] = math.cos(math.radians(
        alpha)) + ((axis[0] ** 2) * (1 - math.cos(math.radians(alpha))))
    rot_matrix[1, 0] = axis[1] * axis[0] * \
        (1 - math.cos(math.radians(alpha))) + \
        a_z * math.sin(math.radians(alpha))
    rot_matrix[2, 0] = axis[2] * axis[0] * (1 - math.cos(math.radians(alpha))) - axis[1] * math.sin(
        math.radians(alpha))

    rot_matrix[0, 1] = axis[0] * axis[1] * \
        (1 - math.cos(math.radians(alpha))) - \
        a_z * math.sin(math.radians(alpha))
    rot_matrix[1, 1] = math.cos(math.radians(
        alpha)) + ((axis[1] ** 2) * (1 - math.cos(math.radians(alpha))))
    rot_matrix[2, 1] = axis[2] * axis[1] * (1 - math.cos(math.radians(alpha))) + axis[0] * math.sin(
        math.radians(alpha))

    rot_matrix[0, 2] = axis[0] * axis[2] * (1 - math.cos(math.radians(alpha))) + axis[1] * math.sin(
        math.radians(alpha))
    rot_matrix[1, 2] = axis[1] * axis[2] * (1 - math.cos(math.radians(alpha))) - axis[0] * math.sin(
        math.radians(alpha))
    rot_matrix[2, 2] = math.cos(math.radians(
        alpha)) + ((axis[2] ** 2) * (1 - math.cos(math.radians(alpha))))

    #print('Rotation Matrix:')
    #print('\t{}\t{}\t{}'.format(rot_matrix[0,0], rot_matrix[1,0], rot_matrix[2,0]))
    #print('\t{}\t{}\t{}'.format(rot_matrix[0,1], rot_matrix[1,1], rot_matrix[2,1]))
    #print('\t{}\t{}\t{}'.format(rot_matrix[0,2], rot_matrix[1,2], rot_matrix[2,2]))
    # self.magnetic_field_unit = rot_matrix.dot(self.magnetic_field_unit)
    return np.around(rot_matrix.dot(inVector),3)


if __name__ == '__main__':
    # ------------------------Testing area--------------------------
    #
    # test the transformation functions
    # buildLUT()
    # print('Lookup table:')
    # for key in LookupTable:
    #     print(f'{key} = {LookupTable[key]}')
    B1 = computeMagneticFieldVector(theta=45, phi=60, magnitude=50)
    print(f'{B1 = }')
    currents = computeCoilCurrents(B1)
    print(f'I1 = {currents[0]}, I2 = {currents[1]}, I3 = {currents[2]}')
    B2 = computeMagField(currents)
    print(f'{B2 = }')

