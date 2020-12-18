""" 
filename: transfomations.py

The following helper functions provide the calculation of magnetic fields and associated currents

Author: Moritz Reinders
        moritzr1@gmail.com

Edited by: Maxwell Guerne-Kieferndorf (QZabre)
           gmaxwell@student.ethz.ch

Date: 09.10.2020
"""

########## Standard library imports ##########
import math
import numpy as np
import pandas as pd

# Lookup table for current/magnetic field values
LookupTable = {}
filepathLUT = r'data_sets\linearization_matrices\20_11_16_180159-86LUT.csv'

def buildLUT(path=filepathLUT):
    global LookupTable
    
    dataLUT = pd.read_csv(path).to_numpy()
    directionsB = dataLUT[:,0:3]
    currentConfigs = dataLUT[:,3:6]
    
    for i in range(len(directionsB)):
        unit = np.around(directionsB[i]/20, 3)
        unitCurrents = np.around(currentConfigs[i]/20, 3)
        
        LookupTable[f'direction {i+1}'] = (unit.tolist(), unitCurrents.tolist())
    

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
    global LookupTable
    buildLUT()

    index = None
    
    magnitude = np.linalg.norm(B_fieldVector)
    unitVector = np.around(B_fieldVector / magnitude, 3)
    print(f'unit vector: ({unitVector[0]:.3f},{unitVector[1]:.3f},{unitVector[2]:.3f})')
    
    for key in LookupTable:
        if LookupTable[key][0] == unitVector.tolist():
            index = key

    if index is not None:
        currVector = np.rint(magnitude * np.array(LookupTable[index][1]))
        
    else:
        actMatrix = np.zeros((3, 3))
        # actMatrix[0, 0] = 0.051343
        # actMatrix[0, 1] = 0
        # actMatrix[0, 2] = -0.051343
        # actMatrix[1, 0] = -0.029643
        # actMatrix[1, 1] = 0.059286
        # actMatrix[1, 2] = -0.029643
        # actMatrix[2, 0] = 0.008820
        # actMatrix[2, 1] = 0.008820
        # actMatrix[2, 2] = -0.008820
        
        # actMatrix = actMatrix * windings

        # These values are still just estimates in the current config. But based on measurements. (06.11.) ~1.5mm above poles
        # dB_x/dI_{1,2,3}
        actMatrix[0, 0] = 44
        actMatrix[0, 1] = -14
        actMatrix[0, 2] = -34
        # dB_y/dI_{1,2,3}
        actMatrix[1, 0] = -35
        actMatrix[1, 1] = 50
        actMatrix[1, 2] = -13
        # dB_z/dI_{1,2,3}
        actMatrix[2, 0] = 29
        actMatrix[2, 1] = -5
        actMatrix[2, 2] = -13 # z-field is positive when I3 is negative and increases when I3 decreases

        
        actMatrix_inverse = np.linalg.inv(actMatrix)
        # print('Inverse actuation matrix: \n', act_matrix_inverse)
        currVector = actMatrix_inverse.dot(B_fieldVector)  # in amps

        currVector = currVector * 1000  # in milliamps

        currVector = np.rint(currVector)  # round to nearest milliamp

    return currVector


def computeMagField(currVector, windings=508):
    """
    Compute magnetic field vector generated by the currents. Actuation matrix as derived from simulations so far

    Args:
    - currVector: Vector of 3 current values, as a np.array, units: [mA]
    x windings: number of windings per coil (508 for the QZabre electromagnet) not needed anymore

    Returns: 
    Vector of 3 B field components (Bx,By,Bz), as a np.array, units: [mT]
    """
    
    actMatrix = np.zeros((3, 3))
    # actMatrix[0, 0] = 0.051343
    # actMatrix[0, 1] = 0
    # actMatrix[0, 2] = -0.051343
    # actMatrix[1, 0] = -0.029643
    # actMatrix[1, 1] = 0.059286
    # actMatrix[1, 2] = -0.029643
    # actMatrix[2, 0] = 0.008820
    # actMatrix[2, 1] = 0.008820
    # actMatrix[2, 2] = -0.008820
    
    # actMatrix = actMatrix * windings
    # dB_x/dI_{1,2,3}
    actMatrix[0, 0] = 44
    actMatrix[0, 1] = -14
    actMatrix[0, 2] = -34
    # dB_y/dI_{1,2,3}
    actMatrix[1, 0] = -35
    actMatrix[1, 1] = 50
    actMatrix[1, 2] = -13
    # dB_z/dI_{1,2,3}
    actMatrix[2, 0] = 29
    actMatrix[2, 1] = -5
    actMatrix[2, 2] = -13 # z-field is positive when I3 is negative and increases when I3 decreases

    currVector = currVector / 1000  # in amps

    B_fieldVector = actMatrix.dot(currVector)  # in millitesla

    return B_fieldVector


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
    buildLUT()
    # print('Lookup table:')
    # for key in LookupTable:
    #     print(f'{key} = {LookupTable[key]}')
    B1 = computeMagneticFieldVector(theta=45, phi=60, magnitude=50)
    print(f'{B1 = }')
    currents = computeCoilCurrents(B1)
    print(f'I1 = {currents[0]}, I2 = {currents[1]}, I3 = {currents[2]}')
    A = rotationMatrix(psi=66,theta=42,alpha=10.1)
    print(A)
