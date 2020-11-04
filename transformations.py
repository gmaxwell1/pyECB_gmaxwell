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
    #print('Requested magnetic field direction: ({},{},{})'.format(unitVector[0],unitVector[1],unitVector[2]))
    #print('Requested magnetic flux density: {} mT'.format(magnitude))

    return unitVector * magnitude


def computeCoilCurrents(B_fieldVector, windings, resistance):
    """
    Compute coil currents (in mA) required to generate the desired magnetic field vector. Actuation matrix derived from simulations so far

    Args:
    - B_fieldVector: np.array containing the B-field, in cartesian coordinates, magnitude units: [mT]
    - windings: number of windings per coil (508 for the QZabre electromagnet)
    - resistance: resistance per coil (ohms)

    Returns: 
    Vector of 3 current values, as a np.array, units: [mA]
    """
    actMatrix = np.zeros((3, 3))
    actMatrix[0, 0] = 0.051343
    actMatrix[0, 1] = 0
    actMatrix[0, 2] = -0.051343
    actMatrix[1, 0] = -0.029643
    actMatrix[1, 1] = 0.059286
    actMatrix[1, 2] = -0.029643
    actMatrix[2, 0] = 0.008820
    actMatrix[2, 1] = 0.008820
    actMatrix[2, 2] = -0.008820
    
    # actMatrix = actMatrix * windings

    # These values are still just estimates in the current config. But based on measurements. (29.10.) ~1-3mm above poles
    # actMatrix[0, 0] = 35
    # actMatrix[0, 1] = -8
    # actMatrix[0, 2] = -26
    # actMatrix[1, 0] = -29
    # actMatrix[1, 1] = 40
    # actMatrix[1, 2] = -10
    # actMatrix[2, 0] = 22
    # actMatrix[2, 1] = 0
    # actMatrix[2, 2] = -8 # z-field is positive when I3 is negative and increases when I3 decreases
    
    actMatrix_inverse = np.linalg.inv(actMatrix)

    # print('Inverse actuation matrix: \n', act_matrix_inverse)

    currVector_amp_turns = actMatrix_inverse.dot(B_fieldVector)  # in amp-turns
    currVector = (currVector_amp_turns / windings) # in amps
    # currVector = actMatrix_inverse.dot(B_fieldVector)  # in amps

    power = ((currVector[0] ** 2) + (currVector[1] **
                                     2) + (currVector[2] ** 2)) * resistance
    #print('Power [W]: ', power)

    currVector = currVector * 1000  # in milliamps
    currVector = np.rint(currVector)  # round to nearest milliamp

    # print('Current vector [mA]: \n', currVector)

    return currVector


def computeMagField(currVector, windings):
    """
    Compute magnetic field vector generated by the currents. Actuation matrix as derived from simulations so far

    Args:
    - currVector: Vector of 3 current values, as a np.array, units: [mA]
    - windings: number of windings per coil (508 for the QZabre electromagnet)

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

    actMatrix[0, 0] = 35
    actMatrix[0, 1] = -8
    actMatrix[0, 2] = -26
    actMatrix[1, 0] = -29
    actMatrix[1, 1] = 40
    actMatrix[1, 2] = -10
    actMatrix[2, 0] = 22
    actMatrix[2, 1] = 0
    actMatrix[2, 2] = -8 # z-field is positive when I3 is negative and increases when I3 decreases

    # print('Inverse actuation matrix: \n', act_matrix_inverse)
    currVector = currVector / 1000  # in amps

    # currVector_amp_turns = (currVector * windings)  # in amp-turns
    B_fieldVector = actMatrix.dot(currVector)  # in millitesla

    # print('Current vector [mA]: \n', currVector)

    return B_fieldVector


def rotationMatrix(inVector, psi, theta, alpha):
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

    a_unit = np.array(([a_x], [a_y], [a_z]))
    #a_unit = a_unit / np.linalg.norm(a_unit)

    #print('Axis of rotation (Cartesian) : \n ', a_unit)

    rot_matrix = np.zeros((3, 3))
    rot_matrix[0, 0] = math.cos(math.radians(
        alpha)) + ((a_unit[0] ** 2) * (1 - math.cos(math.radians(alpha))))
    rot_matrix[1, 0] = a_unit[1] * a_unit[0] * \
        (1 - math.cos(math.radians(alpha))) + \
        a_z * math.sin(math.radians(alpha))
    rot_matrix[2, 0] = a_unit[2] * a_unit[0] * (1 - math.cos(math.radians(alpha))) - a_unit[1] * math.sin(
        math.radians(alpha))

    rot_matrix[0, 1] = a_unit[0] * a_unit[1] * \
        (1 - math.cos(math.radians(alpha))) - \
        a_z * math.sin(math.radians(alpha))
    rot_matrix[1, 1] = math.cos(math.radians(
        alpha)) + ((a_unit[1] ** 2) * (1 - math.cos(math.radians(alpha))))
    rot_matrix[2, 1] = a_unit[2] * a_unit[1] * (1 - math.cos(math.radians(alpha))) + a_unit[0] * math.sin(
        math.radians(alpha))

    rot_matrix[0, 2] = a_unit[0] * a_unit[2] * (1 - math.cos(math.radians(alpha))) + a_unit[1] * math.sin(
        math.radians(alpha))
    rot_matrix[1, 2] = a_unit[1] * a_unit[2] * (1 - math.cos(math.radians(alpha))) - a_unit[0] * math.sin(
        math.radians(alpha))
    rot_matrix[2, 2] = math.cos(math.radians(
        alpha)) + ((a_unit[2] ** 2) * (1 - math.cos(math.radians(alpha))))

    #print('Rotation Matrix:')
    #print('\t{}\t{}\t{}'.format(rot_matrix[0,0], rot_matrix[1,0], rot_matrix[2,0]))
    #print('\t{}\t{}\t{}'.format(rot_matrix[0,1], rot_matrix[1,1], rot_matrix[2,1]))
    #print('\t{}\t{}\t{}'.format(rot_matrix[0,2], rot_matrix[1,2], rot_matrix[2,2]))
    # self.magnetic_field_unit = rot_matrix.dot(self.magnetic_field_unit)
    return rot_matrix.dot(inVector)


if __name__ == '__main__':
    # ------------------------Testing area--------------------------
    #
    # test the transformation functions

    B1 = computeMagneticFieldVector(theta=0, phi=0, magnitude=69)
    currents1 = computeCoilCurrents(B1, 508, 0.47)

    _B1 = computeMagneticFieldVector(theta=90, phi=0, magnitude=20)
    _currents1 = computeCoilCurrents(_B1, 508, 0.47)

    print('B = ({},{},{})^T corresponds to the currents I1 = {}, I2 = {}, I3 = {}'
          .format(B1[0], B1[1], B1[2], currents1[0], currents1[1], currents1[2]))
    print('B = ({},{},{})^T corresponds to the currents I1 = {}, I2 = {}, I3 = {}'
          .format(_B1[0], _B1[1], _B1[2], _currents1[0], _currents1[1], _currents1[2]))
