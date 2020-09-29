""" 
filename: transfomations.py

The following helper functions provide the `low level` calculation of magnetic fields and associated currents

Author: Moritz Reinders
        moritzr1@gmail.com

Edited by: Maxwell Guerne-Kieferndorf (QZabre)
           gmaxwell@student.ethz.ch

Date: 29.09.2020
"""

import math
import numpy as np

# Compute magnetic field vector in cartesian coordinates from given angles in polar coordinates
def computeMagneticFieldVector(theta, phi, magnitude):

    x = math.sin(math.radians(phi)) * math.cos(math.radians(theta))
    y = math.sin(math.radians(phi)) * math.sin(math.radians(theta))
    z = math.cos(math.radians(phi))

    unitVector = np.array(([x], [y], [z]))
    unitVector = unitVector / np.linalg.norm(unitVector)
    print('Requested magnetic field direction: \n ', unitVector)
    print('Requested magnetic flux density [mT]:  \n', magnitude)

    return unitVector * magnitude
     
# Compute coil currents (in mA) required to generate the desired magnetic field vector
# actuation matrix from simulations so far
def computeCoilCurrents(B_fieldVector, windings):

    actMatrix = np.zeros((3, 3))
    actMatrix[0, 0] = 0.051343
    actMatrix[0, 1] = 0
    actMatrix[0, 2] = -0.051343
    actMatrix[1, 0] = -0.029643
    actMatrix[1, 1] = 0.059286
    actMatrix[1, 2] = -0.029643
    actMatrix[2, 0] = 0.008820
    actMatrix[2, 1] = 0.008820
    actMatrix[2, 2] = 0.008820

    actMatrix_inverse = np.linalg.inv(actMatrix)

    # print('Inverse actuation matrix: \n', act_matrix_inverse)

    currVector_amp_turns = actMatrix_inverse.dot(B_fieldVector)  # in amp-turns
    currVector = ((currVector_amp_turns / windings))  # in milliamps

    power = ((current_vector[0] ** 2) + (current_vector[1] ** 2) + (current_vector[2] ** 2)) * self.resistance
    print('Power [W]: ', power)

    currVector = currVector * 1000  # in milliamps
    currVector = np.rint(currVector)  # round to nearest milliamp

    print('Current vector [mA]: \n', currVector)

    return currVector

# Rotate a vector around an axis (axis defined by angles omega and psi), rotation around axis given by eta
def rotationMatrix(unitVector, omega, psi, eta):
    
    a_x = math.sin(math.radians(psi)) * math.cos(math.radians(omega))
    a_y = math.sin(math.radians(psi)) * math.sin(math.radians(omega))
    a_z = math.cos(math.radians(psi))

    a_unit = np.array(([a_x], [a_y], [a_z]))
    a_unit = a_unit / np.linalg.norm(a_unit)

    # print('Axis of rotation (Cartesian) : \n ', a_unit)

    rot_matrix = np.zeros((3, 3))
    rot_matrix[0, 0] = math.cos(math.radians(eta)) + ((a_unit[0] ** 2) * (1 - math.cos(math.radians(eta))))
    rot_matrix[1, 0] = a_unit[1] * a_unit[0] * (1 - math.cos(math.radians(eta))) + a_z * math.sin(math.radians(eta))
    rot_matrix[2, 0] = a_unit[2] * a_unit[0] * (1 - math.cos(math.radians(eta))) - a_unit[1] * math.sin(
        math.radians(eta))

    rot_matrix[0, 1] = a_unit[0] * a_unit[1] * (1 - math.cos(math.radians(eta))) - a_z * math.sin(math.radians(eta))
    rot_matrix[1, 1] = math.cos(math.radians(eta)) + ((a_unit[1] ** 2) * (1 - math.cos(math.radians(eta))))
    rot_matrix[2, 1] = a_unit[2] * a_unit[1] * (1 - math.cos(math.radians(eta))) + a_unit[0] * math.sin(
        math.radians(eta))

    rot_matrix[0, 2] = a_unit[0] * a_unit[2] * (1 - math.cos(math.radians(eta))) + a_unit[1] * math.sin(
        math.radians(eta))
    rot_matrix[1, 2] = a_unit[1] * a_unit[2] * (1 - math.cos(math.radians(eta))) - a_unit[0] * math.sin(
        math.radians(eta))
    rot_matrix[2, 2] = math.cos(math.radians(eta)) + ((a_unit[2] ** 2) * (1 - math.cos(math.radians(eta))))
    # self.magnetic_field_unit = rot_matrix.dot(self.magnetic_field_unit)
    B_fieldVector = rot_matrix.dot(B_fieldVector)