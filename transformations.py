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
def computeMagneticFieldVector(magnitude, theta, phi):

    x = math.sin(math.radians(phi)) * math.cos(math.radians(theta))
    y = math.sin(math.radians(phi)) * math.sin(math.radians(theta))
    z = math.cos(math.radians(phi))

    unitVector = np.array((x, y, z))
    unitVector = unitVector / np.linalg.norm(unitVector)
    #print('Requested magnetic field direction: ({},{},{})'.format(unitVector[0],unitVector[1],unitVector[2]))
    #print('Requested magnetic flux density: {} mT'.format(magnitude))

    return unitVector * magnitude
     
# Compute coil currents (in mA) required to generate the desired magnetic field vector
# actuation matrix from simulations so far
def computeCoilCurrents(B_fieldVector, windings, resistance):

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

    power = ((currVector[0] ** 2) + (currVector[1] ** 2) + (currVector[2] ** 2)) * resistance
    #print('Power [W]: ', power)

    currVector = currVector * 1000  # in milliamps
    currVector = np.rint(currVector)  # round to nearest milliamp

    # print('Current vector [mA]: \n', currVector)

    return currVector

# Rotate a vector around an axis (axis defined by angles omega and psi), rotation around axis given by eta
def rotationMatrix(inVector, omega, psi, eta):
    
    a_x = math.sin(math.radians(psi)) * math.cos(math.radians(omega))
    a_y = math.sin(math.radians(psi)) * math.sin(math.radians(omega))
    a_z = math.cos(math.radians(psi))

    a_unit = np.array(([a_x], [a_y], [a_z]))
    #a_unit = a_unit / np.linalg.norm(a_unit)

    #print('Axis of rotation (Cartesian) : \n ', a_unit)

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
    
    B1 = computeMagneticFieldVector(theta=50, phi=40, magnitude=10)
    currents1 = computeCoilCurrents(B1, 508, 0.416)
    
    print('B = ({},{},{})^T corresponds to the currents I1 = {}, I2 = {}, I3 = {}'
          .format(B1[0],B1[1],B1[2],currents1[0],currents1[1],currents1[2]))

