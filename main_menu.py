"""
filename: main_menu.py

This script is meant to be used as an interface with the ECB 820. The user can choose from various
methods to set currents on the different channels ('coils' in Pantec's terminology) and thus communicate
with the ECB. The standard interface is the command line, but another option is to integrate this into 
a GUI for the best user experience.

Author: Maxwell Guerne-Kieferndorf (QZabre)
        gmaxwell@student.ethz.ch
Based on code by Moritz Reinders

Date: 07.10.2020
"""
########## Standard library imports ##########
import numpy as np
import math
from time import time,sleep
from datetime import datetime
import sys
from os import getcwd, path
from pathlib import Path
import csv

########## local imports ##########
import transformations as tr
from main_comm import *
from measurements import *
from modules.analysis_tools import generate_plots

##########  Current parameters ##########

desCurrents = [0, 0, 0, 0, 0, 0, 0, 0]  # in milliamps
currDirectParam = b'1'

##########  Vector magnet properties ##########

windings = 508  # windings per coil
resistance = 0.47  # resistance per coil


def MainMenu(initialized):
        """
        Main menu for ECB/Vector magnet operation an arbitrary magnitude   

        Args:
        -magnitude: of the B-field, in [mT]
        -theta: angle between desired field direction and z axis
        -phi: azimuthal angle (measured from the x axis)
        -t: time for which the magnetic field should be activated (if not 0)
        -direct: current direct parameter (can usually be left alone)
        """

        # is there a connection?
        if initialized == 0:
                c1 = '0'
                while c1 != 'x':
                        print('--------- Main Menu ---------')
                        print('[x] to exit \n[0]: set currents manually on 3 channels (in mA)')
                        print('[1]: generate magnetic field (specify polar and azimuthal angles, magnitude)')
                        print('[2]: ramp magnetic field (specify polar and azimuthal angles, magnitude range)')
                        print('[3]: sweep multiple current values and make measurment with cube')
                        print('[c]: get currents \n[s]: get ECB status\n[r] roll a die')

                        c1 = input()
                        if c1 == '0':
                                inp1 = input('Channel 1: ')
                                inp2 = input('Channel 2: ')
                                inp3 = input('Channel 3: ')
                                inp4 = input('timer (leave empty -> manual termination) = ')
                                try:
                                        coil1 = int(inp1)
                                except:
                                        print('expected numerical value')
                                        coil1 = 0
                                try:
                                        coil2 = int(inp2)
                                except:
                                        print('expected numerical value')
                                        coil2 = 0
                                try:
                                        coil3 = int(inp3)
                                except:
                                        print('expected numerical value')
                                        coil3 = 0

                                if inp4 == '':
                                        runCurrents(coil1, coil2, coil3, t=0, direct=b'1')
                                else:
                                        try:
                                                timer = int(inp4)
                                        except:
                                                print('expected numerical value')
                                                timer = 0
                                        runCurrents(coil1, coil2, coil3, timer, direct=b'1')

                        elif c1 == '1':
                                inp1 = input('Magnitude in mT = ')
                                inp2 = input('Angle to z axis in deg = ')
                                inp3 = input('Angle to x axis in deg = ')
                                inp4 = input('timer (leave empty -> manual termination) = ')
                                try:
                                        mag = int(inp1)
                                except:
                                        print('expected numerical value')
                                        mag = 0
                                try:
                                        theta = int(inp2)
                                except:
                                        print('expected numerical value')
                                        theta = 0
                                try:
                                        phi = int(inp3)
                                except:
                                        print('expected numerical value')
                                        phi = 0
                                if inp4 == '':
                                        generateMagneticField(mag, theta, phi, t=0, direct=b'1')
                                else:
                                        try:
                                                timer = int(inp4)
                                        except:
                                                print('expected numerical value')
                                                timer = 0
                                        generateMagneticField(mag, theta, phi, timer, direct=b'1')

                        elif c1 == '2':
                                inp1 = input('Angle to z axis in deg = ')
                                inp2 = input('Angle to x axis in deg = ')
                                inp3 = input('starting magnitude in mT = ')
                                inp4 = input('ending magnitude in mT = ')
                                inp5 = input('duration of ramp (at least 10s is recommended): ')
                                inp6 = input('# of steps: ')

                                try:
                                        theta = int(inp1)
                                except:
                                        print('expected numerical value')
                                        theta = 0
                                try:
                                        phi = int(inp2)
                                except:
                                        print('expected numerical value')
                                        phi = 0
                                try:
                                        start_mag = int(inp3)
                                except:
                                        print('expected numerical value')
                                        start_mag = 0
                                try:
                                        end_mag = int(inp4)
                                except:
                                        print('expected numerical value')
                                        end_mag = 0
                                try:
                                        duration = int(inp5)
                                except:
                                        print('expected numerical value')
                                        duration = 0
                                try:
                                        steps = int(inp6)
                                except:
                                        print('expected numerical value')
                                        steps = 0
                                        
                                rampVectorField(theta, phi, start_mag, end_mag, duration, steps)

                        elif c1 == '3':
                                inp1 = input('configuration (z or x-y direction, acceptable inputs: z, xy) = ')
                                inp2 = input('starting current in mA = ')
                                inp3 = input('ending current in mA = ')
                                inp4 = input('# of steps: ')

                                try:
                                        config = inp1
                                except:
                                        print('expected valid input (z or xy), defaulting to z')
                                        config = 'z'
                                try:
                                        start_val = int(inp2)
                                except:
                                        print('expected numerical value')
                                        start_val = 0
                                try:
                                        end_val = int(inp3)
                                except:
                                        print('expected numerical value, defaulting to 1')
                                        end_val = 1
                                try:
                                        steps = int(inp4)
                                except:
                                        print('expected numerical value, defaulting to 1')
                                        duration = 1

                                        
                                sweepCurrents(config, start_val, end_val, steps)

                        elif c1 == 'c':
                                getCurrents()
                        elif c1 == 's':
                                print(getStatus())
                        elif c1 == 'r':
                                print(np.random.randint(1,7))

                c1 = input('Demagnetize coils? [y/n]\t')
                if c1 == 'y':
                        demagnetizeCoils()

        else:
                print('not connected')
                return


def generateMagneticField(magnitude, theta, phi, t=0, direct=b'1'):
        """  
        generate a magnetic field in an arbitrary direction and an arbitrary magnitude

        Args:
        -magnitude: of the B-field, in [mT]
        -theta: angle between desired field direction and z axis
        -phi: azimuthal angle (measured from the x axis)
        -t: time for which the magnetic field should be activated (if not 0)
        -direct: current direct parameter (can usually be left alone)
        """

        B_vector = tr.computeMagneticFieldVector(magnitude, theta, phi)
        I_vector = tr.computeCoilCurrents(B_vector, windings, resistance)

        # make sure that the current is not too high
        if np.amax(I_vector) > ECB_MAX_CURR:
                print("desired current exceeds limit")
                return

        currDirectParam = direct
        # copy the computed current values (mA) into the desCurrents list (first 3 positions)
        # cast to int
        for i in range(len(I_vector)):
                desCurrents[i] = int(I_vector[i])

        # user specified on time
        if t > 0:
                enableCurrents()
                setCurrents(desCurrents, currDirectParam)
                print('Current on coils 1, 2 and 3: [{}, {}, {}]'.format(I_vector[0], I_vector[1], I_vector[2]))
                sleep(t)

                disableCurrents()
        # on until interrupted by user
        elif t == 0:
                enableCurrents()

                setCurrents(desCurrents, currDirectParam)

                # wait until user presses enter
                c1 = '0'
                while c1 != 'q':
                        c1 = input('[q] to disable currents\n[c]: get currents\n[s]: get ECB status\n[t]: get coil temperature\n')
                        if c1 == 'c':
                                getCurrents()
                        elif c1 == 's':
                                print(getStatus())
                        elif c1 == 't':
                                getTemps()


                disableCurrents()
        else:
                return


def runCurrents(*coils, t=0, direct=b'1'):
        """
        run arbitrary currents (less than maximum current) on each channel 

        Args:
        -coils: current values in [mA]
        -t: time for which the magnetic field should be activated (if not 0)
        -direct: current direct parameter (can usually be left alone)
        """
        currDirectParam = direct
        # copy the computed current values (mA) into the desCurrents list (first 3 positions)
        # cast to int
        for i in range(min(len(coils),8)):
                if np.abs(coils[i]) > ECB_MAX_CURR:
                        print("desired current exceeds limit")
                        return
                desCurrents[i] = int(coils[i])

        # user specified time
        if t > 0:
                enableCurrents()
                setCurrents(desCurrents, currDirectParam)

                sleep(t)

                disableCurrents()

        # on until interrupted by user
        elif t == 0:
                enableCurrents()

                setCurrents(desCurrents, currDirectParam)

                # wait until user presses enter
                c1 = '0'
                while c1 != 'q':
                        c1 = input('[q] to disable currents\n[c]: get currents\n[d]: write currents to file\n[s]: get ECB status\n[t]: get coil temperature\n')
                        if c1 == 'c':
                                getCurrents()
                        elif c1 == 's':
                                print(getStatus())
                        elif c1 == 't':
                                getTemps()
#                        elif c1 == 'd':
#                                N = input('No. of measurements: ')
#                                f = input('No. of measurements/second: ')
#                                pollCurrents(int(N),float(f))

                disableCurrents()
        else:
                return


def rampVectorField(theta, phi, start_mag, finish_mag, duration, steps):
        """  
        ramps magnetic field from start_magn to finish_magn in a specified number of steps and over a specified duration

        Args:
        -theta & phi: give the direction of the magnetic field (Polar/azimuthal angles)
        -start/finish_magn: field range
        -duration: time over which to ramp the field 
        -steps: number of steps
        """
        delta_b = (finish_mag - start_mag) / (steps - 1)
        print(delta_b)
        delta_t = duration / steps
        print(delta_t)
        magnitude = start_mag

        for i in range(steps):
                print('Step: ', i+1)
                print('Magnitude: {} mT'.format(magnitude))

                generateMagneticField(magnitude, theta, phi, delta_t, b'1')
                magnitude = magnitude + delta_b


def sweepCurrents(config='z', start_val=0, end_val=1, steps=5):
        """
        sweeps all currents in the (1,1,1) or (1,0,-1) configuration, meaning we have either a z field or an x-y-plane field, measures magnetic field 
        and stores the measured values in various arrays
        
        Args:
        -config: 'z' or 'xy'
        -start/end_val: current to start and end with, mA
        -steps: number of steps
        """

        # initialization of all arrays
        all_curr_steps = np.linspace(start_val, end_val, steps)
        mean_values = np.zeros((steps, 3))
        stdd_values = np.zeros((steps, 3))
        expected_fields = np.zeros((steps, 3))
        directory = ''

        current_direction = np.ndarray(3)
        if config == 'z':
                current_direction[0] = 1
                current_direction[1] = 1
                current_direction[2] = 1
        elif config == 'xy':
                current_direction[0] = 1
                current_direction[1] = 0
                current_direction[2] = -1
        else:
                print('invalid input!')
                return

        enableCurrents()
        #iterate through all possible steps
        for i in range(steps):
                for k in range(3):
                        desCurrents[k] = current_direction[k] * all_curr_steps[i]
                B_expected = tr.computeMagField(current_direction * all_curr_steps[i], windings)
                setCurrents(desCurrents, currDirectParam)
                sleep(0.8)
                # collect measured and expected magnetic field (of the specified sensor in measurements)
                print('measurement nr. ', i)
                mean_data, std_data, directory = measure() # see measurements.py for more details
                mean_values[i] = mean_data
                stdd_values[i] = std_data
                expected_fields[i] = B_expected

                sleep(0.2)
        
        disableCurrents()

        # plotting section
        saveDataPoints((all_curr_steps / 1000), mean_values, stdd_values, expected_fields, directory)
        generate_plots((all_curr_steps / 1000), mean_values, stdd_values, expected_fields, flag_xaxis='I', flags_yaxis='pma', directory=directory, height_per_plot=3)


def demagnetizeCoils(stepSize = 100, amplitude = 1500, dt = 0.5, direct = b'1'):
        """
        Try to minimize residual flux in cores by running square wave shaped (positive-negative) currents on each coil.
        
        Args:
        -stepSize: amount by which to reduce amplitude every period
        -amplitude: Amplitude of square wave
        -dt: half period of pulse
        -direct: see ECB API documentation
        """
        currDirectParam = direct

        enableCurrents()
        print('demagnetizing coils')
        while amplitude > 0:
                desCurrents = [amplitude] * 8
                setCurrents(desCurrents, currDirectParam)
                sleep(dt)
                amplitude = (-1) * amplitude
                desCurrents = [amplitude] * 8
                setCurrents(desCurrents, currDirectParam)
                sleep(dt)
                amplitude = (-1) * (amplitude + stepSize)
        
        disableCurrents()


if __name__ == '__main__':
        ecbInit = openConnection()
        while ecbInit != 0:
                ecbInit = openConnection()

        MainMenu(ecbInit)
        closeConnection()

