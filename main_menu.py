"""
filename: main_menu.py

This code is meant to be used as an interface with the ECB 820. The user can choose from various
methods to set currents on the different channels (coils in Pantec's terminology) and thus communicate
with the ECB. The standard interface is the command line, but another option is to integrate this into 
a GUI for the best user experience.

Author: Maxwell Guerne-Kieferndorf (QZabre)
        gmaxwell@student.ethz.ch

Date: 28.09.2020
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


##########  Current parameters ##########

desCurrents = [0, 0, 0, 0, 0, 0, 0, 0]  # in milliamps
currDirectParam = b'1'

##########  Vector magnet properties ##########

windings = 508  # windings per coil
resistance = 0.45  # resistance per coil


##########  Main menu for ECB/Vector magnet operation ##########
#           an arbitrary magnitude     
#           Args:
#           -magnitude: of the B-field, in [mT]
#           -theta: angle between desired field direction and z axis
#           -phi: azimuthal angle (measured from the x axis)
#           -t: time for which the magnetic field should be activated (if not 0)
#           -direct: current direct parameter (can usually be left alone)

def MainMenu(initialized):

        # is there a connection?
        if initialized:
                c1 = '0'
                while c1 != 'x':
                        print('--------- Main Menu ---------')
                        print('[x] to exit\n[0]: set currents manually on 3 channels (in mA)')
                        print('[1]: generate magnetic fields (specify polar and azimuthal angles, magnitude)')
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
                        elif c1 == 'c':
                                getCurrents()
                        elif c1 == 's':
                                print(getStatus())
                        elif c1 == 'r':
                                print(np.random.randint(1,7))
                
                demagnetizeCoils()

        else:
                print('not connected')
                return


##########  generate a magnetic field in an arbitrary direction and ##########
#           an arbitrary magnitude     
#           Args:
#           -magnitude: of the B-field, in [mT]
#           -theta: angle between desired field direction and z axis
#           -phi: azimuthal angle (measured from the x axis)
#           -t: time for which the magnetic field should be activated (if not 0)
#           -direct: current direct parameter (can usually be left alone)

def generateMagneticField(magnitude, theta, phi, t=0, direct=b'1'):

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
                        elif c1 == 'd':
                                N = input('No. of measurements: ')
                                f = input('No. of measurements/second: ')
                                pollCurrents(int(N),float(f))

                disableCurrents()
        else:
                return


##########  run arbitrary currents (less than maximum current) on each channel ##########
#           Args:
#           -coils: current values in [mA]
#           -t: time for which the magnetic field should be activated (if not 0)
#           -direct: current direct parameter (can usually be left alone)

def runCurrents(*coils, t=0, direct=b'1'):

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
                        elif c1 == 'd':
                                N = input('No. of measurements: ')
                                f = input('No. of measurements/second: ')
                                pollCurrents(int(N),float(f))

                disableCurrents()
        else:
                return


##########  run arbitrary currents (less than maximum current) on each channel ##########
#           Args:
#           -coils: current values in [mA]
#           -t: time for which the magnetic field should be activated (if not 0)
#           -direct: current direct parameter (can usually be left alone)

def demagnetizeCoils():

        currDirectParam = b'1'
        
        stepSize = 100
        amplitude = 1500
        dt = 0.5

        enableCurrents()
        print('demagnitizing coils')
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


##########  gets current values writes the values into a .csv file named according ##########
#           to today's date/time
#           Args:
#           -measurements: an integer specifying how many values are wanted
#           -frequency: how many times per second (approximately!!, not reliable) to get the 
#                       values

def pollCurrents(measurements, frequency=1):
        dt = 1/frequency
        k = 0

        # format the output file/folder
        now = datetime.now().strftime('%H%M%S')
        today = datetime.now().strftime('%Y-%m-%d')
        filename = '{}_{}_{}'.format(today, now, 'meas.csv')
        cwd = getcwd()
        workDir = path.join(cwd, today)
        # create new directory (for current day) if necessary
        Path(workDir).mkdir(parents=True,exist_ok=True)

        with open(path.join(workDir, filename), 'x', newline='') as newFile:
                meas_writer = csv.writer(newFile)
                meas_writer.writerow(['Measurement time', 'Ch1 [mA]', 'Ch2 [mA]', 'Ch3 [mA]',
                                        'Ch4 [mA]', 'Ch5 [mA]', 'Ch6 [mA]', 'Ch7 [mA]', 'Ch8 [mA]'])
                line = []

                startTime = time()

                while k < measurements:
                        elapsedTime = time() - startTime
                        k = k + 1
                        result = getCurrents()
                        if not isinstance(result,list):
                                print(result)
                                return
                        line = result

                        line.insert(0, elapsedTime)
                        meas_writer.writerow(line)

                        # loop_time = time() - startTime - elapsedTime
                        # print(loop_time)
                        # not definitive, just to see how significant of an effect the execution
                        # time has on the sampling period, especially for high frequencies
                        sleep(dt)


if __name__ == '__main__':
        ecbInit = openConnection()
        MainMenu(ecbInit)
        closeConnection()