########## Standard library imports ##########
import numpy as np
import math
from time import time,sleep
from datetime import datetime
from os import getcwd, path
from pathlib import Path
import csv


##########  test method to write a .csv file ##########

def textCSV(str):

    # format the output file/folder
    now = datetime.now().strftime('%H%M%S')
    today = datetime.now().strftime('%Y-%m-%d')
    filename = '{}_{}_{}'.format(today, now, 'meas.csv')
    cwd = getcwd()
    workDir = path.join(cwd, today)
    # create new directory (for current day) if necessary
    Path(workDir).mkdir(parents=True,exist_ok=True)

    with open(path.join(workDir, filename), 'x', newline='') as newFile:
        _writer = csv.writer(newFile)
        _writer.writerow(['Hello, ' + str + '!'])

        elapsedTime = 0
        startTime = time()

        while elapsedTime < 1:
            elapsedTime = time() - startTime

            #loop_time = time() - startTime - elapsedTime

            _writer.writerow([elapsedTime])
            # not definitive, just to see how significant of an effect the execution
            # time has on the sampling period, especially for high frequencies
            sleep(0.1)


##########  TODO: read csv file and save it as a numpy array ##########

def readCSV(filename):
    pass


########## test stuff out ##########

if __name__ == '__main__':
    textCSV('Maxwell')