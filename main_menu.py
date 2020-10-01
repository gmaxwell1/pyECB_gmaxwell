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
from ECB import * # this is ok, since all the function names are pretty much unique and known.
import transformations as tr
from main_comm import generateMagField, getCurrents