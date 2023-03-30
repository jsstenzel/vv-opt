#import argparse
import os
import sys
import shutil
import csv
import fileinput
sys.path.insert(0, "..")

import spectrograph as spec
import matplotlib.pyplot as plt
import observe
from astropy.io import fits
import scipy.signal as ss
import spectrograph as spec
import observe

import numpy as np
import itertools
import multiprocessing as mp
import math
from SALib.sample import saltelli, fast_sampler
from SALib.analyze import sobol, fast
from copy import deepcopy
import scipy.optimize as optimization

os.environ['COATINGS_PATH']="C:/Users/June/Desktop/Repositories/LLAMAS/SIMULATIONS/llamas-etc/COATINGS/"
os.environ['TEST_DATA_PATH']="C:/Users/June/Desktop/Repositories/LLAMAS/TEST DATA/"
