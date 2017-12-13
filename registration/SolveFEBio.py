# -*- coding: utf-8 -*-
import numpy as np                                                                                   
from enthought.mayavi.mlab import *
from fileread import *
from pmeshdeform import *
import ptetmeshdeform as ptet
from scipy.io import *
from scipy.interpolate import Rbf
from scipy.spatial import KDTree
from scipy.stats.mstats import find_repeats
from ptrimeshalt import *
from ptrifeatures import *
from plinedeform import *
import pshapecontext as sh
import pmeshquality as qu
import writeFEBio as wf
import time
import os

os.system('~/Software/MRL/FEBio/febio.lnx -i ./SkullPO_000_10.feb')
os.system('~/Software/MRL/FEBio/febio.lnx -i ./SkullPO_000_20.feb')
os.system('~/Software/MRL/FEBio/febio.lnx -i ./SkullPO_100_10.feb')
os.system('~/Software/MRL/FEBio/febio.lnx -i ./SkullPO_100_20.feb')
