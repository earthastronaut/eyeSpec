#@PydevCodeAnalysisIgnore
## This imports modules eyeSpec is dependent on and performs some checks
# PURPOSE:
#   This file has all the dependencies used throughout eyeSpec
#   This allows for the control of versions
#   
# CATEGORY:
#   eyeSpec core
#
# MODIFICATION HISTORY:
#    13, June 2013: Dylan Gregersen
#    3,  July 2013: Dylan Gregersen
#                   more notes about necessary functions
#    5,  July 2013: Dylan Gregersen
#                   added comments which the code analysis tool can understand


#==============================================================================#
# Modules
import math 
import os 
import sys 
# !! check the version of python to be < 3.0


import time
from copy import deepcopy
import pdb 
import operator
iget = operator.itemgetter
import threading 
import Queue 
import subprocess 
import re 


#==============================================================================#
# astropy
import astropy 
from astropy.io import fits
from astropy import units
from astropy import constants

#==============================================================================#
# wxPython
_found_wx = True

try: 
    import wxversion
except: 
    _found_wx = False

if _found_wx: 
    wxversion.ensureMinimal('2.8')

def _check_for_wx ():
    if not _found_wx: 
        raise ImportError("This package of eyeSpec relies on wxPython and wxversion which weren't found")

#==============================================================================#
# pickle
try: 
    import cPickle as pickle 
except ImportError:
    raise ("module cPickle is required")

#==============================================================================#
# numpy
import numpy as np
nmajor, nminor = [int(n) for n in np.__version__.split('.')[:2]]
if not (nmajor > 1 or (nmajor == 1 and nminor >= 6)):
    raise ImportError(
            'numpy 1.6 or later is required; you have {0}'.format(np.__version__))

import numpy.lib.recfunctions as np_recfunc 


#==============================================================================#
# scipy
import scipy 
nmajor, nminor = [int(n) for n in scipy.__version__.split('.')[:2]]
if not (nmajor > 0 or (nmajor == 0 and nminor >= 10)):
    raise ImportError(
            'scipy 0.10 or later is required; you have {0}'.format(scipy.__version__))

from scipy import sparse 

#==============================================================================#
# Matplotlib
import matplotlib # Must have matplotlib for any eyeSpec internal plotting
mpl = matplotlib

if _found_wx:
    if matplotlib.rcParams['backend'] != 'WXAgg':
        if sys.platform == 'darwin': # on a Mac
            if sys.version.find('64-bit') != -1:
                print "DANGER WILL ROBINSON: on a Mac OSX the WxAgg backend for the 64-bit EPD version of python is currently unavailable. You can TRY to fix this, but I did and (I assume) the EPD folks did and neither of us have found a solution. That said, eventually wxPython will be working with the libsraries of the Mac 64-bit system. Until then I suggest running this on a system which supports wxPython (32-bit on Mac works and is the most thoroughly tested). For now I'm going ahead with switching to the WxAgg backend but it'll probably crash on you. Have a great day!"
                # the other option is to write a different program, similar to below.
                # this basically is an interactive way to output a list of x-points which are the overlaps for different orders, as well as delete parts of orders.
        matplotlib.rcParams['backend'] = 'WXAgg'

from matplotlib.figure import Figure 
from matplotlib.widgets import Button 
import matplotlib.pylab as plt 
#from matplotlib.ticker import LinearLocator, MultipleLocator
from matplotlib.pylab import FormatStrFormatter, savefig 
from matplotlib.path import Path 


#==============================================================================#
# resampling
from .utils import resampling 


