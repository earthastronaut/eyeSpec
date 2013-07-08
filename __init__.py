"""
eyeSpec was designed for performing various data reductions on 1-D spectra especially in the optical wavelengths
"""
#==============================================================================#
# parameters
Params = {'verbose_level':0,
          'goofy_version':False,
          'clobber':True}

#==============================================================================#
# import some modules which are used
import dependencies
from dependencies import np, os, mpl, plt, scipy, pdb, _found_wx

#==============================================================================#
# import basic functions for eyeSpec
from core import (eyeSpec_spec, convert_wavelength_units, inv_var_2_var, var_2_inv_var, verbose, asciiread)

from IO import (wlsoln_coeff_from_header, 
                readin, readin_txt, readin_spectre_files, readin_single_order_files, readin_makee, readin_apogee, readin_hst, readin_spec,
                save, save_spec, load_spec, save_txt_orders, save_txt)
                
from coadd import combine_orders
                
from plotting import plot_spec
plotspec = plot_spec
 
#==============================================================================# 

if _found_wx: from eS_interactive import *

from stellar_atmospheres import *
 
#==============================================================================#
# import SPECTRE emulator
#import SPECTRE_emulator as SP

# Just because:
def _print_eyeSpec ():
    print "================================================================================"
    print "================================================================================"
    print "                                  eeeeeoooo"
    print "                                 e0eeeeeeeoooo"
    print "                                00eeeo   0eooo"
    print "   eeo   eeooo    eeoo   eeo    00eoo    0eooo  eeeeooo      eeo       eeeeoo "
    print " eeooooo  0eoo    0eo  eeooooo    0eeoo         0eeeeeeoo  eeooooo    eeeeeeoo"
    print "0eoo  0oo  0eeo  0eo  0eoo  0oo      0eoo       0eo   0eo 0eoo  0oo  0eo    0o"
    print "0eoo  0eo    0eeoo    0eoo  0eo        0eoo     0eo   0eo 0eoo  0eo  0eo"
    print " 0eeoee       0eeo     0eeoee    eeeoo   0eeoo  0eeeeee    0eeoee    0eo"
    print "  0ee   oo     0eo      0ee   oo 0eeo    00ooo  0eo        0ee   oo 0eo    oo"
    print "   0eeeoo      0eo       0eeeoo   00eoeee0eooo  0eo         0eeeoo   0eoooooo"
    print "    00ee      0eeoo       00ee     0000eeeoo    0eoo          00ee     00eeeo"
    print "================================================================================"
    print "================================================================================"
    
if Params['goofy_version']: 
    _print_eyeSpec()
    os.system("say -v g 'welcome to eyeSpec'")


    
    
    
    
    
