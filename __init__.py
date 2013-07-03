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
from eyeSpec.dependencies import np, mpl, plt, scipy, pdb

#==============================================================================#
# import basic functions for eyeSpec
from eyeSpec.base_classes import eyeSpec_spec

from eyeSpec.base_functions import convert_wavelength_units, inv_var_2_var, var_2_inv_var, find_overlap_pts

from eyeSpec.base_IO import  wlsoln_coeff_from_header, readin, readin_txt

from eyeSpec.extended_IO import save_spec, load_spec, save_spec_orders, save_spec_txt, readin_spectre_files, readin_makee, readin_apogee, readin_hst

from eyeSpec.base_combine import spec_combine #, co_add

from eyeSpec.base_plotting import plot_spec
plotspec = plot_spec

from eyeSpec.MOOG import create_kurucz_atmo_model, moog_synth, moog_ewfind, run_moog_ewfind

from eyeSpec.moog_functions import get_model_name, read_moog_linelist, write_moog_par, write_moog_lines_in, parse_synth_summary_out, parse_synth_standard_out, parse_abfind_summary_out, process_moog_synth_output, crop_data_table
 
#==============================================================================#
# import interactive data editing
from eyeSpec.interactive_classes import Cursor, EventConnections

from eyeSpec.app_edit_var import edit_var

from eyeSpec.app_edit_rv import edit_rv

from eyeSpec.app_edit_ctm import edit_ctm

from eyeSpec.app_edit_data import edit_data

from eyeSpec.app_iplotspec import iplotspec

#==============================================================================#
# import SPECTRE emulator
#import eyeSpec.SPECTRE_emulator as SP


# Just because:
if Params['goofy_version']:
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
    print "  0ee   oo     0eo      0ee   oo 0eeo    00ooo  0eo	    0ee   oo 0eo    oo"
    print "   0eeeoo      0eo       0eeeoo   00eoeee0eooo  0eo	     0eeeoo   0eoooooo"
    print "    00ee      0eeoo       00ee     0000eeeoo    0eoo	      00ee     00eeeo"
    print "================================================================================"
    print "================================================================================"
    from eyeSpec.dependencies import os
    os.system("say -v g 'welcome to eyeSpec'")


if __name__ == "__main__":
    eyeSpec_path = ['/Users/dylangregersen/Desktop/Astrophysics/applications/eyeSpec/eyeSpec']
    
    print "--- external modules ---"
    execfile("dependencies_v1.py")
    
    print "--- classes and functions ---"
    execfile("base_classes_v3.py")
    execfile("base_functions_v1.py")
    
    execfile("base_IO_v2.py")
    execfile("extended_IO_v2.py")
    
    execfile("moog_functions.py")
    execfile("MOOG.py")
    
    
    print "--- 3rd Party ---"
    execfile("resampling.py")
    execfile("piecewise_poly.py")
    execfile("linefinder_dsg.py")
    
    print "--- combining function ---"
    execfile("base_combine_v2.py")
    
    
    print "--- plotting ---"
    execfile("base_plotting_v2.py")
    execfile("interactive_IO_v2.py")
    execfile("interactive_classes_v3.py")
    
    #    print "--- applications ---"
    #    execfile("app_edit_data_v29")
    #    execfile("app_ctm_editor_v10.py")
    #    execfile("app_snr_editor_v11.py")
    #    execfile("app_rv_editor_v4.py")
    #    execfile("app_iplotspec_v1.py")
    
    
    
    
    
    
    
    
    
    
    
    
