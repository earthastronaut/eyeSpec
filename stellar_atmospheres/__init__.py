# the code in this directory deals with stellar atmospheres

# !! todo take the code fro moog and makekurucz and include it in __path__[0]+"/core/"

#===============================================================================#
######################## USER SETUP #############################################

moog_exe = '/Applications/astro/Moog/MOOG'

moog07_exe = "/uufs/astro.utah.edu/common/astro_data/products/moog2007-3/MOOG"

moogsilent_exe = '/Applications/astro/Moog/MOOGSILENT'

makekurucz_exe = '/Applications/astro/makekurandy/makekurucz3.e'

######################## USER SETUP #############################################
#===============================================================================#

# Modules
from kurucz_atmospheres import create_kurucz_atmo_model, run_create_atmo_model

from moog import synth as moog_synth
from moog import ewfind as moog_ewfind
from moog import run_ewfind, ewfind_model

from moog_functions import (get_model_name, read_moog_linelist, write_moog_par, write_moog_lines_in,  #@UnusedImport
                            parse_synth_summary_out, parse_synth_standard_out, parse_abfind_summary_out,  #@UnusedImport
                            process_moog_synth_output, crop_data_table) #@UnusedImport
