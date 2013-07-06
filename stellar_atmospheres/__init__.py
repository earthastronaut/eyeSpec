# the code in this directory deals with stellar atmospheres


from MOOG import create_kurucz_atmo_model, moog_synth, moog_ewfind, run_moog_ewfind
from moog_functions import get_model_name, read_moog_linelist, write_moog_par, write_moog_lines_in, parse_synth_summary_out, parse_synth_standard_out, parse_abfind_summary_out, process_moog_synth_output, crop_data_table