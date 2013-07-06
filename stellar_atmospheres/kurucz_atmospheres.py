# This is used for creating kurucz atmospheres

#=============================================================================#
# import modules
from eyeSpec.dependencies import  np, os, subprocess #@UnresolvedImport
from eyeSpec.core import get_filename #@UnresolvedImport
from eyeSpec.stellar_atmospheres import executables #@UnresolvedImport
from moog_functions import get_model_name

pass
#=============================================================================#
def check_makekurucz (error=True):
    _makekurucz_avail_error = "Must have makekurucz executable declared to use this function not: "+executables['makekurucz']
    check =  os.path.isfile(executables['makekurucz'])
    if not check and error: raise ValueError(_makekurucz_avail_error)
    return check
        
pass
#=============================================================================#
# functions for creating the kurucz atmospheres

def create_kurucz_atmo_model (teff,logg,feh,turb,modtype='ODFNEW',filename='FINALMODEL',verbose=True, clean_up=True):
    """
PURPOSE:
   To use the makekurucz Fortan interpreter to create stellar atmosphere models for MOOG 
   
CATEGORY:
   Stellar Atmospheres

INPUT ARGUMENTS:
    teff : (float) stellar effective temperature
    logg : (float) gravitational acceleration at the surface of the star
    feh  : (float) Normalized solar metallicity [Fe/H]
    vt   : (float) stellar microturbulence
    
INPUT KEYWORD ARGUMENTS:
    modtype : (string) the type of model used, Possible Choices ->
        ODFNEW -> 
        AODFNEW -> Same as the ODFNEW only with an alpha enhancement
        NOVER ->
        KUROLD ->
    filename : (string) Give the output filename. If 'rename' it will use the naming convention
                {teff}p{logg}mp{feh}p{vt}.modtype
    verbose : (bool) If True you get verbose output
    clean_up : (bool) Wll remove files : M1, M2, MOD*

OUTPUTS:
   (bool) Model creation, if the parameters were correct and the model created it returns True

DEPENDENCIES:
   External Modules Required
   =================================================
   subprocess 
   
   External Functions and Classes Required
   =================================================
    get_model_name  
       
NOTES:
   (1) There are limits in the interpolation of Kurucz. If a limit is passed the function will return False
        3500 < teff < 10000
        0 < logg < 5
        -5 < feh < 1.0

    (2) If you give a feh > 1.0 then it will multiply by negative 1 (e.g. 2.13==> -2.13)

    (3) The individual modtypes may have inherent limits as well (e.g. feh in ODFNEW >= -2.5)
    
    (4) This function won't work if the executables['makekurucz'] is not declared
    
    (5) This will also overwrite any file names 'FINALMODEL'
    
EXAMPLE:
   >>> create_kurucz_atmo_model(5000, 4.1, -2.13, 1.1, modtype='aodfnew', filename='rename')
   >>>
   >>> ls
   5000p410m213p110.aodfnew

MODIFICATION HISTORY:
    13, Jun 2013: Dylan Gregersen
    """
    check_makekurucz(True)
    MAKEKURUCZ_EXE = executables['makekurucz']
    
    # impose kurucz limits
    if not (3500 < teff < 10000):
        if verbose: print "Teff is out of kurucz range (3500,10000): "+str(teff)
        return False

    if not (0 < logg < 5):
        if verbose: print "Logg is out of kurucz range (0,5): "+str(logg)
        return False

    if feh > 1.0: 
        if verbose: print "[Fe/H] given larger than 1, multiplying by negative one"
        feh *= -1.0
        
    if not (-5 < feh < 1.0):
        if verbose: print "[Fe/H] is out of kurucz range (-5,1): "+str(feh)
        return False

    modtype=modtype.upper()
    pos_modtypes = ['KUROLD', 'ODFNEW', 'AODFNEW', 'NOVER']
    if modtype not in pos_modtypes:
        if verbose: print "Model type must be in: "+", ".join(pos_modtypes)
        return False

    inp = " ".join([str(teff), str(logg), str(feh), str(turb),"\n"+modtype])+"\n"

    devnull = open('/dev/null', 'w')
    makeKurucz = subprocess.Popen(MAKEKURUCZ_EXE, stdin=subprocess.PIPE,stdout=devnull)
    devnull.close()
    makeKurucz.communicate(input=inp)


    if filename == 'rename': filename = get_model_name(teff,logg,feh,turb,modtype)

    if filename != 'FINALMODEL': os.system("mv FINALMODEL "+str(filename))

    if clean_up: _ = os.popen3('rm -f MOD* M1 M2')
        
    return True

def run_create_atmo_model ():
    """
PURPOSE:
   This interactively runs create_kurucz_atmo_model and prompts the user for the inputs
   
CATEGORY:
   Stellar Atmospheres

DEPENDENCIES:
   External Modules Required
   =================================================
    os, numpy
   
   External Functions and Classes Required
   =================================================
    create_kurucz_atmo_model
    
NOTES:
    
EXAMPLE:
   >>> run_create_atmo_model()

MODIFICATION HISTORY:
    13, Jun 2013: Dylan Gregersen

    """  
    check_makekurucz(True)
    def convert_modtype (modtype):
        if modtype.lower()[0] == 'o': return 'ODFNEW'
        if modtype.lower()[0] == 'a': return 'AODFNEW'
    
    inmodel = raw_input("Please give model parameters: teff  logg  feh  vmicro _OR_ model_filename\n")
    inmodel = inmodel.split()

    if len(inmodel) == 1:
        mod_fname = get_filename("Please give MOOG model file",'r',inmodel[0])
        if os.path.abspath(mod_fname) != os.path.abspath('./FINALMODEL'): os.system("cp "+mod_fname+" ./FINALMODEL")
    else:
        teff,logg,feh,vmicro =  np.array(inmodel[:4],dtype=float)
        modtype = raw_input("Please give model type: "+", ".join(['KUROLD', 'ODFNEW', 'AODFNEW', 'NOVER'])+"\n")
        modtype = convert_modtype(modtype)
        while True:
            check = create_kurucz_atmo_model(teff,logg,feh,vmicro,modtype,filename='FINALMODEL',verbose=True)
            if not check:
                inmodel = raw_input("Please retry model params: teff logg feh vmicro model_type:")
                try:
                    teff,logg,feh,vmicro =  np.array(inmodel[:4],dtype=float)
                    modtype = convert_modtype(inmodel[4])
                except: pass
            else: break