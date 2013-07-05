if __name__ != '__main__':
    import pdb #@UnusedImport
    from dependencies import np, os, deepcopy
    from resampling import get_resampling_matrix, Gaussian_Density

pass
########################################################################################
# various array manipulations
def reduce_output_shape (arr):
    shape = arr.shape
    new_shape = ()
    for i in xrange(len(shape)): 
        if shape[i] != 1: new_shape += (shape[i],)
    return arr.reshape(new_shape)
      
def np_vstack_append (arr,item):
    dtype = arr.dtype.name    
    if len(arr) == 0: arr = np.array([item],dtype=dtype)
    else: arr = np.vstack((arr,item))
    return arr

def np_vstack_delete (arr,index):
    dtype = arr.dtype.name
    N = len(arr)
    if N == 0: return np.array([],dtype=dtype)
    
    if index == 0: outarr = arr[1:]
    elif index == N-1: outarr = arr[:N-1]
    else: outarr = np.concatenate((arr[:index],arr[index+1:]))
    return outarr

pass
########################################################################################
 
def smooth_by_scaler (xy,smoothing_scaler=1,transform_matrix=None):
    """
    takes data points xy.shape = (2,number_points) assumes first column is wavelength and then smoothes Gaussian width which is a linear factor across wavelength

    """    
    if smoothing_scaler == 1: return xy
    if transform_matrix is None:     
        gd = Gaussian_Density(xy[0],smoothing_scaler*xy[0])
        T = get_resampling_matrix(xy[0],xy[0],gd,preserve_normalization=True)
    else: T = transform_matrix
    xy[1] = T*xy[1]
    return xy   
 
def smooth_data_to_resolution (spec_obj,R):
    """
    Smooth some data to a lower resolution
    
    ============    ===========================================================
    Keyword         (type) Description
    ============    ===========================================================
    spec_obj        (eyeSpec_spec) eyeSpec spectrum object
    
    R               (float) the new resolution to smooth to (= delta_lambda/lambda)
    ============    ===========================================================

    """
    print "Currently not available, see code for notes",spec_obj,R
    # check current resolution of spec_obj
    
    # new_wls = wls.copy()
    
    # for i in xrange(1,len(new_wls)):
    #     del_lambda = R*new_wls[i]
    #     new_wls[i] = new_wls[i-1]+del_lambda
    
    # T = get_resampling_matrix(wls,new_wls)
    # new_data = T*data
    
    # new_var = T*var  # << check this
    
    pass 
    
pass
########################################################################################
# plot functions

def figure_adjust_borders(fig, targets):
    "Translate desired pixel sizes into percentages based on figure size."
    # code thanks to user samplebias at http://stackoverflow.com/questions/6066091/python-matplotlib-figure-borders-in-wxpython
    dpi = fig.get_dpi()
    width, height = [float(v * dpi) for v in fig.get_size_inches()]
    conversions = {'top': lambda v: 1.0 - (v / height),
                   'bottom': lambda v: v / height,
                   'right': lambda v: 1.0 - (v / width),
                   'left': lambda v: v / width,
                   'hspace': lambda v: v / height,
                   'wspace': lambda v: v / width}
    opts = dict((k, conversions[k](v)) for k, v in targets.items())
    fig.subplots_adjust(**opts)

def alt_order_colors (i):
    """ 
    This will rotate through colors returning the current color and color scheme index

    can give 1/12 to get the second entry of an array of length 12

    """
    if type(i).__name__ == 'str':
        parts = i.split("/")
        i = parts[0]
        # !! doesn't work like envisioned

    c_i = int(i)
    ColorScheme = ['#003399', # Dark Blue
                   '#996600'] # Pale Orange
    
    # use these to cycle around, e.g. if you past the end go the the beginning
    if c_i < 0: c_i = len(ColorScheme)-1
    elif c_i >= len(ColorScheme): c_i = 0

    return [ColorScheme[c_i],c_i]

pass
########################################################################################
# data operations

def inv_var_2_var (inv_var,fill=1e50):
    inv_var = np.array(inv_var,dtype=float)
    zeros = (inv_var == 0.0)
    bad = (inv_var >= fill)
    inv_var[zeros] = -1.0
    var = 1.0/inv_var
    var[zeros] = fill
    var[bad] = 0.0
    return var

def var_2_inv_var (var,fill=1e50):
    var = np.array(var,dtype=float)
    zeros = (var==0)
    bad = (var>=fill/2.0)
    
    var[zeros] = -1.0
    inv_var = 1.0/var
    
    # set points which are very large to the fill
    inv_var[zeros] = fill
    # set points which are almost zero to zero
    inv_var[bad] = 0.0

    return inv_var

def convert_wavelength_units (obj,from_units="angstrom",to_units="angstrom",take_log=False,reverse_log=False):
    from_units = str(from_units)
    to_units = str(to_units)
    
    name_conventions = {'a':'angstrom',
                        'nm':'nanometer',
                        'um':'micrometer',
                        'ft':'feet'}

    if from_units in name_conventions.keys():
        from_units = name_conventions[from_units]

    if to_units in name_conventions.keys():
        to_units = name_conventions[to_units]


    unit_to_m = {'angstrom':1e-10,
                 'nanometer':1e-9,
                 'micrometer':1e-6,
                 'feet':0.3048}

    if from_units not in unit_to_m.keys(): raise ValueError("From units must be in: "+", ".join(unit_to_m.keys()))
    if to_units not in unit_to_m.keys(): raise ValueError("From units must be in: "+", ".join(unit_to_m.keys()))


    # Convert to a meter
    wl = obj._wl
    wl *= unit_to_m[from_units]

    # Convert to desired from meter
    wl *= 1.0/unit_to_m[to_units]
    
    if take_log:
        wl = np.log10(wl)

    if reverse_log:
        wl = 10**(wl)

    return obj
    
def find_overlap_pts (spec):
    """
    This finds the overlap points simply by finding the midpoint between the max of one order and the min of the next order
    INPUTS:
    spec : (eyeSpec_spec) should already have been check to make sure user gave the proper object
    
    OUTPUS:
    overlap_pts : (list) midpoints shared by orders

    """
    overlap_pts = []
    previous_order = []
    i = 0
    for order in spec:
        if len(previous_order)==0:
            previous_order.append([np.min(order[0]),np.max(order[0])])
            continue
        # check to make sure we're going the right way
        if previous_order[-1][0] > np.min(order[0]):
            print "WARNING: PREVIOUS HAS LARGER VALUE" #!! error
        
        
        for i in range(len(previous_order)):
            i = -1*(i+1) # go in reverse orderer
            if previous_order[i][1] > np.min(order[0]) and i==-1:
                midpt = (previous_order[i][1] + np.min(order[0]))/2
                overlap_pts.append(midpt)

            elif previous_order[i][1] > np.min(order[0]) and i!=-1:
                print 'WARNING: OVERLAP FOUND FOR MORE THAN 2 ORDERS WHEN LOOKING AT RANGE:',np.min(order[0]),'to',np.max(order[0])
                
        previous_order.append([np.min(order[0]),np.max(order[0])])
    
    
    if len(overlap_pts) == 0: #!! AND VERBOSE
        print 'NO OVERLAP FOUND'
    return overlap_pts

pass
########################################################################################
# Misc functions and operations

def get_bounds (prompt,lower_too=False,default=(0,1e20),display_help='No Help Available'):
    while True:
        inbounds = raw_input(prompt)
        if inbounds.lower() in ['help','h','he']:
            print "-"*60
            print display_help
            print "-"*60
            print " "
            continue
        elif inbounds == '': return default
        else:
            spl = inbounds.strip().split()
            try:
                if len(spl)==1 and lower_too: return (float(spl[0]),1e10)
                else: return (float(spl[0]),float(spl[1]))
                
            except: pass
        print "Invalid input. Type 'help' for more info."
       
def yesno (question="Please choose:",default_answer='n',prompt=None):
    """
PURPOSE:
    Prompt for a yes or no question and return True or False respectively
   
CATEGORY:
    User functions

INPUT ARGUMENTS:
    None

INPUT KEYWORD ARGUMENTS:
   question : (string) This becomes the prompt with ('yes','no') appended to the end
   default_answer : (string or boolean) if True or 'y' then enter gives that True, conversely for False and 'n'
   prompt: (string) Over-rules the question and just uses the string for the prompt
   
OUTPUTS:
    (boolean) True if answered yes, False if answered no
       
DEPENDENCIES:
   External Modules Required
   =================================================
   None
   
   External Functions and Classes Required
   =================================================
    None
       
NOTES:
   (1) If default answer is given then the ('yes','no') will become ('yes',['no']) or (['yes'],'no') with the [] giving the default value when enter is hit
    

EXAMPLE:
   >>> if yesno("Please enter yes or no",'n'): print "Yes!"*100

MODIFICATION HISTORY:
    15, July 2013: Dylan Gregersen

    """   
    if type(default_answer) == bool:
        if default_answer: default_answer='y'
        else: default_answer = 'n'
    else: default_answer = default_answer.lower()[0]
    
    question = str(question)

    if prompt is None: 
        if   default_answer == 'n': prompt = question+"('yes',['no'])\n"
        elif default_answer == 'y': prompt = question+"(['yes'],'no')\n"
        else: prompt = question+"('yes','no')\n"
    
    # get the user answer
    while True:
        choice = raw_input(str(prompt)).lower()
        if choice in ('n','no'): return False #'n'
        elif choice in ('y','yes'): return True #'y'
        elif choice == '' and default_answer=='n': return False
        elif choice == '' and default_answer=='y': return True
        else: print "Please answer 'yes','no','y', or 'n'"

def _find_matches (xbounds,obj):
    """
    Find the orders in obj which have values within the xbounds

    INPUTS:
    xbounds : array type of length 2
    obj : of type eyeSpec_spec

    """
    obj._private_check(0) # make sure obj.use_band is appropriate

    wl_arr = deepcopy(obj._wl[obj.use_band])
    
    orders = []
    
    for i in xbounds(len(wl_arr)):
        dmin = np.min(wl_arr[i])
        dmax = np.max(wl_arr[i])
        dmid = (dmin+dmax)/2.
        
        if (dmin > min(xbounds) or dmid > min(xbounds) or dmax > min(xbounds)) and  (dmin < max(xbounds) or dmid < max(xbounds) or dmax < max(xbounds)):
            orders.append(i)
            
    return orders
        
def get_filename (prompt='ENTER FILENAME:', iotype='r', default=None, enter_multi = False, filename=None, find_filename=False):
    """
PURPOSE:
   To interactively get filenames
   
CATEGORY:
   User functions

INPUT ARGUMENTS:
   None

INPUT KEYWORD ARGUMENTS:
   prompt  : (string) The string to appear as the prompt for giving files
   iotype  : (string) either 'r' or 'w' for read or write
   default : (--) This object will be returned as a default value from the function
   enter_multi : (bool) If True it will prompt for entering multiple files
   filename : (string) a filename to check and then prompt for another if it isn't appropriate entry
   fine_filename : (bool) If True then if the file is not found it will prompt for finding it again
   
OUTPUTS:
   if enter_multi: (list) gives a list of the entered values or the default 
   else: (string) gives the filename or the default

DEPENDENCIES:
   External Modules Required
   =================================================
    os, raw_input, deepcopy
   
   External Functions and Classes Required
   =================================================
    yesno
       
NOTES:
   (1) Appropriate entry means that if iotype if 'r' then the file should exist if iotype if 'w' then it shouldn't or if prompts for overwriting the file

   (2) This creates lists of the files using set so you won't get repeats of files 

EXAMPLE:
   >>> fname = get_filename()

MODIFICATION HISTORY:
    5, July 2013: Dylan Gregersen
                   
    """  
    base_prompt = deepcopy(prompt)
    
    check_filename = (filename != None)
      
    files = set([])
      
    i = 0
        
    # Loop     
    while True:
        if enter_multi: prompt = "["+format(i,"2")+"] "+base_prompt
        
        # Get the filename
        if check_filename: fname = str(filename)
        else:            
            fname = raw_input(str(prompt))
            
            # view options
            if fname == '\t': 
                print "PWD>> "+os.path.abspath('.')
                print "-"*60
                os.system('ls')
                print ""
                continue
            # ignore blank lines
            if len(fname.split())==0: 
                if enter_multi: print "To stop entering files enter '.'"
                continue
            
            # spaced objects
            if len(fname.split()) > 1:
                print "Please only enter one file name at a time"
                continue
            
            # abort entering files
            if   fname in ('a','abort'): return default
            elif fname in ('h','help','?'):
                print "This routine allows you to enter files and will check to make sure you gave an appropriate response\n"
                if enter_multi: print "Entering multiple files, to stop enter filename '.'"
                
                print "The current default output is : "+str(default)
                print ""
                print "-"*60
                print "options:"
                print "[a]bort - will break the file enter loop"
                print "[h]elp  - will display this help screen"
                print "tab     - will display 'ls' of the PWD"
                print ".       - when entering multiple files this acts as an EOF"
                print ""
                continue
                
            if enter_multi and fname == '.': break
        
        
        good_file = True
        
        # check if file exists and trying to write over
        if os.path.isfile(fname) and iotype == 'w':
            if check_filename: print "File given '"+fname+"'"
            if yesno("WARNING: File exists, OK to overwrite?",'n'): good_file = True
            else: good_file = False
                    
        elif not os.path.isfile(fname) and iotype == 'r': # file does not exist
            if check_filename: print "File given '"+fname+"'"
            print "File does not exist "
            good_file = False

            
        # you can only check the first filename given
        check_filename=False    
        
        # check 
        if good_file:
            if enter_multi and find_filename and i == 0 : print "[ 0] gave file "+fname
             
            files.add(fname)
            i += 1
            
        # else try again
        elif not good_file and find_filename: continue
        
        if not enter_multi: break

    files = list(files)        
    if   len(files) == 0: return default
    elif enter_multi: return files
    else: return files

def deletevars(varlist):
    for var in varlist:
        try: del var
        except: pass
        


