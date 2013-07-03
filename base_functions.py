if __name__ != '__main__':
    import eyeSpec
    from eyeSpec.dependencies import np, os, time, deepcopy, pdb, pickle
    from eyeSpec.resampling import get_resampling_matrix, Gaussian_Density

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
    print "Currently not available, see code for notes"
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
    previous_ord = []
    i = 0
    for ord in spec:
        if len(previous_ord)==0:
            previous_ord.append([np.min(ord[0]),np.max(ord[0])])
            continue
        # check to make sure we're going the right way
        if previous_ord[-1][0] > np.min(ord[0]):
            print "WARNING: PREVIOUS HAS LARGER VALUE" #!! error
        
        
        for i in range(len(previous_ord)):
            i = -1*(i+1) # go in reverse order
            if previous_ord[i][1] > np.min(ord[0]) and i==-1:
                midpt = (previous_ord[i][1] + np.min(ord[0]))/2
                overlap_pts.append(midpt)

            elif previous_ord[i][1] > np.min(ord[0]) and i!=-1:
                print 'WARNING: OVERLAP FOUND FOR MORE THAN 2 ORDERS WHEN LOOKING AT RANGE:',np.min(ord[0]),'to',np.max(ord[0])
                
        previous_ord.append([np.min(ord[0]),np.max(ord[0])])

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
    question = str(question)
    while True:
        if prompt is None: 
            if default_answer == 'n': prompt = question+"('yes',['no'])\n"
            elif default_answer == 'y': prompt = question+"(['yes'],'no')\n"
            else: prompt = question+"('yes','no')\n"
        choice = raw_input(str(prompt))
        if choice.lower() in ['n','no']: return False #'n'
        elif choice.lower() in ['y','yes']: return True #'y'
        elif choice.lower() == '':
            if default_answer == 'n': return False #'n'
            elif default_answer == 'y': return True #'y'
        else: print "Please answer 'yes','no','y', or 'n'"

def _find_matches (xrange,obj):
    """
    Find the orders in obj which have values within the xrange

    INPUTS:
    xrange : array type of length 2
    obj : of type eyeSpec_spec

    """
    obj._private_check(0) # make sure obj.use_band is appropriate

    wl_arr = deepcopy(obj._wl[obj.use_band])
    
    orders = []
    
    for i in xrange(len(wl_arr)):
        dmin = np.min(wl_arr[i])
        dmax = np.max(wl_arr[i])
        dmid = (dmin+dmax)/2.
        
        if (dmin > min(xrange) or dmid > min(xrange) or dmax > min(xrange)) and  (dmin < max(xrange) or dmid < max(xrange) or dmax < max(xrange)):
            orders.append(i)
            
    return orders
        
def get_filename (propmt='default',iotype='r',default=None,enter_multi = False,filename=None,find_filename=False):
    """
    Prompts user for a file name and then does some checks

    INPUTS:
    iotype : 'r' or 'w' to see whether user wants to read or write information
    enter_multi : (NOT CURRENTLY SUPPORTED) enter multiple file names
    prompt : give a different prompt
    filename : give a name to check

    """
    only_check_file = False
    if filename is not None: only_check_file = True

    def _enter_filename ():
        if prompt == 'default': fname = raw_input("ENTER FILENAME: \n")
        else:fname = raw_input(str(prompt)+" \n")
        return fname.strip()
        
    while True:
        re_enter = False
        break_all = False

        if only_check_file: fname = str(filename)
        else:
            print " "
            fname = _enter_filename()
            # check input
            if len(fname.split())==0:
                # or I could create a file list if more are entered
                print "please enter a file name"
                re_enter = True

            # enter a list of files:
            if len(fname.split()) > 1 and not enter_multi:
                print "plese only enter one file name"
                re_enter = True


        if not re_enter:
            fname = fname.split()
            if fname[0] in ['a','abort']:
                return default

        # check if file exists
        if not re_enter:
            new_f = {}
            for i in range(len(fname)):
                ffile = fname[i]

                if os.path.exists(ffile) and iotype == 'w':
                    while True:
                        clobber = raw_input("WARNING: FILE EXISTS, OK TO OVERWRITE? (y,[n])\n")
                        if clobber in ['yes','y']:
                            clobber = True
                            break
                        elif clobber in ['no','n','']:
                            if find_filename:
                                print "PLEASE ENTER NEW FILENAME"
                                re_enter = True
                            else:
                                clobber = False
                                re_enter = False
                            break
                        else: print "please enter yes or no"

                            
                        
                elif not os.path.exists(ffile): # file does not exist
                    clobber = False
                    if iotype == 'r':
                        #new_f[i] = raw_input("ERROR: File <"+ffile+"> does not exist, enter new file name or abort (a):\n")
                        #if new_f[i] == 'a':
                        #    break_all= True
                        #    break
                        print "ERROR: File does not exist: '"+ffile+"'"
                        if find_filename:
                            print "PLEASE ENTER NEW FILENAME"
                            re_enter = True
                        else:
                            re_enter = False

        # check to make sure
        if break_all or not re_enter: break
    return fname[0]

def deletevars(varlist):
    for var in varlist:
        if var in globals(): del var

def variablename(var):
     import itertools
     # !! put in a try and error statement incase the variable is not in globals
     return [tpl[0] for tpl in itertools.ifilter(lambda x: var is x[1], globals().items())][0]

