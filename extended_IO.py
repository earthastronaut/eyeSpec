
if __name__ != '__main__': 
    from eyeSpec.dependencies import np, os, time, deepcopy, pyfits, pickle, pdb
    from eyeSpec.base_classes import query_fits_header, eyeSpec_spec, eyeSpec_spec
    from eyeSpec.base_IO import readin
    from eyeSpec.base_functions import inv_var_2_var, var_2_inv_var


def save_spec (spec_obj,filename=None,unique_name=True,clobber=False):
    fsuffix = '.spec' #'.pkl'
    
    """ 
    this is used to save a pickled version of the current object
    
    INPUTS:
    ============  ==============================================================
    keyword       (type) Description
    ============  ==============================================================
    spec_obj      (eyeSpec_spec) A eyeSpec spectrum object
    filename      (str) If '.spec' is not added to the end this string will 
                     be used as a prefix for the full filename. 
                     filename = 'foo'  ===> 'foo.spec'
                     Otherwise the filename will be used.
                     If None then the filename will be generic "SpecObjSave"
    unique_name   (bool) if True then it will make sure the file saves by 
                     checking all file names and adding _# if multiples appear
    clobber       (bool) if True it will delete existing files, if False it
                     will raise an error if the file exists,
                     Not applicable if unique_name is True 
    ============  ==============================================================

    """
    # !! I could change this to be a special file
    # Header: contains actual header, info, and notes strings
    # DATA: contains wavelength, data, inv_var information

    unique_name = bool(unique_name)
    clobber = bool(clobber)

    #--------------------------------------------#
    if filename is None:
        filename = 'SpecObjSave'
        #if obj._obj_name is not None: filename += '_'+obj._obj_name
    else: 
        fileprefix = str(filename).strip(fsuffix)
        fileprefix = str(fileprefix)
        if type(fileprefix).__name__ in ['str','string_']: filename = fileprefix
        else: raise TypeError("fileprefix must be of type string not "+type(fileprefix).__name__)
      
    filename += fsuffix
  
    #--------------------------------------------#
    if unique_name:
        # save only unique names by checking if the file already exists
        i = 1
        u_filename = deepcopy(filename)
        while True:
            if os.path.exists(u_filename):
                u_filename = deepcopy(filename[:-4]+"_"+str(i)+filename[-4:])
                i+=1
            else: break
    else:
        if os.path.exists(filename) and not clobber: raise IOError("clobber is set to false and given file exists: "+filename)

    #--------------------------------------------#
    pickle.dump(spec_obj,open(filename,'wb'))

def load_spec (filename): 
    """
    Readin a specified eyeSpec class object from a spectrum file
    """
    # !! can I check if the file given is a pickle file
    out_spec_obj = pickle.load(open(filename,'rb'))
    if out_spec_obj.__class__.__name__ not in ['eyeSpec_spec','eyeSpec_fits']: raise IOError("Loaded an object which is not of class eyeSpec_spec or eyeSpec_fits")
    return out_spec_obj

def save_spec_orders (spec_obj, base_name, band='default',use_cropped=False,order=None,clobber=True,include_varience=True,comment='#',divide_header=True):
    """
    Splits a eyeSpec spectrum object into orders with the file name set by the base name and the order number
    e.g. base_name = 'my_spectrum'
    order 0 = my_spectrum_0.txt
    order 1 = my_spectrum_1.txt
    order 2 = my_spectrum_2.txt
    order 3 = my_spectrum_3.txt
    ....and so on
    

    """
    
    orders = range(spec_obj.shape[1])
    
    for i in orders:
        file = str(base_name)+"_"+str(i)+".txt"
        save_spec_txt(spec_obj,file,
                      band=band,use_cropped=use_cropped,order=[i],
                      clobber=clobber,include_varience=include_varience,comment=comment,divide_header=divide_header)

def save_spec_txt (spec_obj,filename,band='default',use_cropped=False,order=None,clobber=True,include_varience=True,divide_orders=True,comment='#',divide_header=True):
    """
    Outputs the eyeSpec spectrum class into a given file as text data.

    INPUTS:
    =============   ============================================================
    keyword         (type) Description
    =============   ============================================================
    spec_obj        (eyeSpec_spec) spectrum class for eyeSpec
    filename        (str) This gives the filename to save the as
    band            (int,'default') This tells which of the first dimensions of
                      spec_obj to use. 'default' is spec_obj.get_band()
    use_cropped     (bool) If True it will crop off points at the begining and 
                      end of orders which have inverse varience = 0, i.e. have
                      inf errors
    order           (int,array,None) you can specify which orders are output
                      If None then it will output all possible orders
    clobber         (bool) If True then the function will overwrite files of the
                      same file name that already exis

    include_varience (bool) If True the third column which gives the varience 
                      will be included
    divide_orders    (bool) If True it will but a commend line with '#' between
                      each of the orders
    comment          (str) What symbol to use as a comment
    divide_header    (bool,None) If False it will give one long string as the first header line
                                 If True it will divide it up by 80 character lines with a comment of '#:' 
                                 If None then no header will be printed
    =============   ============================================================

    """
    #===============================================#
    # check inputs
    if spec_obj.__class__.__name__ != 'eyeSpec_spec': raise ValueError("spec MUST BE OF CLASS eyeSpec_spec")

    # double check booleans
    clobber = bool(clobber)
    include_varience = bool(include_varience)
    divide_orders = bool(divide_orders)
    if divide_header is not None:
        divide_header = bool(divide_header)

    # if order is not None then convert and check
    if order is not None:
        try: order = np.array(order,dtype=int)
        except: raise ValueError("Order must be convertable to an integer ndarray")
        if order.ndim == 0: order = np.array([order])
        elif order.ndim > 1: 
            print "HeadsUp: Multiple dimensions give for order, taking only the first"
            order = order[0]

    # if band is not default check it and output
    if band != 'default':
        band = spec_obj._check_band_num(band=band)
        spec_obj.set_band(band)

    # check order_delimiter
    if type(comment).__name__ not in ['str','string_']: raise ValueError("Please give comment delimiter as a string")
    
    # filename:
    if type(filename).__name__ not in ['str','string_']: raise ValueError("Filename must be a string")

    if os.path.exists(filename) and not clobber: raise IOError("Clobber set to False and file exists: '"+filename+"'")
    

    #===============================================#
    # edit header to reflect current information 
    # !! not exactly sure yet


    #===============================================#
    # open and write to file
    f = open(filename,'w')


    # !! add check for multiple headers

    if divide_orders:
        aline = comment+"eSorders - eyeSpec multiple order text format"
        if divide_header: aline += "\n"
        else: aline = format(aline,'80')
        f.write(aline)

    #--------------------------------#
    # add header
    if divide_header is not None:
        header_line = ""
        for card in spec_obj.header.ascardlist():
            card = str(card).strip()
            if divide_header:
                f.write(comment+": "+card)
                f.write("\n")
            else: header_line += card
            
        if not divide_header:
            f.write(comment+" ")
            f.write(header_line.strip())
            f.write("\n")

    #--------------------------------#
    # add all data
    # walk through orders
    if order is not None: ran = order
    else: ran = range(spec_obj.shape[1])

    use_cropped = bool(use_cropped)
    spec_obj.set_use_cropped(use_cropped)
    for i in ran:
        # output information
        wl = spec_obj.get_wl(i)
        dat = spec_obj.get_data(i)
        inv_var = spec_obj.get_inv_var(i)

        if divide_orders and len(ran) > 1: 
            f.write(" ".join([comment,"order",str(i),str(len(wl)),"\n"]))
        # go through data points
        for j in range(len(wl)):
            out = [format(wl[j],'>15.10'),
                   format(dat[j],'>15.10')]
            if include_varience:
                if abs(inv_var[j]) < 1e-30: var = 1e30
                else: var = 1.0/inv_var[j]
                out.append(format(var,'>10.10'))
            out.append("\n")
            f.write("  ".join(out))

    spec_obj.set_use_cropped('previous')
    #===============================================#
    f.close()

def readin_spectre_files (filelist,relative_paths=False,verbose=True):
    """
    takes a list of spectre 1D files and creates a single object

    INPUTS:
    =============  =============================================================
    keyword        (type) Description
    =============  =============================================================
    filelist       (string) give the name of a file which contains a list of 
                           spectre files (fits/txt)
                    OR
                   (array) which gives each file name
    relative_paths (bool) if True eyeSpec will look for each file name in the 
                   filelist relative to the current directory. If False it will
                   take the absolute path of the filelist file as the base name
                   to look for each file. Not applicable if filelist is array
    =============  =============================================================
    """
    list_of_files = []
    relative_paths = bool(relative_paths) # !! I think I want to do this every time

    #============================================================#
    # import the files
    #-----------------------------------------------#
    # if given input is a string treat it like it's a file to import
    if type(filelist).__name__ in ['str','string_']:
        # check if file exists
        if not os.path.exists(filelist): raise IOError("Input file not found: '"+filelist+"'")
        dirpath = os.path.dirname(os.path.abspath(filelist))+"/"

        f = open(filelist)
        for file in f:
            file = file.rstrip().split()
            # file = "\ ".join(file)
            file = file[0]
            if not relative_paths:
                bfile = os.path.basename(file)
                file = dirpath+bfile

            if not os.path.exists(file): raise IOError("File doesn't exist: '"+file+"'") # could skip
            else: list_of_files.append(file)

        if len(list_of_files) == 0: raise IOError("No valid files found")
        f.close()
    #-----------------------------------------------#
    # if given input is not a string assume it's a list/array
    else:
        relative_paths = False
        # see if it's array like and can be used to iterate through
        try: list_of_files = list(np.array(filelist,dtype=str))
        except: 
            raise IOError("Input must either be a string giving a file which contains lists of files or array like list of file names")


    #============================================================#
    # now with list_of_files import all the objects
    all_objs = []
    first_file = True
    for file in list_of_files:
        # the code was doing something funny, this is the work around
        if file[0] == "'" and file[-1] == "'": file = file.replace("'","")

        if not os.path.exists(file): raise IOError("not using file '"+file+"' because it doesn't exist") # could skip this
        # the real need for readin is to get the wavelength data
        else:
            new_obj = readin(file)
            if first_file:
                first_file = False
                shape = new_obj.shape

            elif new_obj.shape != shape: raise IOError("the data from file '"+file+"' is not consistent in dimension to the first file") # could skip
            else: all_objs.append(new_obj) 
        
    # each of all_obj is shape  = (1,1,#pts)
    # each is brand new
    if len(all_objs) == 0: raise IOError("NO OBJECTS READ IN")


    #============================================================#
    # now organize the data
    #------------------------------------------#
    # set up a storage for information
    all_wl   = [all_objs[0]._wl[0][0]] # shape (1,#pts)
    all_data = [all_objs[0]._data[0][0]] # shape (1,#pts)

    priv_info   = [deepcopy(all_objs[0].info(0))]
    header_list = [deepcopy(all_objs[0].header)]
    first_shape = deepcopy(all_objs[0].shape)

    # walk through all the objects and store the useful information while also keeping track of sorting information and number of points 
    sortit = [all_wl[0][0]] # first wavelength point of the array
    max_pts = len(all_wl[0])

    for i in range(1,len(all_objs)):
        obj = all_objs[i]

        sortit.append(obj._wl[0][0][0])

        all_wl.append(obj._wl[0][0])
        all_data.append(obj._data[0][0])

        max_pts = max(max_pts,len(obj._wl[0]))

        priv_info.append(deepcopy(obj.info(0)))
        header_list.append(deepcopy(obj.header))


    # sort by the first value of the wavelength array
    sortit = np.array(sortit).argsort()

    #------------------------------------------#
    # initialize the arrays used for output
    ord_num = 0
    new_wl = [] 
    new_data = []
    new_priv_info = {}
    new_hdr_list = []

    #------------------------------------------#
    # walk through the sorted index points
    for j in sortit:
        # check to make sure the length is ok
        if len(all_wl[j]) < max_pts:
            all_wl[j] = np.concatenate((all_wl[j],np.ones(max_pts - len(all_wl[j]))*-.5))
            all_data[j] = np.concatenate((all_data[j],np.ones(max_pts - len(all_data[j]))*-.5))
            if verbose: print "HeadsUp: File #"+str(j)+" has a fewer number of data points, it may not have a wavelength solution"
            
        if first_shape[2] != all_wl[j].shape[0]: raise IOError("NUMBER OF POINTS MUST BE THE SAME IN ORDER TO CONCATENATE THEM")

        new_wl.append(all_wl[j])
        new_data.append(all_data[j])

        # define the private information
        priv_info[j]['id'][1] = (0,ord_num)
        new_priv_info[(0,ord_num)] = priv_info[j]
        ord_num+=1
        
        new_hdr_list.append(header_list[j])

    #============================================================#
    # create output object

    new_obj = eyeSpec_spec(np.array([new_wl],dtype=float),np.array([new_data],dtype=float),header = new_hdr_list[0])
    new_obj._private_info = new_priv_info

    new_obj.header.update('NAXIS',1)
    new_obj.header.update('NAXIS1',len(new_wl[0]))
    new_obj.header.update('NAXIS2',len(new_wl))

    naxis3 = query_fits_header(new_obj.header,'naxis3')
    if naxis3.found: del new_obj.header['naxis3']

    # add this
    new_obj.hdrlist = new_hdr_list

    #============================================================#
    # return new object
    if verbose: print "COMBINED ",len(all_objs)," SPECTRE 1D FILES"
    return new_obj

def combine_objs(obj_list):
    """
    Takes a tuple list of objs and returns one list will all the orders combine

    returns an obj with a single band which combines all the orders for the objects given in the list based on their obj.use_band

    appends the header from the first obj

    """

    obj_list = np.array(obj_list,dtype = object) # all must be of same type
    if obj_list[0].__class__.__name__ != 'eyeSpec_spec': raise TypeError("OBJECTS IN OBJECT LIST MUST BE OF TYPE eyeSpec_spec")

    if obj_list[0].shape[0] != 1: raise ValueError("CAN ONLY COMBINE OBJECTS WITH A SINGLE BAND")

    obj_list[0]._check_band_num() # !! check to make sure obj.use_band is good

    all_data = deepcopy(obj_list[0]._data[obj_list[0].use_band])
    all_wl = deepcopy(obj_list[0]._wl[obj_list[0].use_band])

    obj_out = obj_list[0].copy()
    
    obj_out.combo_headers = [deepcopy(obj_out.header)]

    new_priv_data = {}

    # now add orders to obj_out
    k = deepcopy(obj_out.shape[1])
    for i in range(1,len(obj_list)):
        obj = obj_list[i]

        obj._check_band_num() # check to make sure obj.use_band is good
        if obj.shape[0] != 1: raise ValueError("CAN ONLY COMBINE OBJECTS WITH A SINGLE BAND")

        # add all headers to objec
        obj_out.combo_headers.append(deepcopy(obj_out.header))

        # adopt info and notes for the each object
        for j in range(obj.shape[1]): # for orders
            obj_out._private_info[(0,k)] = deepcopy(obj.info(j))
            obj_out._notes[(0,k)] = deepcopy(obj.notes(j))
            k+=1

        # !! HEY,HEY what happens when you have two spectre files but one has had ends chopped, probably need to do a run through and find the largest lenght and then adjust all the others to match based on something
        if all_data.shape[-1] != obj._data[obj.use_band].shape[-1]: # check # pts
            print "ERROR: OBJECTS NOT OF SAME SHAPE",all_data.shape,obj.shape
            return


        obj_out._wl[0] = np.vstack((obj_out._wl[0],obj._wl[obj.use_band]))
        obj_out._data[0] = np.vstack((obj_out._data[0],obj._data[obj.use_band]))
        obj_out._inv_var[0] = np.vstack((obj_out._inv_var[0],obj._inv_var[obj.use_band]))
        
#        all_wl = np.vstack((all_wl,obj._wl[obj.use_band]))
#        all_data = np.vstack((all_data,obj._data[obj.use_band]))


    all_wl = np.array([all_wl])
    all_data = np.array([all_data])

    obj_out.header['NAXIS2'] = obj_out.shape[1]

    obj_out._initialization_time = time.ctime()
    obj_out._obj_name = None

    obj_out.orders = {}
    for i in range(obj_out.shape[1]):
        obj_out._orders[i] = i

    return new_obj

def readin_apogee (filename,use_row=1):
    """ 
    This takes the pipeline reduced fits from the APOGEE. This should contain several header units each with several image extensions.

    INPUT:
    filename : (fits) APOGEE pipeline reduced data with a 0 header unit similar to the below
    use_order : (int) APOGEE refers to these as rows, default is row1 ("combined spectrum with individual pixel weighting")

    OUTPUT:
    data,tell = (eyeSpec_spec, eyeSpec_spec) returns a two element array the first being the data spectrum (HDU1) the seconde being the telluric spectrum (HDU6) the errors have been converted to inverse varience and incorporated into the spectrum class.


    =================================================================
    Example header 0 header unit:
    
    HISTORY APSTAR: The data are in separate extensions:                      
    HISTORY APSTAR:  HDU0 = Header only                                       
    HISTORY APSTAR:  All image extensions have:                               
    HISTORY APSTAR:    row 1: combined spectrum with individual pixel weighti 
    HISTORY APSTAR:    row 2: combined spectrum with global weighting         
    HISTORY APSTAR:    row 3-nvisis+2: individual resampled visit spectra     
    HISTORY APSTAR:   unless nvists=1, which only have a single row           
    HISTORY APSTAR:  All spectra shifted to rest (vacuum) wavelength scale    
    HISTORY APSTAR:  HDU1 - Flux (10^-17 ergs/s/cm^2/Ang)                     
    HISTORY APSTAR:  HDU2 - Error (10^-17 ergs/s/cm^2/Ang)                    
    HISTORY APSTAR:  HDU3 - Flag mask (bitwise OR combined)                   
    HISTORY APSTAR:  HDU4 - Sky (10^-17 ergs/s/cm^2/Ang)                      
    HISTORY APSTAR:  HDU5 - Sky Error (10^-17 ergs/s/cm^2/Ang)                
    HISTORY APSTAR:  HDU6 - Telluric                                          
    HISTORY APSTAR:  HDU7 - Telluric Error                                    
    HISTORY APSTAR:  HDU8 - LSF coefficients                                 
    HISTORY APSTAR:  HDU9 - RV and CCF structure

    """
    # this is related to the row1
    # can also give it an oid form (band,order)
    
    # use_order = 0 
    use_order = int(use_row-1)
    hdu_header = 0
    hdu_flux = 1
    hdu_err = 2
    hdu_tell = 6
    hdu_tell_er = 7

    non_std_fits=False
    disp_type='log linear'
    preferred_disp='crval'

    def _get_obj (filename,header,use_order,hdu_data,hdu_error):
        x_data = readin(filename,hdu=hdu_data)
        x_data.set_use_cropped(False)
        wl = x_data.get_wl(use_order)
        data = x_data.get_data(use_order)

        x_err = readin(filename,hdu=hdu_error)
        err = x_err.get_data(use_order)
        var = err**2
        inv_var = var_2_inv_var(var)
    
        return eyeSpec_spec(wl,data,inv_var,header)

    headeru = pyfits.open(filename)
    header = headeru[hdu_header].header
    
    data_out = _get_obj(filename,header,use_order,hdu_flux,hdu_err)
    tell_out = _get_obj(filename,header,use_order,hdu_tell,hdu_tell_er)

    data_out.header = header
    tell_out.header = header

    return data_out,tell_out

def readin_makee (filename,varience_filename=None,output_list=False,verbose=False):

    """ 
    Knows how to identify the KOA MAKEE file structure which ships with extracted data
    and apply the eyeSpec function readin to the important directories to obtain a coherent 
    spectrum object from the files


    INPUTS:
    filename : give an individual filename for the star or give the top level Star directory from MAKEE. 
               It will go from TOP_LEVEL/extracted/makee/ and use directories ccd1/ etc to find the appropriate files

    output_list : if it finds multiple chips of data it will return as a list and not a combined object


    """
    non_std_fits=False
    disp_type='default'
    preferred_disp='makee'
    

    def obj_var_2_inv_var (obj,fill=1e50):
        var = deepcopy(obj._data)

        # !! how to treat non values, i.e. negative values
        zeros = (var<=0)
        bad = (var>=fill/2.0)
        infs = (var == np.inf)

        var[zeros] = 1.0/fill
        inv_var = 1.0/var

        # set points which are very large to the fill
        inv_var[zeros] = fill
        # set points which are almost zero to zero
        inv_var[bad] = 0.0
        inv_var[infs] = 0.0

        obj._inv_var = deepcopy(inv_var)
        return inv_var


    filename = str(filename)
    if not os.path.exists(filename): raise ValueError("the given path does not exist")


    is_file = True
    objs = {}
    inv_vars = {}

    if os.path.isdir(filename):
        is_file = False
        if filename[-1:] != '/': filename += "/"
        # !! could make it smarter so it would know from anywhere within the TOP_FILE/extracted/makee/ chain
        full_path = filename+'extracted/makee/'
        if not os.path.exists(full_path): raise ValueError("Must have extracted files:"+full_path)
        
        ccds = os.listdir(full_path)
        

        for ccdir in ccds:
            if not os.path.isdir(full_path+ccdir): continue
            if not os.path.exists(full_path+ccdir+'/fits/'): continue

            if ccdir in objs.keys():
                print "Directory was already incorporated:"+ccdir
                continue

            fitspath = full_path+ccdir+'/fits/'
            fitsfiles = os.listdir(fitspath)

            print ""
            print "="*20+format("GETTING DATA FROM DIRECTORY:"+ccdir,'^40')+"="*20
            for file in fitsfiles:
                fname = fitspath+file

                # !! if file.find("_*.fits")
                # !! I could add in stuff which would go as another band for the current Flux

                if file.find("_Flux.fits") != -1: 
                    print "flux file:"+ccdir+'/fits/'+file
                    objs[ccdir] = readin(fname,preferred_disp=preferred_disp,disp_type=disp_type,non_std_fits=non_std_fits,verbose=verbose)

                elif file.find("_Var.fits") != -1:
                    print "variance file:"+ccdir+'/fits/'+file
                    tmp_obj = readin(fname,preferred_disp=preferred_disp,disp_type=disp_type,non_std_fits=non_std_fits,verbose=verbose)
                    inv_vars[ccdir] = obj_var_2_inv_var(tmp_obj)


    else:
        print "Reading in flux file:"+fname
        objs['file'] = readin(filename,preferred_disp=preferred_disp,disp_type=disp_type,non_std_fits=non_std_fits,verbose=verbose)

        if varience_filename is not None:
            inv_var = readin(varience_filename,preferred_disp=preferred_disp,disp_type=disp_type,non_std_fits=non_std_fits,verbose=verbose)
            inv_vars['file'] = obj_var_2_inv_var(inv_var)


    num_objs = 0

    OUTPUT_list = []
    OUT_header = []
    OUT_wl = []
    OUT_data = []
    OUT_inv_var = []


    for key in objs.keys():
        obj1 = objs[key]
        inv_var1 = inv_vars[key]
        # !! note I'm masking out the inf values
        mask = (inv_var1 == np.inf)
        inv_var1[mask] = 0.0

        if obj1._inv_var.shape != inv_var1.shape:
            print "HeadsUp: object and inverse variance shape are not the same"
        else: obj1._inv_var = deepcopy(inv_var1)
            
        num_objs += 1

        if output_list:
            OUTPUT_list.append(obj1)
        else:
            OUT_header.append(obj1.header)
            if obj1._wl.shape[0] > 1: print "HeadsUp: Multiple bands detected, only using the first"
            OUT_wl.append(obj1._wl[0])
            OUT_data.append(obj1._data[0])
            OUT_inv_var.append(obj1._inv_var[0])

    if output_list:
        if num_objs == 1: return OUTPUT_list[0]
        else: return OUTPUT_list        
    else:
        print ""
        print "="*30+format("COMBINING",'^20')+"="*30
        wl = np.concatenate(OUT_wl)
        data = np.concatenate(OUT_data)
        inv_var = np.concatenate(OUT_inv_var)

        obj = eyeSpec_spec(wl,data,inv_var,OUT_header[0])
        obj.hdrlist = OUT_header
        obj.filename = fitspath
        obj.edit.sort_orders()
        return obj

def readin_hst (filename,get_data=False):
    """
    This function is designed to read in Hubble Space Telescope Archive x1d data
    
    """
    format_error = "Unexpected HST format for fits file. Please use the X1D"
    
    try: hdulist = pyfits.open(filename)
    except: raise ValueError(format_error)
    
    if len(hdulist) != 2: raise ValueError(format_error)

    hdu = hdulist[1] # the data table
    
    wl = hdu.data['WAVELENGTH']
    flux = hdu.data['FLUX']
    var = hdu.data['ERROR']**2
    inv_var = var_2_inv_var(var)
    
    if get_data: return (wl,flux,inv_var)
     
    spec_obj = eyeSpec_spec(wl,flux,inv_var,hdu.header)
    
    # set up private information
    
    spec_obj.filename = filename
    spec_obj._private_info['filename'] = filename
            
    if len(hdulist) > 1: spec_obj.hdrlist = [h.header for h in hdulist]
        
    return spec_obj

    
    
    
    
    
    
    
