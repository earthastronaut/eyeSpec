if __name__ != '__main__':
    import pdb #@UnusedImport
    from base_functions import reduce_output_shape
    from dependencies import np, sys, time, deepcopy, pyfits, scipy
    
 
class Verbose:
    """
    A class to handle reporting.  Set the fileo attribute to any file
    instance to handle the output.  Default is sys.stdout
    
    From matplotlib version '1.1.0'
    """

    levels = ('silent','helpful','debug','debug-annoying')
    vald = dict( [(level, i) for i,level in enumerate(levels)])

    # parse the verbosity from the command line; flags look like
    # --verbose-silent or --verbose-helpful
    _commandLineVerbose = None

    for arg in sys.argv[1:]:
        if not arg.startswith('--verbose-'): continue
        _commandLineVerbose = arg[10:]

    def __init__(self):
        self.set_level('silent')
        self.fileo = sys.stdout

    def set_level(self, level):
        'set the verbosity to one of the Verbose.levels strings'

        if self._commandLineVerbose is not None:
            level = self._commandLineVerbose
        if level not in self.levels:
            raise ValueError('Illegal verbose string "%s".  Legal values are %s'%(level, self.levels))
        self.level = level

    def set_fileo(self, fname):
        std = {
            'sys.stdout': sys.stdout,
            'sys.stderr': sys.stderr,
        }
        if fname in std:
            self.fileo = std[fname]
        else:
            try:
                fileo = file(fname, 'w')
            except IOError:
                raise ValueError('Verbose object could not open log file "%s" for writing.\nCheck your matplotlibrc verbose.fileo setting'%fname)
            else:
                self.fileo = fileo

    def report(self, s, level='helpful'):
        """
        print message s to self.fileo if self.level>=level.  Return
        value indicates whether a message was issued

        """
        if self.ge(level):
            print >>self.fileo, s
            return True
        return False

    def wrap(self, fmt, func, level='helpful', always=True):
        """
        return a callable function that wraps func and reports it
        output through the verbose handler if current verbosity level
        is higher than level

        if always is True, the report will occur on every function
        call; otherwise only on the first time the function is called
        """
        assert callable(func)
        def wrapper(*args, **kwargs):
            ret = func(*args, **kwargs)

            if (always or not wrapper._spoke):
                spoke = self.report(fmt%ret, level)
                if not wrapper._spoke: wrapper._spoke = spoke
            return ret
        wrapper._spoke = False
        wrapper.__doc__ = func.__doc__
        return wrapper

    def ge(self, level):
        'return true if self.level is >= level'
        return self.vald[self.level]>=self.vald[level]
      
class Timer:
    def __init__ (self,dtime=1):
        self._stime = time.time()
        self._ptime = time.time()
        self._dtime = dtime # Seconds
        self._record = []
        
    def check (self,dtime='default',reset=True):
        """
        Returns if current time is larger than the last recorded time
        """
        if dtime is 'default': dtime = self._dtime
        else: dtime = float(dtime)
        truth = (time.time() - self._ptime > dtime)
        if truth and reset: self.reset()
        return truth
        
    def reset (self):
        self._ptime = time.time()
        
    def record (self,info=''):
        self._record.append([info,time.time(),self._ptime(),self._dtime])
    
    def get_record (self):
        return self._record
            
    def get_start_time (self):
        return self._stime
    
    def get_prev_time (self):
        return self._ptime

    def set_delta_time (self,dtime):
        self._dtime = float(dtime)
    
    def get_delta_time (self):
        return self._dtime
    
#### class for dealing with header information ############

class query_fits_header:
    """ 
    This class looks through a pyfits header object and formats queries into useable formats

    EXAMPLE:
    >>> cval = __gethdr__(pyfits_header,'CVAL',no_val=0)
    >>> # if found
    >>> cval.found == True
    >>> cval.val # the value for 'CVAL from the header
    >>> # if not found
    >>> cval.found == False
    >>> cval.val == no_val

    """
    # this class records the header information if it's found and provides a simple way to check 
    def __init__ (self,prihdr,keyword, noval = 0, verbose=False):
        self.keyword = keyword
        if keyword in prihdr.keys():
            self.val = prihdr[keyword]
            self.found = True
        else:
            verbose = bool(verbose)
            self.val = noval
            if verbose: print "keyword:",keyword,"not found"
            self.found = False

    def __repr__ (self):
        return self.val

class EditSpectrumClassData:
    """
    This defines methods used to edit the data in place
    """
    def __init__ (self,parent):
        # !! parent must be of certain class
        self._parent = parent
        # protect data from changes?
        self._protect_data = False
        self.history = [] # I should edit this in


    def __repr__ (self):
        return "Methods used to edit eyeSpec_spec object in place"

    #========================================================================#
    def set_protection (self,truth):
        """
        Protection True will keep any edits from being able to be made
        """
        if type(truth).__name__ != 'bool': raise TypeError("Valuse for truth must be of type boolean")
        self._protect_data = truth

    #========================================================================#
    def get_protection (self):
        """
        Return whether the data can be edited
        """
        return deepcopy(self._protect_data)

    #========================================================================#
    def translate (self,oid=None, # !! change this to *oid ??
                   add_to_data = 0.,
                   add_to_wl = 0.):
        """
        Translate the data by adding to the data and wavelength

        """
        # !! should I only be translating/scaling to points which aren't cropped out?? probably not

        if self._protect_data: raise StandardError("Data Protected, to unprotect run obj.edit.set_protection(False)")

        # make sure values are appropriate
        try: add_to_data = float(add_to_data)
        except: raise TypeError("Value for add_to_data must be convertable to a float")

        try: add_to_wl = float(add_to_wl)
        except: raise TypeError("Value for add_to_wl must be convertable to a float")


        # if oid is not given apply to entire band
        # !! won't break if you change bands to be a list
        if oid is None:
            self._parent._check_band_num()
            self._parent._wl[self._parent._band_index] += add_to_wl
            self._parent._data[self._parent._band_index] += add_to_data

            for j in range(self._parent.shape[1]):
                poid = (self._parent._band_index,j)
                self._parent._private_info[poid]['disp'][0][0] += add_to_wl
                self._parent._private_info[poid]['translate'].append([add_to_wl,add_to_data])
            return
        # make sure oid is appropriate
        oid = self._parent._convert_oid(oid)
        if oid is None: raise TypeError("Please give appropriate oid value")


        self._parent._wl[oid[0]][oid[1]] += add_to_wl
        self._parent._data[oid[0]][oid[1]] += add_to_data

        self._parent._private_info[oid]['disp'][0][0] += add_to_wl
        self._parent._private_info[poid]['translate'].append([add_to_wl,add_to_data])

        self.history.append([time.ctime(),'translated the data :'+str((add_to_wl,add_to_data))])

    #========================================================================#
    def scale (self,oid=None,
                   mult_to_data = 1.0,
                   mult_to_wl = 1.0):
        """
        Scale the data by a constant, scaling the data also causes the inverse varience to be adjusted
        """
        if self._protect_data: raise StandardError("Data Protected, to unprotect run obj.edit.set_protection(False)")

        # make sure values are appropriate
        try: add_to_data = float(add_to_data)
        except: raise TypeError("Value for add_to_data must be convertable to a float")

        try: add_to_wl = float(add_to_wl)
        except: raise TypeError("Value for add_to_wl must be convertable to a float")

        self.history.append([time.ctime(),'scale the data :'+str((mult_to_wl,mult_to_data))])

        # if oid is not given apply to entire band
        # !! won't break if you change bands to be a list
        if oid is None:
            self._parent._check_band_num()
            self._parent._wl[self._parent._band_index] *= mult_to_wl
            self._parent._data[self._parent._band_index] *= mult_to_data

            for j in range(self._parent.shape[1]):
                poid = (self._parent._band_index,j)
                self._parent._private_info[poid]['disp'][0][0] *= mult_to_wl
                self._parent._private_info[oid]['scale'].append([mult_to_wl,mult_to_data])

            return
        # make sure oid is appropriate
        oid = self._parent._convert_oid(oid)
        if oid is None: raise TypeError("Please give appropriate oid value")


        self._parent._wl[oid[0]][oid[1]] *= mult_to_wl
        self._parent._data[oid[0]][oid[1]] *= mult_to_data
        self._parent._inv_var[oid[0]][oid[1]] *= 1./(mult_to_data)**2.

        self._parent._private_info[oid]['disp'][0] *= mult_to_wl
        self._parent._private_info[oid]['scale'].append([mult_to_wl,mult_to_data])

    #========================================================================#
    def apply_function (self,oid=None,
                        fxn_to_data = None,
                        fxn_to_wl = None,
                        fxn_to_inv_var = None):
        """
        Apply a 1-to-1 function to the wavelength, data, or inverse varience
        
        Note that applying functions to the data should change the inverse varience but that must be user provided

        """
        if self._protect_data: raise StandardError("Data Protected, to unprotect run obj.edit.set_protection(False)")

        def _check_function (fxn):
            try: Y = fxn(np.arange(10))
            except: raise ValueError("Function must convert numpy array")
            if Y.shape != (10,): raise ValueError("Function must be 1-to-1")

        _check_function(fxn_to_data)
        _check_function(fxn_to_wl)
        

        hist_info = []
        if fxn_to_data is not None: hist_info.append('data')
        if fxn_to_wl is not None: hist_info.append('wl')
        if fxn_to_inv_var is not None: hist_info.append('inv_var')
        if len(hist_info) > 0:
            self.history.append([time.ctime(),'Apply function to: '+",".join(hist_info)])
        
        print "currently not implemented"
        #if type(fxn_to_data).__name__ != 'function':
        # try: single_val = fxn_to_data(0)
        
        # fxn_to_data = np.vectorize(fxn_to_data)

        # !! make sure you are also applying to inv_var correctly

        return

#         if self._protect_data: raise StandardError("Data Protected, to unprotect run obj.edit.set_protection(False)")

#         # make sure values are appropriate
#         try: add_to_data = float(add_to_data)
#         except: raise TypeError("Value for add_to_data must be convertable to a float")

#         try: add_to_wl = float(add_to_wl)
#         except: raise TypeError("Value for add_to_wl must be convertable to a float")


#         # if oid is not given apply to entire band
#         # !! won't break if you change bands to be a list
#         if oid is None:
#             self._parent._check_band_num()
#             self._parent._wl[self._parent._band_index] *= mult_to_wl
#             self._parent._data[self._parent._band_index] *= mult_to_data

#             for j in range(self._parent.shape[1]):
#                 poid = (self._parent._band_index,j)
#                 self._parent._private_info[poid]['disp'][0][0] *= mult_to_wl
#                 self._parent._private_info[oid]['scale'].append([mult_to_wl,mult_to_data])

#             return
#         # make sure oid is appropriate
#         oid = self._parent._convert_oid(oid)
#         if oid is None: raise TypeError("Please give appropriate oid value")


#         self._parent._wl[oid[0]][oid[1]] *= mult_to_wl
#         self._parent._data[oid[0]][oid[1]] *= mult_to_data

#         self._parent._private_info[oid]['disp'][0] *= mult_to_wl
#         self._parent._private_info[oid]['scale'].append([mult_to_wl,mult_to_data])



    #========================================================================#
    def crop (self,order=None,
              inverse = False,
              replace_data = None,
              # explicit versions for cropping              
              crop_mask = None,
              index_bounds = None,
              wavelength_bounds = None,
              data_bounds = None,
              box_bounds = None,
              # short cut versions for cropping
              cm = None,
              ib = None,
              wb = None,
              db = None,
              bb = None):
        """
        oid must indicate an order
        oid = 10, (10), or (6,10)
           
        keep what's within the bounds, like cropping a photo. If you want to crop out certain points using box_bounds then select inverse=True

        INPUTS:
        ===================   ====================================================
        keyword               (type) Description
        ===================   ====================================================
        Inverse               (bool) reverse what points are kept

        replace_data          (None or Float) if None then it will use the nearby remaining points to interpolate new points, if a value is given then it will replace using that value
                              Note: regardless of the choice the inverse varience for these points is set to zero. So in co-additions and reading out they will be excluded

        cm = crop_mask        (ndarray, bool) must be same length as the obj.shape[2] all True entries will be cropped out
        ib = index_bounds     (array) ['int_left','int_right'] points from that end
        wb = wavelength_bounds  (array) ['float_left','float_right']
        db = data_bounds     (array) ['float_up','float_down']
        bb = box_bounds      (array) ['float_left','float_right','float_down','float_up'])
        ==================   ====================================================

        """
        if self._protect_data: raise StandardError("Data Protected, to unprotect run obj.edit.set_protection(False)")

        if order is None: print "To make a crop selection you must indicate which order" ## !!? error and exit?

        #-----------------------------------------------------------#
        # make it easier to input bounds, though it will prioritize the longer name and will ignore the shorter one if given
        if cm is not None and crop_mask is None: crop_mask = cm
        if ib is not None and index_bounds is None: index_bounds = ib
        if wb is not None and wavelength_bounds is None: wavelength_bounds = wb
        if db is not None and data_bounds is None: data_bounds = db
        if bb is not None and box_bounds is None: box_bounds = bb

        #-----------------------------------------------------------#
        # initialize variables
        # these are the possible inputs, in prioritized order
        kwargs = [index_bounds,wavelength_bounds,data_bounds,box_bounds,crop_mask]
        kwargs_names = np.array(['index_bounds','wavelength_bounds','data_bounds','box_bounds','crop_mask'])
        gave_input = (np.ones(len(kwargs)) < 0.0)
        
        # check the value used to replace data points
        if replace_data is not None:
            try: replace_data = float(replace_data)
            except: raise TypeError("Variable 'replace_data' must be convertable to a float")

        # check the order
        try: oid = int(order)
        except: raise TypeError("Variable 'order' must be an integer")

        # make sure oid (a.k.a. order) is appropriate
        oid = self._parent._convert_oid(oid)
        if oid is None: raise TypeError("Please give appropriate oid value")

        # !! could update this once i have a set_* attribute or something
        # get data to access
        tmp_use_cropped = deepcopy(self._parent.get_use_cropped())
        self._parent.set_use_cropped(False)
        new_wl = self._parent.get_wl(oid)
        new_data = self._parent.get_data(oid)
        new_inv_var = self._parent.get_inv_var(oid)

        #-----------------------------------------------------------#
        # check the inputs, see which have been given values:
        for i in range(len(kwargs)):
            if kwargs[i] is None: continue
            else:
                # everything given must be convertable to some type of numpy array
                try: np.array(kwargs[i])
                except: raise TypeError("The input for "+kwargs_names[i]+" must be convertable to a ndarray")
                gave_input[i] = True

        num_found = len(np.ones(5)[gave_input])
        # check to make sure multiples weren't given, could sort of delete this check because it won't break it'll just use the first one found in priority
        if num_found > 1: raise ValueError("Please provide ONLY ONE type of bounds, you gave "+len(np.arange(5)[gave_input]))

        # make sure bounds were given at all
        elif num_found == 0:
            print "No bounds given" # !! tell user choices
            return

        #-----------------------------------------------------------#
        # get the mask for the data points in an order

        # take the first bounds type given
        bounds_type = kwargs_names[gave_input][0]
        ind = np.arange(len(kwargs))[gave_input][0]
        given_bounds = deepcopy(kwargs[ind])

        #________________________________
        if bounds_type == 'index_bounds':
            # check array
            try: bounds = np.array(given_bounds,dtype=int)
            except:  raise TypeError("The input for index_bounds must be convertable to a ndarray of dtype int")
            if len(bounds) != 2 and bounds.ndim == 1: raise TypeError("The input for index_bounds must be of length 2")

            # use bounds to set mask
            bounds = [np.min(bounds),np.max(bounds)]
            indicies = np.ones(self._parent.shape[2]) # number of pixels in each order
            indicies[:bounds[0]] = 0
            indicies[-bounds[1]:] = 0
            bounds_mask = (indicies == 0.0)

        #________________________________
        elif bounds_type == 'wavelength_bounds':
            # check array
            try: bounds = np.array(given_bounds,dtype=float)
            except:  raise TypeError("The input for wavelength_bounds must be convertable to a ndarray of dtype float")
            if len(bounds) != 2 and bounds.ndim == 1: raise TypeError("The input for wavelength_bounds must be of length 2")

            # use bounds to set mask
            bounds = [np.min(bounds),np.max(bounds)]
            bounds_mask = np.logical_not((self._parent.get_wl(oid) > bounds[0])*(self._parent.get_wl(oid) < bounds[1]))

        #________________________________
        elif bounds_type == 'data_bounds':
            # check array
            try: bounds = np.array(given_bounds,dtype=float)
            except:  raise TypeError("The input for data_bounds must be convertable to a ndarray of dtype float")
            if len(bounds) != 2 and bounds.ndim == 1: raise TypeError("The input for data_bounds must be of length 2")
            
            # use bounds to set mask            
            bounds = [np.min(bounds),np.max(bounds)]
            bounds_mask = np.logical_not((self._parent.get_data(oid) > bounds[0])*(self._parent.get_data(oid) < bounds[1]))

        #________________________________
        elif bounds_type == 'box_bounds':
            # check array
            try: bounds = np.array(given_bounds,dtype=float)
            except:  raise TypeError("The input for box_bounds must be convertable to a ndarray of dtype float")
            if len(bounds) != 4 and bounds.ndim == 1: raise TypeError("The input for box_bounds must be of length 4")
            
            # use bounds to set mask, this picks all the points within the box
            bbounds = [np.min(bounds[:2]),np.max(bounds[:2])]
            bounds_mask_wl = (self._parent.get_wl(oid) > bbounds[0])*(self._parent.get_wl(oid) < bbounds[1])
            bbounds = [np.min(bounds[2:]),np.max(bounds[2:])]
            bounds_mask_dat = (self._parent.get_data(oid) > bbounds[0])*(self._parent.get_data(oid) < bbounds[1])
            bounds_mask = np.logical_not(bounds_mask_wl*bounds_mask_dat) # where both are true

        #________________________________
        elif bounds_type == 'crop_mask':
            bounds = given_bounds
            # check array
            if bounds.dtype.name != 'bool': raise TypeError("crop_mask must be dtype bool")

            if bounds.shape != new_data.shape: raise TypeError("shape of crop_mask must match data of shape:"+str(new_data.shape))            
            # use mask
            bounds_mask = deepcopy(bounds)

        else: raise RuntimeError("bounds_type should have been one of, please debug: "+", ".join(kwargs_names))

        # use number of points found to check if all the points were selected and interpolating (i.e. replace_data is None) because then you'll need atleast two points for the interpolation to work
        num_found = len(np.arange(len(bounds_mask))[bounds_mask])

        # self explanatory
        if inverse: bounds_mask = np.logical_not(bounds_mask)

        # check if bounds mask has any values
        if not np.any(bounds_mask):
            # !! are these even necessary!!??
            if inverse:  print "Found no points to crop given selection" #!! can't be a 'raise' warning/error because it will crash other things, Perhaps shouldn't have a message? #!! verbose
            else: print "No cropping because all points were selected" 


        self.history.append([time.ctime(),'Cropped the data'])

        #-----------------------------------------------------------#
        # now I have bounds, bounds_mask and bounds_type make the changes
        # !! i'm not sure if I should use masked arrays for _wl,_data,_inv_var or not
        # edit data appropriately
        #new_wl[bounds_mask] = unchanged
        new_inv_var[bounds_mask] = 0.0

        if replace_data is not None: new_data[bounds_mask] = replace_data
        else: # interplate data
            if num_found == 0: # if all points were selected for editing unselect two points to interpolate with
                MEAN = np.mean(new_data)
                uhist = np.histogram(new_data,bins=int(np.sqrt(len(new_data))))
                MODE = uhist[1][1:][uhist[0]==np.max(uhist[0])][0] # where the number in the bin (uhist[0]) is equal to the bin with the maximum number, return the corresponding vaue for the data (uhist[1])
                mean1 = np.abs(new_data - MEAN).argmin()
                mean2 = np.abs(new_data[new_data != new_data[mean1]] - MEAN).argmin()
                # make two mid points opposite for interpolation 
                bounds_mask[mean1] = not bounds_mask[mean1]
                bounds_mask[mean2] = not bounds_mask[mean2]

            # wavelength points inside box =>
            edited_wl = new_wl[bounds_mask]
            # date points outside box =>
            wl_pts_outside = deepcopy(new_wl[np.logical_not(bounds_mask)])
            dat_pts_outside = deepcopy(new_data[np.logical_not(bounds_mask)])
            # interpolate to adjust data points inside box
            if len(wl_pts_outside) != 0:
                edited_data = np.array(scipy.interp(edited_wl,wl_pts_outside,dat_pts_outside),dtype=np.float32)
                # use the mask to set the new data, note that both must be of same dtype
                new_data[bounds_mask] = deepcopy(edited_data)
            else: pass # no change to new_data because no points were found


        # set the new data
        self._parent._wl[oid[0]][oid[1]] = new_wl
        self._parent._data[oid[0]][oid[1]] = new_data
        self._parent._inv_var[oid[0]][oid[1]] = new_inv_var
        
        self._parent.set_use_cropped(tmp_use_cropped)


        # record stuff which happened
        i = 0
        # !! should I have multiple crops in one long list or should I create new oids for each crop??!!!!
        self._parent._private_info[oid]['crop'].append([bounds,bounds_type,'inverse = '+str(inverse),'replace_data ='+str(replace_data)])

#         new_oid = 'crop'
#         while True:
#             print "____"
#             if new_oid in self._parent._private_info[oid].keys():
#                 new_oid = "crop_"+str(i)
#                 i+=1
#             else:
#                 # !! for crop_mask should I add this in when 'bounds' currently is? 
#                self._parent._private_info[oid][new_oid] = [bounds,bounds_type,'inverse = '+str(inverse),'replace_data ='+str(replace_data)] # !! could add bound_mask to this mix
#                break

    #======================================================================#
    def apply_rv (self,rv,all_bands = True):
        """ 
        Will apply a radial velocity v_rad [km/s] to the data via z=v_rad/c and 1/(1+z) to all orders in obj.get_band()

        multi = (2. + rv/c)/(2. - rv/c)
        new_wl = old_wl*multi

        """
        self.history.append([time.ctime(),'Apply radial velocity shift:'+str(rv)])

        c = 299792.4580
        # assuming c = 3e5 gives percent difference of:
        # 0.06922855944562016 %

        rv = float(rv)

        # solve  delta_lambda/lambda = rv/c
        # delta_lambda = wl_new - wl_old
        # lambda in [wl_old, wl_new, wl_mid]
        # wl_mid = (wl_old + wl_new)/2.        

        # if you have del_lambda/wl_new = rv/c
        #multi = 1./(1. - rv/c)

        # if you have del_lambda/wl_old = rv/c
        #multi = (1. + rv/c)

        # if you have del_lambda/wl_mid = rv/c
        multi = (2. + rv/c)/(2. - rv/c)


        if all_bands:
            self._parent._wl *= multi
            for i in range(len(self._parent._wl)):
                for j in range(len(self._parent._wl[i])):
                    self._parent._private_info[(i,j)]['rv'].append(rv)
                    
        else: 
            band = self._parent.get_band()
            self._parent._wl[band] *= multi
            self._parent._private_info[(band,j)]['rv'].append(rv)


    #========================================================================#
    def compute_inverse_varience (self,norm=True):
        """
        This computes the inverse varience for the data array in a straight forward way by taking 1/sqrt(N)
        """
        if self._protect_data: raise StandardError("Data Protected, to unprotect run obj.edit.set_protection(False)")
        print "!! not currently available"
#         self._parent.set_use_cropped(False)

#         self._parend.set_use_cropped('previous')
#         self.history.append([time.ctime(),'Compute inverse varience from 1/sqrt(N)'])


    #========================================================================#
    def sort_orders (self):
        """
        Will sort orders based on the begining wavelength
        """
        # get originals
        wl = deepcopy(self._parent.get_wl())
        data = deepcopy(self._parent.get_data())
        inv_var = deepcopy(self._parent.get_inv_var())
        notes = deepcopy(self._parent._notes)
        info = deepcopy(self._parent._private_info)

        # figure out sorting
        sortit = wl.T[0].argsort()

        # create copies
        wl_sorted = deepcopy(wl)
        data_sorted = deepcopy(data)
        inv_var_sorted = deepcopy(inv_var)
        notes_new = deepcopy(notes)
        info_new = deepcopy(info)

        # get the current band
        band = self._parent.get_band()

        # go through and correctly order orders
        for i in range(len(sortit)):
            # i is the new
            # j is the old
            j = sortit[i]

            notes_new[(band,i)] = deepcopy(notes[(band,j)])
            info_new[(band,i)] = deepcopy(info[(band,j)])

            wl_sorted[i] = deepcopy(wl[j])
            data_sorted[i] = deepcopy(data[j])
            inv_var_sorted[i] = deepcopy(inv_var[j])


        # set the new data
        self._parent._wl[band] = wl_sorted
        self._parent._data[band] = data_sorted
        self._parent._inv_var[band] = inv_var_sorted
        self._parent._notes = notes_new
        self._parent._private_info = info_new

        self.history.append([time.ctime(),'Sort the orders by wavelength'])


    #========================================================================#
    def quick_continuum_normalize (self,poly_order=1):
        poly_order = int(poly_order)

        self._parent._check_band_num()
        band = self._parent._band_index
        
        for i in range(self._parent.shape[1]):
            new_dat = deepcopy(self._parent._data[band][i])
            wl = deepcopy(self._parent._wl[band][i])
            Coeff = np.polyfit(wl,new_dat,poly_order)
            Polie = np.polyval(Coeff,wl)
            
            self._parent._data[band][i] = new_dat/Polie
            self._parent._private_info[(band,i)]['continuum divide'].append(Polie)
            # !! instead of 'scale' I may want a 'normalize'
            self.history.append([time.ctime(),'Applied quick continuum normalization'])


    #========================================================================#
#    def del_order (self,order): # !! this won't work if there are multipule bands unless I make the bands a list...., or delete it out of all bands


    #========================================================================#
    def del_band (self,band=None,safety=True,record=True):
        """ 
        this deletes a band from the data 

        INPUTS:
        band : (int) which band to delete
        safety: (bool) if True it will prompt the user to make sure they want to delete
        record: (bool) if True it will record the bands deleted in obj.info('general')
        """
        # check band is appropriate
        band = self._parent._check_band_num(band=band)

        # add safety statment for user
        if safety:
            while True:
                choice = raw_input("Are you sure you want to delete all orders and info of band "+str(band)+"? ('y',['no'])\n")
                if choice.lower() in ['n','no','']: return
                elif choice.lower() in ['y','yes']: break


        # remove band from the array
        self._parent._wl = np.delete(self._parent._wl,band,0)
        self._parent._data = np.delete(self._parent._data,band,0)
        self._parent._inv_var = np.delete(self._parent._inv_var,band,0)

        # previous private info and notes
        prev_info = deepcopy(self._parent._private_info)
        prev_notes = deepcopy(self._parent._notes)
        
        # remove notes and private_info
        for i in range(self._parent.shape[0]): # new shape for bands
            #print ">>",band,i
            if i == band:
                #print "======== del",i
                for j in range(self._parent.shape[1]): # prev shape for orders
                    del self._parent._private_info[(i,j)]
                    del self._parent._notes[(i,j)]
                continue
            else:
                for j in range(self._parent.shape[1]): # prev shape for orders
                    prev_oid = (i,j)
                    if i < band: new_oid = (i,j) 
                    if i > band: new_oid = (i-1,j)
                    #print ">>>>>",prev_oid,"=>",new_oid
                    self._parent._private_info[new_oid] = prev_info[prev_oid]
                    self._parent._notes[new_oid] = prev_notes[prev_oid]


        # !! do i add in self._parent._private_info['bandids']["BANDIDS"+str(band)] = 'deleted'
        if record: self._parent._private_info['general'].append('deleted band '+str(band))
        self._parent.shape = self._parent._data.shape

        self.history.append([time.ctime(),'Deleted band:'+str(band)])

    #========================================================================#
    def convert_wl_to_pts (self,consequative=True):
        """
        Convert all the wavelenths to points
       
        WARNING: Will probably cause artificial stretching to your data

        consequative : (bool) will progress the pixel values so that each begins following the previous
        """
        progressive_pt = 1
        band = self._parent.get_band()
        
        for i in range(len(self._parent._wl[band])):
            wl = self._parent._wl[band][i]

            if consequative: self._parent._wl[band][i] = np.arange(progressive_pt,len(wl)+progressive_pt)
            else: self._parent._wl[band][i] = np.arange(1,len(wl)+1)
            
            progressive_pt += len(wl)

class eyeSpec_spec:
    """
    -----------------------------------------------------------------------
    eyeSpec_spec : Main spectrum class object used in eyeSpec. This 
       contains the wavelength, data, and inverse varience information
       as well as some other meta-data. Also has several inherent methods

    to create:
    >>> spec = eyeSpec_spec(wl,flux,header)

    INPUTS:
    ==============  ===============================================================
    keyword         (type) Description
    ==============  ===============================================================
    wl              (ndarray) wl.shape = (#bands,#orders,#pixels)  wavelength points 
    data            (ndarray) data.shape = (#bands,#orders,#pixels)  data points
    inverse_varience (ndarray) inv_var.shape = (#bands,#orders,#pixels) has the 
                       inverse varience for each point. If not provided then varience 
                       is set to one.
    header          (pyfits.header.Header) pyfits header, if None then it will 
                    create an empty pyfits header
    ==============  ===============================================================
    NOTE: wl, data, inverse_varience must have same shape but doesn't have to be
         three dimensional could just be shape = (#pixels) for each




    =======================================================================
    METHODS IN BRIEF:
    x.copy()    : create a copy of obj x

    x.hdrlist   : list of all the headers which are incorporated, if applicable
    x.header    : dict, can be used to call by header keywords

    x.get_band()  : int, which band to use
    x.set_band(int) : set which band to use
    x.get_wl()   : accessing wl via inputs (#band,#order) or #order for band x._band_index 
    x.get_data() : accessing data via inputs (#band,#order) or #order for band x._band_index 

    x.notes     : dict, allows the user to add arbitrary notes to any order
    x.notes_keys: possible combinations of (#band,#order)

    x.info      : dict, stores useful information such as the dispersion equation and radial velocity which the user should not be asked to change 
    x.info_keys : possible combinations of (#band,#order)
    
    x.add_band : takes an ndarray and concatenates on data and new info
    x.add_order : takes two ndarrays (wl and data) and concatenates them to the data array

    x.edit : methods for editing the data in place

    =======================================================================   
    METHODS FOR ACCESSING DATA:

    >>> obj._band_index = 6
    >>> print "2-D array of wavlength and data for band 6 order 10:",obj[10]
    
    >>> obj.getdata(10)



    Note: if you accidentally overwrite a method (e.g. x.info = 5) the information is not lost because it's stored under x._private_info
    -----------------------------------------------------------------------
    """

    #======================================================================#
    def __init__ (self,wl,data,inverse_varience=None,header=None):

        #--------------------------------------------#
        self._initialize_time = time.ctime() # !! unique

        # these map which band from fits file => index
        self._bands = {0:0} # !! unique
        self._orders = {0:0} # !! unique

        # this if from the fits file and uses self._bands to convert
        self._band_index = 0 # !! unique

        #--------------------------------------------#
        data = np.asarray(data,dtype=float)
        wl = np.asarray(wl,dtype=float)
        if data.shape != wl.shape: raise ValueError("Wavelength array must have same shape as data not: "+str(wl.shape))
        init_shape = data.shape
        
        # reshape so that dim0 is the bands, dim1 is orders, dim2 is the data
        if   data.ndim == 1: newshape = (1,1,len(data))
        elif data.ndim == 2: newshape = (1,) + data.shape
        elif data.ndim == 3: newshape = data.shape
        else: raise ValueError("Data must be less than three dimensions")
        self._data = data.reshape(newshape)   
        self._wl = wl.reshape(newshape)   
        self.shape = newshape 
        self.num_orders = newshape[1]

        # inverse_var for the data
        if inverse_varience is None: self._inv_var = np.ones(self.shape)
        else:
            inv_var = np.asarray(inverse_varience,dtype=float)
            if inv_var.shape != init_shape: raise ValueError("Inverse variance must have same shape as data not: "+str(inv_var.shape))
            self._inv_var = inv_var.reshape(self.shape)

        #--------------------------------------------#
        self._notes = {} # !! unique
        self._private_info = {} # !! unique
        self._reset_notes_and_info()

        #--------------------------------------------#
        if repr(type(header)) != "<class 'pyfits.header.Header'>": print "HeadsUp: creating original header because input header not of type <class 'pyfits.header.Header'> != :"+repr(type(x.header))
        else:
            new_header = pyfits.PrimaryHDU()
            hdulist = pyfits.HDUList([new_header])
            header = hdulist[0].header
            del new_header,hdulist
            header.update('NAXIS1',self.shape[2])
            header.update('NAXIS ',1)
            if self.shape[1] != 1:
                header.update('NAXIS2',self.shape[1])
                header.update('NAXIS ',2)
            if self.shape[0] != 1:
                header.update('NAXIS3',self.shape[0])
                header.update('NAXIS ',3)
                
        self.header = header

        #--------------------------------------------#
        self.filename = 'no_associated_file'
#        try:
#            self.filename = os.path.abspath(hdulist.filename())
#            base = os.path.basename(hdulist.filename()) 
#        except:
#            self.filename = 'no_associated_file'
#            base = '-------' # if hdulist has no associated file

        self._obj_name = None

        self.edit = EditSpectrumClassData(self)


        self._use_cropped = True
        self._prev_use_cropped = True

        self._safe_attributes = ['_safe_attributes'
                                 '_data','_wl','_inv_var',
                                 '_notes','notes','notes_keys',
                                 '_private_info','info','info_keys',
                                 '_band_index','_obj_name','_initialize_time']  
                            
    #======================================================================#
    def __repr__ (self): 
        return "eyeSpec_spec : "+str(self._data.shape)

    #======================================================================#
    def __str__ (self):
        """ print information when directly printed """
        if 'filename' not in self._private_info.keys(): fname = 'unknown'
        else: fname = self._private_info['filename']
        return "".join([format("eyeSpec_spec:",'12'),"(#BANDS,#ORDERS,#PIXELS) = ",
                         str(self._data.shape),
                         "\n",
                         format("FILE NAME:",'12'),fname,
                         "\n",
                         format("CREATED:",'12'),self._initialize_time])

    #======================================================================#
    def __len__ (self):
        """ return the number of bands """
        return deepcopy(self.shape[0])

    #======================================================================#
    def __iter__ (self):
        """
        Will iterate over all orders in band obj.get_band() and return a 2-D array with the first entry as wavelength and second as data

        >>> spec.set_band(6)
        >>> for order in obj:
        >>>     print "wavelength array",order[0]
        >>>     print "data array",order[1]
        >>>     print "inv_var array",order[2]

        """
        orders = []
        for i in xrange(self.num_orders): orders.append(np.vstack((self.get_wl(i),self.get_data(i),self.get_inv_var(i))))
        return iter(orders)


    #======================================================================#
#     def __setattr__ (self,name,val):
#        #if name in dir(self): return
#         print name,val
#         #stop,here=0
#         setattr(self,name,val)

    #======================================================================#
#     def __delattr__ (self,name):
#        if name in self._safe_attributes:
#            print name,"cannot be deleted"
#            return
#        delattr(self,name)

    #======================================================================#
    def _reset_notes_and_info (self,oid=None,skip_notes=False):
        """ this wipes the current information clean """
        if oid is None:
            self._private_info['filename'] = 'unknown'
            self._private_info['general'] = []
            bands = range(len(self._data))
            orders = range(len(self._data[0]))
        else:
            oid = self._convert_oid(oid)
            bands = [oid[0]]
            if oid[1] is None: orders = range(len(self._data[0]))
            else: orders = [oid[1]]
        # add information for each band and order, doesn't add any time to creating

        for i in bands:
            for j in orders:
                if not skip_notes:
                    self._notes[(i,j)] = {}
                    self._notes[(i,j)]['id']=['(band,order)',(i,j)]

                self._private_info[(i,j)] = {}

                self._private_info[(i,j)]['id']=['(band,order)',(i,j)]
                self._private_info[(i,j)]['rv'] = [0]

                self._private_info[(i,j)]['translate'] = []
                self._private_info[(i,j)]['scale'] = []
                self._private_info[(i,j)]['continuum divide'] = []
                self._private_info[(i,j)]['crop'] = []

                self._private_info[(i,j)]['disp'] = [np.zeros(8),'none']

    #======================================================================#
    def _convert_oid (self,*oid):
        """ 
        This takes oid arguments and converts them to a tuple (#band,#order)

        _convert_oid() => (_band_index,0)
        _convert_oid(8) => (_band_index,8)
        _convert_oid(2,3) => (2,3)

        note: it checks to make sure that the band and order are within the correct bounds by using _check_band_num and _check_ord_num

        """
        self._check_band_num() # this checks to make sure that the band is also an integer
        #----------------------------------------#
        # check to see if argument(s) is a tuple
        if len(oid) == 1:
            if oid[0] is None: return None
            elif type(oid[0]) in [tuple,list,np.ndarray]: 
                oid = np.asarray(oid,dtype=int)
                # !! I guess I could take multple tuples to have multiple orders. e.g. (0,1),(0,3),(0,6)
                 
                oid = (oid[0][0],oid[0][1])
            elif type(oid[0]) == float:
                print "!! have it return value"
            
        #----------------------------------------#
        # take the arguments and convert to a tuple (band,order)
        if   len(oid)==0: oid = (self._band_index,None)
        elif len(oid)==1: oid = (self._band_index,oid[0])
        elif len(oid)==2:
            band = self._check_band_num(band=oid[0])
            oid = (band,oid[1])
        else:
            print "HEADS-UP: THIS ONLY USES THE FIRST TWO VALUES" # !! error?
            band = self._check_band_num(band=oid[0])
            oid = (band,oid[1])
            # if there are three maybe take three, (band,order,pixel)


        


        if oid[1] is not None: self._check_order_num(oid[1])
        return oid # tuple with (band,oid)

        # !! next step is to have output oid = (band,order,points)
        # then have things:
        # if oid[3] == int => return value for that index
        # if oid[3] == float => interpolate

        # !! I would need to be smart on the read in so that if input_oid[0] = order and input_oid[1] = this_thing_above I could disentagle that and create (band,order,points)
        # I need to think how this works for objects which want interpolation and those which don't

    #======================================================================#
    def _check_band_num (self,band=None):
        """ checks if current is an appropriate band for the data"""
        # this checks so that you can have specific bands
        #if self._band_index in self._bands.oids():
         
        # !! currently just map self._band_index  to the actual index of the band but I could use _bands dictionary to later have it convert. i.e. there are 7 bands, but you only choose to import 3,5, and 7  then you could have self._band_index = 5 but that corresponds to index 1 <= self._bands[self._band_index(=5)] 
   
        if band is None: band_to_check = int(self._band_index)
        else: band_to_check = int(band) 
        
        if (band_to_check >= self.shape[0]) or (band_to_check < -1*self.shape[0]): raise ValueError("BAND MUST HAVE VALUE < "+str(self.shape[0]))
#        if band_to_check not in self._bands.keys():
#            print "BAND USE MUST BE IN ",self._bands.keys()
#            band_to_check = self._bands[0]

        return band_to_check

    #======================================================================#
    def _check_order_num (self,order_num):
        """
        This checks to make sure the oid is appropriate value
        Should follow a self._check_band_num
        """
        order_num =  int(order_num)        
        if (order_num >= self.shape[1]) or (order_num < -1*self.shape[1]): raise ValueError("GIVEN ORDER IS OUTSIDE OF RANGE ["+str(-1*self.shape[1])+","+str(self.shape[1]-1)+"]") # !! maybe this shouldn't be an error?

    #======================================================================#
    def _find_cropped_bounds (self,*oid):
        """
        Uses the inverse varience to find how many zeros are at the begining and end 
        then returns a two values which give the index bounds to use

        """
        oid = self._convert_oid(*oid)       
        inv_var = self._inv_var[oid[0]][oid[1]]
        if not self._use_cropped: return 0, len(inv_var)
        ind = np.argwhere(inv_var).reshape(-1)
        i = np.min(ind)
        j = np.max(ind)+1
        return i,j
    
    #======================================================================#
    def info (self,*oid):
        """
        Dictionary of information which is unique for each order in each band
        Access is restricted to viewing 

        >>> 
        >>> print "dictionary info for band 6 and order 10", obj.info(6,10)
        >>>
        >>> obj.set_band(6)
        >>> print "dictionary info for band 6 and order 10", obj.info(10)
        >>> print "will also work with", obj.info(10)
        >>>

        """
        if len(oid) == 0: return deepcopy(self._private_info)
        if type(oid[0]) in [np.str,np.str_]:
            val = oid[0].lower()
            if   val == 'filename': return (deepcopy(self._private_info['filename']))
            elif val in ['general','gen']: return (deepcopy(self._private_info['general']))  
            elif val == 'bandids':
                if 'bandids' in self._private_info.keys(): return(deepcopy(self._private_info['bandids']))
                else: raise KeyError("bandids has not been defined") # !! this may annoy someone later
            return
        
        oid = self._convert_oid(*oid)
        # tuple (#bands,#orders) must match obj._private_info.keys()
        if oid != None: return deepcopy(self._private_info[oid])
        else: return None

    #======================================================================#
    def notes (self,*oid):
        """
        Dictionary of information which is unique for each order in each band
        Access is direct to the dictionary containing notes for all orders

        >>> 
        >>> print "dictionary info for band 6 and order 10", obj.notes(6,10)
        >>>
        >>> obj._band_index = 6
        >>> print "dictionary info for band 6 and order 10", obj.notes(10)
        >>> print "will also work with", obj.notes(10)
        >>>

        """
        if len(oid) == 0: return deepcopy(self._notes)
        oid = self._convert_oid(*oid)
        # tuple (#bands,#orders) must match obj._private_info.keys()
        if oid != None: return self._notes[oid]
        else: return None

    #======================================================================#
    def notes_oids (self):
        """=== (band,order) ==="""
        return sorted(self._notes.keys())

    #======================================================================#
    def info_oids (self):
        """=== (band,order) ==="""
        return deepcopy(sorted(self._private_info.keys()))

    #======================================================================#
    def oids (self):
        """=== (band,order) ==="""
        all_oids = sorted(self._private_info.keys())
        r_oids = []
        for oid in all_oids: 
            if type(oid) == tuple: r_oids.append(oid)
        return deepcopy(sorted(r_oids))

    #======================================================================#
    def get_band (self):
        """ this returns the value of the current band index which is used for indexing the data """
        return deepcopy(self._band_index)

    #======================================================================#
    def get_use_cropped (self):
        """
        This tells whether to crop off points at the beginning and end of an order
        whose's inverse varience is 0.
        If True then get_wl and get_data will only return those valid points
        """
        return deepcopy(self._use_cropped)

    #======================================================================#

    def _get_arr (self,arr,*oid):
        if len(oid) == 0:
            self._check_band_num()
            oarr = arr[self._band_index].copy()
            return reduce_output_shape(oarr)
        
        oid = self._convert_oid(*oid)
        # tuple (#bands,#orders) must match obj._private_info.keys()
        if oid is None: return None

        oarr = arr[oid[0]][oid[1]].copy()

        if self._use_cropped:
            i,j = self._find_cropped_bounds(oid)
            return oarr[i:j]
        return oarr
        
    def get_wl (self,*oid):
        """
        Access the fits file wavelength using the value stored in obj._band_index and the oid

        >>> obj.set_band(6)
        >>> print "access data for band 6 in order 10 =", obj.get_wl(10)
        >>> print "access is the same:", obj.get_wl(6,10)

        """
        return self._get_arr(self._wl,*oid)
    
    #======================================================================#
    def get_data (self,*oid):
        """
        Access the fits file data using the value stored in obj._band_index and the oid

        >>> obj.set_band(6)
        >>> print "access data for band 6 in order 10 =", obj.getdata(10)

        """
        # !! it would be cool that if you gave this a float or a numpy array it would use that to interpolate the data and return values for that
        return self._get_arr(self._data,*oid) 

    #======================================================================#
    def get_inv_var (self,*oid):
        """ return the inverse varience in a similar fashion to get_wl and get_data """
        return self._get_arr(self._inv_var,*oid) 

#    #======================================================================#
#    def get_len (self,*odi):
#        """ returns the length of a particualr order """
#        wl = self.get_wl(*oid)
#        return len(wl)


    #======================================================================#
    def get_wlbounds (self,*oid):
        wl = self.get_wl(*oid)
        return (np.min(wl),np.max(wl))
    
    def get_databounds (self,*oid):
        data = self.get_data(*oid)
        return (np.min(data),np.max(data))
    
    def get_min (self,*oid):
        """ returns tuple (min(wl),min(data)) over all the orders of self.get_band() """

        dat = self.get_data(*oid)
        wl = self.get_wl(*oid)

        if len(dat) == 0: return (None,None)
        return (np.min(wl),np.min(dat))

    #======================================================================#    
    def get_max (self,*oid,**kwargs):
        """ returns tuple (max(wl),max(data)) over all the orders of self.get_band() """

        dat = self.get_data(*oid)
        wl = self.get_wl(*oid)

        if len(dat) == 0: return (None,None)
        return (np.max(wl),np.max(dat))
 
    #======================================================================#
    
    def __getitem__ (self,order_index):
        """  return ndarray of wavelength and data for a given order of obj.get_band()"""
        oid = order_index
        dat = self.get_data(oid)
        wl = self.get_wl(oid)
        inv_var = self.get_inv_var(oid)
        return np.vstack((wl,dat,inv_var))

    #======================================================================#
    def set_band (self,band_id):
        """ use this to set the current band self._band_index """
        self._band_index = self._check_band_num(band=band_id)

    #======================================================================#
   
   
    def access_data (self):
        return self._data
    
    def access_wl (self):
        return self._wl

    def access_inv_var (self):
        return self._inv_var
     
    def _set_arr (self, old_arr, arr, oid=None,save_info=False):
        arr = np.asarray(arr)
        oid = self._convert_oid(oid)
        
        # if this is None then replace the entire array
        if oid is None:  
            assert old_arr.shape == arr.shape
            old_arr = arr
            
        # else just replace certain values
        else:
            old_arr[oid[0]][oid[1]] = arr

        if not save_info: self._reset_notes_and_info(oid=oid,skip_notes=True)
        return old_arr
        
    def set_data (self, data, oid = None, save_info=False):
        """
        used for setting the data
        """
        out = self._set_arr(self._data, data, oid, save_info)
        if out is not None: self._data = out
        
    #======================================================================#
    def set_wl (self, wl, oid=None, save_info=False):
        """
        used for setting the data
        """
        out = self._set_arr(self._wl, wl, oid, save_info)
        if out is not None: self._wl = out

    def set_inv_var (self, inv_var, oid=None):
        out = self._set_arr(self._inv_var, inv_var, oid, save_info=True)
        if out is not None: self._inv_var = out
        
    #======================================================================#
    def set_use_cropped (self,truth):
        """ 
        This is used to modify the output of get_data and get_wl, see get_use_cropped
        """
        if truth in ['prev','previous','use previous']:
            truth = deepcopy(self._prev_use_cropped)
        self._prev_use_cropped = deepcopy(self._use_cropped)
        self._use_cropped = bool(truth)
     
    #======================================================================#
#     def set_inv_var (self,oid,inv_var):
#         """
#         returns an editable version of inverse varience

#         if no oid give it returns the section corresponding to the current obj.get_band()

#         """
#         oid = self._convert_oid(oid)

#         if len(oid) == 0:
#             # if no arguments given
#             # use self._band_index  and applies to all the orders
#             self._check_band_num()
#             return self._inv_var[self._band_index]

#         if oid[1] is None: 

#         else:
#             # this allows to set for a specific band or order
#             oid = self._convert_oid(*oid)
#             if oid is None: return None
#             # !! I should probably have a set_data which does this and have get_data return a deepcopy   
#             return self._inv_var[oid[0]][oid[1]]
            
        
    #======================================================================#
    def update_params (self):
        """
        This will update some parameters such as self.shape and perform checks
   
        """
        all_good = [" ","-"*60]
        #----------------------------------------------------#
        # check to make sure self._wl, self._data, self._inv_var are still ndarrays
        def check_arr (arr,all_good,which):
            if type(arr) != np.ndarray: 
                try: arr = np.asarray(arr)
                except: all_good.append("ERROR: could not convert"+which+" array to an ndarray")
            
            if arr.ndim != 3: all_good.append("ERROR: "+which+" array is not of dimension 3, instead:"+str(arr.ndim))
            return arr,all_good
            
        self._wl, all_good      = check_arr(self._wl,all_good,'wavelength')
        self._data, all_good    = check_arr(self._data,all_good,'data')
        self._inv_var, all_good = check_arr(self._inv_var,all_good,'inverse variance')
                
        # update:
        self.shape = self._data.shape

        #----------------------------------------------------#
        # check the shape
        if self.shape != self._wl.shape: all_good.append('ERROR: obj._wl does not match shape:'+str(self.shape))
        if self.shape != self._inv_var.shape: all_good.append('ERROR: obj._inv_var does not match shape:'+str(self.shape))

        #----------------------------------------------------#
        if type(self._private_info).__name__ != 'dict': all_good.append('ERROR: obj._private_info not of type dict\n>>perhaps you should run obj._reset_notes_and_info()')
        if type(self._notes).__name__ != 'dict': all_good.append('ERROR: obj._private_info not of type dict\n>>perhaps you should run obj._reset_notes_and_info()')

        #----------------------------------------------------#
        # check self._bands and self._orders
        self.edit.sort_orders()

        #----------------------------------------------------#
        self._check_band_num()
        
        #----------------------------------------------------#
        if len(all_good) > 2:            
            print "\n".join(all_good)
            raise AttributeError("------ SEE ABOVE ERRORS ------------------")

    #======================================================================#
    def copy (self):
        """ copy the current object, if a value is given it will copy just that particular band over  """
        band = None
        if band is None:
            obj_new = deepcopy(self)
            obj_new._initialize_time = time.ctime()
        else:
            band = self._check_band_num(band=band)
            obj_new = deepcopy(self)
            obj_new._initialize_time = time.ctime()
            print "Currently can only copy all bands" #!! maybe change this, It should work but their's some logic bug
            
        return obj_new


#     #======================================================================#
#     def save_pkl (self,fileprefix=None,unique_name=True,clobber=False):
#         """ 
#         this is used to save a pickled version of the current object
        
#         INPUTS:
#         fileprefix : (str) used as the prefix save for the file, default is "ObjSave" with an added object name if self._obj_name is not None
#                      if value is given it will append '.pkl' i.e. 'foo' => 'foo.pkl'

#         unique_name : (bool) if True then it will make sure the file saves by checking all file names and adding _# if multiples appear
        
#         clobber : (bool) if True it will delete existing files, if False it will raise an error if the file exists, Not applicable if unique_name is True

#         """
#         unique_name = bool(unique_name)
#         clobber = bool(clobber)

#         #--------------------------------------------#
#         if fileprefix is None:
#             filename = 'ObjSave'
#             if self._obj_name is not None: filename += '_'+obj._obj_name
#         else: 
#             if type(fileprefix).__name__ in ['str','string_']: filename = fileprefix
#             else: raise TypeError("fileprefix must be of type string")
        
#         filename += '.pkl'

#         #--------------------------------------------#
#         if unique_name:
#             # save only unique names by checking if the file already exists
#             i = 1
#             while True:
#                 if os.path.exists(filename):
#                     filename = filename[:-4]+"_"+str(i)+filename[-4:]
#                     i+=1
#                 else: break
            
#         else:
#             if os.path.exists(filename) and not clobber: raise IOError("clobber is set to false and given file exists: "+filename)

#         #--------------------------------------------#
#         pickle.dump(self,open(filename,'wb'))


    #======================================================================#
    def __getstate__ (self):
        """ for saving the state,  pickle.dump(obj,open('saveit.pkl','wb')) """
        dump_stuff = self.__dict__
        #del dump_stuff['hdulist'] # if I keep self.hdulist I have to figure out something to do with this
        return dump_stuff
    #======================================================================#
    def __setstate__ (self,state_dict):
        """ for loading the state, new_obj = pickle.load(open('saveit.pkl','rb')) """
        # get the basics back
        self.__dict__ = state_dict
        self._initialization_time = time.ctime()
        # everything except hdulist (because it's a file object)
        # !! if I decide not to keep self.hdulist this is not needed
        #new_header = pyfits.PrimaryHDU(data=self._data,header=self.header)
        #self.hdulist = pyfits.HDUList([new_header])
        #self.__dict__['hdulist'] = hdulist

        #for (name,value) in data_dict.iteritems():
        #    setattr(self,name,value)




