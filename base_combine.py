if __name__ != "__main__":
    from eyeSpec.base_functions import find_overlap_pts, inv_var_2_var, var_2_inv_var
    from eyeSpec.base_classes import query_fits_header, eyeSpec_spec
    from eyeSpec.resampling import get_resampling_matrix
    from eyeSpec.dependencies import np, deepcopy, scipy, pdb


def spec_combine (obj,line_data=None,include_weights=True,include_varience=True,
                  overshoot=4,crop_end_behavior=2,fill_between=True,fill_value=1,
                  linear_interpolate=True):
    """
    This uses the line data to calculate how the orders of obj in obj.get_band() 
    should be combined

    ============================================================================

    INPUTS:
    obj : (eyeSpec_spec) data class, must have obj.shape[1] > 1

    line_data : (float_array,None) Must be convertable to a numpy float array,
    if None then it will calcuate the overlap points based on the mid points 
    between the previous and next orders
    
    include_weights : (bool) If True it will impose a triangle weighting so that
    the next order is added into the previous order based on linear weighting. 
    Weights for previous is one up to the overlap point and then linearly 
    decreases to zero from the overlap point to the end it's data. Weights 
    for the next order are the reverse 

    include_varience: (bool) If True it will use the inverse varience stored in
    the obj when co-adding the points. If False the code applies a equal 
    weighting (depending on scaling weight)

    overshoot : (int) The number of points to be used to overshoot the overlap 
    point which will then be removed to correct for unusual end behavior

    crop_end_behavior : (int) The number of points to be removed from the end 
    of the orders to correct for unusual end behavior

    fill_between : (bool) If True then it will fill in values between orders 
    when they don't overlap. The inverse varience will be set to zero.

    fill_value : (float) If fill_between is True it will use this value to fill 
    in with. Again the inverse varience will be set to zero.

    linear_interpolate : (bool) If True it will just do a linear interpolation 
    and not the resampling. It will set overshoot and crop_end_behavior to 0


    OUTPUTS:    
    obj : (eyeSpec_spec) a single order with the co-added data


    ============================================================================
    Logic explained:

    the previous order representation, ends at point b
    ::::::::::::::::::::b
    
    the next order representation, begins at point a
             a;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    the overlap point is where they intersect
    :::::::::a:::::i::::::b
             a;;;;;i;;;;;;b;;;;;;;;;;;;;;;;;;


    the co-addition is done based on the point i, the 
    two main segments [a,i] and [i,b]

    ===== for segment [a,i]
    Assume that the previous order will have the best quality
    resample the next order onto the the previous, this can induce some bad 
    end behavior which is subsequently controlled by
    inputs crop_end_behavior and overshoot
    (See the code for notes on the use of resampling)
    overshoot = 2
    a;;;;;i;;  ==(match wavelength to)==> a:::::i::

    scaling weights adjust the weight of the next order linearly base on the
    number of points minus the number of cropped end points from 
    input crop_end_behavior and the weight given to the previous order
    is uniformly one this can be controlled by setting the input include_weights
    crop_end_behavior = 2
    overshoot = 3
    a;;;;;i;;  ==(weighted by)==> [0,0,0,.2,.4,.6,.8,1,1,1,1]

    the inverse varience from the obj data is also used to weight the data if 
    the input include_varience is True

    the code then does a weighted average to get the output data
    a;;;;;i;;  +average+ a::::::i::  =  a////////i//

    Finally to deal with the end points I take off the input crop_end_behavior
    from the beginning
    and remove the overshoot points from the end
    a+//////i
    
    ==== for segment [i,b]
    this is done similar to the segment [a,i] except assuming that the next 
    order will have the best quality data

    the previous order is resampled to match the next order. Additionally, if
    input include_weights is True the code will apply a linear weighting from 
    1 to 0 based on the number of points minus the cropped end points (again 
    similar to the [a,i] segment)

    if the input include_varience then the varience will be used when
    averaging the segments of the two orders
    ::i:::::b  +average+ ;;i;;;;;;b  =  \\i\\\\\\\b    ====>   i\\\\\\b-


    ==== put it all together
    the new output is the concatenation of each piece
    [prior,a]+[a,i]+[i,b]+[b,following]
    ::::::a::a+//////i\\\\\\\\b-;;b;;;;;;;;;;;;



    NOTE:
    If no overlap point is given for a region which has overlap then the two 
    will just be spliced together:
    ::::::a::::::_mid_:::::::b
          a;;;;;;_mid_;;;;;;;b;;;;;;;;;;
    co-added gives:
    ::::::a::::::_mid_;;;;;;;b;;;;;;;;;


    NOTE:
    If the two orders don't overlap, depending on the input fill_between
    :::::::a              b;;;;;;;;;;;;;
    if fill_between == True, the code outputs
    :::::::a--fill_value--b;;;;;;;;;;;;
    else it'll just concatenate them and leave the gap in wavelength between
    pixel a and b
    :::::::a_b;;;;;;;;;;;;
    

    """
    #==========================================================================#
    # check input variables
    if obj.__class__.__name__ != 'eyeSpec_spec': raise ValueError("obj MUST BE OF CLASS eyeSpec_spec") 

    if obj.shape[1] == 1: raise ValueError("obj must have more than one order")

    # these are for dealing with bad end behavior
    overshoot = int(overshoot)
    crop_bad_behavior = int(crop_end_behavior)

    # check if want just simple linear interpolation
    linear_interpolate = bool(linear_interpolate)
    if linear_interpolate:
        overshoot = 0
        crop_bad_behavior = 0

    # check the inclusion of weights
    include_weights = bool(include_weights)
    include_varience = bool(include_varience)

    # check the fill variables
    fill_between = bool(fill_between)
    fill_value = float(fill_value) # the fill value must be type floating point

    # if no line_data is given then go with default
    if line_data == None: line_data = find_overlap_pts(obj)

    # check the line data:
    try: line_data = np.array(line_data,dtype=float)
    except: 
        raise TypeError("Variable line_data must be array of single type float")
    if line_data.ndim != 1: 
        raise TypeError("Variable line_data must have only one dimension")

    # set the data to only output the cropped orders
    tmp_use_cropped = obj.get_use_cropped()
    obj.set_use_cropped(True)

    #==========================================================================#
    # do the co-addition work
    # walk through the wavelength orders in increasing order. The code is 
    # building new arrays (A_wl,A_data,A_inv_var). The current order
    # is notated by cur (B_wl,B_data,B_inv_var). The resulting co-add
    # or other is concatenated onto the new arrays until it reaches the end

    #These graphics are used to visualize the possibilities
    # ::::::::::b      <= indicies for A_wl
    #   a;;;;;;;;;;;;; <= indicies for B_wl
    
    # a,a+,b,b-,i represent important indicies in the data used for breaking
    # up the arrays for co-addition

    is_not_normalized = False

    for i in range(obj.shape[1]):
        #----------------------------------------------------------------------#
        print "current order:",i        
        not_normalized_level = 100
        max_number_above = 20        

        # check if it's the first array and initate the new arrays
        if i == 0:
            # initiate the new arrays 
            A_wl = obj.get_wl(i)
            A_data = obj.get_data(i)
            A_inv_var = obj.get_inv_var(i)
            C = np.ones(len(A_data))
            if len(C[(A_data > not_normalized_level)]) > max_number_above: is_not_normalized = True
            continue
        else:
            # get current arrays
            B_wl = obj.get_wl(i)
            B_data = obj.get_data(i)
            B_inv_var = obj.get_inv_var(i)
            C = np.ones(len(B_data))
            if len(C[(B_data > not_normalized_level)]) > max_number_above: is_not_normalized = True

        #----------------------------------------------------------------------#
        # what if: => 
        # :::::::::::::::b
        #                a::::::::::::::::::
        # if they are equal make the fill arrays len = 0
        # then concatenated them together while dropping the first point of the
        # current order, assumed to be of poorer quality
        if np.min(B_wl) == np.max(A_wl):
            A_wl = np.concatenate((A_wl,B_wl[1:]))
            A_data = np.concatenate((A_data,B_data[1:]))
            A_inv_var = np.concatenate((A_inv_var,B_inv_var[1:]))
            continue

        #----------------------------------------------------------------------#
        # what if: => have some filler value between or just put the two arrays 
        # together and have a large gap in wavelength
        # ::::::::::b    mid
        #                mid     a::::::::::::::::::
        # what if: ignores and treats like prevous
        # ::::::::::|:::::b    mid
        #                      mid    a::::::::::::::::::
        if np.min(B_wl) >= np.max(A_wl):
            print ("no overlap found with next order, I'm filling it with a linear spacing.\n     This probably isn't the right thing to do because the del_wl spacing in the previous order may be larger than in the next order")

            # when filling in with artificial data points I have to make a 
            # decision on what resolution is (# pixels per wavelength). 
            # I don't necessarily assume this is linear. Instead I extrapolate
            # the new and current orders on their wavelength solution. However,
            # I ran into problems extrapolating the current order backwards so
            # instead the new order is just extrapolated forward to the 
            # beginning of the current. There will likely be a discontinuity at
            # this point.

            # discontinuous at end point
            mid_pt = np.min(B_wl)

            def __poly_fit (X,coeff):
                return coeff[7] + X*(coeff[6] + X*(coeff[5] + X*(coeff[4] + X*(coeff[3] + X*(coeff[2] + X*(coeff[1] + X*(coeff[0])))))))

            # from a => mid_pt = b
            pts = np.arange(1,len(A_wl)+1) 
            coef_new = scipy.polyfit(pts,A_wl,7)
            A_extend = __poly_fit(np.arange(len(A_wl),len(A_wl)+1e4),coef_new)
            mid_i = np.abs(A_extend-mid_pt).argmin()
            if mid_i == len(A_extend)-1: 
                print "Whoops, I didn't give the extension enough points"
            wl_fill = A_b2mid = A_extend[:mid_i]
            
            # attempt to extrapolate the wavelength solution 
            # for end of B back to mid
            # combining:
            # ::::: a :::: mid ;;;; b ;;;;;;;
            # from mid_pt => b
#             pts = np.arange(1,len(B_wl)+1)
#             coef_cur = scipy.polyfit(pts,A_wl,7)
#             
#             B_extend = __poly_fit(np.arange(0,1e4),coef_cur)
#             if mid_i == 0: print "Whoops, I didn't give the extension enough points"
#             mid_i = np.abs(B_extend-mid_pt).argmin()
#             B_mid2a = B_extend[mid_i:]
#             wl_fill = np.concatenate((A_b2mid,B_mid2a))

            data_fill = np.ones(len(wl_fill))*fill_value
            inv_var_fill = np.zeros(len(wl_fill))
            
            if fill_between:
                A_wl = np.concatenate((A_wl,wl_fill,B_wl))
                A_data = np.concatenate((A_data,data_fill,B_data))
                A_inv_var = np.concatenate((A_inv_var,inv_var_fill,B_inv_var))
            else:
                A_wl = np.concatenate((A_wl,B_wl))
                A_data = np.concatenate((A_data,B_data))
                A_inv_var = np.concatenate((A_inv_var,B_inv_var))
            continue


        #----------------------------------------------------------------------#
        all_mask = (line_data > np.min(B_wl))*(line_data < np.max(A_wl))
        overlap_pt = line_data[all_mask]

        # what if: data is scaled larger in B_
        # vvvvvvvvvvavvvvvvb
        #           a^^^^^^b^^^^^^^^^^^^^

        found_overlap_pt = True
        if overlap_pt.shape[0] == 0:
            found_overlap_pt = False
            print "HeadsUp: there are no overlap points found, splicing at the middle point"
            overlap_pt = (np.max(A_wl) + np.min(B_wl))/2.
        else: overlap_pt = overlap_pt[0]

        A_i_correct = np.abs(A_wl - overlap_pt).argmin()
        B_i_correct = np.abs(B_wl - overlap_pt).argmin()

        # to correct for weird end behavior when resampling if the user
        # is just interpolating then the overshoot will just be zero
        # :::::::::::::a::a+::::::i-::::i+:::::b-::b
        #              a::a+::::::i-::::i+:::::b-::b:::::::::
        # the i+ and i- fully recover the original
        # the a+ and b- loose those data points from the ends
        
        # !! overshoot = 4
        if A_i_correct+overshoot >= len(A_wl): print "For this order carefully examine the end point, there may be some unusual behavior" 
        if B_i_correct-overshoot <= 0: print "For this order carefully examine the end point, there may be some unusual behavior" 

        A_i_plus = np.clip(A_i_correct+overshoot,0,len(A_wl)-1)
        B_i_plus = np.clip(B_i_correct+overshoot,0,len(B_wl)-1)

        A_i_minus = np.clip(A_i_correct-overshoot,0,len(A_wl)-1)
        B_i_minus = np.clip(B_i_correct-overshoot,0,len(B_wl)-1)

        #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # what if: line is given outside of being between a region
        # :::|:::a:::::::::b
        #        a:::::::::b:::::::::
        # what if: no line is give for a regions
        # ::::::a:::::::::b
        #       a:::::::::b:::::::::
        # then return a splice of the two at the mid point
        # ::::::a::::i;;;;b;;;;;;;;

        # If you want to have the code find it's own middle and then do the 
        # interpolating and co-adding
        # instead of splicing then uncomment this line ==>
        # found_overlap_pt = True
        
        if not found_overlap_pt: 
            A_wl = np.concatenate((A_wl[:A_i_correct],B_wl[B_i_correct:]))
            A_data = np.concatenate((A_data[:A_i_correct],B_data[B_i_correct:]))
            A_inv_var = np.concatenate((A_inv_var[:A_i_correct],B_inv_var[B_i_correct:]))
            continue



        #+++++++++++++++++++++ determine index points ++++++++++++++++#
        # what if: ideal
        # ::::::a:::|::::::b
        #       a;;;|;;;;;;b;;;;;;;;;;;
        # then result:  :::::: a+ ////// i \\\\\\ b- ;;;;;;;;

        #-------------------- possible cross-correlation -------------#
        # !! do I need to cross correlate the segments to match the features up?
        
        A_a = np.abs(A_wl-np.min(B_wl)).argmin()
        B_a = 0 

        A_b = len(A_wl)-1
        B_b = np.abs(B_wl - np.max(A_wl)).argmin()

        # set up an end point cropping to deal with a few pixels 
        # at the end which the re-sampling doesn't do well with
        # what if: ideal
        # ::::::a::a+:::|:::::b-:::b
        #       a;;a+;;;|;;;;;b-;;;b;;;;;;;;;;;       
        # crop_bad_behavior = 2
        
        A_a_plus = np.clip(A_a+crop_bad_behavior,0,len(A_wl)-1)
        B_a_plus = np.clip(crop_bad_behavior,0,len(B_wl)-1)
        
        A_b_minus = np.clip(A_b-crop_bad_behavior,0,len(A_wl)-1)
        B_b_minus = np.clip(B_b-crop_bad_behavior,0,len(B_wl)-1)

        # if there are multiple it will take the first on but this could give:
        # what if: won't happen in wavelength space
        # ::::a;;;;;;;;;;;;b
        #          a:::::::b::::::::::::

        # these segments are used for co-adding
        # set up an end point cropping to deal with a few pixels 
        # at the end which the re-sampling doesn't do well with
        # what if: ideal
        # ::::::a::a+::::i-::i::i+:::::b-:::b
        #       a;;a+;;;;i-;;i;;i+;;;;;b-;;;b;;;;;;;;;

        #======================================================================#
        def _calc_new_segment (ai,bi,weights_seg):
            """
            NOTE: This function assumes certain global variables
            from the spec_combine function:
            linear_interpolate
            include_varience
            A_wl,A_data,A_inv_var
            B_wl,B_data,B_inv_var
            
            and functions:
            spectrum_interp
            spectrum_combine

            """
            
            # get the data for the segment of interest
            A_seg,B_seg,B_co_var = spectrum_interp([A_wl, A_data, A_inv_var],
                                                   [B_wl, B_data, B_inv_var],
                                                   ai = ai,
                                                   bi = bi,
                                                   linear_interpolate = linear_interpolate)
            # A_wl_seg, A_data_seg, A_inv_var_seg = AB_seg[0]
            # B_wl_seg, B_data_seg, B_inv_var_seg = AB_seg[1]

            # decide what input to give spectrum_combine
            if include_varience:
                A_in_seg = A_seg[:3]
                B_in_seg = B_seg[:3]
            else:
                A_in_seg = A_seg[:2]
                B_in_seg = B_seg[:2]
                               
            Z_wl_seg, Z_data_seg, Z_inv_var_seg = spectrum_combine(A_in_seg, B_in_seg,
                                                                   weights=weights_seg,
                                                                   B_co_varience = B_co_var)

            return Z_wl_seg, Z_data_seg, Z_inv_var_seg

        #======================================================================#
        #++++++++++++++++++++++++++ do co-addition ++++++++++++++++++++++++++++#
        # this applies the various functions just defined in order to combine
        # the overlap of order A and order B


        # artificial triange weights range from 0 to 1. For cropping bad end 
        # some filler has been applied. 
        # Order A => ::::a:::::a+:::::::::i-:::::::::::i+:::::::::::b-:::::b
        #                _ 111 __ 1==>0.5 __ 0.5===0.5 __ 0.5====>0 __ 000 _
        #
        # Order B =>     a;;;;;a+;;;;;;;;;i-;;;;;;;;;;;i+;;;;;;;;;;;b-;;;;;b;;;;
        #                _ 000 __ 0==>0.5 __ 0.5===0.5 __ 0.5====>1 __ 111 _

        # Determine the filler values for the overshooting
        zeros_fix_bad_behavior = np.zeros(crop_bad_behavior) 
        ones_fix_bad_behavior = np.ones(crop_bad_behavior)
        fill_overshoot = 0.5*np.ones(overshoot)

        #=========================== for segment [a, i+)  ====================#
        # assume that the A_data array has the best quality of data for this 
        # region so the B_data is interpolated/resampled onto it. For this 
        # segment this means that likely the high-resolution end of the current
        # order is interpolated onto the low-resolution end of the new data.
        
        # determine artificial triangle weighting
        if include_weights:
            num_pts = len(np.arange(A_a,A_i_plus)) - crop_bad_behavior-overshoot

            if num_pts < 0: raise ValueError("The inputs for crop_end_behavior and overshoot resulted in negative points here")

            A_in_ai_wei = np.concatenate((ones_fix_bad_behavior,
                                          np.linspace(1,.5,num_pts),
                                          fill_overshoot))  # 1 => .5
            B_in_ai_wei = np.concatenate((zeros_fix_bad_behavior,
                                          np.linspace(0,.5,num_pts),
                                          fill_overshoot))  # 0 => .5
            weights_ai = [A_in_ai_wei, B_in_ai_wei]
        else: weights_ai = None

        # determine new array:
        Z_wl_ai, Z_data_ai, Z_inv_var_ai = _calc_new_segment([A_a,A_i_plus],
                                                             [B_a,B_i_plus],
                                                             weights_ai)

        #=========================== for segment [i-,b] =======================#
        # for this segment I assume that the B_data has the best quality
        # so interpolate/resample the A_data onto the B_data. This means
        # interpolating the low_resolution end of the A_data onto the 
        # high-resolution end of the B_data.

        # determine artificial triangle weighting
        if include_weights:
            num_pts = len(np.arange(A_i_minus,A_b)) - crop_bad_behavior-overshoot
            if num_pts < 0: raise ValueError("The inputs for crop_end_behavior and overshoot resulted in negative points here")

            A_in_ib_wei = np.concatenate((fill_overshoot,
                                          np.linspace(0.5,0.0,num_pts),
                                          zeros_fix_bad_behavior)) # .5 => 0

            B_in_ib_wei = np.concatenate((fill_overshoot,
                                          np.linspace(0.5,1.0,num_pts),
                                          ones_fix_bad_behavior))  # .5 => 1
            weights_ib = [A_in_ib_wei, B_in_ib_wei]
        else: weights_ib = None

        # determine new array:
        Z_wl_ib, Z_data_ib, Z_inv_var_ib = _calc_new_segment([A_i_minus,A_b], # including i
                                                             [B_i_minus,B_b+1],  # including b and i
                                                             weights_ib)

        #======================================================================#
        #+++++++++++++ correct for end behavior and concatenated  +++++++++++++#
        # correct for overshoot and bad end behavior
        # a//a+/////////i///i+  ==>  a+//////////i
        Z_wl_ai = Z_wl_ai[crop_bad_behavior:-overshoot-1]
        Z_data_ai = Z_data_ai[crop_bad_behavior:-overshoot-1]
        Z_inv_var_ai = Z_inv_var_ai[crop_bad_behavior:-overshoot-1]

        # i-\\i\\\\\\\\\b-\\b   ==>  i\\\\\\\b-
        Z_wl_ib = Z_wl_ib[overshoot:-crop_bad_behavior-1]
        Z_data_ib = Z_data_ib[overshoot:-crop_bad_behavior-1]
        Z_inv_var_ib = Z_inv_var_ib[overshoot:-crop_bad_behavior-1]
        
        # put all the segments together:
        # :::::: a+ ////// i \\\\\\ b- ;;;;;;;;
        A_wl = np.concatenate((A_wl[:A_a_plus],Z_wl_ai,Z_wl_ib,B_wl[B_b_minus:]))
        A_data = np.concatenate((A_data[:A_a_plus],Z_data_ai,Z_data_ib,B_data[B_b_minus:]))
        A_inv_var = np.concatenate((A_inv_var[:A_a_plus],Z_inv_var_ai,Z_inv_var_ib,B_inv_var[B_b_minus:]))



    # what happens if: many orders overlap
    #:::::::::a:::::|::::|:::b
    #         a:::::|:c::::::b:::|:::::d
    #                 c::|:::::::|:::::d::::::::
    # then it will deal with them in order of appearing

    #==========================================================================#
    if is_not_normalized: print "HeadsUP: For 1-D spectra without the blaze removed, the co-addition may not work well"

    # return previous use cropped preference
    obj.set_use_cropped(tmp_use_cropped)

    # create the new data object
    A_obj = eyeSpec_spec(A_wl,A_data,inverse_varience=A_inv_var,header=deepcopy(obj.header))

    # set the correct lengths in the header files
    naxis2 = query_fits_header(A_obj.header,'naxis2')
    naxis3 = query_fits_header(A_obj.header,'naxis3')

    A_obj.header.update('NAXIS',1)
    A_obj.header.update('NAXIS1',len(A_wl))
    
    if naxis2.found: del A_obj.header['naxis2']
    if naxis3.found: del A_obj.header['naxis3']

    if 'combo_hdr' in dir(obj): A_obj.combo_hdr = deepcopy(obj.combo_hdr)    
    return A_obj

def co_add (spec_list,use_orders=None,no_interp=False,linear_interpolate=False):
    spec_list = np.array(spec_list)
    if len(spec_list) <= 1: raise ValueError("spec list must have two or more elements")
    if spec_list[0].__class__.find("eyeSpec_spec") == -1: raise ValueError("Elements of spec list must be of class eyeSpec_spec")

    num_orders,num_points = spec_list[0].shape[1:]

    Z_orders = []
    Z_headers = []

    if use_orders is None: ran = range(num_orders)
    else:
        use_orders = np.array(use_orders,dtype=int)
        ran = use_orders
    
    # go through all the orders
    for i in ran:
        # for each order go through all the spectra
        C_order = []
        for j in range(len(spec_list)):
            spec = spec_list[j]
            spec.set_use_cropped(False)
            if spec.shape[1:] != (num_orders,num_points):
                print "Spectrum number "+str(j)+" in list is not of the correct shape, i.e. (#orders,#points)="+str(num_orders)+", "+str(num_points)
                del spec_list[j]
                continue

            # record the headers once
            if i == 0: Z_headers.append(deepcopy(spec.header))

            spec_order = [spec.get_wl(i),spec.get_data(i),spec.get_inv_var(i)]
            # if it's the first spectra take that regularly
            if j == 0: 
                C_order = deepcopy(spec_order)
                continue

            # if it's not the first then co_add to get the next
            C_order = co_add_order(C_order,spec_order,no_interp=no_interp,linear_interpolate=linear_interpolate)
            spec.set_use_cropped('use previous')
        Z_orders.append(C_order)
    Z_orders = np.array(Z_orders)

    spec_Z = eyeSpec_spec(Z_orders[0],Z_orders[1],Z_orders[2],Z_headers[0])
    spec_Z.combo_hdr = Z_headers
    return spec_Z

def co_add_order (ord_A,ord_B,no_interp=False,ai=None,bi=None,wl_bounds=None,linear_interpolate=False,weights=None):
    if no_interp:
        A_seg = A
        B_seg = B
        B_co_var = None
    else:
        A_seg, B_seg, B_co_var = spectrum_interp(A,B,ai,bi,wl_bounds,linear_interpolate)
        
    Z_wl, Z_data, Z_inv_var = spectrum_combine(A_seg,B_seg,weights,B_co_var)
    return (Z_wl, Z_data, Z_inv_var)



#======================================================================#

def spectrum_interp (A,B,ai=None,bi=None,wl_bounds=None,linear_interpolate=False):
    """
Assumes spectrum A has higher signal-to-noise and interpolates B onto A's wavelengths


++++++++++++++++ get resampling for specific regions ++++++++++++
 To match the data of the two orders instead of using straight linear 
 interpolation instead a resampling (re-binning) is performed. This was
 prompted for other science for resampling high resolution data to lower
 resolution for comparison. However, here it provides a more accurate
 matching and correct determination of the errors. 

 Using straight interpolation will often cause a small gain in flux for
 absorption spectra. Features will be cut in the very bottom of their 
 core, resulting in a gain in flux. By summing flux over the same 
 region of a linearize (linearly interpolated) and non-linearlized
 data set, there is a difference on the order of 0.3% either positive
 or negative. Because this is absorption spectra the non-linearlized 
 minus linearlized often results in a -0.3% difference because flux
 was gained when linearlizing.

 The resampling (re-binning) is done in a flux conserving way. The 
 flux for an input is spread over the bins for the output. The code
 used to do this (get_resampling_matrix) returns a transformation 
 matrix which converts the flux in the current pixel to flux in the
 new bin. 
 flux_new[i] =  ...+C2*flux[i-1] + C3*flux[i] + C4*flux[i+1]+...
        
 The C coeffients are encapsulated in the transformation matrix. We
 perform a normalization by the coefficients to return something on 
 the same scale as the new binned data. On high-resolution data this 
 is approximately the same as linear interpolation with the small
 correction.

 This also improves the handling of the varience for each flux.
 (flux varience = flux_er2). If we want the varience on the new
 flux bin then we perform the expectation of the error calculation:
 flux_A_er2[i] =  <(...+C2*flux[i-1] + C3*flux[i] +
                         C4*flux[i+1]+...)**2>
                =  C2**2*<flux[i-1]>**2+C3*<flux[i]>**2+
                          C4**2*<flux[i+1]>**2+ CROSS_TERMS+....

 but <flux[i]>**2 = flux_er2[i]

 flux_A_er2[i] = C2**2*flux_er2[i-1]+C3**2*flux_er2[i]+
                           C4**2*flux_er2[i+1]+ CROSS_TERMS+....
        
 and we assume the correlation between different errors (i.e. the
 cross terms) to be zero. In other words, the errors are all 
 independent from one another so we drop the CROSS_TERMS and have:
  flux_A_er2[i] = C2**2*flux_er2[i-1]+C3**2*flux_er2[i]+
                          C4**2*flux_er2[i+1]+....
        
 but this is just the element wise square of the transform matrix we
 calculated before (all the C coefficients) multiplied by the varience
 array for the original flux.

 The caveat is that the data needs to be normalized first for the 
 resampling to be the best.

    """
    #------------------------------------------------------------------#
    # check input variables
    linear_interpolate = bool(linear_interpolate)

    A = np.array(A,dtype=float)
    B = np.array(B,dtype=float)
    if A.shape[0] != 3 or A.ndim != 2: raise ValueError("The spectrum A needs to be of form [wavelengths,data,inv_var] = A")
    if B.shape[0] != 3 or B.ndim != 2: raise ValueError("The spectrum B needs to be of form [wavelengths,data,inv_var] = B")

    A_wl, A_data, A_inv_var = A
    B_wl, B_data, B_inv_var = B
    
    if wl_bounds is not None:
        wl_bounds = np.array(wl_bounds,dtype=float)
        if i.ndim != 1 or i.shape[0] != 2: 
            raise ValueError("Wavelength bounds must be in form wl_bounds = [start_wl, end_wl]")
        
        ai,bi = [],[]
        ai.append(np.abs(A_wl - wl_bounds[0]).argmin())
        ai.append(np.abs(A_wl - wl_bounds[1]).argmin())

        bi.append(np.abs(B_wl - wl_bounds[0]).argmin())
        bi.append(np.abs(B_wl - wl_bounds[1]).argmin())
        # !! might be good to induce my own overshoot and cropping
        # for behavior of the resampling at the ends

    def _check_indicies (i,wl_spec,error_msg):
        if i is None: return np.array([0,len(wl_spec)],dtype=int)
        # !! could also include a way to give wavelength bounds
        i = np.array(ai,dtype = int)
        if i.ndim != 1 or i.shape[0] != 2: raise ValueError(error_msg)
        return i

    ai = _check_indicies(ai,A_wl,"The bounds for spectrum A must be in form ai = [start_i,end_i]")
    bi = _check_indicies(bi,B_wl,"The bounds for spectrum B must be in form ai = [start_i,end_i]")

    A_wl_seg, A_data_seg, A_inv_var_seg = A_wl[ai[0]:ai[1]],  A_data[ai[0]:ai[1]],  A_inv_var[ai[0]:ai[1]] 
    B_wl_seg, B_data_seg, B_inv_var_seg = B_wl[bi[0]:bi[1]],  B_data[bi[0]:bi[1]],  B_inv_var[bi[0]:bi[1]]  

    #------------------------------------------------------------------#
    # put B onto A
    if linear_interpolate: 
        #>>>>>>>>>>>>>> Interpolates B onto A
        B_data_interp = scipy.interp(A_wl_seg, B_wl, B_data)
        # because you assume the varience pixel-to-pixel 
        # from the original data to be independent you can 
        # just interpolate
        B_inv_var_interp = scipy.interp(A_wl_seg, B_wl, B_inv_var)

        # Note: Can't get much more out of including the co-varience from
        # interpolation 
        B_co_var = None

    else:
        # check the normalization. I'm not sure how necessary this is
        not_normalized_level = 100
        max_number_above = 10
        C = np.ones(len(A_data))
        if len(C[(A_data > not_normalized_level)]) > max_number_above: print "HeadsUp: Spectrum A data does not seem to be normalized"
        if len(C[(B_data > not_normalized_level)]) > max_number_above: print "HeadsUp: Spectrum B data does not seem to be normalized"

        #>>>>>>>>>>>>>> resamples B onto A
        T = get_resampling_matrix(B_wl_seg,A_wl_seg)
        # normalize the values by one, assuming the data is 
        # continuum normalized
        norm = T*np.ones(len(B_data_seg))
        norm[(norm == 0)] = 1e-100
        norm_diag = scipy.sparse.dia_matrix((1.0/norm,0),shape=(len(norm),len(norm)))
        T = norm_diag*T
        # For each data point in B_data_seg this is making it a linear
        # combination of all other data points
        # B_data_interp[i] = t0*B_data_seg[0] + t1*B_data_seg[1] +...
        # most values will be zero
        B_data_interp = T*B_data_seg

        B_var_seg = inv_var_2_var(B_inv_var_seg)
        T2 = T.copy()
        T2.data **= 2
        # Doing a similar thing as with the data only the coefficients
        # for the varience must be squared
        # B_var_interp[i] = t0**2*B_var_seg[0]+t1**2*B_var_seg[1]+....
        B_var_interp = (T2*B_var_seg)
        B_inv_var_interp = var_2_inv_var(B_var_interp)
        
        # calcuate the co-varience of the varience between pixels
        # induced by the resampling. This causes a 'beating' effect
        # with the maximum effect where the distance between 
        # pixels of A and B are at a maximum offset
        B_inv_var_seg = var_2_inv_var(B_var_seg)
        B_var_seg_diag = scipy.sparse.dia_matrix((B_inv_var_seg,0),shape=(len(B_var_seg),len(B_var_seg)))
        T_tran = T.transpose(True)
        B_co_var = T*B_var_seg_diag*T_tran



    A_seg = [A_wl_seg, A_data_seg, A_inv_var_seg]
    B_seg = [A_wl_seg, B_data_interp, B_inv_var_interp]
    return A_seg, B_seg, B_co_var


#======================================================================#

def spectrum_combine (A,B,weights=None,B_co_varience=None):
    """
    This assumes that A and B have the same number of points and will 
    co-add A and B

    INPUTS:
    A,B : (array) generally give by [wavelenght,data,inv_var] where each
          is the same length
          note: Exclude Varience by only giving A = [wavelength,data]

    weights : (array) given by [weights_A, weights_B] where each is 
              the same length as the corresponding data A and B
              note: Exclude extra weights by giving weights = None

    """
    #------------------------------------------------------------------#
    # check input variables
    B_co_var = None
    if B_co_varience is not None:
        B_co_var = scipy.sparse.csr.csr_matrix(B_co_varience)


    def check_input (SPEC,error_msg):
        SPEC = np.array(SPEC,dtype=float)
        # if the shape is not (#_spec,#_pts_per)
        if len(SPEC.shape) != 2: raise ValueError(error_msg)
        # if #_spec = 3 then [wavelengths,data,inv_var]
        if SPEC.shape[0] == 3: SPEC_wl, SPEC_data, SPEC_inv_var = SPEC   
        # if #_spec = 2 then [wavelengths,data]
        elif SPEC.shape[0] == 2: 
            SPEC_wl, SPEC_data = SPEC
            # create unweighted inverse varience
            SPEC_inv_var = np.ones(len(SPEC_wl))
        else: raise ValueError(error_msg)
        return SPEC_wl, SPEC_data, SPEC_inv_var

    errorA = "The spectrum A needs to be of form A = [wavelengths,data,inv_var] or [wavelengths,data]"
    errorB = "The spectrum B needs to be of form B = [wavelengths,data,inv_var] or [wavelengths,data]"
    A_wl, A_data, A_inv_var = check_input(A,errorA)
    B_wl, B_data, B_inv_var = check_input(B,errorB)

    if len(A_wl) != len(B_wl): raise ValueError("Spectrum A must have the same number of data points as spectrum B")

    if weights is not None:
        weights = np.array(weights,dtype=float)
        
        if weights.ndim != 2: raise ValueError("The extra weights must be in the form weights = [weights_A,weights_B]")
        A_wei = weights[0]
        B_wei = weights[1]
        if len(A_wei) != len(A_wl): raise ValueError("Weights for spectrum A and B must have same length")
    else:
        A_wei = np.ones(len(A_wl))
        B_wei = np.ones(len(A_wl))

    # lengths are all equal: for A,B  _wl, _data, _inv_var, _wei
    # if inv_var or weights are not desired they will be set to one

    #------------------------------------------------------------------#
    #>>>>>>> Determine new wavelenghts
    # !! could add something here which does a cross correlation 
    # and adjustment
    Z_wl = deepcopy(A_wl)

    #>>>>>>> Determine new data as the weighted average of previous
    # first calculate the sum of your weights
    var_sums = A_wei*A_inv_var + B_wei*B_inv_var
    xA = A_wei*A_inv_var*A_data
    # Note: If there is co-varience in the B varience spectrum 
    # technically it should be included in the data calculation
    # but gives very little gain
    xB = B_wei*B_inv_var*B_data
    var_sums[(var_sums==0.0)] = 1e-100
    # determine the new data from the sum of the previous points
    # times their respective artificial weights (*_wei) and 
    # inverse varience. Normalize by the sum of the weights
    Z_data = (xA+xB)/var_sums
    
    #>>>>>>> Determine new inverse varience
    # calculate new inverse varience for that point
    
    # get the varience for A and B
    A_var = inv_var_2_var(A_inv_var)
    B_var = inv_var_2_var(B_inv_var)
    
    wei_sums = A_wei+B_wei
    # normalize just by the artificial weights
    A_wei /= wei_sums
    B_wei /= wei_sums
    

    # Z_var = <(a*A_err+b*B_err)> 
    #       = a**2*<A_err>**2 + b**2*<B_err>**2 + 2*a*b*<A_err*B_err>
    #       = a**2*A_var + b**2*B_var + 2*a*b*AB_co_var

    # we assume that the varience in A and B are not co-variant
    # and hence ignore the AB_co_var term

    # the B_co_var is included because the B varience is dependant 
    # on it's nearby points because of the resampling or interpolation

    if B_co_var is None:
        # excluding any co-varience in the B spectrum is correct
        # if the points don't depend on the neighboring points
        Z_var = (A_wei**2*A_var + B_wei**2*B_var) 
    else:
        # When the B spectrum is resampled or interpolated this induces
        # co-varience between pixels. The effect is most strongly seen
        # in the resampling routine. In that case the co-varience is 
        # accounted for here
        cross_terms = B_wei*B_co_var*B_wei
        Z_var = A_wei**2*A_var + cross_terms

    # inverse varience 
    Z_inv_var = var_2_inv_var(Z_var)

    #>>>>>>> Return 
    return Z_wl, Z_data, Z_inv_var
