################################################################################
# These are the basic input functions for reading fits and text files
#
#
#
################################################################################
# import modules
if __name__ != '__main__':
    import eyeSpec
    from eyeSpec.base_functions import var_2_inv_var, np_vstack_append
    from eyeSpec.base_classes import query_fits_header, eyeSpec_spec
    from eyeSpec.dependencies import np, os, time, deepcopy, pdb
    from eyeSpec.dependencies import scipy, pyfits
    verbose = eyeSpec.Params['verbose_level']
    
pass
################################################################################
# initiate some parameters


pass
#################################################################################
# classes to store the coefficients and wavelength solution as well as to solve
class WavelengthSolutionFunctions:
    """
    A class with defined wavelenth solutions. These are functions which take points/pixels and converts them to wavelengths
    """

    def __init__ (self):
        ft = {'no solution'    : ["wl = np.arange(start_pt,start_pt+len(pts))",self.wl_soln_no_solution_prog],
              'none'           : ["same as 'no solution'",self.wl_soln_no_solution_prog],
              None             : ["same as 'no solution'",self.wl_soln_no_solution_prog],
              'pts'            : ["wl = pts",self.wl_soln_pts],
              'ordinary poly'  : ["wl = fxn(pts)  where fxn is an ordinary 7th degree polynomial",self.wl_soln_ordinary_poly],
              'poly'           : ["same as 'ordinary poly'",self.wl_soln_ordinary_poly],
              'chebyshev poly' : ["wl = fxn(pts) where fxn is a Chebyshev polynomial",self.wl_soln_chebyshev_poly],
              'spline3'        : ["wl = fxn(pts) where fxn is a cubic spline function",self.wl_soln_cubic_spline],
              'legrandre poly' : ["wl = fxn(pts) where fxn is a Legrandre polynomial",self.wl_soln_legrandre_poly],
              'linear'         : ["wl = fxn(pts) where fxn is a polynomial of order 1, uses coeff[0],coeff[1], defauts to coeff[1] = 1 (i.e. no dispersion)",self.wl_soln_linear],
              'log linear'     : ["wl = 10**(fxn(pts)) where fxn is the same polynomial described in 'linear'",self.wl_soln_log_linear_2_linear]}
        self._function_types = ft

    def __call__ (self, pts, coeff, function_type='no solution', default_function=True):
        """
        convert points to wavelength in Angstroms via the wavelength solution
        
        INPUT:
        pts : array, contains the points used to convert element-by-element to wavelength
        coeff : array len(coeff) > 2, contains the dispersion coefficients to use in the polynomial to the wavelength solution
        function_type :  can specify a specific solution to be applied
                    'no solution'    : wl = np.arange(disp[0],disp[0]+len(pts)), can use disp[0] to advance where the starting point is
                    'none' or None   : same as 'no solution'
                    'pts'            : wl = pts
                    'ordinary poly'  : wl = fxn(pts)  where fxn is an ordinary 7th degree polynomial
                    'poly'           : same as 'ordinary poly'
                    'chebyshev poly' : wl = fxn(pts) where fxn is a Chebyshev polynomial
                    'spline3'        : wl = fxn(pts) where fxn is a cubic spline function
                    'legrandre poly' : wl = fxn(pts) where fxn is a Legrandre polynomial
                    'linear'         : wl = fxn(pts) where fxn is a polynomial of order 1, uses coeff[0],coeff[1], defauts to coeff[1] = 1 (i.e. no dispersion)
                    'log linear'     : wl = 10**(fxn(pts)) where fxn is the same polynomial described in 'linear'
        
        default_function : (bool) If True it will use a default function if not of the above are found
        
        """
        ft = function_type
        if ft == 'no solution': coeff = [0,1]
                
        pts,coeff,no_solution = self._check_pts_coeff(pts,coeff)        
        if no_solution: ft = 'no solution'
        
        if ft not in self._function_types:
            if default_function: ft = 'no solution'
            else:raise ValueError("Unknown function type:"+str(function_type))
        
        if ft is 'pts': return self.wl_soln_pts(pts)
        
        func = self._function_types[ft][1]
        return func(pts,coeff)

    def get_function_types (self):
        return (self._function_types.keys())

    def get_func_name (self,name):
        if name not in self._function_types: raise ValueError("Unknown function type:"+name)
        else: return name

    def _check_pts_coeff (self,pts,coeff):
        no_solution = False
        try: pts = np.array(pts)
        except: no_solution = True
        try: coeff = np.array(coeff)
        except: no_solution = True
        
        if not no_solution:
            # check length of dispersion equation
            if len(coeff) < 2:
                print "WARNING: ARRAY coeff MUST BE len(coeff) >= 2, RUNNING coeff = array([0.,1.])"
                coeff = np.array([0.,1.])
                
            if len(coeff) < 8:
                coeff = np.concatenate((coeff,np.zeros(8-len(coeff))))
                
            # check the first order coefficient, if zero there is no solution
                if coeff[1] == 0: coeff_type = 'no solution'

        return pts,coeff,no_solution

    def wl_soln_no_solution_prog (self,pts,coeff):
        start_pt = 1

        try: coeff = np.array(coeff,dtype=int)
        except: coeff = None

        if str(type(coeff[0])).find('int') != -1:
            start_pt = int(coeff[0])

        return np.arange(start_pt,start_pt+len(pts))

    def wl_soln_no_solution (self,pts):
        return self.wl_soln_no_solution_prog(pts,None)

    def wl_soln_pts (self,pts):
        return pts

    def wl_soln_log_linear_2_linear (self,pts,coeff):
        return 10**(coeff[0] + coeff[1]*pts)
    
    def wl_soln_linear (self,pts,coeff):
        return coeff[0] + coeff[1]*pts
    
    def wl_soln_legrandre_poly (self,pts,coeff):
        xpt = (2.*pts - (len(pts)+1))/(len(pts)-1)
        wl = coeff[0] + coeff[1]*xpt
        return wl

    def wl_soln_cubic_spline (self,pts,coeff):
        print "!! DOUBLE CHECK SPLINE SOLUTION"
        # This matches what spectre gives but it seems like it give redundant solutions, i.e. all wl are the same
        s= (pts - 1)/(len(pts)-1) * coeff[7]
        J = np.array(s,dtype=int)
        a = (J+1) - s
        b = s - J
        z0 = a**3
        z1 = 1. + 3.*a*(1.+a*b)
        z2 = 1. + 3.*b*(1.+a*b)
        z3 = b**3
        c  = [coeff[x] for x in J]
        c1 = [coeff[x+1] for x in J]
        c2 = [coeff[x+2] for x in J]
        c3 = [coeff[x+3] for x in J]
        wl = c*z0 + c1*z1 + c2*z2 + c3*z3
        return wl 

    def wl_soln_chebyshev_poly (self,pts,coeff):
        if len(coeff) < 4: 
            raise ValueError("Number of coefficients insufficent for Chebyshev")
        #c20    p = (point - c(6))/c(7)
        #c      xpt = (2.*p-(c(9)+c(8)))/(c(9)-c(8))
        # !! is this right?
        xpt = (2.*pts - (len(pts)+1.))/(len(pts)-1.) 
        
        # xpt = (2.*point-real(npt+1))/real(npt-1)

        wl =  coeff[0] + xpt*coeff[1] + coeff[2]*(2.0*xpt**2.0-1.0) + coeff[3]*xpt*(4.0*xpt**2.0-3.0)+coeff[4]*(8.0*xpt**4.0-8.0*xpt**2.0+1.0)
        return wl

    def wl_soln_ordinary_poly (self,pts,coeff):
        # maximum degree polynomial will be determined
        # as degree 7
        degree = 7
        coeff = coeff[:degree]
        coeff = list(coeff)
        # must reverse my coeff
        coeff.reverse()
        return scipy.polyval(coeff,pts)
    
wlsolvefxn = WavelengthSolutionFunctions()

class WavelengthSolutionCoefficients:
    def __init__ (self, equ_type='none' , extra='none'):
        """
        equ_type correpsonds to the WavelengthSolutionFunctions
        coeffs correponds to the function but for a polynomial
             coeffs = c
             y = c[0]+c[1]*x+c[2]*x**2+c[3]*c**3+....
        extra     extra info from the extraction process
        """
        self.equ_type = self.set_equation_type(equ_type)
        self.coeffs = []
        self.extra = str(extra)
    
    def __len__ (self):
        return len(self.coeffs)
    
    def add_coeffs (self,coeff,ordi=None):
        if ordi is None: self.coeffs.append(coeff)
        else: self.coeff[ordi] = coeff
        
    def get_coeffs (self,ordi=None):
        index = ordi - 1
        if index not in xrange(len(self.coeffs)): index = len(self.coeffs)-1
        return self.coeffs[index]
    
    def get_equation_type (self):
        return deepcopy(self.equ_type)
    
    def set_equation_type (self,equ_type):
        self.equ_type = wlsolvefxn.get_func_name(equ_type)
    
    def get_coeffsicients (self, ordi, reverse=False):
        output = deepcopy(self.coeffs)
        
    def get_extra_info (self):
        return deepcopy(self.extra)

WSC = WavelengthSolutionCoefficients

pass
#################################################################################
# misc functions which could be used
def pts_2_phys_pixels (pts,bzero=1,bscale=1):
    """
    convert points to wavelength in Angstroms via the wavelength solution
    INPUT:
    pts : array, contains the points used to convert element-by-element to wavelength

    bzero : float, from the header file which gives the starting point for the physical pixels
    bscale : float, from the header file which gives the scaling for the physical pixels
            pts = bzero + pts*bscale


    bzero = query_fits_header(prihdr,'BZERO',noval=1) # for scaled integer data, here is the zero point
    bscale = query_fits_header(prihdr,'BSCALE',noval=0) # for scaled integer data, here is the multiplier
   
    I'm pretty sure pyfits uses these when it reads in the data so I don't need to
   
    """
    if bzero != 0 or bscale !=1:
        print "Whoops, I don't know exactly what to do with bzero!=1 or bscale!=0 :<",bzero,"><",bscale,">"
        #pts +=1

    pts = bzero + pts*bscale

    # should return [1,2,3,......,#pts]
    return pts

def check_for_txt_format (filename,**np_kwargs):
    try: txt_data = np.loadtxt(filename,unpack=True,dtype=float,**np_kwargs)
    except: return False, None
    return True, txt_data

def _check_header_type (pyfits_header):
    if repr(type(pyfits_header)) != "<class 'pyfits.header.Header'>": 
        raise ValueError("pyfits_header input not of type <class 'pyfits.header.Header'> :"+repr(type(pyfits_header)))

pass
#################################################################################
# functions to extract wavelength solution coefficients and equation type from
# a pyfits header

def coeff_from_ctype1 (pyfits_header):
    _check_header_type(pyfits_header)
    wlcoeff = WSC()
    #==========================================================================#

    ctype1 = query_fits_header(pyfits_header,'CTYPE1',noval='')
    crval1 = query_fits_header(pyfits_header,'CRVAL1',noval=0) # for linear dispersions, the starting wavelength
    crpix1 = query_fits_header(pyfits_header,'CRPIX1',noval=0) # for linear dispersions, the pixle to which CRVAL1 refers
    cdelt1 = query_fits_header(pyfits_header,'CDELT1',noval=0) # for linear dispersion, here is the dispersion
    
    start_pix = 1.0 - crpix1.val # this is because I start the pixel counting at 1 later

    if (ctype1.found and ctype1.val == 'LINEAR'):
        if crval1.found and crpix1.found and cdelt1.found:
            coeff = np.array([crval1.val + start_pix*cdelt1.val, cdelt1.val])
            wlcoeff.extra = 'used header to get crval1 and cdelt1, to apply wl = crval1 + cdelt1*pts'
            wlcoeff.set_equation_type('linear')            
            wlcoeff.add_coeffs(coeff)

    if (ctype1.found and ctype1.val == 'LOG-LINEAR'):
        if crval1.found and crpix1.found and cdelt1.found:
            coeff = [crval1.val + start_pix*cdelt1.val,cdelt1.val]
            wlcoeff.extra = 'LOG-LINEAR, used header to get crval1 and cdelt1, to apply wl = 10**(crval1 + cdelt1*pts)'
            wlcoeff.set_equation_type('log linear')            
            wlcoeff.add_coeffs(coeff)
              
    return wlcoeff
         
def coeff_from_crvl (pyfits_header):
    _check_header_type(pyfits_header)
    wlcoeff = WSC()

    def get_for_order (ordi):
        ordi = format(ordi,'02')
        #==========================================================================#
        linintrp = query_fits_header(pyfits_header,'LININTRP') # string with infor about the type of linear interpretation
        crvl1_ = query_fits_header(pyfits_header,'CRVL1_'+ordi,noval=1) # for linear dispersions, the starting wavelength
        cdlt1_ = query_fits_header(pyfits_header,'CDLT1_'+ordi,noval=0) # for linear dispersions, the pixle to which 
        
        if crvl1_.found and cdlt1_.found:
            if linintrp.found and linintrp.val.find('linear') == -1: print "WARNING: KEYWORD LININTRP HAS NO REFERENCE TO 'linear' BUT PERFORMING A LINEAR DISPERSION"
            coeff = [crvl1_.val, cdlt1_.val]
            wlcoeff.extra = 'used header to get crvl1_ and cdlt1_ depending on order, to apply wl = crvl1_? + cdlt1_?*pts'
            wlcoeff.set_equation_type('linear')
            wlcoeff.add_coeffs(coeff)
            return False
        else: 
            if ordi == 0: return False
            return True
            
    for i in xrange(100):
        if get_for_order(i): break
        
    if len(wlcoeff) > 0: print "!! FIRST TIME WITH CRVL AND CRDEL, CHECK THE OUTPUT OF WAVELENGTH VS FLUX"
        
    return wlcoeff

def coeff_from_wcs (pyfits_header, apply_rv=False):
    _check_header_type(pyfits_header)
    wlcoeff = WSC()

    #==========================================================================#
    wat0_001 = query_fits_header(pyfits_header,'WAT0_001')
    wat1_001 = query_fits_header(pyfits_header,'WAT1_001') 
    wat2_001 = query_fits_header(pyfits_header,'WAT2_001') 
    wat3_001 = query_fits_header(pyfits_header,'WAT3_001') 

    if not wat0_001.found: return wlcoeff
    # found the wat0_001 keyword
    
    wat0 = wat0_001.val.split("=")
    if wat0[0].lower() != 'system': raise ValueError("IRAF WCS, Unknown keyword in WAT0_001 = "+wat0[0])

    wcs_system = wat0[1]
    if wcs_system not in ['physical','multispec']: 
        # print "IRAF WCS: Unknown system given in WAT0_001 = "+wcs_system
        return wlcoeff
    
    # can't be 'world' or 'equispec' because those are for 2D data not 1D
    if wat3_001.found: print "HeadsUp: CHECK OUT WAVELENGTH DISPERSION, FOUND KEYWORD 'WAT3_001' WHICH I DON'T KNOW HOW TO DEAL WITH" 
    # similarly, wat3 corresponds to higher dimensions than expected
    
    if not wat1_001.found: return wlcoeff
        
    #========== now it is on system mutlispec
    wat1_001_dict = {}
    for val in pyfits_header['wat1_001'].split():
        val = val.split("=")
        wat1_001_dict[val[0]] = val[1]

    # WAT1_001= 'wtype=linear label=Wavelength units=Angstroms' 
    # WAT1_001= 'wtype=multispec label=Wavelength units=angstroms'
    
    # check wtype:
    if wat1_001_dict['wtype'] == 'linear': raise ValueError("IRAF WCS keyword WAT1_001 has wtype = linear, I don't know how to deal with this")
    elif wat1_001_dict['wtype'] != 'multispec': raise ValueError("Unknown value for WAT1_001 wtype:"+wat1_001_dict['wtype'])
   
    # check wavelength
    if wat1_001_dict['label'].lower() != 'wavelength': raise ValueError("IRAF WCS for WAT1_001 keyword I expected label=wavelength but got label="+str(wat1_001_dict['label']))
    
    # check units
    if wat1_001_dict['units'].lower() != 'angstroms': raise ValueError("ERROR: DON'T KNOW HOW TO HANDLE IRAF WAT1_001 FOR UNITS GIVEN"+wat1_001_dict['units'])
 
 
    #======== now has 'wtype=multispec label=Wavelength units=angstroms'        
    
    #        ordi = format(ordi,'02')
    #        wat2_0 = query_fits_header(pyfits_header,'WAT2_0'+ordi,'02'))
    #        if not wat2_0.found: return True
    #    

    def unpack_WCS_multispec (root_keyword):     
        # creates order_disp which is a dictionary that has the dispersion for a given order, this would be inefficient to make the whole thing every time when you only need one value every time it's run. But with the way IRAF layed out these WAT2 stuff and how my code is written it's sort of necessary
        wat_str = ''
        for key in pyfits_header.keys():
            if key.find(root_keyword) == 0: wat_str += format(pyfits_header[key],'68')
        
        cut_str =  wat_str[:30].split()
        if cut_str[0].lower() != 'wtype=multispec':
            raise ValueError("IRAF WCS: For root keyword "+root_keyword+" didn't get wtype=multispec :"+wat_str[:30])
            

        # separate lines by spec 
        sep_spec = wat_str[15:].split("spec")
        order_coeff = {}
        def _check_string_for_splits (dat,k,ret = 'prefix'):
            """
            This is used to break up values like:
            4.21234D+50.234    =>  4.21234E+5  and  0.234
            -2.12345E-40.234  => -2.12345E-4  and  0.234
            2.3451650.456     => -2.345165    and  0.456
            2.123-32.123      =>  2.123       and  -32.123
            2.124E-4-32.123   =>  2.123E-4    and  -32.123

            can't deal with:
            +2.345   and   +23.4      => +2.345+23.4
            2.123    and   +32.12     =>  32.123+32.123
            
            sdat is the line split by '.'

            ret = 'prefix' or 'suffix'
            """
            sdat = dat.split('.')
            #>> deal with the first case, assume sdat[0] is fine
            # check for scientific notation
            sdat[k].upper()
            if sdat[k].find('D') != -1: sdat[k] = sdat[k].replace('D','E') # change D to E
            if sdat[k].find('+') != -1 and sdat[k].find('E+') == -1: print "I didn't code anything for '+', check the output. ORDER:"+cur_order+" problem line:"+dat

            checkfor = ['E-','E','E+']
            found_checkfor = False
            for cval in checkfor:
                if sdat[k].upper().find(cval) != -1: 
                    found_checkfor = True
                    ssdat = sdat[k].split(cval) # assume ssdat[0] if fine (e.g. 2345E-40 => ['2345','40])
                    if cval == 'E': cval = 'E+'
                    ss_dat = ssdat[1].split('-')
                    if len(ss_dat) == 1: # (e.g. '-23.2345E-40.23' =1> '2345E-40' => ['2345','40] => [40])
                        if ss_dat[0][-1] != '0' and ret == 'prefix' and k!=0: raise ValueError("ORDER:"+cur_order+" I expected a zero, bad input:"+dat)
                        else:
                            suffix = ssdat[0]+cval+ss_dat[0][:-1] # (e.g. 2345E-40 => 2345E-4)  !! what if I wanted 2345E-10
                            prefix = '0'
                    elif len(ss_dat) == 2: # (e.g. 2345E-4-34 => ['2345','4-34] => ['4','34'])
                        suffix = ssdat[0]+cval+ss_dat[0]
                        prefix = "-"+ss_dat[1]
                    else: # (e.g. 2345E-4-3-24 => ['2345','4-3-24] => ['4','3','24'])
                        raise ValueError("ORDER:"+cur_order+" unexpected value has multiple dashes/negative signs: "+dat)
                    break
            # if it's not scientific notation
            if not found_checkfor:
                ss_dat = sdat[k].split('-')
                if len(ss_dat) == 1: # (e.g. '234540' => ['234540'])
                    if ss_dat[0][-1] != '0' and ret == 'prefix' and k!=0: raise ValueError("ORDER:"+cur_order+" I expected a zero, bad input:"+ss_dat[0])
                    else:
                        suffix = sdat[k] # (e.g. 234540 => 234540) 
                        prefix = '0'
                elif len(ss_dat) == 2: # (e.g. '2345-34' => ['2345','34])
                    suffix = ss_dat[0]
                    prefix = "-"+ss_dat[1]
                else: # (e.g. 2345-3-24 => ['2345','3','24'])
                    raise ValueError("ORDER:"+cur_order+" unexpected value has multiple dashes/negative signs: "+dat)
            if ret == 'prefix': return prefix
            elif ret == 'suffix': return suffix
            else: raise TypeError("Whoops, I expected either 'prefix' or 'suffix'")

        for val in sep_spec:
            # split by the equals order_num = stuff
            sep_spec2 = val.split("=")
            
            is_problem_child = False
            debug_problem_children = False
            if len(sep_spec2)>1:
                cur_order = str(sep_spec2[0])
                sep_spec_val = sep_spec2[1].replace('"','').split()
                # go through and make sure no values have extra decimal (e.g. -0.0040.234) which happened for one fits file when the wavelength solution was given twice
                for i in range(len(sep_spec_val)):
                    dat = deepcopy(sep_spec_val[i])

                    sdat = dat.split(".")
                    new_entries = []
                    if len(sdat) > 2:
                        if debug_problem_children: print "=== "+cur_order+" ===== problem child ===>",dat
                        is_problem_child = True
                        #prefix,suffix = _check_string_for_splits (dat,1) 
                        #new_entries.append(sdat[0]+"."+suffix)

                        for j in range(1,len(sdat)): # len>2
                            prefix = _check_string_for_splits (dat,j-1,ret='prefix')
                            suffix = _check_string_for_splits (dat,j,ret='suffix')
                            new_entries.append(prefix+"."+suffix)

                    # now stick in new_entries
                    if len(new_entries) != 0:
                        new_entries.reverse()
                        del sep_spec_val[i]
                        for new_entry in new_entries:
                            sep_spec_val.insert(i,new_entry) 

                if is_problem_child and debug_problem_children:
                    for i in range(len(sep_spec_val)):
                        print i,">>",sep_spec_val[i]
                    skip = raw_input("")
                    if skip == 'a': stop,here=0

                order_coeff[int(sep_spec2[0])] = np.array(sep_spec_val,dtype=float)
        return order_coeff, wat_str
    
    order_coeff, wat2_str = unpack_WCS_multispec('WAT2_')
    # REFERNCE: http://iraf.net/irafdocs/specwcs.php
    # The dispersion functions are specified by attribute strings with the identifier specN where N is the physical image line. The attribute strings contain a series of numeric fields. The fields are indicated symbolically as follows. 
    #     specN = ap beam dtype w1 dw nw z aplow aphigh [functions_i]
    # 
    # order_coeff[?][0] = aperture number
    # order_coeff[?][1] = beam number
    # order_coeff[?][2] = dispersion type, dcflag = -1 no disp, 0 linear, 1 log-linear, 2 nonlinear
    # order_coeff[?][3] = c0, first physical pixel
    # order_coeff[?][4] = c1, average disperasion interval
    # order_coeff[?][5] = npts, number valid pixels !! could have problem if different orders have different lengths
    # order_coeff[?][6] = rv,z, applies to all dispersions coordinates by multiplying 1/(1+z)
    # order_coeff[?][7] = aplow, lower limit of aperture
    # order_coeff[?][8] = aphigh, upper limit of aperture
    # order_coeff[?][9] = N/A for linear or log-linear
    # ----------------
    #            function_i =  wt_i w0_i ftype_i [parameters] [coefficients]
    # order_coeff[?][9]  = wieght wt_i 
    # order_coeff[?][10] = zeropoint offset w0_i
    # order_coeff[?][11] = type dispersion fxn = 1 cheby, 2 legrandre, 3 cubic spline3, 
    #                        4 linear spline, 5 pixel coordinate array, 6 sampled coordinate array
    # order_coeff[?][12+] = [parameters...]
    # order_coeff[?][12++] = [coefficients...]

    if 1 not in order_coeff: raise ValueError("IRAF WCS: I didn't find a spec1:"+wat2_str[:40])
    dcflag = order_coeff[1][2]
    equ_type = 'none'
    
    LTV_dimn = [1,1]
  
    def get_order_wlsoln (ordi):
        if ordi not in order_coeff: 
            if i == 0: return False
            return True
        # w = sum from i=1 to nfunc {wt_i * (w0_i + W_i(p)) / (1 + z)}
        if order_coeff[ordi][2] != dcflag: raise ValueError("Encountered orders with different functions dcflag1="+str(dcflag)+"  dcflag"+str(ordi)+"="+str(order_coeff[ordi][2]))
        z = order_coeff[ordi][6]
        if z != 0 and apply_WCS_rv: print "HeadsUp: Found a radial velcoity shift in IRAF WCS, applying Vrad = ",1./(1.+z)," to order "+str(ordi)
        else: z = 0
    
        if dcflag < 0:
            equ_type='none'
            coeff = [0,1]
        
        elif dcflag == 0:
            equ_type = 'linear'
            coeff = [order_coeff[ordi][3],
                     order_coeff[ordi][4]]
            
        elif dcflag == 1 or dcflag == 2:
            polytype = order_coeff[ordi][11]
            ltv = query_fits_header(pyfits_header,'LTV'+str(LTV_dimn[0]),noval=0) # IRAF auxiliary wavelenth solution parameters 
            ltm = query_fits_header(pyfits_header,'LTM'+str(LTV_dimn[0])+'_'+str(LTV_dimn[0]),noval=0) # IRAF auxiliary wavelenth solution parameters 
        
            if (ltv.found or ltm.found): print ("IRAF WCS: found WAT keywords with system=multispec and found dcflag = "+str(dcflag)+" for order "+str(ordi)+" and LTV and LTM keywords which I don't know what to do with")
        
            if polytype == 1:
                # order_coeff[?][11] = type dispersion fxn = 1 cheby
                # order_coeff[?][12] = order
                # order_coeff[?][13] = xxmin
                # order_coeff[?][14] = xxmax
                # order_coeff[?][15:15+order] = [coefficients...]
        
                # !! the points used to calculate the wavelength solution: np.arange(xxmin,xxmax)
                xxmin = order_coeff[ordi][13]
                if xxmin != 1: print "WARNING: From WCS I got a xxmin of",xxmin,"but hardwired in is a value of 1"
        
                equ_type = 'chebyshev poly'
                poly_order = order_coeff[ordi][12]
                xxmin,xxmax = order_coeff[ordi][13:15]
                coeff = order_coeff[ordi][15:15+poly_order] # I think
        
                #coeff[5] = ltv.val # from SPECTRE == c(6)
                #coeff[6] = ltm.val # ==c(7)          
                #c*****a Chebyshev polynomial solution
                #c20    p = (point - c(6))/c(7)
                #c      xpt = (2.*p-(c(9)+c(8)))/(c(9)-c(8))


        elif polytype == 2:
            equ_type = 'legrandre poly'
            print "WARNING: NOT SURE IF I'M GETTING THE LEGRANDRE COEFFICIENTS FROM DATA CORRECTLY"
            coeff = order_coeff[ordi][14:] # I think

            #coeff[5] = ltv.val  # from SPECTRE
            #coeff[6] = ltm.val                                

        elif polytype == 3: 
            equ_type = 'spline3'
            print "WARNING: NOT SURE IF I'M GETTING THE CUBIC SPLINE COEFFICIENTS FROM DATA CORRECTLY"
            coeff = order_coeff[15:]
            
        elif polytype == 4:
            print "WARNING: I CAN'T CURRENTLY COMPUTE A LINEAR SPLINE, ASSUMING LINEAR DISPERSION"
            equ_type = 'linear'
            coeff = [0,1]
            coeff = [order_coeff[ordi][3],
                     order_coeff[ordi][4]]
        
        elif polytype > 5:
            print "ERROR: NOT SET UP TO HANDLE THIS TYPE OF DISPERSION",polytype # !! error exit

        # apply the radial velocity shift given in the WCS
        coeff = np.array(coeff)*(1.0/(1.0+z))
        
        wlcoeff.set_equation_type(equ_type)
        wlcoeff.add_coeffs(coeff)
        return False

    for i in xrange(100):
        if get_order_wlsoln(i): break
        
    wlcoeff.extra = ['used header to get parameters and coefficients, function: '+equ_type+', to apply wl = function(pts)',order_coeff]
    return wlcoeff

def coeff_from_w0 (pyfits_header):
    _check_header_type(pyfits_header)
    wlcoeff = WSC()

    w0 = query_fits_header(pyfits_header,'W0',noval=0) # for linear dispersions, starting wavelength
    wpc = query_fits_header(pyfits_header,'WPC',noval=0) # for linear dispersion, here is th dispersion
    if w0.found and wpc.found:
        coeff = [w0,wpc]
        wlcoeff.extra = 'used header to get W0 and WPC, to apply wl = W0 + WPC*pts'
        wlcoeff.set_equation_type('linear')
        wlcoeff.add_coeffs(coeff)
        print "!! FIRST TIME WITH W0 AND WPC, CHECK THE OUTPUT OF WAVELENGTH VS FLUX"
    return wlcoeff

def coeff_basic_linear (pyfits_header,wlsol=True):
    _check_header_type(pyfits_header)
    wlcoeff = WSC()
    #==========================================================================#
    if wlsol:
        equ_type = 'linear'
        extra_data = 'no header info used, just a wl = pts'
    else: 
        equ_type = 'none'
        extra_data = 'following checks found that the first order coefficient is zero, setting basic linear wl = pts'

    wlcoeff.extra = extra_data
    wlcoeff.set_equation_type(equ_type)
    wlcoeff.add_coeffs([0,1])
    return wlcoeff

def coeff_from_makee_wv (pyfits_header):
    _check_header_type(pyfits_header)
    #==========================================================================#
    wlcoeff = WSC()

    def get_for_order (order_id):
        WV_0_ = query_fits_header(pyfits_header,'WV_0_'+format(int(order_id),'02')) # first 4 coefficients
        WV_4_ = query_fits_header(pyfits_header,'WV_4_'+format(int(order_id),'02')) # second 4 coefficients
        
        if not WV_0_.found: return True
        
        coeff = np.zeros(8)
        coeff[:4] = np.array(WV_0_.val.split(),dtype=float)
        if WV_4_.found: coeff[4:8] = np.array(WV_4_.val.split(),dtype=float)
        else: raise IOError("Expected to get "+'WV_4_'+format(int(order_id),'02')+" along with keyword "+'WV_0_'+format(int(order_id),'02'))

        # !! could put a check to make sure WV_4_ follows WV_0_ in the header
        wlcoeff.extra = 'used header to get MAKEE keywords WV_0_? and WV_4_? depending on order, to apply polynomial coefficients given by WV_0_? and WV_4_?'
        wlcoeff.set_equation_type('poly')
        wlcoeff.add_coeffs(coeff)
        return False
            

    for i in xrange(1,100):
        if get_for_order(i): break
    
    #==========================================================================#   
    return wlcoeff

def coeff_from_makee_c0 (pyfits_header):
    _check_header_type(pyfits_header)
    wlcoeff = WSC()
    #==========================================================================#
    def get_from_order (order_id):
        CO_0_ = query_fits_header(pyfits_header,'CO_0_'+format(int(order_id),'02')) # first 4 coefficients
        CO_4_ = query_fits_header(pyfits_header,'CO_4_'+format(int(order_id),'02')) # second 4 coefficients
            
        if CO_0_.found or CO_4_.found: 
            print "WARNING: KEYWORDS",'CO_0_'+format(int(order_id),'02'),"AND",'CO_4_'+format(int(order_id),'02'),"FOUND BUT I DON'T KNOW WHAT TO DO WITH THEM"
            
        if CO_0_.found and CO_4_.found:
            coeff = np.zeros(8)
            coeff[:4]  = np.array(CO_0_.val.split(),dtype=float)
            coeff[4:8] = np.array(CO_4_.val.split(),dtype=float)
            wlcoeff.add_coeffs(coeff)
            return False    
        return True
    
    for i in xrange(1,100):
        if get_from_order(i): break
    
    if len(wlcoeff) > 0: 
        wlcoeff.set_equation_type('none')
        wlcoeff.extra = 'coefficients from C0_0_? and C0_4_? makee pipeline keywords'  
         
    return wlcoeff

def coeff_from_SPECTRE (pyfits_header,get_hist=False):
    _check_header_type(pyfits_header)
    wlcoeff = WSC()
    #==========================================================================#
    
    if not query_fits_header(pyfits_header,'HISTORY',noval=0).found: return wlcoeff # old SPECTRE-stype dispersion information
    history_lines = []
    spectre_history = {}
    kount = 0
    for i in range(len(pyfits_header.ascardlist())):
        line=str(pyfits_header.ascardlist()[i]).rstrip()
        if line[0:8] == 'HISTORY ':
            history_lines.append(line)
            
            if line[23:28] == 'DISP=':
                print "!! NEED TO RESOLVE HOW TO READ THIS TYPE OF DISPERSION FROM SPECTRE"                    
#                   read (head(k+1:k+80),1022) (disp(i),i=1,4)
# 1022              format(28x,1p4e13.5)

            if line[19:26] == 'D1,2,3:':
                sline = line.split(":")
                day = int(sline[0].split()[1])
                month = int(sline[1])
                year = int(sline[2].split()[0])
                timetag = time.mktime((year,month,day,0,0,kount,0,0,0))
                kount += 1
        
                coeff = np.zeros(9)
                line = line.replace('D','e')
                coeff[:3] = np.array([line[26:44],line[44:62],line[62:80]],dtype=float)

                line2= str(pyfits_header.ascardlist()[i+1]).rstrip().replace('D','e')
                if line2[0:8] != 'HISTORY ': raise IOError('EXPECTED NEXT LINE TO HAVE TAG HISTORY')
                coeff[3:6] = np.array([line2[26:44],line2[44:62],line2[62:80]],dtype=float)

                line3= str(pyfits_header.ascardlist()[i+2]).rstrip().replace('D','e')
                if line3[0:8] != 'HISTORY ': raise IOError('EXPECTED NEXT LINE TO HAVE TAG HISTORY')
                disp_info = np.array([line3[26:44],line3[44:62],line3[62:80]],dtype=float)

                # cheby poly may need disp_info[0]
                # disp_info[0] == c(7)
                # disp_info[1] == c(8)
                # from SPECTRE: 
                # c20    p = (point - c(6))/c(7)
                # c      xpt = (2.*p-(c(9)+c(8)))/(c(9)-c(8))
                if disp_info[2] == 1: disp_type = 'chebyshev poly'
                elif disp_info[2] == 2: disp_type = 'legrendre poly'
                elif disp_info[2] == 3:
                    print "WARNING: check the output, I think I may need to use disp_info[1] "
                    # from SPECTRE: s = (point-1.)/(real(npt)-1.)*c(8)
                    # c(8) == disp_info[1] ==> true
                    disp_type = 'spline3'
                else: disp_type = 'poly' # but likely the orders > 1 have coefficients zero
                extra_data= ['used header to get SPECRE HISTORY tags, function:'+disp_type+', to apply wl=function(pts)',[line,line2,line3]]
                
                spectre_history[timetag] = (line,extra_data,disp_type,coeff)
    
    if get_hist: return spectre_history,history_lines
                
    # if you found history lines use one with the most recent data                
    if len(spectre_history) > 0:
        most_recent = sorted(spectre_history.keys())[-1]
            
        extra_data, disp_type, coeff = spectre_history[most_recent][1:]
        
        wlcoeff.extra = extra_data
        wlcoeff.set_equation_type(disp_type)
        wlcoeff.add_coeffs(coeff)
            
            
    # NOTE: spectre_history has all the history tags for spectre keyed by time stamp                          
    return wlcoeff



################################################################

pass
#################################################################################
# Functions to do all the coefficient types above and then resolve which
# are more important than others

def wlsoln_coeff_from_header (pyfits_header, apply_WCS_rv=False, preferred=None, output='sel'):
    """
    This uses the header and tries out all possible functions for getting the wavelength solution coefficients
    
    returns a dictionary with keywords in ['linear','ctype','crvl','wcs','w0','co_0','wv_0','spectre']
    """
    # coefficient choices
    cc = {}
    #========================================================================#
    # linear dispersion
    cc['linear'] = coeff_basic_linear(pyfits_header)

    #========================================================================#
    # using keywords ctype, crval, crpix, cdelt
    cc['ctype1'] = coeff_from_ctype1(pyfits_header)

    #========================================================================#
    # linear dispersion using keywords linintrp, crvl1_?, cdlt1_?
    # from IRAF, order by order  !! do I need to look up what the 1_ means?
    # some of these are doubled by WAT0_001 stuff
    cc['crvl'] = coeff_from_crvl(pyfits_header)
    # if preferred_disp == 'any' or preferred_disp == 'linear' or preferred_disp == 'crvl' or preferred_disp == 'makee linear':
  
    #========================================================================#
    # IRAF WCS keywords WAT?_001 
    #if preferred_disp == 'any' or preferred_disp == 'IRAF_WCS':
    cc['wcs'] = coeff_from_wcs(pyfits_header,apply_WCS_rv)

    #========================================================================#
    # linear dispersion for keywords w0 and wpc
    cc['w0'] = coeff_from_w0(pyfits_header)
    #if preferred_disp == 'any' or preferred_disp == 'linear' or preferred_disp == 'w0':

    #========================================================================#
    # MAKEE type dispersion using keywords co_0_? and co_4_?
    # I'm not sure what type of coefficients these are !!
    #cc['co_0'] = coeff_from_makee_c0(pyfits_header)
#     if preferred_disp == 'any' or preferred_disp == 'makee' or preferred_disp == 'co_0':

    #========================================================================#
    # MAKEE coeffificients using keywords wv_0_? and wv_4_?
    cc['wv_0'] = coeff_from_makee_wv(pyfits_header)
    #if preferred_disp == 'any' or preferred_disp == 'makee' or preferred_disp == 'wv_0':

    #========================================================================#
    # spectre type dispersion
    cc['spectre'] = coeff_from_SPECTRE(pyfits_header)
    #if preferred_disp == 'any' or preferred_disp == 'spectre':

    #========================================================================#
    #========================================================================#
    
    if output == 'all': return cc
    return resolve_wlsoln_coeffs(cc,preferred)
    
def resolve_wlsoln_coeffs (wlcoeffs,preferred=None):
    res_order = ['spectre',
                 'wv_0',
                 #'co_0',
                 'w0',
                 'wcs',
                 'crvl',
                 'ctype1',
                 'linear',
                 'no solution']
    
    found_soln = False
    
    if preferred is not None: 
        if preferred not in res_order: raise ValueError("Given preference for wavelength solution is not in : "+", ".join(res_order))
        elif preferred not in wlcoeffs: 
            print "HeadsUp: Preferred dispersion not in options contining through list: " +", ".join(res_order)
        else: return wlcoeff[preferred]
    
    
    for coeff_type in res_order:
        if coeff_type in wlcoeffs:
            if len(wlcoeffs[coeff_type]) == 0: continue
            return wlcoeffs[coeff_type]
        else: print "key '"+str(coeff_type)+"' not in list:"+str(tuple(res_order))
            
    return None


pass
#################################################################################
# these functions are specific versions for reading in text files

def readin_multi_order_txt (filename,comment='#',skiprows=0,delimiter=None,progressive_pixel = True):
    
    comment = str(comment)
    skiprows = int(skiprows)
    delimiter = str(delimiter)
    
    f = open(filename)
    lines = f.readlines()
    f.close()
    
    eSmo = False # eyeSpec multiple order text file
    txt_data = {}
    which_order = 0
    pixel_val = 1
    
    print_once = [True,True]
    
    for i in range(len(lines)):
        line = lines[i].rstrip()
        # check to see if its a multiple order file
        if i == 0 and line.find('eSorders') != -1:
            eSmo = True
            continue
        
        # skip number of rows 
        if i < skiprows: continue 
        
        # ignore comments
        ci = line.find(comment)
        if ci > -1:
            oi = line.find("order")
            if eSmo and oi != -1:
                sline = line.split()
                which_order = int(sline[2]) 
            elif oi == 0: continue
            else: 
                line = line[:ci]
                while ci > -1:
                    ci = line.find(comment)
                    line = line[:ci]
                        
        # now should just be data
        sline = line.split(delimiter)
        var = 1e30
        if len(sline) == 0: continue
        elif len(sline) == 1:
            if print_once[0]:
                print "only found one data column assuming wl is index"
                print_once[0] = False
            if progressive_pixel:
                pixel_val += 1
            else:   
                if which_order not in txt_data:
                    pixel_val = 1
                else: 
                    pixel_val = len(txt_data[which_order])
        
            wl = float(pixel_val)
            data = float(sline[0])
        elif len(sline) == 2:
            wl = float(sline[0])
            data = float(sline[1])
            
        elif len(sline) >= 3:
            if len(sline) != 3 and print_once[1]:
                print "number of columns is greater than 3, only using the first three for wavelength,data,varience"
                print_once[1] = False
            wl = float(sline[0])
            data = float(sline[1])
            var = float(sline[2])
            
        inv_var = var_2_inv_var([var])[0]
        if which_order in txt_data:
            txt_data[which_order].append([wl,data,inv_var])
        else:
            txt_data[which_order] = [[wl,data,inv_var]]
        
        
    all_orders = np.array(txt_data.keys()).sort()
    np_txt_data = []
    for ori in all_orders:
        
        # !! need some way to deal with multiple order text files
        
        # check the input txt_data
        if txt_data.shape[0] == 1:
            print "HeadsUp: NO WAVELENGTH DATA FOUND, USING FIRST COLUMN AS DATA"
            txt_data = np.vstack((np.arange(1,len(txt_data)+1),txt_data))
        #elif txt_data.shape[0] == 2: perfect!
        elif txt_data.shape[0] > 2: 
            print "HeadsUp: MORE COLUMNS FOUND IN TEXT FILE THAN I KNOW WHAT TO DO WITH, TAKING THE 1ST TO BE WAVELENGTH AND 2ND TO BE DATA"
            txt_data = txt_data[:2]
            # !! I have a previous note about this but I could take a third column to be inverse varience, or varience
            # !! could create a system which has three columns (wavelenght,data,inv_varience)

        # create new header
        new_header = pyfits.PrimaryHDU()
        hdulist = pyfits.HDUList([new_header])
        prihdr = hdulist[0].header

        # add/adjust fits header
        prihdr['NAXIS'] = 1 # !! should I have a way to recognize multi order text files?
        prihdr.update('FILENAME',os.path.basename(filename))# !!,comment=' ')

        wl = np.array([[txt_data[0]]])
        data = np.array([[txt_data[1]]])
        # inv_var = np.array([[txt_data[2]]])

        spec_obj = eyeSpec_spec(wl,data,header=prihdr)#!!,inverse_varience = inv_var)
        spec_obj.filename = os.path.abspath(filename)

        # set up private information
        spec_obj._private_info['filename'] = filename

        pts = np.arange(len(data[0][0]))+1. 
        
        # find  wavelength solution
        # !! I should have something here for when their's just a single column which is data, i.e. no wavelength solution
        if preferred_disp == 'linear':
            poly_order = 1
            disp_type = 'linear'
        elif preferred_disp == 'any' or preferred_disp in ['poly','ordinary poly']:
            poly_order = 6
            disp_type = 'poly'
          
        else:
            print "WARNING: CURRENTLY CAN ONLY APPLY LINEAR AND POLYNOMIAL FITS TO TEXT DATA, USING POLYNOMIAL"
            # !! when I create a fancy get_disp for a wl array then I can put that in here
            poly_order = 6
            disp_type = 'poly'

        disp = scipy.polyfit(pts,wl[0][0],poly_order)
        disp = np.flipud(disp) # switch the order which the dispersion is given [0] = zero-ith order
        disp = np.concatenate((disp,np.array([0.0]))) # array length = 8

        spec_obj._private_info[(0,0)]['disp']=[disp,disp_type]
        spec_obj._private_info[(0,0)]['rv'] = [0]

def readin_txt (filename,txt_data=None,get_data=False,**np_kwargs):
    if txt_data is None:
        if 'unpack' in np_kwargs: del np_kwargs['unpack']
        if 'dtype' in np_kwargs: del np_kwargs['dtype']
        txt_data = np.loadtxt(filename,unpack=True,dtype=float,**np_kwargs)
    
    # check the input txt_data
    if txt_data.ndim == 1:
        print "HeadsUp: NO WAVELENGTH DATA FOUND, USING FIRST COLUMN AS DATA"
        data = txt_data 
        wls = np.arange(len(data))+1
        var = np.ones(len(data))
    elif txt_data.ndim == 2:
        wls = txt_data[0]
        data = txt_data[1]
        var = np.ones(len(data))
    elif txt_data.ndim == 3: wls,data,var = txt_data
    elif txt_data.shape[0] > 2: 
        print "HeadsUp: Found more than 3 columns in text file '"+filename+"' taking the first three to be wavelength, data, variance respectively"
        wls,data,var = txt_data[:3]
      
    inv_var = var_2_inv_var(var)
          
    # create new header
    new_header = pyfits.PrimaryHDU()
    hdulist = pyfits.HDUList([new_header])
    prihdr = hdulist[0].header
    
    # add/adjust fits header
    prihdr['NAXIS'] = 1 
    prihdr.update('FILENAME',os.path.basename(filename))# !!,comment=' ')
    prihdr.update('NAXIS1',len(wls)) 
    
    if get_data:
        return (wls,data,inv_var)
    
        
    spec_obj = eyeSpec_spec(wls,data,inv_var,header=prihdr)
    spec_obj.filename = os.path.abspath(filename)
    
    # set up private information
    spec_obj._private_info['filename'] = filename
    
    return spec_obj

#    pts = np.arange(len(data[0][0]))+1. 
#    
#    # find  wavelength solution
#    # !! I should have something here for when their's just a single column which is data, i.e. no wavelength solution
#    if preferred_disp == 'linear':
#        poly_order = 1
#        disp_type = 'linear'
#    elif preferred_disp == 'any' or preferred_disp in ['poly','ordinary poly']:
#        poly_order = 6
#        disp_type = 'poly'
#      
#    else:
#        print "WARNING: CURRENTLY CAN ONLY APPLY LINEAR AND POLYNOMIAL FITS TO TEXT DATA, USING POLYNOMIAL"
#        # !! when I create a fancy get_disp for a wl array then I can put that in here
#        poly_order = 6
#        disp_type = 'poly'
#    
#    disp = scipy.polyfit(pts,wl[0][0],poly_order)
#    disp = np.flipud(disp) # switch the order which the dispersion is given [0] = zero-ith order
#    disp = np.concatenate((disp,np.array([0.0]))) # array length = 8
#    
#    spec_obj._private_info[(0,0)]['disp']=[disp,disp_type]
#    spec_obj._private_info[(0,0)]['rv'] = [0]
    
    






############################################################################
# readin is the main function for input


pass
#################################################################################
# This is the general read in function for fits and text files

def readin (filename, hdu=0,  non_std_fits=False,
            text_comments='#', text_skiprows=0, get_data=False, verbose=False,
            apply_WCS_rv=False):
    """
    Read data from a fits or text file and compute the wavelength
    
    INPUTS:
    ==============   ===========================================================
    Keyword          (type) Description
    ==============   ===========================================================    
    filename         (string) Location relative to current working directory
                       file must be fits or a txt file of wl,data pairs (can 
                       also include a column for inverse varience)
    hdu              (int) Which header to use to get the data from
    text_comments    (string) Information to the right of this symbol will be ignored
    text_skiprows    (int) Number of rows to skip in text document
    get_data         (bool) If True then return tuple of (wl,data,inv_var) else return eyeSpec_spec
    apply_WCS_rv   (bool) If True and if the wavelength solution is stored in 
                       IRAF WCS format with a non-zero radial velocity then it is
                       applied to the data
                      
    ==============   ===========================================================
    
    
    
    """
    multi_order_txt = False
    use_naxis2='all'
    use_naxis3='all'
    
    
    preferred_wlsoln=None # !! need to fix this
    # !! should also be able to input wavelength solution?
    
    if preferred_wlsoln is not None: preferred_wlsoln = wlsolvefxn.get_func_name(preferred_wlsoln)
    
    #### check if file exists   ####### #############
    if not os.path.exists(filename): raise IOError("File does not exist:'"+filename+"'")


    #### check if file is text#############  
    np_kwargs = {'comments':text_comments,
                 'skiprows':text_skiprows}
    is_text_file, txt_data = check_for_txt_format(filename,**np_kwargs)

    #### if it is a text file ######################
    if is_text_file:
        spec_obj = readin_txt(filename,txt_data,get_data)        
        return spec_obj      

    #### now check how it behaves as a fits file
    if non_std_fits: hdulist = pyfits.open(filename)
    else:
        # give standard pyfits readin a try
        try: hdulist = pyfits.open(filename)
        except: raise IOError("PYFITS DOES NOT LIKE THE FILE YOU GAVE ('"+filename+"'), TO SEE WHAT ERROR IT GIVES TRY: hdulist = pyfits.open('"+filename+"')")


    #### open up fits file ##############################
    hdulist = pyfits.open(filename)

    # select which header unit ot use
    if len(hdulist) > 1: 
        hdu = int(hdu)
        hdu = np.clip(hdu,0,len(hdulist)-1)
    else: hdu = 0

    # specify the current header unit
    header_unit = hdulist[hdu]
    prihdr = header_unit.header

    # can display some useful information 
    if verbose: 
        print "="*60
        print (hdulist.info(),'\n')
        if len(hdulist) > 1:
            print "="*20+" USING HEADER: "+"="*20
            print repr(hdulist[hdu])

    ##### fill in the data class
    # not get header info of relevance
    simple = query_fits_header(prihdr,'SIMPLE',noval=False)
    xtension = query_fits_header(prihdr,'XTENSION')
    if simple.found:
         if not simple.val: print "HeadsUp: Header Keyword SIMPLE is False, you may encounter unexpected behavior"
    else:
        if not xtension.found: print "HeadsUp: No extension keyword found in headers, you may encounter unexpected behavior"
            
            
    #### read in important information from header, if present
    ibits = query_fits_header(prihdr,'BITPIX') # how many bits per pixel in the data? Not currently necessary, numpy will adapt
    
    naxis  = query_fits_header(prihdr,'NAXIS' ,noval=0) # how many dimenstions?
    naxis1 = query_fits_header(prihdr,'NAXIS1',noval=0) # number of points per order
    naxis2 = query_fits_header(prihdr,'NAXIS2',noval=0) # number of orders
    naxis3 = query_fits_header(prihdr,'NAXIS3',noval=0) # number of different spectra

    apformat = query_fits_header(prihdr,'APFORMAT')
    if apformat.found: print "WARNING: I'M NOT SURE HOW TO DEAL WITH APFORMAT VALUES" # !! though I think it's just the spec files

    if not naxis.found: raise IOError("ERROR: Keyword NAXIS not found")

    bzero = query_fits_header(prihdr,"BZERO",noval=0)
    bscale = query_fits_header(prihdr,"BSCALE",noval=1)

    ###### read in data ##############################################
    data = header_unit.data

    if data is None:
        wl, data, inv_var = np.zeros(3).reshape((3,1))
        if get_data: return (wl,data,inv_var)
        else: return eyeSpec_spec(wl,data,inv_var,header_unit.header)
    else:
        # check that data matches up with at least one of the dimensions
        if data.ndim != naxis.val: raise ValueError("Dimension of data "+str(data.ndim)+" does not match keyword naxis "+str(naxis.val))
        
        statement = 'Dimension does not match data.shape = '+str(data.shape)+" fits file (naxis1, naxis2, naxis3) "+str(tuple([naxis1.val,naxis2.val,naxis3.val]))
        if   data.ndim == 1: 
            assert data.shape  ==  (naxis1.val,) , statement
            data = data.reshape((1,1,)+data.shape)
            
        elif data.ndim == 2: 
            assert data.shape == (naxis2.val, naxis1.val), statement
            data = data.reshape((1,)+data.shape) 
               
        elif data.ndim == 3: 
            assert data.shape == (naxis3.val, naxis2.val, naxis1.val), statement
    
    ##### Determine the which data is useful   
    # which orders to read in  
    nband = np.arange(data.shape[0])+1
    nord = np.arange(data.shape[1])+1

        
    ##### Calculate the wavelengths for the data
    # set up wavelength and inverse_variance
    wl = np.ones(data.shape)
    
    # get the wavelength coefficients
    wlcoeff = wlsoln_coeff_from_header(header_unit.header, apply_WCS_rv, preferred_wlsoln)
    
    # the same wavelength solution is applied to all bands so just pick the first and broadcast
    band = 0
    priv_info = {}
    
    # go through all the orders
    do_progress = True
    progressive_pt = 1 # this will advance and be used when there is no wavelength solution
    for i in xrange(len(nord)):
        order_i = nord[i]

        # get the coefficients and function type    
        equ_type = wlcoeff.get_equation_type()
        if equ_type in ['none',None,'no solution'] and do_progress: 
            coeff = [progressive_pt,1]
            equ_type = 'pts'
        else: coeff = wlcoeff.get_coeffs(order_i)
        
        # pts[0] = 1 :: this was definitely the right thing to do for SPECTRE's 1-D output but may not be for other equations, may need pts[0]=0,  this may be for bzero,bscale
        pts = np.arange(len(wl[0][i]))+1    
        # apply function
        wl[0][i] = wlsolvefxn(pts, coeff, equ_type)    
    
        progressive_pt += len(pts)
  
        for j in xrange(len(nband)): 
            band_j = nband[j]
            if (band_j,order_i) not in priv_info: priv_info[(band_j,order_i)] = {} 
            # record the private information
            priv_info[(band_j,order_i)]['disp']= [coeff, equ_type]
            priv_info[(band_j,order_i)]['rv'] = [0]                
            priv_info[(band_j,order_i)]['disp extr'] = deepcopy(wlcoeff.extra)
        
    # now propogate the solution to the other bands
    stdwl = wl[0]
    for i in xrange(1,len(nband)): wl[i] = stdwl      
    
    inv_var = np.ones(data.shape)
    #=================================================================#
    # return the data .OR. go on and create the spec_obj
    if get_data: return (wl, data, inv_var)

    #=================================================================#    
    spec_obj = eyeSpec_spec(wl,data,inv_var,header_unit.header)
    # set up private information
    priv_info['filename'] = filename
    spec_obj.filename = filename
    
    bands = np.array(np.arange(1,len(data)+1),dtype=str)
    band_info = {}
    i = -1
    for key in prihdr.keys():
        if key[:6] != 'BANDID': continue
        if i < len(bands):
            i+=1
            bands[i] = prihdr[key]
            band_info[key] = prihdr[key]
        else: raise IOError("MORE BANDID KEYWORDS IN HEADER THAN FIRST DIMENSION OF DATA") 

    # add band info if available:
    if len(band_info) != 0: priv_info['bandids'] = band_info
    else: priv_info['bandids'] = None
    
    # match up the private info created during read in to the spec_obj
    for key in priv_info: spec_obj._private_info[key] = priv_info[key]
    
    # map fits value => acutal index
    #    spec_obj._bands = {}
    #    spec_obj._orders = {}
    #    for i in range(len(nspec)): spec_obj._bands[nspec[i]] = i
    #    for i in range(len(nord)): spec_obj._orders[nord[i]] = i
    # 
        
    if 7 in nband: spec_obj.set_band(6) # this is where Magellian data stores it's object data, i.e. BANDID7 which is index 6

    if len(hdulist) > 1: spec_obj.hdrlist = [h.header for h in hdulist]
        
    return spec_obj
    


pass
#################################################################################
# extended IO functions










