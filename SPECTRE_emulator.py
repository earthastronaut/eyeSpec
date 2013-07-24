# these are equavlient to SPECTRE's commands
# this should be a module level eyeSpec.sp.
# from eyeSpec.sp import *    <== would act like spectre only with rd() instead of rd

#from Readin import readin


data_arrays = {}


def rd ():
    """
    Read in new fits file
    
    INPUT:
    select_orders : do you want to just input one or a selected amount of orders?

    """
    fname =  getfilename()
    # !! add error satements
    if fname == 'no_file':
        # aborting from fname
        return 
    

    while True:
        arr = raw_input("INPUT TO WHICH ARRAY? (x,y,z)\n")

        # to add options I need them to be explicit
        if arr == 'x' or arr == '':
            global x
            x = readin(fname,verbose=True)
            x._obj_name = variablename(x)
            break
        elif arr == 'y':
            global y
            y = readin(fname,verbose=True)
            y._obj_name = variablename(y)
            break
        elif arr == 'z':
            global z
            z = readin(fname,verbose=True)
            z._obj_name = variablename(z)
            break
        # can add more types just by adding based on the format above
        #         elif arr == 'u':
        #             global u
        #             u = readin(fname,verbose=True)
        #             u._obj_name = variablename(u)
        #             break

        elif arr in ['abort','a','q']:
            return
        else:
            #exec("global "+arr) # < this doesn't work for some reason
            #exec(arr+" = readin(fname)")

            print "PLEASE ENTER x,y,z\n"
        
    stat = data_set._additem(arr)
    #if stat != 0 :
    #    return


#    x = readin(fname)
#                  select_orders=select_orders,
#                  select_band=select_band)

def d0 (obj=None):
    if obj is None:
        if _check_for_x(): p_obj = x
        else: print "Array x is not defined"
    else:
        # !! check to make sure obj is of right type
        p_obj = obj
    plotall(p_obj,fig_num=100)
    plt.show()

def dx ():
    if not _check_for_xyz('x'): return
    
def dy (obj=None):
    # !! so similar check to d0, this would add the x data
    # !! would check if fig_num = 100 is defined, if so add it to that plot
    specoplot(y)

def dz (z):
    specoplot(z)


# def du (ord,filename = 'manual file'):
#     """
#     dumps a single order to a file
#     """

#     if filename == 'manual file':
#         fname = getfilename(iotype='w') 
#     elif type(filename).__name__!='str':
#         print "ERROR: VARIABLE filename MUST BE A string" # error exit
#         return
    
#     if len(filename.strip())!=1:
#         print "ERROR: VARIABLE filename MUST BE A SINGLE VALUE" # !! error exti
#         return

#     fname = filename # !! need to put in something for getfilename to check a given input filename

def du_all (obj,filename = 'manual file',no_comments=False):
    """
    dumps an object to a single text file
    INPUTS:

    """

    if filename == 'manual file':
        fname = getfilename(iotype='w') 
    elif type(filename).__name__!='str':
        print "ERROR: VARIABLE filename MUST BE A string" # error exit
        return
    
    if len(filename.split())!=1:
        print "ERROR: VARIABLE filename MUST BE A SINGLE VALUE" # !! error exti
        return

    fname = filename # !! need to put in something for getfilename to check a given input filename

    # !! if verbose
    # print "Using band", obj.get_band()

    f = open(fname,'w')
    if not no_comments:
        f.write("#  wavelength  flux\n")
        f.write("#  Order 1\n")
    kount = 1
    for ord in obj:
        data_out = np.dstack((ord[0],ord[1]))[0]
        for pair in data_out:
            f.write("   ".join([format(pair[0],'15'),format(pair[1],'15')]))
            f.write("\n")
        kount += 1
        f.write("# Order "+str(kount)+"\n")

    f.close()


###########################################################################

def hd (input='manual file',number_lines=10,verbose=True,prompt_lines=True,return_header=False):
    """
    D header information for a file

    INPUTS:
    file : takes a filename to display, the default prompts the user for the file
    number_lines : the number of lines shown between prompts
    verbose print : information to the screen
    prompt_lines : break up the header by "SEE MORE?" prompt
    return_header : return a list of the header lines
    
    OUTPUTS:
    Header_list (opt) : a list of the header lines as strings

    """
    print "!! CHECK: do you like the ability to scroll through or is it unnecessary?"
#    input = []
#    for val in args:
#        input.append(val)
#        
#    print input
#    print "---"
#    print args

    multiple_headers = False

    if input == 'manual file':
        file = getfilename()
        header = fits.getheader(file)

    elif type(input).__name__ == 'str' and input!='manual file':
        file = input
        header = fits.getheader(file)

    elif input.__class__.__name__ == 'eyeSpec_spec':
        file = "OBJ of file:"+str(input.filename)
        header = deepcopy(input.header)
        if 'combo_hdr' in dir(input):
            print "!! multiple headers are combined"
            multiple_headers = True
    else:
        print "MUST GIVE A FILE NAME OR OBJ"
        return

    print "__"*15


    pauseit=number_lines

    for i in range(len(header.ascardlist())):
        # !! can change this to a while loop and be able to go backwards
        if not verbose: break
        print "   ",str(header.ascardlist()[i]).rstrip()
        
        if i == pauseit and prompt_lines:
            while True:
                stopit = raw_input("\n SEE MORE? ([y],n)")
                if stopit in ['y','yes','']:
                    pauseit += number_lines
                    break
                elif stopit == 'n':
                    break
                else:
                    print "INVALID ENTRY"

            if stopit == 'n':
                break
            print "__"*10

    print "__"*15
    print ("HEADER FOR:"+file+"\n")

    if return_header:
        header_list=[]
        for line in header.ascardlist():
            header_list.append(format(str(line).rstrip(),'80')) # formated for standard fits header file

        return header_list


#x = []
#y = []
#z = []



#def cd (obj=None):


def he(function='display list'):
    """
    Display help file for commands
    INPUT:
    function : specify a specific function to return help on

    """
    print "list of commands not currently available"
    if function != 'display list':
        print help(function)


def vs():
    """
    Define radial velocity
    """
    global rv
    
    while True:
        rv = raw_input("GIVE RADIAL VELOCITY:\n")
        try:
            rv = float(rv)
            break
        except:
            print "VRAD MUST BE A FLOATING POINT"


    arr = raw_input("APPLY TO WHICH ARRAY?\N")
        
    


def xy ():
    """
    switch the x and y data sets
    """
    tmp = x.copy()
    x = y.copy()
    y = tmp

    return

#yx = xy()
#xy = xy()
    
