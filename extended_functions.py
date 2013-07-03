# these are other function wrappers instead of only the sp ones


def jo (Load=None):
    """
    This prompts the user for some files, or one file and then joins those orders
    along the way you have options to Save, 

    INPUTS:
    Load : (SaveProgress) this is of class SaveProgress and should have the attribute self.jo 

    """
    # keep objs along the way and then delete at the end, same with clean-up files
    
    # read in data (could be from multiple files
    # Setoverlap regions => line_data, new_obj
    # join_order => newNEW_obj 
    # saves the newNEW_obj to a fits file
    pass


def combine (obj=None):

    # if obj: skip file get
    # print "This will walk through the process of combining several orders into one. Please give a file with multiple orders or a list of single order files

    # filename = getfilename(iotype = 'r')

    # is filename a single file? 
    # check the file, is it a fits => load as obj
    # is it a list of files? => if spectre files then use load_spectre_files
    #   if not spectre files, load them and check that they are 1-D, combine them
    pass

    # line_data = setoverlap(obj)

    # new_obj = spec_combine(obj,line_data)

    # save_to_file(new_obj)
 
    

