if __name__ == '__main__':
    check = ['np','os','deepcopy','wx','math','FigureCanvas','NavigationToolbar2Wx','Figure','savefig']
    for val in check:
        if val not in globals(): raise StandardError("Not in globals:"+val)
else:
    from eyeSpec.dependencies import os, sys, time, deepcopy, pdb
    from eyeSpec.dependencies import np, math
    from eyeSpec.dependencies import plt, Figure, FormatStrFormatter, savefig
    from eyeSpec.dependencies import wx, FigureCanvas, NavigationToolbar2Wx

 
################################################################################
class AChoice:
    def __init__ (self,save_message,open_message,wildcard,defaultfile,save_method,open_method,progress_routine):
        self.smessage = save_message
        self.omessage = open_message
        self.wildcard = wildcard
        self.dfile = defaultfile
        self.smethod = save_method
        self.omethod = open_method
        self.sprog_routine = progress_routine
        
    def can_save (self):
        if self.smethod is None: return False
        else: return True

    def can_open (self):
        if self.omethod is None: return False
        else: return True

class SaveOpenChoices:
    def __init__ (self):
        self._choices = {}
        self._set_save_choices = []
        self._set_open_choices = []

    def __repr__ (self):
        return 'SaveOpenChoices'

    def __iadd__ (self,save_open_choices):
        return self.__add__(save_open_choices)

    def __add__ (self,save_open_choices,clobber=False):
        if repr(save_open_choices) != repr(self): raise ValueError("Not correct concatenation")
        soc = save_open_choices
        
        # add the choices
        for ch in soc._choices:
            if ch in self._choices and not clobber:
                print "Keyword already used:"+ch
                continue
            self._choices[ch] = soc._choices[ch]
    
        # add the save choices
        for ch in soc._set_save_choices:
            if ch in self._set_save_choices and not clobber:
                print 'save keyword already used :'+ch
            self._set_save_choices.append(ch)
        
        # add the open choices
        for ch in soc._set_open_choices:
            if ch in self._set_open_choices and not clobber:
                print 'open keyword already used :'+ch
            self._set_open_choices.append(ch)   
             
        return self
    
    def add (self,keyword,
             wildcard='All files (*)|*',
             default_file='none',
             filetype=None,
             save_method=None,
             open_method=None,
             progress_routine = 'undefined'):
        
        if keyword == 'SAVE PROGRESS': 
            filetype = 'PROGRESS SAVE'
            
        
        if filetype is None: filetype = keyword
        smessage = 'Save '+filetype+' File as...'    
        omessage = 'Open '+filetype+' File:'
                
        wcard = filetype+' Files '+wildcard
        defaultfile = default_file
                        
        smeth = save_method
        ometh = open_method
    
        self._choices[keyword] = AChoice(smessage,
                                         omessage,
                                         wcard,
                                         defaultfile,
                                         smeth,
                                         ometh,
                                         progress_routine)
        
        if smeth is not None: self._set_save_choices.append(keyword)
        if ometh is not None: self._set_open_choices.append(keyword)

    def set_choices (self,which,choices):
        ch_list = []
        for ch in choices:
            if ch not in self._choices:
                print "Keyword not defined:"+ch
                continue
            ch_list.append(ch)

        which = which.lower()
        if which == 'save':
            self._set_save_choices = ch_list
        elif which == 'open':
            self._set_open_choices = ch_list
        
    def get_choices (self,which):
        which = which.lower()
        if which == 'save':
            return self._set_save_choices
        elif which == 'open':
            return self._set_open_choices
        else:
            print "which must be open or save"
        
    def check_choices (self,which='both'):
        def checkit (which,set_choices):
            for kw in set_choices:
                if kw not in self._choices:
                    print which+" keyword not in choices :"+kw
                    continue
                copen = (which == 'open' and not self._choices[kw].can_open())
                csave = (which == 'save' and not self._choices[kw].can_save())
                if copen or csave: 
                    print "can't "+which+" keyword :"+kw                                            
     
        if which in ['save','both']:
            checkit('save',self._set_save_choices)
        if which in ['open','both']:
            checkit('open',self._set_open_choices)
        
    def _check_keyword (self,keyword):
        keyword = str(keyword)
        if keyword not in self._choices:
            print "keyword not defined : "+keyword
            return True
        return False    

    def get_save (self,keyword):
        if self._check_keyword(keyword): return None
        
        if not self._choices[keyword].can_save():
            print "can't save: "+keyword
            return None
        
        ch = self._choices[keyword]
        out = {'message':ch.smessage,
               'wildcard':ch.wildcard,
               'defaultfile':ch.dfile,
               'method':ch.smethod,
               'routine':ch.sprog_routine}
        return out
    
    def get_open (self,keyword):
        if self._check_keyword(keyword): return None
        if not self._choices[keyword].can_save():
            print "can't open: "+keyword
            return None
                
        ch = self._choices[keyword]
        out = {'message':ch.omessage,
               'wildcard':ch.wildcard,
               'defaultfile':ch.dfile,
               'method':ch.omethod,
               'routine':ch.sprog_routine}
        return out        


################################################################################

class ProgressSave:
    """
    This class is used for saving progress for eyeSpec. From given routines it will create/use a directory in which all the information to restore a given session is kept


    #==========================================================================#
    To create a new directory to save progress
    sprog = ProgressSave()
    # will automatically craete a directory sprog_? where ? is the next largest number for the sprog_* directories

    sprog = ProgressSave("new_sprog/")

    To open up an existing directory
    sprog = ProgressSave('sprog_4',mode='existing')


    """

    def __init__ (self,sprog_dir=None,mode='new'):
        """
        INPUTS: see help(SaveProgress.set_current_dir)
        
        """
        self.set_current_dir(sprog_dir=sprog_dir,mode=mode)
        
        if 'record_contents' not in dir(self):
            self._create_blank_record_contents()
            self.set_current_save_num (0,save_type='automatic')


        all_save_nums = self.record_contents.keys()
        if len(all_save_nums) == 0: self.set_current_save_num (0,save_type='automatic')
        else: self._current_save_num = np.max(all_save_nums)
        

    def get_current_save_num (self):
        return deepcopy(self._current_save_num)

    def set_current_save_num (self,save_num=None,save_type='manual',overwrite=False):
        """
        If the value for save_num is in self.record_contents.keys() then it will set to access that particular entry depending on overwrite, else it will create a new entry and set to access that
        """
        


        # set a new number based on the largest value
        all_save_nums = self.record_contents.keys()


        is_existing = False
        if save_num is None:
            if len(all_save_nums) == 0: self._current_save_num = 0 
            else: self._current_save_num = np.max(all_save_nums)+1
        else:
            save_num = int(save_num) # must be an integer

            if len(all_save_nums) == 0 or save_num not in all_save_nums: self._current_save_num = save_num

            elif save_num in all_save_nums and not overwrite:
                #self._current_save_num = np.max(all_save_nums)+1
                #print "Given save number already exists, creating a new one:"+self._current_save_num
                self._current_save_num = save_num
                is_existing = True
                
            elif save_num in all_save_nums and overwrite:
                self._current_save_num = save_num
            
        # now that self._current_save_num is set
        if not is_existing:
            self._create_blank_content(save_type=save_type,overwrite=overwrite) 
        

    def delete_save (self,save_num):
        all_save_nums = self.record_contents.keys()
        if len(all_save_nums) == 0: return

        save_num = int(save_num) # must be integer
        if save_num not in self.record_contents.keys(): return

        contents = self.record_contents[save_num]

        # if there are files associated then delete them
        for value in contents.values():
            cdir = self.get_current_dir()
            if os.path.exists(cdir+value): os.system("rm "+cdir+value)
        
        # delete the entry
        del self.record_contents[save_num]

        if save_num == self.get_current_save_num(): self.set_current_save_num()

    def _check_for_sprog_dir (self):
        if not os.path.exists(self._sprog_dir): raise ValueError("Directory no longer exists! '"+self._sprog_dir+"'")
        # !! could create a new directory....

    def add_file_to_dir (self,filename,save_num = None,keyword='',overwrite=False):
        filename = str(filename)
        if not os.path.exists(filename): raise IOError("Whoops, file name does not exist:"+filename)

        if os.path.exists(self._sprog_dir+filename) and not overwrite: raise IOError("Whoops, file name already exists:"+self._sprog_dir+filename)

        self.add_record_info(save_num,keyword=keyword,value=filename,overwrite=overwrite)
        os.system('mv '+filename+' '+self._sprog_dir)


    def focus_on_last_record (self,routine=None):
        
        if len(self.record_contents.keys()) == 0: return False

        if routine is None:
            all_save_nums = self.record_contents.keys()
        else:
            routine = str(routine)
            if routine not in self.record_by_routine.keys():
                print "No entry for routine:"+routine
                return False
            all_save_nums = self.record_by_routine[routine]    
        self.set_current_save_num(np.max(all_save_nums))
        return True


    def keys (self):
        if self._current_save_num is None: return deecopy(self.record_contents.keys())
        else: return deepcopy(self.record_contents[self._current_save_num].keys())


    def add_record_info (self, save_num, keyword = '', value='',add_comment='',overwrite = False):
        """
        A more restrictive way then using getitem/setitme which will take save_num = self._current_save_num  and overwrite = True



        """

        # now add the appropriate information if desired
        if save_num is None: pass #!! return  #!! I don't know if this is the thing I want to do
        if save_num is 'current': save_num = deepcopy(self._current_save_num)
        save_num = int(save_num)

        if save_num not in self.record_contents.keys(): 
            print "desired save_num is not an option so no action is done"
            return
            
        keyword = str(keyword).strip().upper()
        value = str(value)
        add_comment = str(add_comment)
        if keyword == '': keyword = 'EXTRA_1'
        
        if len(keyword) > 10:
            print "HeadsUp: only keywords of len < 10 can be used, I'm cropping the current given value"
            keyword = keyword[:10]

        while True:
            if keyword in self.record_contents[save_num].keys():
                if overwrite:
                    self.record_contents[save_num][keyword]=value
                    break
                else:
                    if keyword.find('EXTRA_')==0: 
                        val = int(keyword.split('_')[1]) + 1
                        if val >= 999: raise ValueError("Extra value can not exceed 999")
                        keyword = 'EXTRA_'+str(val)
                    else:
                        raise TypeError("This keyword ('"+keyword+"') already has content in the desired save_num:"+str(save_num))

            else:
                if len(add_comment)!=0: value += "# "+add_comment
                self.record_contents[save_num][keyword] = value
                break


    def write_record_file (self,overwrite=False):
        filename = self._sprog_dir+'record.sprog'
        if os.path.exists(filename) and not overwrite: raise IOError("Whoops, file name already exists:"+filename)
        # !! instead of an error I could go ahead and open that record file

        
        f = open(filename,'w')
        lines = []
        lines.append("#"+"="*78+"#")
        lines.append("START OF HEADER")
        lines.append("#"+"="*78+"#")
        lines.append("#>>> eyeSpec Progress Save Record File")
        lines.append("#>>> This file is formatted so that the eyeSpec code can read this directory and restore your previous session. This file is specifically formated so that eyeSpec can read this directory. You can put whatever you would like in the lines below this, but before the END OF HEADER. ")
        
        # header information goes here
        for line in self.record_header:
            lines.append(line)

        lines.append("#"+"="*78+"#")
        lines.append("END OF HEADER")
        lines.append("#"+"="*78+"#")        

        
        save_order = ['CREATOR','DATE_TIME','SAVE_TYPE','ROUTINE']

        for i in range(1,len(self.record_contents)+1):
            lines.append("")
            lines.append("#"+"-"*78+"#")
            
            # write in reverse order
            save_num = self.record_contents.keys()[-i]
            save_contents = self.record_contents[save_num]
            
            # save the number
            lines.append(format("SAVE_NUM",'<10')+": "+str(save_num))

            # save standard information
            for keyword in save_order:
                if keyword not in save_contents.keys(): raise TypeError("Could not find necessary keyword:"+keyword) 
                value = save_contents[keyword]
                line = format(keyword,'<10')+": "+str(value)
                # !! depending on keyword I could add other comment information
                # line += " # can be manual from the user or automatic"
                lines.append(line)

            # find the other keywords
            extra_keywords = []
            for key in save_contents.keys():
                if key not in save_order: extra_keywords.append(key)

            # save the extra information for the routine
            for keyword in extra_keywords:
                value = save_contents[keyword]
                line = format(keyword,'<10')+": "+str(value)
                lines.append(line)



        f.write('\n'.join(lines))
        f.close()


    def _create_blank_content (self,save_num=None,save_type = 'automatic',overwrite=False):
        if save_num is None: save_num = self._current_save_num
        else: save_num = int(save_num) # must be an integer

        all_save_nums = self.record_contents.keys()
        if save_num in all_save_nums and not overwrite:
            raise ValueError("This save number already exists:"+save_num)

    
        self.record_contents[save_num] = {}
        
        self.record_contents[save_num]['CREATOR'] = os.getlogin()+"@"+socket.gethostname()
        self.record_contents[save_num]['DATE_TIME'] = time.ctime()
        self.record_contents[save_num]['SAVE_TYPE'] = save_type
        self.record_contents[save_num]['ROUTINE'] = 'none'

        if 'none' not in self.record_by_routine.keys(): self.record_by_routine['none'] = []
        self.record_by_routine['none'].append(save_num)


    def _create_blank_record_contents (self,save_type = 'automatic',overwrite=False):

        if 'record_contents' in dir(self) and not overwrite:
            print "self.record_contents already exists"
            return

        if save_type not in ['automatic','manual']: raise TypeError("save_type must be either 'automatic' or 'manual'")

        self.record_contents = {}
        self.record_header = []
        self.record_by_routine = {}

        
        # need to put some information into these
        # self._create_blank_content(0,save_type=save_type,overwrite=overwrite)
        
    def add_header_info (self,line):
        line = str(line)
        self.record_header.append(line)


    def get_current_dir (self):
        return deepcopy(self._sprog_dir)

    def set_current_dir (self,sprog_dir=None,mode='transfer'):
        """
        INPUTS:
        sprog_dir : (string) This is the directory the save progress focuses on. If None then it will create a new directory sprog_? where ? is the largest interger which doesn't conflict with any other sprog_* directories

        mode : (string) 
               'existing' will focus on an existing directory or 
               'transfer' (default) will maintain the current self.record_ information which will then be moved to the sprog_dir directory
               'new' will create new self.record_ information for this sprog_dir

        """


        if mode not in ['existing','transfer','new']: return

        # if no save progress directory given create new
        if sprog_dir is None:
            ls_list = os.listdir('.')
            sprog_dirs = []
            sprog_nums = []
            # look through all the current sprog_* directories
            for line in ls_list:
                if line.find('sprog_') == 0:
                    sprog_dirs.append(line)
                    sprog_nums.append(int(line.split('_')[1]))
            # take the maximum number and add one
            cur_num = np.max(sprog_nums) + 1
            self._sprog_dir = 'sprog_'+str(cur_num)+'/'
            os.system('mkdir '+self._sprog_dir)

            # if mode not 'new' then it will take the existing records
            if mode == 'new': self._create_blank_record_contents()

        else:
            sprog_dir = str(sprog_dir) # must be a string value
            if sprog_dir[-1] != '/': sprog_dir += '/'
            # does the give directory exist
            if os.path.exists(sprog_dir):
                self._sprog_dir = sprog_dir
                rfile = self._sprog_dir+"record.sprog"
                if mode == 'new':
                    if os.path.exists(rfile): print "HeadsUp: a 'record.sprog' file already exists for dir "+sprog_dir+" new information is being created that when saved will overwrite the existing"
                    self._create_blank_record_contents()
                else: 
                    if os.path.exists(rfile) and mode != 'transfer': self.open_record_file()
                    elif os.path.exists(rfile) and mode == 'transfer': pass
                    elif not os.path.exists(rfile): self._create_blank_record_contents()

            # create new directory            
            else:
                self._sprog_dir = sprog_dir
                os.system('mkdir '+self._sprog_dir)
                if mode == 'new': self._create_blank_record_contents()
        

    def __getitem__ (self,key):
        save_num = self._current_save_num
        key = str(key).upper()
        return self.record_contents[save_num][key]

    def __setitem__ (self,key,val):
        save_num = self._current_save_num
        key = str(key).upper()
        self.record_contents[save_num][key] = val
        


    def check_dir_contents (self,save_num = None):
        # !! this goes through the record_contents and makes sure all the files which are suppose to be there are
        # if save_num is None it will check all
        # if save_num is specified it will check for just those associated with that particular save_num

        pass

    def open_record_file (self,overwrite=False):

        if 'record_contents' in dir(self) and not overwrite:
            print "self.record_contents already exists"
            return

        filename = self._sprog_dir+'record.sprog'
        if not os.path.exists(filename): raise IOError("Whoops, file name does not exist:"+filename)
    
        f = open(filename)
        

        self.record_contents = {}
        self.record_header = []
        self.record_by_routine = {}

        end_header = False
        for line in f:
            line = line.rstrip()

            # check if you're still in the header section
            if line.find('END OF HEADER') == 0:
                end_header = True
                continue

            if not end_header:
                if line.find('START OF HEADER') == 0:continue
                if line.find('#>>>') == 0: continue
                if line.find("#"+"="*30) == 0: continue
                # if len(line) == 0: continue

                self.record_header.append(line)
                continue
        
            # Ignore everything to the right of a comment
            no_c_line = line.split('#')[0]

            # if the line is blank then skip
            if len(no_c_line) == 0: continue
            
            # Get the important section of the line
            split_i = no_c_line.find(':')
            if split_i == -1:
                print "HeadsUp: Formatting error for line, skipping:'"+no_c_line+"'"
                continue
            keyword = no_c_line[:split_i]
            keyword = keyword.strip().upper()
            value = no_c_line[split_i+1:]
            value = value.strip()
            
            # if there's no keyword, no value
            if len(keyword) == 0: continue
            if len(value) == 0: continue

            # Is this a new save number?
            if keyword.find('SAVE_NUM') == 0: 
                save_num = int(value) # value must be an integer
                if save_num in self.record_contents.keys(): raise TypeError("Format of the record file is not correct, each SAVE_NUM must be unique. Found "+value+" twice")
                self.record_contents[save_num] = {}
                continue

            else:
                if keyword.find('ROUTINE') == 0: 
                    if value not in self.record_by_routine.keys(): self.record_by_routine[value] = [save_num]
                    else: self.record_by_routine[value].append(save_num)
                self.record_contents[save_num][keyword] = value

class InputOutput:
    def __init__ (self): pass
    def Output_To_File (self, save_choices, progress_message=True):
        """
        Create and show the Save FileDialog
        """
#         # progress save is formatted differently: 
#         save_choices['PROGRESS SAVE'] = 'snr_progress_save_function' #==> snr_progress_save_function(sprog)
# #         save_choices['PROGRESS SAVE'] = {}
# #         save_choices['PROGRESS SAVE']['ROUTINE'] = 'SNR Edit'
# #         save_choices['PROGRESS SAVE']['SNR_SEL'] = ['something here'] # = savefunction   which will take filename and do approgriate stuff
# #         save_choices['PROGRESS SAVE']['ORIG_OBJ'] = ['something here']
# #         save_choices['PROGRESS SAVE']['EDIT_OBJ'] = ['something here']
# #         save_choices['PROGRESS SAVE']['HISTORY'] = ['something here']
# #         save_choices['PROGRESS SAVE']['VAR_PARAMS'] = ['something here']
        
        
        ################################
        # Decide what thing to save
        dlg = wx.SingleChoiceDialog(None,"Choose one to save:","Save Options",save_choices.get_choices('save'))
        if dlg.ShowModal() == wx.ID_OK: which_to_save = dlg.GetStringSelection()
        else: return
        dlg.Destroy()

        # defaults:
        message = "Save file as..."
        wildcard = 'All files (*)|*'
        defaultFile = ""
        def _default_function (path):
            print "You choose :"+str(path)
        save_function = _default_function
            

        # from input
        specifics = save_choices.get_save(which_to_save)        
        if specifics is None: raise ValueError("Problem with save option "+which_to_save)
        
        
        if specifics['method'] is None: 
            print "SAVE: "+str(which_to_save)+" currently unavailable"
            return
        
        if which_to_save == 'PROGRESS SAVE':
            dlg = wx.DirDialog(None,"Choose a Save Progress Directory:",
                               style = wx.DD_DEFAULT_STYLE)
            if dlg.ShowModal() == wx.ID_OK:
                sprog_dir = dlg.GetPath()
            dlg.Destroy()
            
            if progress_message: print "=== SAVING PROGRESS INTO DIRECTORY:"+sprog_dir
            save_choices['PROGRESS SAVE'](sprog_dir,'existing')
            if progress_message: print "=== PROGRESS SAVED"

            # need to have a button: NEW PROG
            # which will create a new sprog_? directory

            # ==> get path for the sprog directory, e.g. sprog_?

            # sprog = ProgressSave(sprog_?,mode='existing')
            # sprog = ProgressSave(sprog_?,mode='new')

            #save_function = specifics # snr_progress_save_function
            # save_function(sprog)

        else:
            message = specifics['message']
            wildcard = specifics['wildcard']
            defaultFile = specifics['defaultfile']
            save_function = specifics['method']

            dlg = wx.FileDialog(None, message=message,
                                defaultDir=os.getcwd(), wildcard=wildcard,
                                defaultFile=defaultFile, style=wx.SAVE | wx.OVERWRITE_PROMPT)
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                save_function(path)
            dlg.Destroy()
 
    def Input_From_File (self, open_choices):
                ################################
        # Decide what thing to save
        dlg = wx.SingleChoiceDialog(None,"Choose one to open:","Open Options",open_choices.get_choices('open'))
        if dlg.ShowModal() == wx.ID_OK: which_to_open = dlg.GetStringSelection()
        else: return

        dlg.Destroy()

        # defaults:
        message = "Open file:"
        wildcard = 'All files (*)|*'
        defaultFile = ""
        def _default_function (path):
            print "You choose :"+str(path)
        open_function = _default_function
            
        ########################
        #if which_to_open != save_choices['CHOICES'][0]: print "whoot==",which_to_open
        #which_to_open = save_choices['CHOICES'][0] # !! this is just for initial checking

        specifics = open_choices.get_open(which_to_open)        
        if specifics is None: raise ValueError("Problem with save option "+which_to_open)

        if specifics['method'] is None: 
            print "OPEN: "+str(which_to_open)+" currently unavailable"
            return
        if which_to_open == 'PROGRESS SAVE':
            dlg = wx.DirDialog(None,"Choose a Save Progress Directory:",
                               style = wx.DD_DIR_MUST_EXIST)
            if dlg.ShowModal() == wx.ID_OK:
                sprog_dir = dlg.GetPath()
            dlg.Destroy()
            sprog = ProgressSave(sprog_dir,'existing')

            routine = specifics['routine']
            
            routine_found = sprog.focus_on_last_record(routine)
            if routine_found:
                all_save_nums = sprog.record_by_routine[routine]
                if len(all_save_nums) > 1:
                    all_save_nums = np.array(all_save_nums)
                    sort_i = all_save_nums.argsort()
                    all_save_nums = np.array(all_save_nums[np.flipud(sort_i)],dtype='a15')

                    choices = []
                    for i in range(len(all_save_nums)):
                        j = all_save_nums[i]
                        choices.append(j+" -- "+sprog.record_contents[int(j)]['DATE_TIME'])                        

                    dlg = wx.SingleChoiceDialog(None,"Multiple Saves, Choose One:","Progress Save Choices",list(choices))
                    if dlg.ShowModal() == wx.ID_OK: 
                        which_save_num = dlg.GetStringSelection().split("--")[0].strip()
                        sprog.set_current_save_num(which_save_num)
                    else: sprog = None
                    dlg.Destroy()
                    
            return sprog

        else:
            message = specifics['message']
            wildcard = specifics['wildcard']
            defaultFile = specifics['defaultfile']
            open_function = specifics['method']
            
            dlg = wx.FileDialog(None, message=message,
                                defaultDir=os.getcwd(), wildcard=wildcard,
                                defaultFile=defaultFile, style=wx.OPEN)
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                open_function(path)
            dlg.Destroy()
            return None

            

