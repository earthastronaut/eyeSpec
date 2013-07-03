if __name__ != '__main__':
    from eyeSpec import __path__ as path_2_eyeSpec
    from eyeSpec.interactive_IO import ProgressSave, InputOutput
    from eyeSpec.interactive_classes import (EventConnections, History, KeyboardConfiguration, RandomPanel, SysOutListener, State, EditSubplotDialog, eyeSpecBaseDataPanel,
                                            eyeSpecBaseEventManager, eyeSpecBaseApp, eyeSpecBaseFrame, eyeSpecBaseFrame,eyeSpecBaseMainPanel, eyeSpecBaseLineEditor, eyeSpecBaseDataPlot)
    from eyeSpec.extended_IO    import save_spec, load_spec
    from eyeSpec.base_functions import find_overlap_pts,alt_order_colors
    from eyeSpec.app_edit_data  import InteractiveDataEditor
    from eyeSpec.dependencies   import (np, os, sys, time, iget, deepcopy, pdb, scipy, math,
                                        plt, FormatStrFormatter, savefig,
                                        wx, FigureCanvas, NavigationToolbar2Wx, Figure, Button, Path)
      
if wx is None:  raise StandardError("This will not work without wxpython")

# this takes a list of known wavelengths in stationary frame, plots the data against those, then allows for a simple click to radial velocity shift the data to match that point.

# it will walk through each of the stationary points until you quit or it's done


# have the option to pull up a new plot of the solar spectrum for that particular region!!! Awesome

# !! add the solar and arcturus data
# create a directory with the solar data broken into easy to read pieces
# can add arcturus, alpha cen,

# !! could use this for "cross-correlating" and lining up two spectra

def _rv_what_happened ():
    # General
    f = open('TMP_WHAT_JUST_HAPPENED.txt','w')
    f.write("This file is intended to help in the case that the code crashes and you're left with only TMP_* files\n")
    f.write("If this isn't helpful please email Dylan Gregersen <dylan.gregersen@utah.edu>\n")
    f.write("\n")
    
    # What was saved? Why?
    f.write("This file was created from the set radial velocity routine\n")
    f.write("The original spectrum information is in the file: TMP_OBJ_SAVE_ORIG.pkl\n")
    f.write("The final radial velocity was stored in: TMP_FINAL_VRAD.txt\n")
    f.write("The edited spectrum was not saved")
    f.write("\n")

    # How to recover with these files.
    f.write("To restore what you were working on:\n")
    f.write(">>> from eyeSpec import load_spec")
    f.write(">>> import numpy as np")
    f.write(">>> orig_obj = load_spec('TMP_OBJ_SAVE_ORIG.pkl')\n")
    f.write(">>> rv_stat = np.loadtxt('TMP_FINAL_VRAD.txt')\n")
    f.write(">>> edit_obj = orig_obj.copy()\n")
    f.write(">>> edit_obj.edit.apply_rv(rv_stat[0])\n")

    f.close()
    
class OPlotSun:
    def __init__ (self,parent,data_path,xbounds,**mpl_kwargs):
        # where does the data live?
        self.data_path = data_path
        #
        self.wl = np.empty(1)
        self.data = np.empty(1)
        self.spectra, = self.parent.ax.plot(wl,data,**mpl_kwargs)
        pass

    def plot_range (self,xbounds):
        # check to see if self.spectra.get_xdata() is in xbounds
        # if not in bounds, go looking for next range
        # if found next range then read in and append
        pass

    def set_visible (self,truth):
        pass
    
    def get_visible (self):
        pass

    def toggle_visible (self):
        pass

class OplotArcturus:
    def __init__ (self,parent,data_path,xbounds,**mpl_kwargs):
        # where does the data live?
        self.data_path = data_path
        #
        self.wl = np.empty(1)
        self.data = np.empty(1)
        self.spectra, = self.parent.ax.plot(wl,data,**mpl_kwargs)
        pass

    def plot_range (self,xbounds):
        # check to see if self.spectra.get_xdata() is in xbounds
        # if not in bounds, go looking for next range
        # if found next range then read in and append
        pass

    def set_visible (self,truth):
        pass
    
    def get_visible (self):
        pass

    def toggle_visible (self):
        pass

class LineInfo:
    def __init__ (self):

        self.line_dir = "/Users/dylangregersen/Desktop/Astrophysics/data/arcturus/DATA_LINES/"

        self.optical = self._read_file('INFO_LINES.txt')
        self.uv = self._read_file('INFO_LINES_UV.txt')
        self.arcturus_emission = self._read_file('INFO_LINES_appi.txt')
        self.solar_emission = self._read_file('INFO_LINES_appii.txt')
        self.solar_absorb = self._read_file('INFO_LINES_appiii.txt')
        self.arcturus_absorb = self._read_file('INFO_LINES_appiv.txt')


        #------------------------------------#
        # concatenations of the above
        self.all_lines = self.uv+self.optical

        self.arcturus_UV = self.arcturus_emission + self.arcturus_absorb
        self.arcturus_UV = sorted(self.arcturus_UV,key=iget(0,1))

        self.solar_UV = self.solar_emission + self.solar_absorb
        self.solar_UV = sorted(self.solar_UV,key=iget(0,1))
        #------------------------------------#


    def search (self,val,nearby=0,which='all',ret='value',tol='closest'):
        try:val = float(val)
        except: raise TypeError("Value must be convertable to floating point")

        nearby = int(nearby)

        if which == 'all': choose_lines = deepcopy(self.all_lines)
        else: return None

        alines = np.array(map(iget(0),choose_lines))

        i = np.abs(alines-val).argmin()
        mindist = np.min(np.abs(alines-val))

        if ret == 'index': return i

        if type(tol).__name__ == 'float':
            if mindist < tol: return choose_lines[i]
            else: return None
        else: 

            if nearby != 0:
                iless = np.clip(i-nearby,0,len(choose_lines)-1)
                iplus = np.clip(i+nearby,0,len(choose_lines)-1)
                return choose_lines[iless:iplus]
            else: return choose_lines[i]



    def _read_file (self,filename):
        output = []
        filename = self.line_dir+filename
        if not os.path.exists(filename):
            print "Whoops, didn't find file:"+filename
            return None
        
        f = open(filename)
        for line in f:
            line = line.rstrip()
            if line[:1] != '#': output.append([float(line[:12]),line[12:].strip()])
        return output

    def annotate_on_plot (self,ax,xrange,which_lines='all',**mpl_kwargs):
        if len(xrange) != 2: return
        self.ax = ax
        self.annotations = []
        
        label_color='r'

        if which_lines == 'all': info_lines = self.all_lines


        yran = self.ax.axis()[3]-self.ax.axis()[2]
        xran = self.ax.axis()[1]-self.ax.axis()[0]

        
        text_y = 0.9*yran+self.ax.axis()[2]
        data_y = 0.1*yran+self.ax.axis()[2]
        text_x_offset = 0.05*xran


        rolling_i = 0
        for i in range(len(info_lines)):
            if np.min(xrange) < info_lines[i][0] < np.max(xrange):
                if rolling_i==0: 
                    label_hi=text_y
                    rolling_i += 1
                elif rolling_i == 1:
                    label_hi += .02*yran
                    rolling_i += 1
                elif rolling_i == 2:
                    label_hi += .02*yran
                    rolling_i = 0


                annote = self.ax.annotate(info_lines[i][1].strip(),
                                          xy=(float(info_lines[i][0]),
                                              data_y),
                                          xytext=(float(info_lines[i][0])+text_x_offset,
                                                  label_hi),
                                          
                                          ha="right",va="center",rotation=-15,
                                          arrowprops=dict(arrowstyle='->',
                                                          connectionstyle="angle,angleA=-15,angleB=90,rad=10"),
                                          color=label_color,**mpl_kwargs)
                
                self.annotations.append(annote)


    def on_update (self):
        for i in range(len(self.annotations)):
            # adjust the y and x offsets
            pass

class OplotStandards (EventConnections):

    def __init__ (self,parent):
        self.parent = parent


        self.STDspectra = [] # spectra,lines,


        #---------------------------------------------------------#
        # these attributes are standard to my editor classes
        # keyboard configuration
        self.key_cfg = KeyboardConfiguration()
        self.keyboard_cfg = {'display':['1'],
                             'info':'',
                             '1':'Toggle on/off solar data'}
                             
        self.init_connection_callbacks(self)

    def key_press_callback (self,event):
        # restore to normal
        # shift up
        # shift down
        # add constant
        # mult constant

        # 1 toggle solar
        # 2 toggle arcturus
        # 3 ...

        pass



    def update_callback (self,event):
        """  Here I want to have it check to see if the solar data is visible and if so does it fall in the correct range"""

pass
###############################################################################

class RVRecord:
    def __init__ (self):
        self._record_rv = []
        self._reci = None
    
    def get_current_state (self):
        curstate = State("RV_record")
        curstate['record_rv'] = self._record_rv
        curstate['reci'] = self._reci
        
    def apply_state (self,prevstate):
        prevstate.check_id("RV_record")
        self._record_rv = prevstate['record_rv']
        self._reci = prevstate['reci']
        
    def __repr__ (self): 
        return 'RVRecord'
     
    def __getitem__ (self,index):
        return self._record_rv[index]
        
    def get_text_repr (self,index):
        rec = self[index]
        wl = rec[0]
        info = rec[2]
        return format(info,'>10')+" ("+format(wl,'^10.3f')+")"
        
    def add (self,rv,wl,info,verbose=True):
        self._record_rv.append([wl,rv,info])
        self._reci = len(self._record_rv) - 1        
        
        txt = self.get_text_repr(self._reci)
        if verbose: print "RECORDED==> "+format(rv,'11.3f')+" [km/s]  at  "+txt

    def __delitem__ (self,index):
        self.delete(index)

    def delete (self,index=None,verbose=True):
        if index is None: 
            if self._reci is None: return
            index = self._reci
        
        N = len(self._record_rv)
        if index < 0 or index >= N: raise ValueError("Index out of bounds "+str(index)+" not in "+str((0,N)))
        
        if verbose:
            rv = self._record_rv[index][1]
            txt = self.get_text_repr(index)
            print "DELETED==> "+format(index,'3')+"  "+format(rv,'8.3f')+" [km/s]  at  "+txt
        del self._record_rv[index]
        
        if index == N-1: self._reci = N-2
        if self._reci < 0: self._reci = None
               
    def display_record (self): 
        lines = []
        for i in xrange(len(self._record_rv)):
            rv = self._record_rv[i][1]
            txt = self.get_text_repr(i)
            outxt = ''
            if i == self._reci: outxt += "selected >>> "
            outxt += " "+txt+" : "+format(rv,'<10.2f')
            lines.append(outxt)
        
        print "-"*60
        print " Recorded Radial Velocities ===>"
        print "\n".join(lines)
        print " -"*30
        self.calc_stats()
        print " "
          
    def get_record (self): 
        return deepcopy(self._record_rv)
    
    def scan_through_rv (self,direction):
        if len(self._record_rv)==0: return
        
        if self._reci is None: 
            self._reci = 0
            return self._record_rv[self._reci][0]
            
        ind = 0
        if   direction == '+': ind = +1
        elif direction == '-': ind = -1
        
        self._reci = np.clip(self._reci+ind,0,len(self._record_rv)-1)
        return self._record_rv[self._reci][0]
        
    def calc_stats (self,verbose = True):
        N = len(self._record_rv)
        if N == 0:
            print "No lines have been recorded"
            return (0, 0, 0, 0)

        rv = np.array(map(iget(1),self._record_rv))

        rv_mean = np.mean(rv)
        rv_std = np.std(rv)
        rv_std_mean = rv_std/np.sqrt(N)
        if verbose: print "Radial Velocity = "+format(rv_mean,'10.3f')+" +- "+format(rv_std_mean,'5.3f')+"  STD = "+format(rv_std,'5.3f')+"  N = "+str(N)
        return (rv_mean, rv_std, rv_std_mean, N)

    def get_record_i (self):
        return deepcopy(self._reci)    

    def save_list (self,filename,clobber=True):
        if os.path.exists(filename) and not clobber: raise ValueError("File already exists: '"+filename+"'")
        
        lines = ["==> RV RECORD"]
        for i in xrange((self._record_rv)):
            wl, rv, info = self._record_rv[i][:3]
            line = "  ".join((format(wl,'>10.3f'),format(rv,'>10.3f'),info))
            lines.append(line)
            
        f = open(filename,'w')
        f.write("\n".join(lines)+"\n")
        f.close()
        
    def open_list (self,filename):
        
        f = open(filename)
        lines = f.readlines()
        f.close()
        if lines[0].find("==> RV RECORD") != 0: raise IOError("File format unexpected")
        
        self._record_rv = []
        for line in lines:
            line = line.rstrip()
            if line.strip() == '' or line.strip()[0] == '#': continue
            sline = line.split()            
            wl = float(sline[0])
            rv = float(sline[1])
            info = " ".join(sline[2:])
            self._record_rv.append([wl,rv,info])
            
        self._reci = len(self._record_rv)-1
        return self._reci
 
 
class RVLineLists:
    """
This holds the linelist information if you want to change
    """
    
    _datapath = os.path.abspath(path_2_eyeSpec[0])+"/eS-data/linelists/rv_linelists/"
    
    def __init__ (self,custom_linelist=None):
        
        self._llists = {}
        
        self._load_sparse_linelist()
        self._load_moderate_linelist()
        self._load_dense_linelist()

        self._llists['custom'] = []
  
  
  
    def add_linelist (self,which):
        which = str(which)
        if which in self.get_choices(): 
            print "HeadsUp: option already exists "+which
        self._llists[which] = []
      
    def get_datapath (self):
        return deepcopy(self._datapath)
      
    def read_rv_linelist (self,fname):
        fpath = self.get_datapath()+fname
        if not os.path.exists(fpath):
            print "Whoops, file does not exists: '"+fpath+"'"
            return []
        
        f = open(fpath)
        lines = f.readlines()
        f.close()
        llist = []
        
        for line in lines:
            line = line.rstrip()
            sline = line.strip().split()
            
            if len(sline) == 0 or sline[0][0] == '#': continue
            
            wl = float(sline[0])
            llist.append((wl," ".join(sline[1:])))

        f.close()
        return llist
        
    def _load_sparse_linelist (self):
        self._llists['sparse'] = self.read_rv_linelist("sparse_linelist.txt")
    
    def _load_moderate_linelist (self):
        self._llists['moderate'] = self.read_rv_linelist("moderate_linelist.txt")
    
    def _load_dense_linelist (self):
        self._llists['dense'] = self.read_rv_linelist("dense_linelist.txt")
    
    def get_choices (self):
        return self._llists.keys()
 
    def get_lines_info (self,which):
        if which not in self.get_choices(): return    
        lines = map(iget(0),self._llists[which])
        info = map(iget(1),self._llists[which])
        return lines, info
      
    def get_lines (self,which):
        if which not in self.get_choices(): return    
        lines = map(iget(0),self._llists[which])
        return lines

    def get_info (self,which):
        if which not in self.get_choices(): return
        info = map(iget(1),self._llists[which])
        return info
    
    def save_linelist (self,fname,which,clobber=True):
        if os.path.exists(fname) and not clobber: raise ValueError("File exists:'"+fname+"'")
        if which not in self.get_choices(): raise ValueError("Unknown option:"+which)
        f = open(fname,'w')
        llist = self._llists[which]
        for i in xrange(len(llist)):
            wl = format(llist[i][1],">10.3f")
            f.write(wl+" "+llist[i][0]+"\n")
        f.close()
                       
class InteractiveRVEditor:
    
    c = 299792.4580
    def __init__ (self,ax,spec_obj,rv_record,linelist=None,history=None):
        
        if history is None: self.history = HistoryAdvanced()
        else: self.history = history
        
        self.ax = ax
        self.dp = self.dataplot = eyeSpecBaseDataPlot(self.ax,spec_obj,'center')
               
        self._xrange = 10            
        self._shifted_once = False

        # record cursor click points
        self._recorded_pts = []
        self._recorded_lim = 20 # mostly only care about the -1 indexed point

        self._sum_rv = 0
        self.rv_record = rv_record

        self._rv_max = 1000
        self._apply_rv = False

        self._scale_rv = 1e-1
        self._start_data = None

        self._use_rv_type = 'mid' # 'old' or 'new'

        self._pressed = {'c':False}        
        
        self._linelists = RVLineLists(linelist)
        
        lines,info = self._linelists.get_lines_info('sparse')
        self.ile = eyeSpecBaseLineEditor(self.ax,lines,info,lock_order=True,lock_line_positions=True)        
        self.ile.text.set_visible(False)
        
        self._changed_line = False
        
        self.rv_text = self.ax.text(0.0, 1.05, 'Current Radial Velocity: 0.00',
                                      transform=self.ax.transAxes, va='top')

        check = self.scan_through_lines('+')
        if check is None: print "!! need a different line list"
            
        self.update_fixed_line()    
     
    def record_rv_history (self):
        curstate = State('rv')
        curstate['rv'] =  self.get_current_rv()

        self.history.add(curstate,self.undo_redo,id=id(self),info='changed the radial velocity')
        return curstate
        
    def undo_redo_rv (self,prevstate):
        prevstate.check_id('rv')
        curstate = State('rv')
        curstate['rv'] =  self.get_current_rv()
        
        rv = prevstate['rv']
        self.set_rv(rv,hist=False)
        
        return curstate
              
    def change_line_list (self,which):
        if which not in self._linelists.get_choices(): raise ValueError("Which must be one of: "+", ".join(choices))
        lines, info = self._linelists.get_lines_info(which)
        self.ile.change_linelist(lines,info)
        
    def _auto_scale_line_select (self):
        if not self.ile.sel_line.is_selected(): return
        
        xpt = self.ile.sel_line.get_xdata()
        xmin,xmax = self.ax.axis()[:2]        
        if (xmin < xpt < xmax): return
        halfx = self._xrange/2.0
        self.ax.set_xlim(xpt-halfx,xpt+halfx)      
    
    def _smart_auto_scale_line_select (self):
        # use windows of 2? angstroms find where the median is > 0.98
        
        # on both sides of the feature.
        
        # min(5A, right_wl) and same 
        
        # set the window based on that window
        pass
    
    
    def apply_rv (self,rv,hist=True):   
        if hist: self.record_rv_history()
        
        # rv_diff = rv - self._sum_rv
        self._sum_rv += rv

        self.dp.spec_obj.edit.apply_rv(rv)
        self.dp.plot_data.update_plot_data(self.dp.spec_obj)
        # self._sum_rv = rv

    def set_rv (self,rv,hist=True):
        if hist: self.record_rv_history()
        rv_diff = rv - self._sum_rv
    
        self.dp.spec_obj.edit.apply_rv(rv_diff)
        self.dp.plot_data.update_plot_data(self.dp.spec_obj)
        self._sum_rv = rv        

    def apply_stat_rv (self):
        rv_stat = self.rv_record.calc_stats()
        rv = rv_stat[0]
        self.set_rv(rv)
     

    def interactive_set_rv (self):
        dlg = wx.TextEntryDialog(None,"Please enter radial velocity","Text Entry","100.0",style=wx.OK|wx.CANCEL)
        if dlg.ShowModal() == wx.ID_OK:
            val = dlg.GetValue()
            val = val.strip()
            if val == '':
                print "No value given"
                dlg.Destroy()
                return

            try: rv = float(val)
            except:
                print "Invalid radial velocity"
                return
            self.set_rv(rv,hist=True)
        dlg.Destroy()
 
    def interactive_input_line (self):
        dlg = wx.TextEntryDialog(None,"Please enter wavelength and information:","Text Entry","6562.89 H alpha",style=wx.OK|wx.CANCEL)
        if dlg.ShowModal() == wx.ID_OK:
            val = dlg.GetValue()
            val = val.split()
            if len(val) == 0: 
                print "No value given"
                dlg.Destroy()
                return
            elif len(val) == 1: val.append("Unknown")
            
            try: wl = float(val[0])
            except:
                print "first value must be given as a float and then a space"
                return
            
            info = ""
            for i in range(1,len(val)): info+=val[i]+" "
            self.wavelengths_list.append([wl,info])
            self._wl_i = self._walk_backwards_check() 
            
            vymin = self.ax.axis()[2]
            vymax = self.ax.axis()[3]
            self.show_wl_lines._paths.append(Path([[wl,vymin],[wl,vymax]]))

        dlg.Destroy() 
        

    def record_rv (self):
        if not self.ile.sel_line.is_selected(): return
        wl = self.ile.sel_line.get_xdata()
        ind = self.ile.plotlines.get_index_by_wl(wl)
        info = self.ile.plotlines.get_line_info(ind)              
        self.rv_record.add(self._sum_rv, wl, info, True)

    def calc_rv (self,wl_fixed,wl_to_shift):
        delta = wl_fixed - wl_to_shift # difference in wavelength
        xmid = np.abs(wl_to_shift + wl_fixed)/2. # calculate the mid point
        
        if   self._use_rv_type == 'old': rv = delta*self.c/wl_to_shift # calculated based on old wavelength
        elif self._use_rv_type == 'new': rv = delta*self.c/wl_fixed # calculate based on new wavelength
        elif self._use_rv_type == 'mid': rv = delta*self.c/xmid # calculate based on the centeral
        
        return round(rv,20)

    def calc_wl_scaler (self,rv):
        # solve  delta_lambda/lambda = rv/self.c
        # delta_lambda = wl_new - wl_old
        # lambda in [wl_old, wl_new, wl_mid]
        # wl_mid = (wl_old + wl_new)/2.        
            
        # if you have del_lambda/wl_new = rv/self.c
        if self._use_rv_type == 'new': multi = 1./(1. - rv/self.c)
            
        # if you have del_lambda/wl_old = rv/self.c
        elif self._use_rv_type == 'old':  multi = (1. + rv/self.c)

        # if you have del_lambda/wl_mid = rv/self.c
        else: multi = (2. + rv/self.c)/(2. - rv/self.c)

        return multi 
     

    def scan_through_lines (self,direction):
        # looks for lines within the data
        dmin,dmax = self.dp.spec_obj.get_wlbounds()
        buffer = 10
        
        lmin,lmax = self.ile.plotlines.get_bounds()
                
        no_lines = False
        crash = 0
        while True:
            crash += 1
            if crash > 2.0*len(self.ile.plotlines): raise StandardError("Whoops, must be a logic problem")
            
            xpt = self.ile.scan_through_lines(direction,verbose=False)
            if xpt is None: continue
            
            # check if xpt is in a data range
            if dmin-buffer <= xpt <= dmax+buffer: break 
        
            # if you reach the end of the line list
            if xpt in [lmin,lmax]:    
                if xpt == lmin and direction == '+': continue
                if xpt == lmax and direction == '-': continue
                
                if no_lines:
                    print "No lines found which match the data"
                    return None

                # change direction
                if direction == '+':   direction = '-'
                else: direction = '+'
                no_lines = True
                
        self._auto_scale_line_select()             
        return xpt
    
    # !! history may be unnecessary

    def scan_through_rv_record (self,direction):
        self.rv_record.scan_through_rv(direction)
        self.rv_record.display_record()
        
    def get_current_rv (self):
        return deepcopy(self._sum_rv)
        
    def get_current_state (self):
        curstate = State("interactive_rv_editor")
        curstate['_sum_rv'] = self._sum_rv
        curstate['dataplot'] = self.dataplot.get_current_state()
        curstate['ile'] = self.ild.get_current_state()
        curstate['rv_rec'] = self.rv_record.get_current_state()
        
    def apply_state (self,prevstate): pass
    def undo_redo (self,prevstate): pass
    
    def delete_selected_rv (self):
        self.rv_record.delete()
    
    def get_rv_record (self):
        return self.rv_record.get_record()
 
    def key_move_data (self,event):
        if event.key == 'right': pass
        elif event.key == 'left': pass
        return
    
    def key_pressed (self,event,which):
        
        #    def key_release_callback (self,event):
        #        self.key_pressed(event,'release')
        #
        #    def key_press_callback (self,event):
        #        if self.key_pressed(event,'press'): return
                
        # which = 'press' 'release'
        # check the pressed keys
        for key in self._pressed:
            # if any are already pressed then return
            # and you are not clicking the current key
            if self._pressed[key] and key != event.key:
                return True
    
        if event.key in self._pressed:
            key = event.key
            if which == 'press': self._pressed[key] = True
            if which == 'release': self._pressed[key] = False
    
        bool_pressed = np.array(self._pressed.values())
        return (np.any(bool_pressed))

    def btn_press_rv (self,event): # !! perhaps could get rid of this 
        if event.inaxes is None: return False
        if event.button != 1: return False
        if self.dp.is_toolbar_button_on(): return False

        # record mouse clicks
        self._recorded_pts.append([event.xdata,event.ydata])
        if len(self._recorded_pts) > self._recorded_lim: del self._recorded_pts[0]
        
        
        c1 = False
        if self._pressed['c']: 
            epsx = self.dp.get_epsilon('x')
            c1 = self.ile.btn_press_lines(event,epsx)
            self._changed_line = c1
        
        return c1

    def mot_notify_rv (self,event): pass # possibly could drag data

    def btn_release_rv (self,event): pass
    
    def btn_click_rv (self,event):
        if event.button != 1: return False
        if event.inaxes is None: return False
        if self.dp.is_toolbar_button_on(): return False
        
        if self._changed_line: 
            self._changed_line = False
            return False
        
        if self.ile.sel_line.is_selected() and not self._pressed['c']:
            # shift data left and right
            wl_to_shift = event.xdata #self._recorded_pts[-1][0]
            wl_fixed = self.ile.sel_line.get_xdata()
            
            rv = self.calc_rv (wl_fixed,wl_to_shift)

            if np.abs(rv) < self._rv_max:
                self.apply_rv(rv) 
                return True 
            else: print "Radial Velocity is > "+str(self._rv_max)+" [km/s]"
        return False

    def user_input_line (self):
            dlg = wx.TextEntryDialog(None,"Please enter wavelength and information:","Text Entry","6562.89 H alpha",style=wx.OK|wx.CANCEL)
            if dlg.ShowModal() == wx.ID_OK:
                val = dlg.GetValue()
                val = val.split()
                if len(val) == 0: 
                    print "No value given"
                    dlg.Destroy()
                    return
                elif len(val) == 1: val.append("Unknown")

                try: wl = float(val[0])
                except:
                    print "first value must be given as a float and then a space"
                    return

                info = ""
                for i in range(1,len(val)): info+=val[i]+" "
                self.wavelengths_list.append([wl,info])
                self._wl_i = self._walk_backwards_check() 
            
                vymin = self.ax.axis()[2]
                vymax = self.ax.axis()[3]
                self.show_wl_lines._paths.append(Path([[wl,vymin],[wl,vymax]]))
                


            dlg.Destroy()
            
            # val_input = raw_input(
            # val_input = '2345.2'  => [2345.2,'mystry line']
            # val_input = '2345.2 H alpha cool' => [2345.2,'H alpha cool']

            # !! should I add the line to the wavelengths to check? perhaps insert at the current self._wl_i ?
    
    def update_fixed_line (self):
        N = len(self.ile.plotlines) -1
        if not self.ile.sel_line.is_selected(): wl_fixed_txt = 'Fixed Line (None/'+str(N)+")  : "
        else:
            wl = self.ile.sel_line.get_xdata()
            ind = self.ile.plotlines.get_index_by_wl(wl)
            info = self.ile.plotlines.get_line_info(ind) 

            wl_fixed_txt = "Fixed Line ("+str(ind)+"/"+str(N)+") at "+format(wl,'>6.3f')+"  :  "+info
            

        if 'wl_fixed_txt' not in dir(self):
            self.wl_fixed_txt = self.ax.text(0.40,1.05,wl_fixed_txt, transform=self.ax.transAxes, va='top')
            # originally at 0.42 and 1.05
            return True
            
        if wl_fixed_txt == self.wl_fixed_txt.get_text(): return False   
        self.wl_fixed_txt.set_text(wl_fixed_txt)
        return True

    def update (self):
        c1,c2 = False,False    
        new_txt = 'Current Radial Velocity: '+format(self._sum_rv,'5.3f')
        if new_txt == self.rv_text.get_text(): c1 = False
        else:
            self.rv_text.set_text(new_txt)
            c1 = True
            
        c2 = self.update_fixed_line()
        return (c1 or c2)
            
class EditRVManager (eyeSpecBaseEventManager):
    def __init__ (self,rv_panel,rv_record,spec_obj,linelist=None):
        eyeSpecBaseEventManager.__init__(self)
        self.ppanel = rv_panel
        self.ax = rv_panel.ax
        self._dragging = False
        
        self.rve = InteractiveRVEditor(self.ax,spec_obj,rv_record,linelist)
        
        self.init_connection_callbacks(self)
        
        self.key_cfg.add_key('d',"Same as 'backspace")
        self.key_cfg.add_key('esc',"same as 'q'")
        self.key_cfg.add_key('g','Toggle grid on/off')
        self.key_cfg.add_key('h','Display this screen')
        
        
        self.key_cfg.add_key('p',"Toggle pan/zoom tool (note: editor won't work while engaged")
        
        self.key_cfg.add_key('q','Close and return')
        
        
        self.key_cfg.add_key(";",'Toggle data between scatter and line plot options')
        
        self.key_cfg.add_key('`','Toggle auto scaling options')
        
        self.key_cfg.add_key('[','Undo data edit')
        self.key_cfg.add_key(']','Redo data edit')
        
        self.key_cfg.add_key('backspace','Delete selected radial velocity record')
        self.key_cfg.add_key("up/down",'Scroll through radial velocity records')
        
        self.key_cfg.add_key('b','Scan backwards through line list')
        self.key_cfg.add_key(' ','Scan forwards through line list')
        
        self.key_cfg.add_key('v','Apply the statistical radial velocity based on recorded radial velocities')
        #self.key_cfg.set_display_order(['backspace','d','esc','g','h','p','q','z','`',';','[',']',
        self.key_cfg.add_key('z',"Toggle zoom rect tool (note: editor won't work while engaged")
        
        
        
        self.key_cfg.set_display_order(['backspace','esc','g','h','p','q','v','b',' ','z','`',';','[',']','up/down'])
        
        self.key_cfg.check_display()
        self.key_cfg.add_mpl_key_convert('up','up/down')
        self.key_cfg.add_mpl_key_convert('down','up/down')
        
        
        self.key_cfg.add_key("space",'add keys')
        self.key_cfg.set_display_order(['space'])
        self.key_cfg.check_display()
        self.key_cfg.add_mpl_key_convert(' ','space')        

    def key_press_callback (self,event):
        if self.rve.key_pressed(event, 'press'): return False
        
    def key_release_callback (self,event): 
        if self.rve.key_pressed(event,'release'): return False
        if event.key == 'backspace': event.key = 'd'
        
        if event.key == 'h': self.display_help('Radial Velocity Editor')
        
        elif event.key == 'g': self.rve.dp.toggle_grid()
            
        elif event.key == 'p': self.rve.dp.toggle_pan_tool()
            
        elif event.key in ('q','esc') and 'Close' in dir(self.ppanel.pframe): self.ppanel.pframe.Close()
            
        elif event.key == 'z': self.rve.dp.toggle_zoom_tool()  
                  
        elif event.key == '`': self.rve.dp.set_auto_scaling() 
            
        elif event.key == ';': self.rve.dp.plot_data.toggle_data_line_scatter()
          
        elif event.key == 'alt': self.rve.record_rv()
          
        elif event.key == 'v': self.rve.apply_stat_rv()
          
        elif event.key == 'i': self.rve.interactive_set_rv()
        
        elif event.key == 'b':self.rve.scan_through_lines('-')
        
        elif event.key == ' ':self.rve.scan_through_lines('+')
        
        elif event.key == 'd': self.rve.delete_selected_rv()
        
        
        elif event.key == 'l':
            dlg = EditSubplotDialog(self.ax)
            dlg.ShowModal()
            dlg.Destroy()  
        
        # elif event.key == : display stats on rv records
        # if event.key == '.': scan forward
        # if event.key == ',': scan backwards
        # elif event.key == : change line list  
        elif event.key == '[': self.rve.history.undo()   
        elif event.key == ']': self.rve.history.redo()        
        elif event.key == 'up': self.rve.scan_through_rv_record("-")
        elif event.key == 'down': self.rve.scan_through_rv_record("+")
        
        self.update()
        
    def button_press_callback (self,event):
        self._dragging = False
        if self.rve.btn_press_rv(event): self.update()

    def motion_notify_callback (self,event):
        self._dragging = True
        if self.rve.dp.is_toolbar_button_on(): return
        if self.rve.dp.bounds_changed(True): self.update()
        
    def button_release_callback (self,event):
        if not self._dragging: self.button_click_callback(event)
        self._dragging = True
        self.update()
          
    def button_click_callback (self,event): 
        self.rve.btn_click_rv(event)
        self.update()
        
    def plot_change_update (self): pass
    
    def update (self):
        c1 = self.rve.update()
        c2 = self.rve.ile.update()
        c3 = self.rve.dp.update()
        if c1 or c2 or c3: self.ax.figure.canvas.draw()
        
class EditRVPanel (eyeSpecBaseDataPanel):
    def __init__ (self,rv_main_panel,rv_frame,rv_record):
        eyeSpecBaseDataPanel.__init__(self,rv_main_panel,rv_frame,False)
        self.spec_obj = self.pframe.inputs[0]
        self.linelist = self.pframe.inputs[1]
        self.rvManager = EditRVManager(self,rv_record,self.spec_obj,self.linelist)
        self.rvManager.disconnect()
        
        del self.canvas.callbacks.callbacks['motion_notify_event'][self.statusbar_cid]
        self.canvas.mpl_connect('motion_notify_event',self.RVUpdateStatusBar)
        
    def RVUpdateStatusBar (self,event):
        scale_txt = "Auto Scale "
        scale_opt,info = self.rvManager.rve.dp.get_auto_scaling_opt()
        if scale_opt == 0: scale_txt += 'X,Y'
        elif scale_opt == 1: scale_txt += 'X'
        elif scale_opt == 2: scale_txt += 'Y'
        elif scale_opt == 3: scale_txt += 'None'
        
        cur_vrad = self.rvManager.rve.get_current_rv()
        rvstring = " |  RV="+format(cur_vrad,'<6.2f')


        st = format(scale_txt,'16') + rvstring
        self.UpdateStatusBar(event,st)        

    def OnStart (self,event):
        
        del self.canvas.callbacks.callbacks['key_press_event'][self._onkeystart_cid]
        del self.canvas.callbacks.callbacks['button_press_event'][self._onbutstart_cid]
        
        print ""
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "QUESTIONS TO USER:"
        questions = ["o  please let me know if you have any clever thoughts on how to make this better",
                     "o  do you like the way it auto scales?"]
        print ("\n".join(questions))

        print ""

        self.rvManager.rve.ile.select_line_by_index(0)
        self.rvManager.rve.scan_through_lines('+')
        self.rvManager.rve.scan_through_lines('-')
        xcenter = self.rvManager.rve.ile.get_selected_x()

        xran = 50.0
        self.ax.set_xlim(xcenter-xran/2.0,xcenter+xran/2.0)

        self.rvManager.connect()
        self.rvManager.update()
  
class EditRVMainPanel (eyeSpecBaseMainPanel):
    def __init__ (self,parent_frame):
        eyeSpecBaseMainPanel.__init__(self, parent_frame,'split_top_left')        
        self.rv_record = RVRecord()
        self.rvpanel = EditRVPanel(self.Split1,self.pframe,self.rv_record)
         
        self.randpanel = RandomPanel(self.Split1,self.pframe)
        self.split_top_left(self.rvpanel,self.randpanel)
        



class EditRVFrame (eyeSpecBaseFrame):
    def __init__ (self,parent_window,inputs):
        self.inputs = inputs
        self.spec_obj = inputs[0]
        title = 'Edit Raidal Velocity: '+os.path.basename(self.spec_obj.filename)
        eyeSpecBaseFrame.__init__(self, parent_window, title)    
        self.panel = EditRVMainPanel(self)    

    def OnFinish (self):
        self.Backup()
        spec_obj = self.panel.rvpanel.rvManager.rve.dp.spec_obj
        rv_record = self.panel.rvpanel.rvManager.rve.get_rv_record()
        return spec_obj, rv_record
       
    def Backup (self):
        print "Closing: If this hangs up look at the file TMP_WHAT_JUST_HAPPENED.txt"
        # you need these to recover the data after app.MainLoop() closes
        rv_stat = self.panel.rvpanel.rvManager.rve.rv_record.calc_stats()

        f = open("TMP_FINAL_VRAD.txt",'w')
        f.write("# this is the final radial velocity (mean,std_of_mean,std,N)\n")
        _rv_what_happened()
        f.write("  ".join([format(rv_stat[0],'30.15f'), # mean
                           format(rv_stat[2],'30.15f'), # stdev of the mean
                           format(rv_stat[1],'30.15f'), # stdev
                           str(int(rv_stat[3])), # N
                           "\n"]))
        f.close()
        time.sleep(1)     

pass
###############################################################################

def edit_rv (spec_obj,linelist=None,clean_up=True):
    """

PURPOSE:
    This takes an eyeSpec spectrum object and allows for visual radial velocity setting
   
CATEGORY:
    Spectral Reductions

INPUT ARGUMENTS:
    spec_obj : (eyeSpec_spec) An eyeSpec spectrum object

INPUT KEYWORD ARGUMENTS:
    linelist : (array or None) will take lines to compare against [[wl,info],[wl,info],...]
    clean_up : (boolean) If true then it will remove temporary files it creates, files 'TMP_*'

OUTPUTS:
    (eyeSpec_spec) Edited eyeSpec spectrum object
    
DEPENDENCIES:
   External Modules Required
   =================================================
    Numpy, os, time, operator, copy, scipy, math,
    matplotlib, wxpython
   
   External Functions and Classes Required
   =================================================
    SysOutListener, eyeSpecBaseApp, eyeSpecBaseFrame, eyeSpecBaseMainPanel, eyeSpecDataPanel, eyeSpecBaseEventManager
    EventConnections, InteractiveDataEditor
    ProgressSave, InputOutput, History, KeyboardConfiguration,
    save_spec, load_spec, find_overlap_pts, alt_order_colors
    
   
NOTES:
   (1) if you know the radial velocity (e.g. 100 km/s) you want you can just apply it
    spec_obj.edit.apply_rv(100)

EXAMPLE:
   >>> spec = readin("stardata.fits")
   >>> edited_spec = edit_rv(spec)

MODIFICATION HISTORY:
    13, Jun 2013: Dylan Gregersen

    """

    if linelist is not None: print "Can't currently declare a linelist"

    if spec_obj.__class__.__name__ != 'eyeSpec_spec': raise ValueError("spec_obj MUST BE OF CLASS eyeSpec_spec")
        
    edit_obj = spec_obj.copy()
    save_spec(spec_obj,filename='TMP_OBJ_SAVE_ORIG',clobber=True)
    
    # set_title = 'Set Radial Velocity for: '+os.path.basename(edit_obj.filename)

    # def _run_app (edit_obj,line_data):
    ##########################################
    # run application
    # _app_run_rv(edit_obj)
    app = eyeSpecBaseApp(EditRVFrame,[edit_obj,None])
    sys.stdout = SysOutListener()
    try: app.MainLoop()
    finally:
        app.ExitMainLoop()
        final_spec,rv_record = app.Finish()
        del app


    ##########################################    
    
    line =  "-"*10+format("Set Radial Velocity Complete",'^32')+"-"*10
    print "-"*len(line)
    print line

    # load data after app.MainLoop() has exited
    if os.path.exists("TMP_FINAL_VRAD.txt"): 
        rv_stat = np.loadtxt("TMP_FINAL_VRAD.txt")
        print " "
        print "FINAL RADIAL VELOCITY = "+format(rv_stat[0],'8.3f')+" +- "+format(rv_stat[1],'5.3f')+"   STD = "+format(rv_stat[2],'5.3f')+"   N = "+str(int(rv_stat[3]))

        # !! could also have function which determines the heliocentric radial velocity
        edit_obj.edit.apply_rv(rv_stat[0])

    # clean up temporary files
    if clean_up:
        if os.path.exists('TMP_WHAT_JUST_HAPPENED.txt'): os.system('rm TMP_WHAT_JUST_HAPPENED.txt')
        if os.path.exists('TMP_OBJ_SAVE_ORIG.pkl'): os.system('rm TMP_OBJ_SAVE_ORIG.pkl')        
        if os.path.exists("TMP_FINAL_VRAD.txt"): os.system('rm TMP_FINAL_VRAD.txt')
    
    return edit_obj.copy(), rv_record, rv_stat
    





