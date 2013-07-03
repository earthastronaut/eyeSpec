
if __name__ != '__main__':
    from eyeSpec.interactive_IO import ProgressSave, InputOutput, SaveOpenChoices
    from eyeSpec.interactive_classes import EventConnections, Cursor, History, KeyboardConfiguration, PlotLines, eyeSpecBaseEventManager, eyeSpecBaseDataPanel, eyeSpecBaseMainPanel, eyeSpecBaseFrame
    from eyeSpec.extended_IO import save_spec, load_spec
    from eyeSpec.base_functions import find_overlap_pts, inv_var_2_var, var_2_inv_var,  alt_order_colors, np_vstack_append
    from eyeSpec.base_classes import Timer
    from eyeSpec.app_edit_data import InteractiveDataEditor
    from eyeSpec.linefinder_dsg import continuum_finder
    from eyeSpec.dependencies import (np, os, sys, time, iget, deepcopy, pdb, scipy, math,
                                      plt, FormatStrFormatter, savefig,
                                      pyfits, pickle,
                                      wx, FigureCanvas, NavigationToolbar2Wx, Figure, Button, Path)

    




# this takes a list of known wavelengths in stationary frame, plots the data against those, then allows for a simple click to radial velocity shift the data to match that point.

# it will walk through each of the stationary points until you quit or it's done
# have the option to pull up a new plot of the solar spectrum for that particular region!!! Awesome

# !! add the solar and arcturus data
# create a directory with the solar data broken into easy to read pieces
# can add arcturus, alpha cen,


def _snr_what_happened ():
    # general
    f = open('TMP_WHAT_JUST_HAPPENED.txt','w')
    f.write("This file is intended to help in the case that the code crashes and you're left with only TMP_* files\n")
    f.write("If this isn't helpful please email Dylan Gregersen <dylan.gregersen@utah.edu>\n")

    f.write("\n")

    # What was saved? Why?
    f.write("This file was created from the set SNR routine\n")
    f.write("The original spectrum information is in the file: TMP_OBJ_SAVE_ORIG.pkl\n")
    f.write("The edited spectrum information is in the file: TMP_OBJ_SAVE_EDIT.pkl\n")
    f.write("\n")
    
    # How to recover with these files.
    f.write("To restore what you were working on:\n")
    f.write(">>> from eyeSpec import load_spec")
    f.write(">>> orig_obj = load_spec('TMP_OBJ_SAVE_ORIG.pkl')\n")
    f.write(">>> spec_obj = load_spec('TMP_OBJ_SAVE_EDIT.pkl')\n")
    
    f.close()

class PlotInverseVariance (EventConnections):
    def __ini__ (self,parent):
        # create a small window with the plot data
        pass
    
class RegionGuesserDialog (wx.Dialog, EventConnections):
    
    def __init__(self, guess_snr, *args, **kw):
        super(RegionGuesserDialog, self).__init__(*args, **kw) 
        
        self.guess_snr = guess_snr

        self.InitUI()
        self.SetSize((400, 100))
        # !! make this smarter, based off the size of the parent
        self.SetPosition((860,30))

        self.SetTitle("Parameters for SNR regions guess")

        self.Bind(wx.EVT_WINDOW_DESTROY,self.OnClose)
        self.Bind(wx.EVT_KEY_DOWN,self.OnKeyDown)

        self.key_cfg = KeyboardConfiguration()
        self.key_cfg.add_key('Enter','Apply guessing routine')
        self.key_cfg.add_key('h','Display this help screen')
        self.key_cfg.add_key('q','Quit sub-window')

        self.key_cfg.set_display_order(['Enter','h','q'])
        self.key_cfg.check_display()
        self.key_cfg.add_mpl_key_convert('enter','Epace')

        
    def InitUI(self):

        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        sb = wx.StaticBox(pnl, label='Parameters')
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)        
#         sbs.Add(wx.RadioButton(pnl, label='Null', 
#             style=wx.RB_GROUP))
#         sbs.Add(wx.RadioButton(pnl, label='Null1'))
#         sbs.Add(wx.RadioButton(pnl, label='Null2'))
        
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)        
        #hbox1.Add(wx.RadioButton(pnl, label='GUESS SNR:'))
        self.txtc = (wx.TextCtrl(pnl))
        hbox1.Add(txtc, flag=wx.LEFT, border=5)
        sbs.Add(hbox1)
        
        pnl.SetSizer(sbs)
       
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        okButton = wx.Button(self, label='Ok')
        closeButton = wx.Button(self, label='Close')
        hbox2.Add(okButton)
        hbox2.Add(closeButton, flag=wx.LEFT, border=5)

        vbox.Add(pnl, proportion=1, 
            flag=wx.ALL|wx.EXPAND, border=5)
        vbox.Add(hbox2, 
            flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)

        self.SetSizer(vbox)
        
        okButton.Bind(wx.EVT_BUTTON, self.set_initial_parameters)
        closeButton.Bind(wx.EVT_BUTTON, self.OnClose)
        


    def OnKeyDown (self,event):
        key = self.key_cfg.convert_wx_code_2_key(event.GetKeyCode())
        if key is None: return

        if not self.key_cfg.check_key_press_callback(key):
            return

        if key == 'enter': print "!! run the guesser"
        if key == 'h': self.display_help()
        if key == 'q': self.Close()


    def display_help (self):
        print "="*60
        editor_name = 'Parameters for SNR regions guesses'
        print ("\n<<<< "+format(editor_name,'^26')+" >>>>")
        # !! put in info about how this works
        print ("-"*60)
        self.key_cfg.display_short_info()

    def set_initial_parameters (self,event):
        self.Close()

    def OnClose(self,event):
        self.Close()
                            
class HighlightSNRSelection:
    def __init__ (self,ax):
        self.ax = ax
        self.selection_box = ax.add_patch(plt.Rectangle((0,0),1,1,facecolor='none',lw=2.5,edgecolor='y',alpha=.3,visible=False,zorder=-10))
        self.snri = None
        self._snr_edit_region = 'none'
        self._starting_bounds = np.zeros(4)
        #         none
        #         top_left     top      top_right
        #         left        center        right
        #         bot_left     bot      bot_right

    def get_bounds (self):
        get_bounds = [self.selection_box.get_bbox().xmin, self.selection_box.get_bbox().xmax,
                      self.selection_box.get_bbox().ymin, self.selection_box.get_bbox().ymax]
        return get_bounds

    def set_selection_box (self,snr_bounds):
        if snr_bounds is None: return
        if self.snri is None: return
        
        xb1,xb2,yb1,yb2 = snr_bounds 
        xmin,xmax,ymin,ymax = self.ax.axis()
        xran = xmax-xmin
        yran = ymax-ymin

        scaley = 0.02
        scalex = 0.01

        deltax = scalex*xran
        deltay = scaley*yran

        nx = xb1-deltax
        ny = yb1-deltay
        nwid = abs(xb2-xb1)+2.0*deltax
        nhei = abs(yb2-yb1)+2.0*deltay
        
        self.selection_box.set_bounds(nx,ny,nwid,nhei)

    def get_snr_index (self):
        return self.snri
        
    def set_highlight_visible (self,snr_index,snr_bounds = None):
        visible = True
        self.snri = snr_index
        if self.snri is None or snr_bounds is None: visible = False
        else: self.set_selection_box(snr_bounds)
        self.selection_box.set_visible(visible)
        
    def deselect (self):
        self.set_highlight_visible(None)
        
    def get_visible (self):
        visible = self.selection_box.get_visible()
        if self.snri is None: 
            if visible: self.selection_box.set_visible(False)
            visible = False
            self._starting_bounds = np.zeros(4)
        return visible
    
    def get_edit_region (self):
        return self._snr_edit_region
        
    def get_starting_bounds (self):
        return deepcopy(self._starting_bounds)
   
    def set_edit_region (self,xy,inner_bounds,force_set = None):
        x,y=xy

        if self.get_snr_index() is None: 
            self._snr_edit_region = 'none'
            return  'none'

        # get bounds for the selected box
        bb = inner_bounds
        self._starting_bounds = bb

        if force_set is not None:
            self._snr_edit_region = force_set
            return force_set

        # get bounds for the highlighted box
        hbb = self.get_bounds()
        
        self._snr_edit_region = 'check'

        # check is it outside the highlighted box
        if xy not in self: self._snr_edit_region = 'none'

        # check if it's inside the selected box
        elif (bb[0] <= x <= bb[1]) and (bb[2] <= y <= bb[3]): self._snr_edit_region = 'center'
        
        # check for each of the regions
        elif (bb[0] < x < bb[1]) and (y >= bb[3]): self._snr_edit_region = 'top'
        elif (bb[0] < x < bb[1]) and (y <= bb[2]): self._snr_edit_region = 'bot'

        elif (x >= bb[1]) and (bb[2] < y < bb[3]):  self._snr_edit_region = 'right'
        elif (x <= bb[0]) and (bb[2] < y < bb[3]):  self._snr_edit_region = 'left'

        elif (x >= bb[1]) and (y >= bb[3]): self._snr_edit_region = 'top_right'
        elif (x <= bb[0]) and (y >= bb[3]): self._snr_edit_region = 'top_left'

        elif (x >= bb[1]) and (y <= bb[2]): self._snr_edit_region = 'bot_right'
        elif (x <= bb[0]) and (y <= bb[2]): self._snr_edit_region = 'bot_left'
        
        if self._snr_edit_region == 'check': raise StandardError("Whoops, this should have been set by one of the above. Recheck the logic")
        return self._snr_edit_region   
   
    def __contains__ (self,xy):
        xpt,ypt = xy
        hbb = self.get_bounds()
        return ((xpt > hbb[0]) and (xpt < hbb[1]) and (ypt > hbb[2]) and (ypt < hbb[3]))
                     
class SNRSelection:

    def __init__ (self,ax,inputs,visible=True,color='k'): 
        self._timer_stats = Timer(0.3)       
        self.ax = ax

        try: order, selection_limits, i = inputs
        except: 
            prev_state = inputs
            self.apply_prev_state(prev_state, init=ax)
            return            
                  
        zorder = -1
        self.order_i = i # associated order
        #----------------------------------------------------------------------#
        # calculate statistics on the selected data
        bbounds,dbounds = self.selection_limits_2_values(selection_limits)
        self.calc_stats(order,dbounds,set_stats=True)
               
        #----------------------------------------------------------------------#
        # create spaces for the plot data
        self.plotit(ax, selection_limits, zorder, color, visible)
               
    def __repr__ (self):
        return "SNR_Region"
    
    def __contains__ (self,val):
        xpt,ypt = val
        xmin,xmax,ymin,ymax = self.get_bounds()
        check = (xmin <= xpt) and (xpt <= xmax) and (ymin <= ypt) and (ypt <= ymax)
        return check

    def set_window_visible (self,xbounds=None,truth=True):
        if xbounds is None:
            xbounds= self.ax.axis()[:2]
        xmin,xmax = xbounds
        bb = self.get_bounds()
        checks = np.array([bb[0]<xmax,
                           bb[1]>xmin])
        if np.all(checks): self.set_visible(truth)
        else: self.set_visible(truth)

    def get_bounds (self):
        """ return the bounds of the box which defines the SNR region [xmin,xmax,ymin,ymax] """
        bbounds = [self.selection_box.get_bbox().xmin,self.selection_box.get_bbox().xmax,   
                   self.selection_box.get_bbox().ymin,self.selection_box.get_bbox().ymax]
        return bbounds

    def selection_limits_2_values (self,selection_limits):
        sl = selection_limits
        xmin,xmax = min(sl[:2]),max(sl[:2])
        ymin,ymax = min(sl[2:]),max(sl[2:])
        wid = xmax - xmin
        hei = ymax - ymin
        
        bbounds = [xmin,ymin,wid,hei]
        dbounds = [xmin,xmax,ymin,ymax]
        return [bbounds,dbounds]

    def update_plot_objects (self,selection_limits):
        sl = selection_limits
        self.update_box(sl)
        self.update_text(sl)
        self.update_mean_line(sl)

    def update_box (self,selection_limits):
        if self.selection_box is None: return
        xmin,ymin,wid,hei = self.selection_limits_2_values(selection_limits)[0]
        self.selection_box.set_bounds(xmin,ymin,wid,hei)
        
    def update_text (self,selection_limits):
        if self.selection_box is None: return
        if self.txt is None: return
        
        bbounds,dbounds = self.selection_limits_2_values(selection_limits)
        
        xmin,xmax,ymin,ymax = dbounds
        wid,hei = bbounds[2:]
        
        
        yval = ymax-0.01*abs(hei)
        ybottom,ytop = self.ax.axis()[2:]
        yran = ytop - ybottom
        
        yval = min([yval, ytop-0.1*yran]) 
        
        self.txt.set_x(xmin+0.05*abs(wid))
        self.txt.set_y(yval)
        self.txt.set_text(format(self.snr,'3.3f'))

    def update_mean_line (self,selection_limits=None):
        if self.selection_box is None: return
        if self.meanline is None: return
        self.meanline.set_ydata([self.stats[0],self.stats[0]])
        if selection_limits is not None:
            bbounds,dbounds = self.selection_limits_2_values(selection_limits)
            xmin,xmax = dbounds[:2] 
            self.update_mean_line_x(xmin, xmax)  
    
    def update_mean_line_x (self,xmin,xmax):
        self.meanline.set_xdata([xmin,xmax])    
    
    def plotit (self,ax,selection_limits,zorder=0,color='b',visible=True,**mpl_kwargs):
        """  add the SNR region selection box with the mean line and SNR text """
        bbounds,dbounds = self.selection_limits_2_values(selection_limits)
        x,y,wid,hei = bbounds
        
        alpha = 0.35
        self.selection_box = ax.add_patch(plt.Rectangle((x,y),wid,hei,
                                                        zorder=zorder,facecolor=color,alpha=alpha,visible=visible,
                                                        **mpl_kwargs))
        
        xmin = min([x,x+wid])
        xmax = max([x,x+wid])
        ymax = max([y,y+hei])
        
        ymean = self.stats[0]
        self.meanline, = ax.plot([xmin,xmax],[ymean,ymean],linestyle='-',lw=2,
                                 zorder=zorder+1,color=color,alpha=alpha,visible=visible,
                                 **mpl_kwargs)

        self.txt = ax.text(xmin+0.05*abs(wid),
                           ymax-0.01*abs(hei),
                           format(self.snr,'3.3f'),
                           color=color,fontsize='small',
                           horizontalalignment='left',verticalalignment='top',
                           zorder=zorder+2)

    def get_zorder (self):
        """
        This has zorder from the selection box +2
        """
        return self.selection_box.get_zorder()

    def set_zorder (self,zorder):
        zorder = int(zorder)
        self.selection_box.set_zorder(zorder)
        self.meanline.set_zorder(zorder+1)
        self.txt.set_zorder(zorder+2)

    def get_color (self):
        return self.meanline.get_color()
    
    def set_color (self,color):
        if self.selection_box is not None: self.selection_box.set_color(color)
        if self.meanline is not None: self.meanline.set_color(color)
        if self.txt is not None: self.txt.set_color(color)        

    def set_visible (self,truth):
        """ Define if the plot objects are visible  """
        truth=bool(truth)
        self.selection_box.set_visible(truth)
        self.meanline.set_visible(truth)
        self.txt.set_visible(truth)

    def get_visible (self):
        if self.selection_box is None: return False
        else: return self.selection_box.get_visible()

    def get_stats (self):
        return deepcopy(self.stats)

    def set_stats (self,stats):
        self.stats = stats
        self.snr = stats[4]
      
    def calc_mean_std (self,data,inv_var=None):
        N = len(data)
        if N == 0: return np.zeros(3)
        
        if inv_var is None: meantype ='simple'
        else: meantype = 'weighted'
        
        if meantype == 'weighted':
            inv_var_sum = np.sum(inv_var)
            if inv_var_sum == 0: inv_var_sum = 1e-30
            sel_mean = np.sum(inv_var*data)/inv_var_sum
            sel_std = np.sqrt(np.sum(inv_var*(data-sel_mean)**2)/inv_var_sum)

        elif meantype == 'simple':
            sel_mean = np.mean(data)
            sel_std = np.std(data)
        
        return sel_mean, sel_std, N
        
    def calc_stats (self,order,data_bounds=None, meantype='simple',mode='full',set_stats=False):
        if meantype not in ['weighted','simple']: raise TypeError("Whoops, wrong input")
        if mode not in ['mean','full']: raise TypeError("Whoops, wrong input")

        # get the data
        wl, data, inv_var = order

        if data_bounds is None:
            sel_wl = wl
            sel_data = data
            sel_inv_var = inv_var
        else:
            # if data bounds are given then crop based on those
            xmin,xmax,ymin,ymax = data_bounds
            xmask = (wl > xmin)*(wl < xmax)
            ymask = (ymin < data)*(data < ymax)
            mask = xmask*ymask

            sel_wl = wl[mask]
            sel_data = data[mask]
            sel_inv_var = inv_var[mask]

        # crop out any zeros from the remainder
        zeromask = (sel_inv_var != 0)
        if not np.any(zeromask):
            sel_wl = sel_wl[zeromask]
            sel_data = sel_data[zeromask]
            sel_inv_var = sel_inv_var[zeromask]
        
        # if no points are left then return
        N = len(sel_data)
        if N == 0: 
            if set_stats: self.set_stats(np.zeros(5))
            return np.zeros(5)
        
        # calculate mean, sigma, and N
        if meantype == 'weighted':
            sel_mean,sel_std,sel_N = self.calc_mean_std(sel_data,sel_inv_var)
        if meantype == 'simple':
            sel_mean,sel_std,sel_N = self.calc_mean_std(sel_data)

        # calculate the standard deviation of the mean
        sel_std_mean = sel_std/np.sqrt(sel_N)
        
        # calculate the signal to noise ratio (SNR)
        if sel_std != 0: snr = sel_mean/sel_std
        else: snr = np.NaN
 
        # return all calculated stats
        output = (sel_mean,sel_std_mean,sel_std,sel_N,snr)
        if set_stats:self.set_stats(output)
        else: return output

    def get_stats_dict (self,order,data_bounds=None, meantype='weighted',mode='full'):
        sel_mean,sel_std_mean,sel_std,sel_N,snr = self.calc_stats(order,data_bounds,meantype,mode)
        out_stats_dict = {'mean':sel_mean,
                          'std_mean':sel_std_mean,
                          'std':sel_std,
                          'N':sel_N,
                          'SNR':snr}
                    
        return out_stats_dict
              
    def get_current_state (self):
        the_state = {'state':'snr region selection',
                     'stats':self.stats,
                     'snr':self.snr,
                     'bounds':self.get_bounds(),
                     'vis':self.selection_box.get_visible(),
                     'zorder':self.selection_box.get_zorder(),
                     'ordi':self.order_i,
                     'color':self.get_color()}
        return the_state
        
    def apply_prev_state (self,prev_state,init=None):
        assert prev_state['state']=='snr region selection'
        self.stats = prev_state['stats']
        self.snr = prev_state['snr']
        bbounds = prev_state['bounds']
        zorder = prev_state['zorder']
        color = prev_state['color']
        vis = prev_state['vis']
        
        if init is not None:
            ax = init
            self.plotit(ax,bbounds,zorder,visible=vis,color=color)
        else:
            self.update_plot_objects(bbounds)
            self.set_zorder(zorder)
            self.set_visible(vis)
            self.set_color(color)
            
        self.order_i = prev_state['ordi']
                            
class OrderSNRSelections:
    def __init__ (self,ax,ordi,xbounds):
        self.list = []
        self._list_visible = np.array([],dtype=bool)
        self.all_bounds = np.array([],dtype=float) # a numpy array [[xmin,xmax,ymin,ymax],[xmin,xmax....]]
        
        self.order_index = ordi
        self.ax = ax
        self._xbounds = xbounds

         
    def append (self,item):
        if repr(item) != 'SNR_Region': 
            print "Input with wrong repr:",repr(item)
            return
        
        bounds = item.get_bounds()
        self.all_bounds = np_vstack_append(self.all_bounds,bounds)
                
        vis = item.get_visible()
        self._list_visible = np.append(self._list_visible,vis)
        
        self.list.append(item)
    
    def __str__ (self):
        return "OrderForRegions : "+str(tuple(self.all_bounds))
    
    def __repr__ (self):
        return "OrderForRegions ("+str(len(self))+")"
    
    def __contains__ (self,snri):
        return (snri in xrange(len(self.list)))
            
    def __len__ (self):
        return len(self.list)
    
    def __getitem__ (self,key):
        if key is None: print "Whoops, input snr index is none"
        return self.list[key]

    def __setitem__ (self,index,item):
        self.list[index] = item
        bounds = item.get_bounds()
        self.all_bounds[index] = bounds
        self._list_visible[index] = item.get_visible()
        
    def __iter__ (self):
        return iter(self.list)

    def __delitem__ (self,index):
        data_bounds = self.list[index].get_bounds()
        abounds = self.all_bounds[index]
        if not np.all(data_bounds == abounds): pdb.set_trace()
        
        del self.list[index]
        
        self.all_bounds = np_vstack_delete(self.all_bounds,index)
        self._list_visible = np.delete(self._list_visible,index)
        
    def any_visible (self):
        return np.any(self._list_visible)
 
    def set_window_visible (self,xbounds,truth=True):
        xmin,xmax = xbounds
        if len(self.all_bounds) == 0: return False
        dmins,dmaxs = self.all_bounds.T[:2]        
        
        inmask = (xmin < dmaxs)*(dmins < xmax)
        outmask = np.logical_not(inmask)
        
        vis = (self._list_visible == True)        
        not_truth = (self._list_visible != truth)

        ran = np.arange(len(self._list_visible))

        ran_2_off = ran[outmask*vis]
        ran_2_change = ran[inmask*not_truth]
        
        #!!print "order regions: ",len(self._list_visible),len(ran_2_off),len(ran_2_change)
        changed = False
        # which are visible and out, deselect
        for i in ran_2_off:
            self.list[i].set_visible(False)
            self._list_visible[i] = False
            changed = True

        # for those in which don't equal the truth
        #for i in ran_off:
        for i in ran_2_change:
            self.list[i].set_visible(truth)
            self._list_visible[i] = truth
            changed = True
        
        return changed
              
    def get_xbounds (self):
        return self._xbounds
              
    def get_order_index (self):
        return deepcopy(self.order_index)
    
    def iter_bounds (self):
        return iter(self.all_bounds)

    def update_snr_box_bounds (self,snri,selection_limits):
        if snri not in self: return
        self.list[snri].update_box(selection_limits)
        if len(self.all_bounds) == 0: self.all_bounds = np.array([selection_limits])
        else: self.all_bounds[snri] = selection_limits
        
    def hide_snr_regions (self,force=False):
        if force: ran = xrange(len(self.list))
        else: ran = np.argwhere(self._list_visible == True)[0]
        for i in ran:     
            self._list_visible[i] = False
            self.list[i].set_visible(False)

    def check_for_overlap (self,new_bounds,exclude=None):
        if len(self.all_bounds) == 0: return False
        
        to_exclude = (exclude != None)
        
        xb1,xb2,yb1,yb2 = new_bounds

        for i in xrange(len(self.all_bounds)):
            if to_exclude and i in exclude: continue
            bb = self.all_bounds[i]
            if (bb[0] < xb2) and (bb[1] > xb1) and (bb[2] < yb2) and (bb[3] > yb1): 
                print "HeadsUp: You can't overlap two selection regions"
                return True
 
        return False

    def sort_by_wavelength (self,snri=None):
        if len(self.list) == 0: return None
        
        region_mins = self.all_bounds.T[0]
        if len(self.all_bounds) == 1: return 0
        
        sortit = np.argsort(region_mins)
        # if no sorting needed then don't
        if np.all(sortit == np.arange(len(region_mins))): return snri
        
        if snri is not None: snri = np.argmin(np.abs(sortit - int(snri)))
        new_list = []
        
        for i in sortit: new_list.append(self.list[i])
        self.all_bounds = self.all_bounds[sortit]
        self._list_visible = self._list_visible[sortit]

        self.list = new_list
        return snri

    def get_current_state (self):
        region_states = []
        for i in xrange(len(self.list)):
            cur_state = self.list[i].get_current_state()
            region_states.append(cur_state)
            
        prev_state = {'state':'order snr regions',
                      'region_states':region_states,
                      'ordi':self.order_index,
                      '_list_visible':self._list_visible}        
        return prev_state
        
    def apply_prev_state (self,prev_state):
        assert prev_state['state'] == 'order snr regions'
        self.hide_snr_regions(True)
        
        region_states = prev_state['region_states']
        self.list = []
        self.all_bounds = np.array([],dtype=float)
        self._list_visible = np.array([],dtype=bool)
                
        for pstate in region_states:
            snr_region = SNRSelection(self.ax,pstate)
            snr_region.set_window_visible(self.ax.axis()[:2])
            self.append(snr_region)            
            
        self.order_index = prev_state['ordi']
        self.sort_by_wavelength() 

        

################################################################################

class SpectrumSNRSelections (dict):
    def __init__ (self,ax):
        #dict.__init__(*args,**kwargs)
        self._order_bounds = np.array([],dtype=float)
        self._order_visible = np.array([],dtype=bool)
        self.ax = ax
        
        
    def __setitem__ (self,ordi,oss):
        ordi = int(ordi) # should be in the spec_obj range
        if repr(oss).find("OrderForRegions") == -1: raise ValueError("Trying to set incorrect value :"+repr(oss))
        xbounds = oss.get_xbounds()
        self._order_bounds = np_vstack_append(self._order_bounds,xbounds)
        vis = oss.any_visible()
        
        self._order_visible = np.append(self._order_visible,vis)
        super(SpectrumSNRSelections,self).__setitem__(ordi,oss)

    def update_visible_list (self):
        for i in xrange(len(self)):
            self._order_visible[i] = self.__getitem__(i).any_visible()

    def set_window_visible (self,xbounds,truth=True):
        xmin,xmax = xbounds
        
        # which orders lie in the bounds
        xmins, xmaxs = self._order_bounds.T
        inmask = (xmin < xmaxs)*(xmins < xmax)
              
        changed = False
        # lie outside the borders, but has some which are visible
        if len(self._order_visible) > 1:   
            outmask = np.logical_not(inmask)
                    
            vis = (self._order_visible == True)        
            not_truth = (self._order_visible != truth)

            ran = np.arange(len(self._order_visible))

            ran_2_off = ran[outmask*vis]
            ran_2_change = ran[inmask]
            for i in ran_2_off:
                self.__getitem__(i).set_window_visible(xbounds,False)
                self._order_visible[i] = False    
                changed = True
            
            # now change the ones in the x bounds
            for i in ran_2_change:
                self.__getitem__(i).set_window_visible(xbounds,truth)
                self._order_visible[i] = self.__getitem__(i).any_visible()
                changed = True
        else:
            changed = self.__getitem__(0).set_window_visible(xbounds,truth)
            self._order_visible[0] = self.__getitem__(0).any_visible()
            
        return changed

    def get_current_state (self):
        all_ord_sel = []
        for ordi in self:
            all_ord_sel.append([ordi,self.__getitem__(ordi).get_current_state()])
            
        cur_state = {'state':'all snr regions',
                     'all':all_ord_sel,
                     '_order_bounds':self._order_bounds,
                     '_order_visible':self._order_visible}
        return cur_state
            
    def apply_prev_state (self,prev_state):
        assert prev_state['state'] == 'all snr regions'
        
        self._order_bounds = np.array([],dtype=float)
        self._order_visible = np.array([],dtype=bool)
        
        xbounds = self.ax.axis()[:2]
        xmin,xmax = xbounds
        
        for pstate in prev_state['all']:
            ordi = pstate[0]
            if ordi not in self:
                print "HeadsUp: couldn't restore order:",ordi
                continue
            self.__getitem__(pstate[0]).apply_prev_state(pstate[1])
            self.__getitem__(pstate[0]).set_window_visible(xbounds)
            
            dbounds = self.__getitem__(pstate[0]).get_xbounds()
            vis = self.__getitem__(pstate[0]).any_visible()
            
            self._order_bounds = np_vstack_append(self._order_bounds,dbounds)
            self._order_visible= np.append(self._order_visible,vis)
                        
class InteractiveVarianceIO (InputOutput):
    """
    This contains the Input and Output methods for the InteractiveVarianceEditor Class
    """
    def __init__ (self,parent):
        self.parent = parent
        self.p = parent

        # self.p.history.save_history(filename,clobber=True)
        # self.p.history.open_history(filename)
        T = time.gmtime()
        curdate = format(T.tm_mon,'02')+format(T.tm_mday,'02')+str(T.tm_year) 


        self.iochoices = SaveOpenChoices()
        self.iochoices.add('Signal-to-Noise Regions',
                      '(*_snr.txt)|*_snr.txt| All files (*)|*',
                      'REGIONS'+curdate+"_snr.txt", 
                      'SNR Regions',
                      self.save_snr_regions,
                      self.open_snr_regions)
        
        self.iochoices.add('Signal-to-Noise Regions - SMART',
                           '(*_continuum.dat)|*_continuum.dat',
                           'REGIONS'+curdate+"_continuum.dat",
                           'SNR Regions SMART',
                           self.save_snr_for_smart)

        self.iochoices.add('SAVE PROGRESS',
                           save_method = None,
                           open_method = None)

    def _current_params_var (self):
        lines = []
        lines.append("# These are editor parameters for the SNR Editor. If you change the order of these then eyeSpec won't be able to read back in")
        snri = self.p.region_highlight.get_snr_index()
        lines.append("_snr_i = "+str(snri))
        dat_lines = self.p._current_params_data()
        return lines+dat_lines

    def _apply_params_var (self,lines):
        # !! GAH, this I can't specialize input for SNR editor
        self._apply_params_data(lines)


    def save_current_params_var (self,filename):
        print "!! not currently working"
#        f = open(filename,'w')
#        lines = self._current_params_var()
#        f.write("\n".join(lines))
#        f.close()

    def save_snr_regions (self,filename):
        filename = str(filename)
        f = open(filename,'w')
        f.write("# SPECTRA FILE:"+self.p.spec_obj.filename+" \n")
        lines = self.p.get_snr_regions(True)
        f.write("\n".join(lines))
        f.close()

    def save_snr_for_smart (self,filename,clobber=True):
        """
        This saves the continuum regions in the same way as Bayard Stringer's code SMART likes it

        fileprefix =>   fileprefix+"_continuum.dat"
        """
        filename = str(filename)

        if filename[-14:] != '_continuum.dat':
            filename = filename+"_continuum.dat"

        if os.path.exists(filename) and not clobber: 
            print "Error: File exists: '"+filename+"'"
            return
        
        f = open(filename,'w')
        f.write("# FILE:"+self.p.spec_obj.filename+" divided by orders\n")

        output = [format('xmid','<10'),
                  format('mean','<10'),
                  format('stdev','<10')]
        
        f.write("# "+" ".join(output)+"\n")
        
        for ord in self.p._snr_selections.keys():
            f.write("# order "+str(ord)+" "+str(len(self.p._snr_selections[ord]))+"\n")
            for snr_sel in self.p._snr_selections[ord]:
                bb = snr_sel.get_bounds()
                stats = snr_sel.get_stats()
                output = [format((bb[0]+bb[1])/2.,'<10.3f'),
                          format(stats[0],'<10.3f'),
                          format(stats[2],'<10.3f')]
                f.write("   ".join(output))
                f.write("\n")
        

        f.close()
               
    def save_snr_progress_save (self,sprog_dir,mode,add_to=False):
        # Define the save progress directory based on previous or new
        if add_to:
            if not sprog_dir.__class__.__name__ == 'ProgressSave': raise TypeError("Input not of correct type ProgressSave")
            sprog=sprog_dir
        else: sprog = ProgressSave(sprog_dir,mode)
        if not add_to: sprog.set_current_save_num()

        # check the save progress number
        save_num = sprog.get_current_save_num()

        # if starting here then define the routine
        if save_num is None: raise StandardError("Whoops, save_num shouldn't be None")
        if not add_to:
            sprog['ROUTINE'] = "SNR Editor"
            sprog['SAVE_TYPE'] = 'manual'

        #----------------------------------------------------------------------#
        # DEFINE THE SPECIFIC PIECES TO GO INTO THE SAVE PROGRESS =>
        # save SNR Regions
        filename = 'SNR_REGIONS_'+str(save_num)+'_snr.txt'
        self.save_snr_regions(filename)
        sprog.add_file_to_dir(filename,save_num,'SNR_SEL',overwrite=True)
        # save original Spectrum
        filename = 'ORIG_SPEC_'+str(save_num)+'.pkl'
        self.p.dp.save_orig_obj(filename)
        sprog.add_file_to_dir(filename,save_num,'ORIG_SPEC',overwrite=True)
        # save edited spectrum
        filename = 'EDIT_SPEC_'+str(save_num)+'.pkl'
        self.p.dp.save_spec_obj(filename)
        sprog.add_file_to_dir(filename,save_num,'EDIT_SPEC',overwrite=True)
        # save history
        filename = 'SNR_HISTORY_'+str(save_num)+'_hist.pkl'
        self.p.history.save_history(filename)
        sprog.add_file_to_dir(filename,save_num,'HISTORY',overwrite=True)
        # save current editor parameters
        filename = 'SNR_PARAMS_'+str(save_num)+'_snr.par'
        self.save_current_parameters(filename)
        sprog.add_file_to_dir(filename,save_num,'SNR_PARAMS',overwrite=True)
        #----------------------------------------------------------------------#        

        # write record file if you're not adding to another
        if not add_to: sprog.write_record_file(overwrite=True)
  
    def open_current_params_var (self,filename):
        print "!! not currently working"
#        filename = str(filename)
#        if os.path.exists(filename): raise IOError("File does not exist:"+filename)
#        
#        f = open(filename,'r')
#        lines = f.readlines()
#        self._apply_params_var(lines)
#        f = open(filename)
 
    def open_snr_regions (self,path):
        if not os.path.exists(path): raise IOError("File does not exist:"+path)

        self.p.history.add('edit all regions',self.p._snr_selections.get_curent_state(),
                           info='open new file of snr regions :'+str(path))

        print "=== READING IN SNR REGIONS"
        try: snr_info = np.loadtxt(path)
        except: raise IOError("File must be a text file")
        
        for row in snr_info:
            if row.ndim != 1 : 
                print "Dimension of row is not correct"
                break
            if len(row) < 4:
                print "Length of row is not correct"
                break

            ordi = int(row[0])
            xmin,xmax = row[1],row[2]
            ymin,ymax = row[3],row[4]
            
            # if the order number is within the correct range
            if 0 <= ordi < self.p.spec_obj.shape[1]:
                snr_sel = self.p.convert_to_snr_region(ordi=ordi,bounds=[xmin,xmax,ymin,ymax])
                if snr_sel is not None:
                    self.p._snr_selections[ordi].append(snr_sel)

        for ordi in self.p._snr_selections:
            self.p._snr_selections[ordi].sort_by_wavelength()

    def open_snr_progress_save (self,sprog):
        save_num = sprog.get_current_save_num()

        sprog_dir = sprog.get_current_dir()
        # this gives what I need to restore progress. The Keyword and then the default file, if the file was given a different name then it will take that from the record.sprog file
        needed_objects = {'SNR_SEL'   : sprog_dir+'SNR_REGIONS_'+str(save_num)+'.txt',
                          'ORIG_SPEC' : sprog_dir+'ORIG_SPEC_'+str(save_num)+'.pkl',
                          'EDIT_SPEC' : sprog_dir+'EDIT_SPEC_'+str(save_num)+'.pkl',
                          'HISTORY'   : sprog_dir+'SNR_HISTORY_'+str(save_num)+'_hist.pkl',
                          'SNR_PARAMS': sprog_dir+'SNR_PARAMS_'+str(save_num)+'.par'}

        for extra_keyword in needed_objects.keys():
            if extra_keyword not in sprog.keys(): raise IOError("SNR Editor progress save expected a keyword not found: "+extra_keyword)
            else: 
                value = sprog[extra_keyword] 
                needed_objects[extra_keyword] = sprog_dir+value

        # does the file actually exist where you think it is?
        for file in needed_objects.values():
            if not os.path.exists(file): raise IOError("SNR Editor for save number: "+str(save_num)+" needs file:"+file)
        
        filename = needed_objects['SNR_SEL']
        self.open_all_snr_regions(filename)

        filename = needed_objects['ORIG_SPEC']
        self.p.dp.open_orig_obj(filename)

        filename = needed_objects['EDIT_SPEC']
        self.p.dp.open_spec_obj(filename)

        filename = needed_objects['HISTORY']
        self.p.history.open_history(filename)

        filename = needed_objects['SNR_PARAMS']
        self.open_current_parameters(filename)

class InteractiveVarianceEvents:
    """
    This contains the Event methods for the InteractiveVarianceEditor Class
    """
    def __init__ (self,parent):
        self.p = parent
        self._started_in_snr_region = False
        self._prev_region_state = None
        self._is_button_down = False
        
        
        self.ive_evt_timer = Timer(0.3)

        
# #=========== NEEDED PIECES =============================#
#

    def is_xy_in_snr_region (self,event):
        """
        Go through all the snr regions for the current order and see if the x,y event is in that region
        """
        ordi = self.p.get_order_index()

        if ordi is None: return -1
        # if you click check to see if you clicked on a selection box
        xy = (event.xdata,event.ydata)
        for i in xrange(len(self.p._snr_selections[ordi])):   
            snr_sel = self.p._snr_selections[ordi][i]
            if xy in snr_sel: return i
            
        return -1

    def _editing_snr_region (self):
        self._started_in_snr_region = True
        # !! this might slow it down some
        self._prev_region_state = self.p.hist.get_current_state_snr_regions()

    def _snr_region_edit (self,event):
        """
        Return True if need to update
        """
        ordi = self.p.get_order_index()
        if ordi is None: return False
        if len(self.p._snr_selections[ordi]) == 0: return False
        
        rh = self.p.region_highlight
        
        snri = rh.get_snr_index()
        if snri is None: return False

        snr_edit_region = rh.get_edit_region()
        if snr_edit_region is "none": return False

        ordi = ordi
        start_lims = rh.get_starting_bounds()
        start_pt = self.p.get_recorded_xy(-1)

        # set the fixed points
        xc1,xc2,yc1,yc2 = rh.get_starting_bounds()


        if snr_edit_region is "top_left":
            xc1 = start_lims[0] + (event.xdata - start_pt[0])
            yc2 = start_lims[3] + (event.ydata - start_pt[1])

        elif snr_edit_region is "top": 
            yc2 = start_lims[3] + (event.ydata - start_pt[1])
            
        elif snr_edit_region is "top_right": 
            xc2 = start_lims[1] + (event.xdata - start_pt[0])
            yc2 = start_lims[3] + (event.ydata - start_pt[1])
            
        elif snr_edit_region is "left":
            xc1 = start_lims[0] + (event.xdata - start_pt[0])
            
        elif snr_edit_region is "center":
            xc1 = start_lims[0] + (event.xdata - start_pt[0])
            xc2 = start_lims[1] + (event.xdata - start_pt[0])
            yc1 = start_lims[2] + (event.ydata - start_pt[1])
            yc2 = start_lims[3] + (event.ydata - start_pt[1])
            
        elif snr_edit_region is "right": 
            xc2 = start_lims[1] + (event.xdata - start_pt[0])
            
        elif snr_edit_region is "bot_left": 
            xc1 = start_lims[0] + (event.xdata - start_pt[0])
            yc1 = start_lims[2] + (event.ydata - start_pt[1])
            
        elif snr_edit_region is "bot":
            yc1 = start_lims[2] + (event.ydata - start_pt[1])
            
        elif snr_edit_region is "bot_right":
            xc2 = start_lims[1] + (event.xdata - start_pt[0])
            yc1 = start_lims[2] + (event.ydata - start_pt[1])
 
        xn1 = min((xc1,xc2))
        xn2 = max((xc1,xc2))
        yn1 = min((yc1,yc2))
        yn2 = max((yc1,yc2))
        
        # update boxes
        self.p._snr_selections[ordi].update_snr_box_bounds(snri,[xn1,xn2,yn1,yn2])
        self.p._snr_selections[ordi][snri].update_mean_line_x(xn1,xn2)
        self.p.region_highlight.set_selection_box([xn1,xn2,yn1,yn2])
        
        # update stats of the box
        self._snr_region_edit_update_text_mean(ordi,snri,[xn1,xn2,yn1,yn2])
        return True

    def _snr_region_edit_update_text_mean (self,ordi,snri,bounds,force_update=False):
        if not force_update:
            if not self.ive_evt_timer.check(): return False
                    
        spec_obj = self.p.dataplot.spec_obj
        order = [spec_obj.get_wl(ordi),
                 spec_obj.get_data(ordi),
                 spec_obj.get_inv_var(ordi)]
        
        snr_sel = self.p._snr_selections[ordi][snri]
        snr_sel.calc_stats(order,bounds,set_stats=True)
        snr_sel.update_text(bounds)
        if force_update: snr_sel.update_mean_line(bounds)
        else: snr_sel.update_mean_line()
        return True

# #=========== Event call-backs =============================#
#

    def btn_press_snr_regions (self,event):
        """
        For used during button_press_callback
        
        returns changed,kind
        changed is True if need to be updated
        kind is 'deselect' if the region was deselected and 'select' if selected, and None if other
        """
        
        ordi = self.p.get_order_index()
        if ordi is None: return False,None
        self._started_in_snr_region = False
        order_regions = self.p._snr_selections[ordi]
        if len(order_regions) == 0: return False,None
        self._is_button_down = True

        # !! could stick a bunch of the 'selected'/'deseleted stuff into 
        # self._snr_selections[ordi].is_selected()

        changed = (False,None)
        self._prev_region_state = None
        # is an SNR region selected?
        rh = self.p.region_highlight
        if rh.get_visible():
            snri = rh.get_snr_index()
            snr_bounds = order_regions[snri].get_bounds()  
            xy = (event.xdata,event.ydata)
              
            if xy in rh:
                self.p.region_highlight.set_edit_region(xy,snr_bounds)                
                self._editing_snr_region()
                return (False,'select')
            else:
                self.p.region_highlight.deselect()
                self._started_in_snr_region = False
                changed = (True,'deselect')
    
        # if you didn't click in the current active region or 
        # if no region is active
        
        # check the xy if it's in an snr_region
        index = self.is_xy_in_snr_region(event)
        
        # didn't click in any region
        if index == -1: 
            self._started_in_snr_region = False
            return changed
        
        snr_bounds = order_regions[index].get_bounds()
        # found a snr selection
        self._editing_snr_region()
        self.p.region_highlight.set_highlight_visible(index,snr_bounds)
        self.p.region_highlight.set_edit_region((0,0),snr_bounds,'center')
        return (True,'select')
    
    def mot_notify_snr_regions (self,event):
        if event.inaxes is None: return False
        changed = False
        if self.p.region_highlight.get_visible() and self._is_button_down:
            changed = self._snr_region_edit(event)
        return changed 
      
    def btn_click_snr_regions (self,event):
        """
        For modifying a snr region
        given during a button_release_callback
        """
        self._is_button_down = False
        if event.inaxes is None: return False
        
        if self.p._ordi is None:
            if self.p.region_highlight.get_visible():
                self.p.region_highlight.deselect()
                return True
            return False            
        
        # if you right click
        if event.button == 3 and not self._started_in_snr_region:
            changed = self.p._check_and_convert_to_snrsel()
            return changed     
            
    def btn_release_snr_regions (self,event):
        if not event.inaxes: return False
        ordi = self.p.get_order_index()
        if ordi is None: return False
        if event.button != 1: return False
        self._is_button_down = False
        

        # combinations of start/end-in/out of SNR region and dragging
        # determined you must be dragging for the following     
        if not self._started_in_snr_region: 
            return False
       
       
        # if no highlight box is visible then skip this
        rh = self.p.region_highlight
        if not rh.get_visible(): return False
        snri = rh.get_snr_index()
        
        # update the plot objects
        bounds = self.p._snr_selections[ordi][snri].get_bounds()
        spec_obj = self.p.dataplot.spec_obj
        order = [spec_obj.get_wl(ordi),
                 spec_obj.get_data(ordi),
                 spec_obj.get_inv_var(ordi)]
        
        stats = self.p._snr_selections[ordi][snri].calc_stats(order,data_bounds=bounds)
        found_merger =  self.p._snr_selections[ordi].check_for_overlap(bounds,exclude=[snri])
        if self._prev_region_state is not None:
                self.p.history.add('edit region',self._prev_region_state,
                                   info='moved a region:'+str((ordi,snri)))
        
        if stats[0] == 0 or found_merger: 
            # no points remained, so delete the box
            # or found an overlap with the new bounds
            # then delete
            self.p.delete_snr_region()
            return True
        
        # if the box wasn't deleted
        self.p._snr_selections[ordi].update_snr_box_bounds(snri,bounds)
        self._snr_region_edit_update_text_mean(ordi,snri,bounds,force_update=True)
        self.p.region_highlight.set_selection_box(bounds)
        return True

class InteractiveVarianceHistory:
    """
    This contains the History undo/redo methods for the InteractiveVarianceEditor Class
    """
    def __init__ (self,parent):
        self.p = parent
        
# #=============  ALL SNR REGIONS ===============================#
#

    def undo_redo_all_snr_regions (self,prev_state):
        cur_state = self.p._snr_selections.get_current_state()
        self.p._snr_selections.apply_prev_state(prev_state)
        return cur_state

# #=============  ONE SNR REGIONS ===============================#
#

    def get_current_state_snr_regions (self):
        ordi = self.p.get_order_index()
        if ordi is None: return None
        cur_state = self.p._snr_selections[ordi].get_current_state()
        if cur_state is not None: 
            cur_state['snri'] = self.p.region_highlight.get_snr_index()
            cur_state['ordi'] = ordi
        return cur_state

    def apply_prev_state_snr_regions (self,prev_state):
        if prev_state is None: return None
        self.p.region_highlight.deselect()
        self.p.deselect_all()
        
        ordi = prev_state['ordi']
        snri = prev_state['snri']
        if ordi not in self.p._snr_selections: raise ValueError("Whoops, you'll need to add some code here to deal with that")
        self.p._snr_selections[ordi].apply_prev_state(prev_state)
        self.p.select_order_by_index(ordi)
        if snri is not None:
            snr_bounds = self.p._snr_selections[ordi][snri].get_bounds()
            self.p.region_highlight.set_highlight_visible(snri, snr_bounds)
                
    def undo_redo_edit_region (self,prev_state):
        cur_state = self.get_current_state_snr_regions()
        self.apply_prev_state_snr_regions(prev_state)
        return cur_state

class InteractiveVarianceEditor (object, InteractiveDataEditor):
    """
    This allows for setting regions of continuum determines the signal-to-noise ratio (SNR)
    and does a calculation to determine the inverse variance for the data based on 
    these regions    
    """

    def __init__ (self,ax,spec_obj,snr_regions=None):
        InteractiveDataEditor.__init__(self,ax,spec_obj,auto_scale_focus='full')

        # set the starting auto scaling to be scale y axis only
        self.dp.set_auto_scaling(3)

        # !! need to be able to load snr_regions
        self._snr_selections = SpectrumSNRSelections(ax)
        for i in xrange(self.spec_obj.shape[1]):
            xmin,ymin = self.spec_obj.get_min(i)
            xmax,ymax = self.spec_obj.get_max(i)
            self._snr_selections[i] = OrderSNRSelections(self.ax,i,[xmin,xmax])

        self.hist = InteractiveVarianceHistory(self)
        self.evts = InteractiveVarianceEvents(self)
        self.io = InteractiveVarianceIO(self)

        self.region_highlight = HighlightSNRSelection(self.ax)
        
        self.history.set_options('edit region',self.hist.undo_redo_edit_region)
        self.history.set_options('edit all regions',self.hist.undo_redo_all_snr_regions)

    def __repr__ (selfs):
        return 'InteractiveVarianceEditor'

    def guess_continuum_regions (self,order_list=None):
        """
        Create a guess of continuum regions
        """
        if order_list is not None: order_list = list(order_list)
        else: order_list = xrange(self.spec_obj.shape[1])
    
        #dlg = RegionGuesserDialog(guess_snr,None,-1,'Add continuum to plot:')
        #dlg.ShowModal()
       
        dlg = wx.TextEntryDialog(None, 'Enter SNR guess:',defaultValue = '100')
        if dlg.ShowModal() == wx.ID_OK:
            try: guess_snr = float(dlg.Value)
            except: 
                print "Could not convert given string to floating point value:"+dlg.Value
                return
        else: guess_snr = 0.8
   
        if guess_snr == 0: guess_snr = 1e-10
        dlg.Destroy()

        # !! give min number of points in a region
        min_number_pts = 10.0
        
        # !! how many sigma wide for features
        sigma_width = 4.0
        
        # !! how many sigma below
        sigma_cutoff = 3.5
        sig = 1.0/guess_snr
        
        
        # !! density of regions,  number_of_regions/angstrom
        density_regions = None
        
        
        cutoff_depth = 1.0 - sigma_cutoff*sig
        if cutoff_depth < 0.1:
            print "Cutoff too conservative 0.1 > :"+str(cutoff_depth)
            return
        
        self.history.add('edit all regions',self._snr_selections.get_current_state(),
                         info='gave a guess at all the SNR regions :'+str(guess_snr))
        
        for ordi in order_list:
            self._order_add_snr(ordi,cutoff_depth,sigma_width,min_number_pts)
        
        self._snr_selections.set_window_visible(self.ax.axis()[:2])
        
            
    def _order_add_snr (self,ordi,cutoff_depth,sigma_width,min_number_pts):
        print "-"*40
        print "GUESS CONTINUUM FOR ORDER:",ordi
        try: flux = self.spec_obj.get_data(ordi)
        except:
            print "Error in getting data, perhaps the given order number does not exist:"+str(ordi)
            return
        
        wl = self.spec_obj.get_wl(ordi)
        order = np.vstack((wl,flux))

        # !! edit this function with Brian to take a SNR which corresponds to the level which the cut is subsequently done
        pairs = continuum_finder(order,cutoff_depth,sigma_width)
             
        if len(pairs) == 0: print "None"
        for pair in pairs:
            xmin,xmax = pair
            mask = (xmin < wl)*(wl < xmax)
            if not np.any(mask): return

            ymin = np.min(flux[mask])
            ymax = np.max(flux[mask])
            mwl = wl[mask]
            if len(mwl) <= min_number_pts: return
            
            snr_sel = self.convert_to_snr_region(ordi=ordi,bounds=[xmin,xmax,ymin,ymax])
            if snr_sel is not None:
                self._snr_selections[ordi].append(snr_sel)
                stat = snr_sel.stats
                
                xbounds = snr_sel.get_bounds()[:2]
                outbounds = [format(xbounds[0],'4.3f'),format(xbounds[1],'4.3f')]
                outputval = [format('('+outbounds[0]+','+outbounds[1]+')','^20'),
                             ' : ',
                             format(snr_sel.snr,'<10.3f'),
                             format(stat[0],'<10.3f'),
                             format(stat[2],'<10.4f'),
                             format(stat[3],'<10.4f')]
                
                print "RECORDED == "+" ".join(outputval)        
        
    def get_snr_regions (self,as_string=False):
        lines = []
        first_line = fl = ['order','xmin','xmax','ymin','ymax','snr','mean','stdev','N']
        if as_string: 
            output = [format(fl[0],'<6'),
                      format(fl[1],'<10'),
                      format(fl[2],'<10'),
                      format(fl[3],'<10'),
                      format(fl[4],'<10'),
                      format(fl[5],'<10'),
                      format(fl[6],'<10'),
                      format(fl[7],'<10'),
                      format(fl[8],'<10')]
            first_line = "#"+"  ".join(output)
            
        lines.append(first_line)
        for ordi in self._snr_selections.keys():
            for snr_sel in self._snr_selections[ordi]:
                bb = snr_sel.get_bounds()
                stats = snr_sel.get_stats()
                snr = snr_sel.snr
                if as_string:
                    output = ["  "+format(ordi,'<6'),
                              format(bb[0],'<10.3f'),
                              format(bb[1],'<10.3f'),
                              format(bb[2],'<10.3f'),
                              format(bb[3],'<10.3f'),
                              format(snr_sel.snr,'<10.3f'),
                              format(stats[0],'<10.3f'),
                              format(stats[2],'<10.3f'),
                              format(int(stats[3]),'<10')]
                    output = '  '.join(output)
                else: output = [ordi]+list(bb)+[snr,stats[0],stats[2],stats[3]]
                lines.append(output)  
                
        return lines

#################################################

    def delete_snr_region (self,ordi=None,snri=None):
        """
        Delete an SNR region based on the order and SNR index
        """
        if ordi is None:
            if self._ordi is None: return
            else: ordi = self._ordi
        else: raise StandardError("Whoops, ordi not currently implemented")

        used_active_region = False
        if snri is None:
            if self.region_highlight.get_visible():
                used_active_region = True
                snri = self.region_highlight.get_snr_index()
            else: return
        else: raise StandardError("Whoops, snri not currently implemented")
        
        check_no_other_box = (not self.selection_box.get_visible())
        check_ordi = (ordi is not None)

        if (ordi in self._snr_selections.keys()):
            if len(self._snr_selections[ordi]) == 0: return
        else: return

        if check_no_other_box and check_ordi:
            if snri is None: snri = len(self._snr_selections[self._ordi])-1
                
            snr_sel = self._snr_selections[self._ordi][snri]
            stat = snr_sel.stats
            xbounds = snr_sel.get_bounds()[:2]
            outbounds = [format(xbounds[0],'4.3f'),format(xbounds[1],'4.3f')]
            outputval = [format('('+outbounds[0]+','+outbounds[1]+')','^20'),
                                ' : ',
                                format(snr_sel.snr,'<10.3f'),
                                format(stat[0],'<10.3f'),
                                format(stat[2],'<10.4f'),
                                format(stat[3],'<10.4f')]
            print "DELETED == "+" ".join(outputval)

            self.history.add('edit region',self.hist.get_current_state_snr_regions(),
                             info='deleted a region in bounds'+str(tuple(xbounds)))

            snr_sel = self._snr_selections[self._ordi][snri]
            snr_sel.set_visible(False)
            if used_active_region: self.region_highlight.deselect()
            
            del self._snr_selections[self._ordi][snri]
            
            snri -= 1
            if snri <= 0: snri = 0
            if len(self._snr_selections[self._ordi]) == 0: return
            
            snr_bounds = self._snr_selections[self._ordi][snri].get_bounds()
            self.region_highlight.set_highlight_visible(snri,snr_bounds)
        
    def scan_through_snr_selections (self,direction):
        if self._ordi is None: return
        if self._ordi not in self._snr_selections: return
        if len(self._snr_selections[self._ordi]) == 0: return
        
        snri = self.region_highlight.get_snr_index()
        if snri is None: return
        
        if direction == '+': inc = +1
        elif direction == '-': inc = -1
        else: inc = 0
        
        snri += inc
        if snri < 0: snri = 0
        N = len(self._snr_selections[self._ordi])
        if snri >= N: snri = N-1
        
        snr_bounds = self._snr_selections[self._ordi][snri].get_bounds()
        self.region_highlight.set_highlight_visible(snri,snr_bounds)
        self.print_snr_selections(self._ordi)
        
    def print_snr_selections (self,ordi):
        if ordi is None:return
        if ordi not in self._snr_selections.keys(): return

        order_selections = self._snr_selections[ordi]

        print "==="*10+format('ORDER:'+str(ordi),'^10')+"==="*10    
        outputval = [format('(xmin,xmax)','^20'),
                     ' : ',
                     format('SNR','<10'),
                     format('mean','<10'),
                     format('stdev','<10'),
                     format('N','<10')]
                     
        
        print "  "+" ".join(outputval)
        if len(order_selections) == 0:
            print " NO VALUES RECORDED. USE ALT KEY (UNSHIFTED) TO RECORD "
            print "==="*20
            return
        
        for i in xrange(len(order_selections)):
            snr_sel = order_selections[i]
        
            stat = snr_sel.stats

            xbounds = snr_sel.get_bounds()[:2]
            outbounds = [format(xbounds[0],'4.3f'),format(xbounds[1],'4.3f')]
            outputval = [format('('+outbounds[0]+','+outbounds[1]+')','^20'),
                         ' : ',
                         format(snr_sel.snr,'<10.3f'),
                         format(stat[0],'<10.3f'),
                         format(stat[2],'<10.4f'),
                         format(stat[3],'<10.4f')]
                     

            print_out = str(i)+"  "+" ".join(outputval)
            snri = self.region_highlight.get_snr_index()
            if i == snri and snri is not None: print_out = ">> "+print_out+" << selected"
            print print_out
        print "==="*20

    def _check_for_snr_selection_overlap (self,ordi,snr_sel_new,exclude=[]):
        # !! instead you could make this merge overlapping SNR regions into one
        if ordi is None: return snr_sel_new
        if ordi not in self._snr_selections.keys(): return snr_sel_new
        order_selections = self._snr_selections[ordi]
        new_bounds = snr_sel_new.get_bounds()
        xb1,xb2,yb1,yb2 = new_bounds 
        found_merger = order_selections.check_for_overlap(new_bounds)
        if found_merger: "HeadsUp: You can't overlap two selection regions"
        return found_merger

    def select_snr_region_by_index (self,snri,always_find=False):
        if self._ordi is None: return
        snri = int(snri) # must be an integer
        snr_selections = self._snr_selections[self._ordi]

        # check to make sure the snri is valid
        out_of_bounds = False
        if snri >= len(snr_selections): out_of_bounds = 'end'
        elif snri < 0: out_of_bounds = 'begin'

        if out_of_bounds: 
            if always_find:
                if out_of_bounds == 'end': snri = len(snr_selections) - 1
                if out_of_bounds == 'begin': snri = 0
            else: 
                print "HeadsUp: You gave invalid selection for SNR Region in order:", self._ordi
                return

        snr_bounds = self._snr_selections[self._ordi][snri].get_bounds()
        # apply that to current 
        self.region_highlight.set_highlight_visible(snri,snr_bounds)

    def _check_and_convert_to_snrsel (self,add_hist=True):
        if self._ordi is None: return False
        if not self.selection_box.get_visible(): return False
        
        # turn off selection box
        self.selection_box.set_visible(False)
        
        # convert region
        snr_sel = self.convert_to_snr_region()
        if snr_sel is None: return False
        
        if add_hist: self.history.add('edit region',self.hist.get_current_state_snr_regions(),
                                      info='added a region to order'+str(self._ordi))
        
        # add snr region to the list
        self._snr_selections[self._ordi].append(snr_sel)
        
        # get from snr selection
        stat = snr_sel.get_stats()
        bounds = snr_sel.get_bounds()
                        
        # select the newest region        
        self.region_highlight.deselect()
        snri = len(self._snr_selections[self._ordi])-1
        snri = self._snr_selections[self._ordi].sort_by_wavelength(snri)
        self.region_highlight.set_highlight_visible(snri,bounds)                
        
        # print out what was recored
        xbounds = bounds[:2]
        outbounds = [format(xbounds[0],'4.3f'),format(xbounds[1],'4.3f')]
        outputval = [format('('+outbounds[0]+','+outbounds[1]+')','^20'),
                     ' : ',
                     format(snr_sel.snr,'<10.3f'),
                     format(stat[0],'<10.3f'),
                     format(stat[2],'<10.4f'),
                     format(stat[3],'<10.4f')]
        print "RECORDED == "+" ".join(outputval)
        
        
        self._snr_selections.update_visible_list()
        return True
        
    def convert_to_snr_region (self,ordi = None,bounds=None):        
        if ordi is None: 
            if self._ordi is None: raise ValueError("whoops")
            i = self._ordi
        else: i = int(ordi)

        if bounds is None: bbounds = self.selection_box.get_bounds()
        else: bbounds = bounds

        # check bounds for overlap with previous
        if self._snr_selections[i].check_for_overlap(bbounds): return None

        # get the data
        wl = self.spec_obj.get_wl(i)
        data = self.spec_obj.get_data(i)
        inv_var = self.spec_obj.get_inv_var(i)      
        order = [wl,data,inv_var]  
        
        # get the color
        fcolor = self.dataplot.plot_data[i].get_color()
        
        # create signal to noise (snr) selection
        snr_sel = SNRSelection(self.ax,[order,bbounds,i],color=fcolor)
        return snr_sel

    def calc_order_inv_var (self,ordi,function='linear',get_checks=False):
        if ordi is None: return None
        if ordi not in self._snr_selections.keys(): return None
      
        # get_checks = False # !! change this to give checks if wanted
      
        obj = self.spec_obj

        # use the entire wl array
        obj.set_use_cropped(False)
        
        # get the current data
        wl = obj.get_wl(ordi)
        data = obj.get_data(ordi)
        inv_var = obj.get_inv_var(ordi)

        order = [wl,data,inv_var]

        # mask out points which have inv_var == 0
        mask_for_deleted_points = (inv_var == 0.0)
        mask_not_deleted_points = np.logical_not(mask_for_deleted_points)

        mask_for_set_points = (np.ones(len(inv_var)) < 0)
        
        # create new array which will eventually be the new inv_var
        new_fix_sigmas = np.ones(inv_var.shape)*99.9
        new_fix_means = np.ones(inv_var.shape)

        kount = 0
        order_selections = self._snr_selections[ordi]
        # go through all snr regions in the particular order
        for i in xrange(len(order_selections)):
            snr_sel = order_selections[i]
            # get the wavelength points with the snr region
            mask_wl = (wl >= snr_sel.get_bounds()[0])*(wl <= snr_sel.get_bounds()[1])
            
            # create mask of points within the bounds and != 0
            mask_set = mask_wl*mask_not_deleted_points
            mask_for_set_points[mask_set] = True

            # check to see if any points are in the mask
            if np.any(mask_set): kount += 1
            else: continue
                        
            # use this mask to set all those points equal to derived error
            mmean,sigma = snr_sel.stats[0],snr_sel.stats[2]
            
            # >>>>>> FILL IN POINTS WHICH HAVE SIGMA SET <<<<< #
            new_fix_sigmas[mask_set] = sigma
            new_fix_means[mask_set] = mmean

        # now get sigmas for all the points inbetween these regions
        mask_no = (new_fix_sigmas == 1.0) # all points which didn't get sigmas
        mask_yes = np.logical_not(mask_no)

        # >>>>>> FILL IN POINTS BETWEEN SIGMA REGIONS USING SET <<<<< #
        # use the points which were set using sigmas to interpolate the points
        # which are not deleted but were not set before
        mfsp = mask_for_set_points
        mask_not_set_not_deleted = mnsnd = mask_not_deleted_points*np.logical_not(mfsp)

        # sigma => variance
        new_fix_var = new_fix_sigmas**2
   
        if np.any(mnsnd):
            if np.any(mfsp):
                new_fix_var[mnsnd] = scipy.interp(wl[mnsnd],wl[mfsp],new_fix_var[mfsp])
                new_fix_means[mnsnd] = scipy.interp(wl[mnsnd],wl[mfsp],new_fix_means[mfsp])
            # put in the default value if no points were set
            else: 
                new_fix_var[mnsnd] = 1.0
                new_fix_means[mnsnd] = data[mnsnd]
        if get_checks:      
            fix_plot_means = new_fix_means
            fix_plot_var_pls = new_fix_means+new_fix_var
            fix_plot_var_min = new_fix_means-new_fix_var
        
        # based on possion noise, scale the variance based on an offset 
        # with the mean assuming fixed SNR
        new_var = (data/new_fix_means)*new_fix_var

        # new_var += 0.0001 # add variance from the read noise

        # Now all the points which are not deleted (inv_var == 0) have been
        new_inv_var = inv_var_2_var(new_var)
        
        # delete points in the new inv_var
        new_inv_var[mask_for_deleted_points] = 0.0
        obj.set_use_cropped('previous')

        print "Calculated inverse variance based on",kount,"regions for order:",ordi
        
        if get_checks:
            return new_inv_var,[wl,fix_plot_means,fix_plot_var_pls,fix_plot_var_min,new_var,data]
        else:
            return new_inv_var

    def apply_new_inv_var (self):
        spec_obj = self.dataplot.spec_obj
        
        band = spec_obj.get_band()
        global stuff
        comp_wl = []
        comp_data = []
        fix_m = []
        fix_vp = []
        fix_vm = []
        new_v = []

        for i in xrange(spec_obj.shape[1]):
            new_inv_var = self.calc_order_inv_var(i,get_checks=False)
#            comp_wl.append(checks[0])
#            fix_m.append(checks[1])
#            fix_vp.append(checks[2])
#            fix_vm.append(checks[3])
#            new_v.append(checks[4])
#            comp_data.append(checks[5])
            
            spec_obj._inv_var[band][i] = new_inv_var
        self.dataplot.spec_obj = spec_obj

#        p_wl = np.concatenate(comp_wl)
#        p_d = np.concatenate(comp_data)
#        p_fm = np.concatenate(fix_m)
#        p_vp= np.concatenate(fix_vp)
#        p_vm = np.concatenate(fix_vm)
#        p_nv = np.concatenate(new_v)
#
#        stuff = np.array([p_wl,p_d,p_fm,p_vp,p_vm,p_nv])
#        np.savetxt("TESTING_STUFF.txt",stuff
#         ax = plt.figure(10).add_subplot(111)
#         s = stuff
#         ax.plot(s[0],s[1],color='b')
#         ax.plot(s[0],s[2],color='r',lw=2,alpha=.4)
#         ax.plot(s[0],np.sqrt(s[3]),color='y',lw=2,alpha=.4,linestyle='--')
#         ax.plot(s[0],np.sqrt(s[4]),color='y',lw=2,alpha=.4,linestyle='--')
#         ax.plot(s[0],10*s[5]+1,color='c')
        
    def get_final_obj (self):
        self.apply_new_inv_var()
        return self.dataplot.spec_obj
                
    def update (self):
        c1 = False
        if self.dataplot.bounds_changed():
            c1 = self._snr_selections.set_window_visible(self.ax.axis()[:2])
            
        c2 = super(InteractiveVarianceEditor,self).update()
        
        return (c1 or c2)
             
class EditVarianceManager (eyeSpecBaseEventManager):
    def __init__ (self, edit_var_panel, spec_obj, linelist = None):
        eyeSpecBaseEventManager.__init__(self)
        self.ppanel = edit_var_panel
        self.ax = edit_var_panel.ax
                
        # add editors
        self.ive = InteractiveVarianceEditor(self.ax,spec_obj)

        if linelist is not None:
            info_lines = None
            if len(linelist) > 1: info_lines = linelist[1]
            # !! the updating of the text on the plots takes a very long time. Currently this ability is supressed
            self.ile = PlotLines(self.ax,linelist[0])#,info_lines) 
        else:
            self.ile = None

        self.init_connection_callbacks(self) 
        self.key_cfg = KeyboardConfiguration()       
        self.key_cfg.add_key('d','Does same as backspace')
        self.key_cfg.add_key('h','Display this screen')
        self.key_cfg.add_key('g','Toggle grid on/off')
        self.key_cfg.add_key('p',"Toggle pan/zoom tool (note: editor won't work while engaged")
        self.key_cfg.add_key('z',"Toggle zoom rect tool (note: editor won't work while engaged")
        self.key_cfg.add_key('q','Close and return')
        self.key_cfg.add_key('`','Toggle auto scaling options')
        self.key_cfg.add_key('[','Undo data edit')
        self.key_cfg.add_key(']','Redo data edit')
        self.key_cfg.add_key('< and >','Scan through orders')
        self.key_cfg.add_key('backspace','Delete selected data or Signal-to-Noise region')
        self.key_cfg.add_key('esc',"same as 'q'")
        self.key_cfg.add_key(";",'Toggle data between scatter and line plot options')        
        #!! add more
        
        self.key_cfg.set_display_order(['backspace','d','esc','g','h','p','q','z','`',';','[',']',
                                        '< and >'])
        
        self.key_cfg.check_display()
        self.key_cfg.add_mpl_key_convert('<','< and >')
        self.key_cfg.add_mpl_key_convert('>','< and >')

        T = time.gmtime()
        curdate = format(T.tm_mon,'02')+format(T.tm_mday,'02')+str(T.tm_year) 

        self.iochoices = SaveOpenChoices()
        self.iochoices += self.ive.io.iochoices # save regions, ProgressSave
        self.iochoices += self.ive.dp.iochoices # save input/edited data
        
        self.iochoices.set_choices('save',['Signal-to-Noise Regions',
                                           'Signal-to-Noise Regions - SMART',
                                           'Input Spectrum',
                                           'Edited Spectrum',
                                           'SAVE PROGRESS'])

        self.iochoices.set_choices('open',['Signal-to-Noise Regions',
                                           'Input Spectrum',
                                           'Edited Spectrum',
                                           'SAVE PROGRESS'])
        self.iochoices.check_choices()

#        # progress save is formatted differently: 
#        self._save_choices['PROGRESS SAVE'] = None #self._save_snr_progress_save #==> snr_progress__save_function(sprog)        
#
#        self._open_choices['PROGRESS SAVE'] = 'SNR Editor'

    def delete_item (self):
        if self.ive.region_highlight.get_visible():
            self.ive.delete_snr_region()
        elif self.ile is not None:
            print "!! have a line deleter"
        else:
            self.ive.delete_data_selection()    
       
    def key_press_callback (self, event):        
        keypress = self.ive.key_held_callback(event, 'press')
        if keypress: return           
    
    def key_release_callback (self, event):
        self._pressing = False
        self._key_call(event)
 
    def _key_call (self, event):       
        if event.key == 'backspace': event.key = 'd'
        if event.key == 'esc': event.key = 'q'

        keypress = self.ive.key_held_callback(event, 'release')
        if keypress: return
        
        if event.key == '/':
            i = self.ive.get_order_index()
            if i is None: return
            osnr = self.ive._snr_selections[i]
            vis = osnr._list_visible
            print "======= for order",i
            for i in xrange(len(vis)):
                print "====",i,vis[i],osnr[i].get_visible(),osnr[i].get_bounds()[:2]
                            
        if event.key == 'd':
            self.delete_item()
            
        elif event.key == 'h':
            self.display_help('Variance Editor')
            
        elif event.key == 'g':
            self.ive.dp.toggle_grid()
            
        elif event.key == 'p':
            self.ive.dp.toggle_pan_tool()
            
        elif event.key == 'q':
            if 'Close' in dir(self.ppanel.pframe): self.ppanel.pframe.Close()
            
        elif event.key == 'z':
            self.ive.dp.toggle_zoom_tool()  
                  
        elif event.key == '`':
            self.ive.dp.set_auto_scaling() 
            
        elif event.key == ';':
            self.ive.dp.plot_data.toggle_data_line_scatter()
            
        elif event.key == '[':
            self.ive.history.undo()
            
        elif event.key == ']':
            self.ive.history.redo()
            
        elif event.key == ',': 
            self.ive.scan_through_orders("-")
            
        elif event.key == '.': 
            self.ive.scan_through_orders("+")
            
        elif event.key == 'c':
            self.ive._check_and_convert_to_snrsel()
        
        elif event.key == 'm':
            self.ive.guess_continuum_regions()
        
        elif event.key == ' ':
            self.ive.dp.scan_through_walking('+')
        
        elif event.key == 'b':
            self.ive.dp.scan_through_walking('-')
        
        elif event.key == 's':
            self.Output_To_File(event,self.iochoices)
        
        elif event.key == 'o':
            sprog = self.Input_From_File(event,self.iochoices)
            if sprog is not None: print "!! pass to self.ive.io.open_save_progress"
        
        elif event.key == 'down':
            self.ive.scan_through_snr_selections('+')
        
        elif event.key == 'up':
            self.ive.scan_through_snr_selections("-")
            
        self.update()
        
    def key_held_callback (self, event):
        pp = self._pressing
        crash = 0
        while pp:
            time.sleep(.5)
            tc = ((time.time() - self._pressing_time) > self._pressing_dtime)   
            if tc:
                print "do command again"
                self._key_call(event)
            crash += 1
            if crash > 10: break   

    def button_press_callback (self, event):
        """ Data Edit Manager  """
        if event.inaxes is None: return
        if self.ive.dp.is_toolbar_button_on():
            self.ive.dp.skip_auto_scale()
            self.update()
            return
        
        self._dragging = False
        self._clicked = True
        
        if event.button == 1: 
            self.ive.add_to_recorded_pts(event.xdata, event.ydata)

        
        c1 = self.ive.btn_press_selection_box(event)
        c2 = False
        
        kind = None            
        # if you're not clicking in a prexisting selection box
        if (event.xdata,event.ydata) not in self.ive.selection_box:
            c2,kind = self.ive.evts.btn_press_snr_regions(event)  
             
        # if the snr box changed and wasn't just deselected
        # then it is selected and the selection box should be 
        if c2 and kind != 'deselect':        
            if self.ive.selection_box.get_visible(): self.ive.deselect_selection()
            else: self.ive.selection_box.clear_box_variables() 
            
        
        if (c1 or c2): self.update()
        
    def motion_notify_callback (self, event):
        """ Data Edit Manager  """
        if event.inaxes is None: return
        if self.ive.dp.is_toolbar_button_on():
            self.ive.dp.skip_auto_scale()
            self.update()
            return
        if event.button != 1: return
        
        self._dragging = True
        self._clicked = False
        
        changed = False
        if self.ive.region_highlight.get_visible():
            changed = self.ive.evts.mot_notify_snr_regions(event)
            
        if not changed:
            changed = self.ive.mot_notify_selection_box(event)
        
        if changed: self.update()
                  
    def button_release_callback (self, event):
        """ Data Edit Manager  """
        if self.ive.dp.is_toolbar_button_on():
            self.ive.dp.skip_auto_scale()
            self.update()
            return
        
        # did you just click 
        if self._clicked:      
            c1 = self.ive.evts.btn_click_snr_regions(event)
            if not c1: c1 = self.ive.btn_click_data_editor(event)
        # or were you dragging           
        else: 
            c1 = self.ive.evts.btn_release_snr_regions(event)
        
        self._dragging = False  
        self._clicked = False      
        if c1: self.update()
        
    def update (self):
        self.ive.update()
        if self.ile is not None:
            self.ile.update()
        self.ax.figure.canvas.draw()
    
class SNR_MasterControl (EventConnections):#LinelistContinuumEditor (SNREditor,LineEditor):
    """
    This is the master class for SNREditor and LineEditor and decides which thing to use
    
    Controls:
    History
    key_press_callback
    button_press_callback
    motion_notify_callback
    button_release_callback

    """
    def __init__ (self,parent,line_editor,snr_editor): 
        self.parent = parent
        self.ax = self.parent.ax

        self._snr_editor = snr_editor
        self._line_editor = line_editor

        self._highlight_smart_line = self.ax.plot([0,0],[1,1],lw=10,color='r',alpha=.4,zorder=4,visible=False)

        self._use_which = None # 'snr_editor' 'line_editor'

        self._snr_editor.set_auto_scaling(1)        
        self._line_editor.auto_scale_x = False

        self._snr_editor.disconnect()
        self._line_editor.disconnect()

        self._line_editor.axes_update = False


        self._is_smart_engaged = False
        self._shifted = False


        self._SMART_output = [] # line, order, xmin, ymax, ymin, ymax, snr, mean, stdev, #pts 
        self._smart_i = None

        
        self._save_choices_meta = {}
        self._save_choices_meta['CHOICES'] = ['SMART output',
                                              'SNR Regions',
                                              'Line List Data',
                                              #'SNR Regions - SMART',
                                              'Original Spectrum',
                                              'Edited Spectrum',
                                              'History',
                                              'Editor Parameters', #!! add function to save all parameters
                                              'PROGRESS SAVE']

        self._save_choices_meta['SNR Regions'] = self._snr_editor._save_choices['SNR Regions']
        self._save_choices_meta['Line List'] = self._line_editor._save_choices['Line List']
        #self._save_choices_meta['SMART Output'] = ["Save"] # !! need a prefix and to save the *_continuum.dat and *_output.dat
        self._save_choices_meta['PROGRESS SAVE'] = self._save_progress


        #----------------------------------------------------------------------#
        #----------------------------------------------------------------------#
        self.init_connection_callbacks(self)

        #----------------------------------------------------------------------#
        #----------------------------------------------------------------------#
        self.keyboard_cfg = {'display':['s','o','x','h'],
                             's':'Save',
                             'o':'Open',
                             'x':'Toggle scrolling type for space and b',
                             'h':'display this help menu'}


        


        self.bindings = {}
        self.bindings['key_char'] = [None,wx.EVT_CHAR,self.key_char_callback]
#        self.bindings['mouse_any'] = [None,wx.EVT_MOUSE_EVENTS,self.mouse_any_callback]


        out = self.parent.Bind(self.bindings['key_char'][1],self.bindings['key_char'][2])
        self.bindings['key_char'][0] = out



    def _save_progress (self,sprog_dir,mode,add_to=False):
        # Define the save progress directory based on previous or new
        if add_to:
            if not sprog_dir.__class__.__name__ == 'ProgressSave': raise TypeError("Input not of correct type ProgressSave")
            sprog=sprog_dir
        else: sprog = ProgressSave(sprog_dir,mode)
        save_num = sprog.get_current_save_num()
        # check the save progress number
        if save_num is None: raise StandardError("Whoops, save_num shouldn't be None")

        # if starting here then define the routine
        if not add_to:
            sprog['ROUTINE'] = "SNR Editor Plus"
            sprog['SAVE_TYPE'] = 'manual'
            
        #----------------------------------------------------------------------#
        # DEFINE THE SPECIFIC PIECES TO GO INTO THE SAVE PROGRESS =>
        self._snr_editor._save_snr_progress_save(sprog,mode,True)
        self._line_editor._save_lines_progress(sprog,mode,True)

        filename = 'SNR_MASTER_PARAMS_'+str(save_num)+'.par'
        self._save_master_current_parameters(filename)
        sprog.add_file_to_dir(filename,save_num,'MSNR_PARAM',overwrite=True)
        #----------------------------------------------------------------------#

        # write record file if you're not adding to another
        if not add_to: sprog.write_record_file(overwrite=True)


    def _save_master_current_parameters (self,filename,clobber=True):
        filename = str(filename)
        if os.path.exists(filename) and not clobber: raise IOError("File already exists:"+filename)
        
        f = open(filename,'w')
        lines = []
        lines.append("# These are the parameters for the SNR_MasterControl. If you change the order of these then eyeSpec who't be able to read back in")
        lines.append("_use_which = "+str(self._use_which))
        lines.append("_is_smart_engaged = "+str(self._is_smart_engaged))
        lines.append("_shifted ="+str(self._shifted))
        lines.append("_smart_i = "+str(self._smart_i))
        f.write("\n".join(lines))
        f.close()

    def _open_master_current_parameters (self):
        filename = str(filename)
        if os.path.exists(filename): raise IOError("File does not exist:"+filename)
        
        f = open(filename)        
        # read comment line:
        f.readline()
        def _get_val():
            return f.readline().split('=')[1].strip()
        
        self._use_which = str(_get_val())
        self._is_smart_engaged = bool(_get_val())
        self._shifted = bool(_get_val())
        self._smart_i = int(_get_val())
        

    def _change_to_line_editor (self):
        """ returns boolean whose value depends on whether objects in the SNR editor are selected
        if they are not it returns True """
        if self._is_smart_engaged: return False
        # if SNR region is not selected and no selection box
        if not self._snr_editor.region_highlight.get_visible() and not self._snr_editor._selection_box.get_visible(): return True
        else: return False

    def _change_to_snr_editor (self):
        """ returns boolean whose value depends on whether objects in the line editor are selected
        if they are not it returns True """
        if self._is_smart_engaged: return True

        if self._line_editor._mindex is None or self._line_editor._kindex is None:
            self._line_editor._deselect_all()
            return True
        else: return False

    def _check_nearby_snr_regions (self,wl,max_dist,choose='highest SNR'):
        """  
        this finds the closests snr regions within max_dist 
        if found it will return a tuple (maximum_distance,value_i)

        """
        
        choose = str(choose)
        if choose not in ['closest','highest SNR']: raise TypeError("The variable choose must be either 'closest' or 'highest SNR'")
        wl = float(wl)
        max_dist = abs(float(max_dist))


        ordi = self._snr_editor._ordi
        if ordi is None: raise StandardError("Whoop, this shouldn't be None")

        # check the SNR regions for a newest neighbor
        snr_selections = self._snr_editor._snr_selections[ordi]
        
        if len(snr_selections) == 0:
            return False
#             self.ax.set_xlim([wl-max_dist,
#                               wl+max_dist])
#             self.update()
#             return 0

        # look for the nearest SNR region within max_dist of wl
        closest_dist = 1e20
        prev_closest = deepcopy(closest_dist)
        min_i = None

        current_snr = -1.0
        prev_max_snr = -1.0
        max_snr_i = None

        xmin,xmax = self.ax.axis()[:2]
        for i in xrange(len(snr_selections)):
            snr_sel = snr_selections[i]
            current_snr = deepcopy(snr_sel.snr)

            snr_xmin,snr_xmax = snr_sel.get_bounds()[:2]
            test_dat = np.linspace(snr_xmin,snr_xmax,100)

            min_dist = np.min(np.abs(test_dat - wl))
            closest_dist = min([closest_dist,min_dist])
            
            # check base on limit
            if closest_dist > max_dist: continue

            # if it finds a closer value
            if closest_dist < prev_closest:
                min_i = i
                prev_closest = deepcopy(closest_dist)
               
            # find maximum snr
            if current_snr > prev_max_snr:
                max_snr_i = i
                prev_max_snr = deepcopy(current_snr)
                
        # if you don't find any return False
        if min_i == None and max_snr_i == None: return False
        
        if choose == 'closest': choose_i = min_i
        if choose == 'highest SNR': choose_i = max_snr_i

        if choose_i is None: return False
        else:
            snr_sel = snr_selections[choose_i]
            snr_xmin,snr_xmax = snr_sel.get_bounds()[:2]
            test_dat = np.linspace(snr_xmin,snr_xmax,100)
            largest_dist = np.max(np.abs(test_dat - wl))
            
            xran = 1.5*largest_dist
            
            if largest_dist == 0: largest_dist = 1e-10
            return (largest_dist,choose_i)
           

    def _check_over_data (self,wl):
        wl = float(wl)
        
        # select the order that's nearby
        ordi = self._snr_editor._ordi

        # find all the orders which overlap the current wl
        overlapping_orders = []
        for i in xrange(len(self._snr_editor.plot_data)):
            xdat = self._snr_editor.plot_data[i].get_xdata()
            if np.min(xdat) < wl < np.max(xdat): overlapping_orders.append(i)

        if len(overlapping_orders) == 0: return False
           
        if not (ordi in overlapping_orders):
            if len(overlapping_orders) > 1: print "Found multiple overlapping orders"
            ordi = overlapping_orders[0]
  
        # highligh selected order
        self._snr_editor._ordi = ordi
        self._snr_editor._update_order_box(True)

        return True

    def _scan_smart (self,direction):
        max_dist = 10.0 # [\dot{A}] 
        start_direction = deepcopy(direction)

        # use line editor to move
        while True:
            self._line_editor.scan_through_lines(direction)        
            wl = self._line_editor.return_current()
            tocheck = self._check_over_data(wl)

            if tocheck: break

            # check to see if you reached an end
            kindex = self._line_editor._kindex
            if (kindex == 0 and direction == '-'):
                if start_direction == '+': 
                    print "No Data Found!!"
                    return
                print "Reached the begining of the line list which matches the data"
                # try going back the other direction
                direction = '+'

            if (kindex == len(self._line_editor.line_data)-1 and direction == '+'):
                if start_direction == '-': 
                    print "No Data Found!!"
                    return
                print "Reached the end of the line list which matches the data"
                # try going back the other direction
                direction = "-"


        # check to see if there is a nearby snr region within max_dist
        tocheck = self._check_nearby_snr_regions(wl,max_dist)

        # update the window accordingly         
        if tocheck:
            xmin,xmax,ymin,ymax = self.ax.axis()
            xran = xmax-xmin
            largest_dist = float(tocheck[0])
            deltax = largest_dist+0.1*xran
            self.ax.set_xlim([wl-deltax,
                              wl+deltax])
            self._snr_editor.select_snr_region_by_index(tocheck[1])
        else:
            self.ax.set_xlim([wl-max_dist,
                              wl+max_dist])

        self.update()

    def key_char_callback (self,event):
        self._shifted = event.ShiftDown()

    def key_press_callback (self,event):

        possible_keys = self.keyboard_cfg['display']
        if self._is_smart_engaged: possible_keys += ['b',' ','d']

        if event.key == 'backspace': event.key = 'd'
        if event.key == 'shift': self._shifted = True
        if event.key in possible_keys:
            if event.key == 'x':
                self._is_smart_engaged = not self._is_smart_engaged
                if self._is_smart_engaged: print "Smart scanning with space and 'b' engaged"
                else: print "Smart scanning with space and 'b' disengaged"

            # if smart scanning is engaged then it will access these if 
            # statements, else it will pass the key onto the other editors
            if event.key == ' ': self._scan_smart("+")
            if event.key == 'b': self._scan_smart("-")
            
            if event.key == 'd':
                if self._shifted: print "!! delete line"
                else: print "!! delete snr region"

            self.update()
        else:
            if self._change_to_snr_editor(): self._snr_editor.key_press_callback(event)
            elif self._change_to_line_editor(): self._line_editor.key_press_callback(event)

    def key_release_callback (self,event):
        if event.key == 'shift': self._shifted = False
        
    def button_press_callback (self,event):
        if self._change_to_line_editor(): self._line_editor.button_press_callback(event)

        # if nothing was selected from the lines then move onto the data
        if self._change_to_snr_editor():
            self._snr_editor.button_press_callback(event)
            self._use_which = 'snr_editor'
        elif self._change_to_line_editor:
            self._use_which = 'line_editor'
            
            #if self._line_editor.sel_line.get_visible() and self._snr_editor._ordi is not None:
            #    self._snr_editor._update_deselect_all()
            #    self.update()

    def motion_notify_callback (self,event):
        if self._use_which is None: return
        if self._use_which is 'snr_editor': 
            self._snr_editor.motion_notify_callback(event)
        #if self._use_which is 'line_editor': self._line_editor.motion_notify_callback(event)

    def button_release_callback (self,event):
        if self._change_to_snr_editor(): self._snr_editor.button_release_callback(event)
        elif self._change_to_line_editor(): self._line_editor.button_release_callback(event)
            #if self._line_editor.sel_line.get_visible() and self._snr_editor._ordi is not None:
            #    self._snr_editor._update_deselect_all()

        self.update()

    def _update_status_bar (self):
        if self._is_smart_engaged: self.parent.add_status_bar_text = format("Smart Scanning Engaged",">30")
        else: self.parent.add_status_bar_text = ''

    def update (self):
        if self._is_smart_engaged: self._line_editor.sel_line.set_color('b')
        else:  self._line_editor.sel_line.set_color('y')

        self._line_editor.update()
        self._snr_editor.update()

class EditVariancePanel (eyeSpecBaseDataPanel):
    def __init__ (self,edit_data_main_panel,edit_data_frame):
        parent_panel = edit_data_main_panel
        parent_frame = edit_data_frame
        
        eyeSpecBaseDataPanel.__init__(self, parent_panel, parent_frame)
        self.spec_obj = edit_data_frame.spec_obj
        
        self.varManager = EditVarianceManager(self,self.spec_obj,linelist=edit_data_frame.linelist)
        self.varManager.disconnect()
        
        del self.canvas.callbacks.callbacks['motion_notify_event'][self.statusbar_cid]
        self.canvas.mpl_connect('motion_notify_event',self.VarianceUpdateStatusBar)
        
    def VarianceUpdateStatusBar (self,event):
        scale_txt = " Auto Scale "
        scale_opt,info = self.varManager.ive.dp.get_auto_scaling_opt()
        if scale_opt == 0: scale_txt += 'X,Y'
        elif scale_opt == 1: scale_txt += 'X'
        elif scale_opt == 2: scale_txt += 'Y'
        elif scale_opt == 3: scale_txt += 'None'
                        
        current_order = self.varManager.ive.get_order_index()
        ord_txt = "Order: "
        if current_order is None: ord_txt += 'None'
        else: ord_txt += str(current_order)+"/" + str(self.spec_obj.shape[1] - 1)

        st = format(scale_txt,'16')+" |  "+format(ord_txt,'17')
        self.UpdateStatusBar(event,st)
 
    def OnStart (self,event):
        
#        super(eyeSpecBaseDataPanel,self).OnStart(event)
        # !! couldn't get super working???
        
        del self.canvas.callbacks.callbacks['key_press_event'][self._onkeystart_cid]
        del self.canvas.callbacks.callbacks['button_press_event'][self._onbutstart_cid]

        print ""
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "QUESTIONS TO USER:"
        questions = ["o  right now none of the data editing is saved (i.e. if you remove a point) is there any you want to retain? perhaps deleting ends of orders?",
                     "o  how much would you like to be able to undo/redo deleting SNR selections",
                     "o  One the start are there any things you want to happen?",
                     "o  how much would you like to be able to adjust SNR selections after they've initially been created?",
                     "o  should the data regions be full boxes or just bounded in wavelength (including all the data points within a given wavelength range)?",
                     "o  should the scanning forward/back be done with space/'b' or Tab/Shift_Tab"]
        print ("\n".join(questions))
        print "="*60
        print ""

        ive = self.varManager.ive
        ive.select_order_by_index(0)
        
        xy = ive.dataplot.plot_data.plot_data[0].get_xydata().T
        ymin,ymax = 0,np.max(xy[1])*1.2
        xran = 10.0
        xmin = np.min(xy[0])-0.05*xran
        xmax = np.min(xy[0])+xran
        self.varManager.ive.ax.axis([xmin,xmax,ymin,ymax])
        
        self.varManager.connect()
        self.varManager.update()

class EditVarianceMainPanel (eyeSpecBaseMainPanel):
    def __init__ (self,parent_frame):
        eyeSpecBaseMainPanel.__init__(self, parent_frame,'split_top_left')        
        self.datapanel = EditVariancePanel(self.Split1,self.pframe)
        self.randpanel = RandomPanel(self.Split1,self.pframe)
        self.split_top_left(self.datapanel,self.randpanel)

class EditVarianceFrame (eyeSpecBaseFrame):
    def __init__ (self,parent_window,inputs):
        self.spec_obj   = inputs[0]
        self.snr_regions = inputs[1]
        self.save_regions= inputs[2]
        self.linelist    = inputs[3]

        title = 'Edit Variance: '+os.path.basename(self.spec_obj.filename)
        eyeSpecBaseFrame.__init__(self, parent_window, title)
        self.panel = EditVarianceMainPanel(self)

    def OnFinish (self):
        self.Backup()
        spec_obj = self.panel.datapanel.varManager.ive.get_final_obj()
        snr_regions =self.panel.datapanel.varManager.ive.get_snr_regions()
        return (spec_obj,snr_regions)

    def Backup (self):
        print '!! save data, snr region files'
        time.sleep(.5)

#
#class AppEditVar (wx.App):
#    def __init__ (self,input):
#        self.edit_spec = input[0]
#        self.snr_regions = input[1]
#        self.save_regions = input[2]
#        self.ll_in = input[3]
#
#        redirect=False
#        filename=None
#        wx.App.__init__(self,redirect,filename)
#        
#    def Finish (self,event):
#        print "save the data"
#        self.Exit()
#
#    def OnInit (self):
#        set_title = 'Edit variance for '+os.path.basename(self.edit_spec.filename)
#        self.window = wx.Window(None)
#        self.window.SetId(10)
#        self.window.Bind(wx.EVT_WINDOW_DESTROY,self.Finish)
#
#        self.frame = CanvasFrame_SNR(None,-1,set_title,
#                                     self.edit_spec,snr_regions=self.snr_regions,save_regions=self.save_regions,linelist=self.ll_in)
#
#        self.frame.SetSize((1000,650))
#        self.frame.Show(True)
#        return True
# 
############################################################################

def edit_var (spec_obj,clean_up=True,save_regions=True,linelist=None):
    """
    This takes an eyeSpec spectrum object and by setting regions of SNR it will 
    calculate the inverse variance for all the points in the spectrum


    INPUTS:
    =============  ===============================================================
    keyword        (type) Description
    =============  ===============================================================
    spec_obj       (eyeSpec_spec) eyeSpec spectrum class object
    save_regions   (bool) If True then the SNR regions will be saved automatically
                   with either choice you have the option for manual saves
    clean_up       (bool) If True then it will remove temporary files it creates
    =============  ===============================================================


    LOGIC:
    ==============================================================================
    This code calculates inverse variance assuming Possion noise. The regions the 
    user sets define the mean and standard deviation thus the signal-to-noise-
    ratio (SNR). The standard deviation is squared to become the variance. 
    The mean and variance is interpolated between each region across the order.
    Then use equation (1) below to get the variance for each point.
    
    For example:
    > set a SNR box which gives mean = 0.98 and STD = 0.016
    > the variance is var = 0.000256, consequently the SNR = 61.25
    > then for the nearby data points the variance would be calculated
    
    new_var[i] = (0.98/data[i])*(0.000256)                (1)
    
    > the mean and var are interpolated so that will change but
      equation (1) is still the same

    NOTE 2: For very deep lines the read noise becomes significant where the 
    "true" variance for each data point should be:    

    true_var[i] =  new_var[i] + read_noise_var            (2)
    
    the eyeSpec_spec class allows to to edit the data by adding in another 
    variance to the data such as the read noise. But you can't get this just
    from the data so it was not incorporated here. Again this only affects where
    read noise is significant in the cores of large, saturated lines.
  

    """
    snr_regions=None
    load_progress=None


    # !! load_progress gives the directory name for the progress file


        
    if spec_obj is not None: 
        if spec_obj.__class__.__name__ != 'eyeSpec_spec': raise ValueError("spec MUST BE OF CLASS eyeSpec_spec")
        edit_spec = spec_obj.copy()
        save_spec(spec_obj,filename='TMP_OBJ_SAVE_ORIG',clobber=True)

    # def _run_app (edit_spec,line_data):
    set_title = 'Set Signal-to-Noise Ratio for: '+os.path.basename(edit_spec.filename)
    ##########################################
    # run application
    # _app_run_rv(spec_obj)
    app = eyeSpecBaseApp(EditVarianceFrame,[edit_spec,snr_regions,save_regions,linelist])
    sys.stdout = SysOutListener()    
    try: app.MainLoop()
    finally:
        app.ExitMainLoop()
        spec_final,snr_regions = app.Finish()
        del app

    ##########################################
    print "-"*60
    print "-"*20+format("Set SNR Complete",'^26')+"-"*20

    output_spec = None

    # load data after app.MainLoop() has exited
    if os.path.exists('TMP_OBJ_SAVE_EDIT.pkl'): print "get backup file"# load_spec('TMP_OBJ_SAVE_EDIT.pkl')

    # clean up temporary files
    try: 
        if bool(clean_up):
            if os.path.exists('TMP_WHAT_JUST_HAPPENED.txt'): os.system('rm TMP_WHAT_JUST_HAPPENED.txt')
            if os.path.exists('TMP_OBJ_SAVE_ORIG.pkl'): os.system('rm TMP_OBJ_SAVE_ORIG.pkl')
            if os.path.exists('TMP_OBJ_SAVE_EDIT.pkl'): os.system('rm TMP_OBJ_SAVE_EDIT.pkl')
    finally:
        return spec_final,snr_regions
    
     

