# classes and functions associated with normalizing the spectrum data

if __name__ != '__main__':
    from eyeSpec.interactive_IO import ProgressSave, InputOutput
    from eyeSpec.interactive_classes import (EventConnections, Cursor, History, KeyboardConfiguration, SysOutListener, RandomPanel,
                                             eyeSpecBaseEventManager, eyeSpecBaseDataPanel, eyeSpecBaseMainPanel, eyeSpecBaseFrame, eyeSpecBaseDataPlot)
    from eyeSpec.IO import save_spec, load_spec
    from eyeSpec.base_functions import find_overlap_pts, alt_order_colors
    from eyeSpec import __path__ as path_2_eyeSpec
    # now import basic dependencies from other modules
    from eyeSpec.dependencies import (np, os, sys, time, iget, deepcopy, pdb,
                                      scipy, math, plt, FormatStrFormatter, savefig,
                                      wx, FigureCanvas, NavigationToolbar2Wx, Figure, Button, Path, resampling)


def _ctm_what_happened ():
    # General
    f = open('TMP_WHAT_JUST_HAPPENED.txt','w')
    f.write("This file is intended to help in the case that the code crashes and you're left with only TMP_* files\n")
    f.write("If this isn't helpful please email Dylan Gregersen <dylan.gregersen@utah.edu>\n")
    f.write("\n")

    # What was saved? Why?
    f.write("This file was created from the set continuum routine\n")
    f.write("The original spectrum information is in the file: TMP_OBJ_SAVE_ORIG.pkl\n")
    f.write("The edited spectrum information is in the file: TMP_OBJ_SAVE_EDIT.pkl\n")
    f.write("\n")

    # How to recover with these files.
    f.write("To restore what you were working on:\n")
    f.write(">>> from eyeSpec import load_spec")
    f.write(">>> orig_obj = load_spec('TMP_OBJ_SAVE_ORIG.pkl')\n")
    f.write(">>> spec_obj = load_spec('TMP_OBJ_SAVE_EDIT.pkl')\n")

    f.close()

pass
#############################################################################

class SynthesisBasedGuess (wx.Dialog, EventConnections):
    def __init__ (self, linelist_file, tol=0.99):     
        self.set_tolerance(tol)
        self.set_linelist(linelist_file)
        self._model = [0,0,0,0]
        self._model_type = None
        self._model_file = 'None'
        
        self._syn_data = {}
       
    def get_tolerance (self):
        return self._tol
    
    def set_tolerance (self,tol):
        self._tol = float(tol) 

    def set_linelist (self,linelist_file):
        self._llist_file = linelist_file
        #self._llist = np.empty(0)
        # read in data
       
    def check_for_moog (self):
        pass
        
    def guess_continuum_points (self,data,ordi):
        
        # !! if your system doesn't have moog then raise issue
        
        if self._model_type == None: self.user_get_model()
        
        if ordi not in self._syn_data: self._get_syn_data_for_order(data,ordi)
       
        syndat = self._syn_data[ordi]
        starts = []
        stops = []
        
        min_num = 10

        wls = syndat[0]
        flux_mask = syndat[1] < self._tol       
        kount = 0
        starti = 0
        vals = np.arange(len(syndat[0]))[flux_mask]
        previ = vals[0]-1
        for i in vals:
            if kount == 0: starti = 0
            kount += 1
            
            if previ < 0 :continue
            if i != previ+1:
                if kount >= min_num: regions.append([starti,i])
                kount = 0
                           
            previ = i
            
        if kount >= min_num: regions.append([starti,i])
    
        ctm_pts = []
        for inds in regions:
            wl_seg = wls[inds[0]:inds[1]]
            xpt = (wl_seg[-1]+wl_seg[0])/2.0
            ypt = np.median(syndat[1][inds[0]:inds[1]])
            ctm_pts.append((xpt,ypt))
        
        
        ctm_pts = np.asarray(ctm_pts).T
        return ctm_pts
        
    def _get_syn_data_for_order (self,data,ordi):
        wlmin,wlmax = np.min(data[0]),np.max(data[0])
        xy_syn = self.get_syn_data(wlmin,wlmax)
        self._syn_data[ordi] = xy_syn
        
    def get_syn_data (self,wlmin,wlmax,cleanup=True):
        
        # get linelist_file
        # create linelist subset? or just put in self._llist_file
        
        # get_model <=== self._model_file

        # create batch.par <== linelist_file, 

        # run MOOG for synthesis
        
        # get data from synthesis and return
        return [],[]
    
    def user_get_model (self):
        # get model parameters  OR model filename
        # get model type
        # define
        self._model
        self._model_type
        self._model_file
        pass
       
class PlotContinuumDialog (wx.Dialog, EventConnections):
    
    def __init__ (self,parent,*args,**kwargs):
        super(PlotContinuumDialog,self).__init__(*args,**kwargs)
        self.parent = self.p = parent
        
        hei = self.InitUI()
        self.SetSize((400,hei))

        # !! make this smarter, based off the size of the parent
        self.SetPosition((860,30))

        self.SetTitle("Select continuum to use:")
        self.Bind(wx.EVT_WINDOW_DESTROY,self.OnClose)
        self.Bind(wx.EVT_KEY_UP,self.OnKeyUp)

        # if "I need a way to unbind the keys associated with the plot :P so the only binding is to the keys here!" 
        


        self.key_cfg = KeyboardConfiguration()
        self.key_cfg.add_key('a','Toggle select all then deselect all')
        self.key_cfg.add_key('Space','Average checked continuum and apply to current order')
        self.key_cfg.add_key('h','Display this help screen')
        self.key_cfg.add_key('q','Quit sub-window')

        self.key_cfg.set_display_order(['a','Space','h','q'])
        self.key_cfg.check_display()
        self.key_cfg.add_mpl_key_convert(' ','Space')

    def InitUI (self):
        # Get the numbers for all the orders with continuum
        cur_ctm_i = self.parent.order_ctm.get_ctm_list()
                 
        ordi = self.parent.dp.plot_data.get_current_i()
        if ordi is None: raise StandardError("This shouldn't be None if I am always selected whenever I plot")

        prev_i = -1
        options = []
        self.options = {}
        for i in cur_ctm_i:
            xmin,xmax = self.parent.spec_obj.get_min(i)[0],self.parent.spec_obj.get_max(i)[0]
            if xmin is None:
                opt = 'DELETED: '+str(i)+" can't be used"
                options.append(opt)
                self.options[i] = opt
                continue

            # check for normalization in current data
            mask = (self.parent.spec_obj.get_data(i) > 100)
            norm = 'True'
            if np.sum(mask) > 50 : norm = 'False'

            opt = []
            
            # get the current state of the order
            if i == self.parent.dp.plot_data.get_current_i():
                opt.append(format('CURRENT: '+str(i),'<10')) # order #
            else:
                opt.append(format('ORDER: '+str(i),'<10')) # order #


            opt.append("("+format(xmin,'>10.3f')+",") # min wavelength for order
            opt.append(format(xmax,'<10.3f')+")") # max wavelength for order
            opt.append(" NORM: "+format(norm,'<6')) # is it normalized?
            options.append(" ".join(opt))
            self.options[i] = " ".join(opt)

        label = "Orders with continuum ("+str(len(options))+"/"+str(self.parent.spec_obj.shape[1])+"):"

        pnl = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)
        
        sb = wx.StaticBox(pnl, label = label)
        sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)
        
        self._checkboxes = []
        self._check_evts = {}

        for i in range(len(options)):
            val = options[i]
            checkbox = wx.CheckBox(self, -1, val)

            if i == ordi:             
                if ordi in self.parent.order_ctm.keys():
                    checkbox.SetValue(self.parent.order_ctm[ordi].get_visible())
            elif ordi in self.parent._selected_order_ctm.keys():
                if i in self.parent._selected_order_ctm[ordi].keys():
                    checkbox.SetValue(True)

            self._checkboxes.append(checkbox)
            sbs.Add(self._checkboxes[i],0,wx.ADJUST_MINSIZE,0)
            self.Bind(wx.EVT_CHECKBOX,self.OnClick, self._checkboxes[i],id=100+i)
            self._check_evts[val] = [i,checkbox.GetValue()]
            
        pnl.SetSizer(sbs)


        # initiate button panel
        hbox = wx.BoxSizer(wx.HORIZONTAL)

        # buttons:
        # 1: change overplot colors
        # 2: average all selected, if only one is selected then adopt that one
        # 3: interpolate between outside most selected (must select >2)
        # 4: select all (if all selected then deselect)

        # !! multiply two continua together into a new one and normalize
        # 
        # btn0 = wx.Button(self,label='Combine')
        # btn0.Bind(wx.EVT_BUTTON, self.Combine_Continuum_Apply)
        # hbox.Add(btn0)


        btn1 = wx.Button(self,label='Average')
        btn1.Bind(wx.EVT_BUTTON, self.Average_Continuum_Apply)
        hbox.Add(btn1)
        
#         btn2 = wx.Button(self,label='Interp ')
#         btn2.Bind(wx.EVT_BUTTON, self.Interpolate_Outer_Apply)
#         hbox.Add(btn2)

        self._are_all_selected = False
        btn3 = wx.Button(self,label='Select All')
        btn3.Bind(wx.EVT_BUTTON, self.Select_All)
        hbox.Add(btn3)

        btn4 = wx.Button(self,label='Select None')
        btn4.Bind(wx.EVT_BUTTON, self.Deselect_All)
        hbox.Add(btn4)

        btn5 = wx.Button(self,label='Close') # !! maybe not necessary, just close window or when scanning through orders
        btn5.Bind(wx.EVT_BUTTON, self.OnClose)
        hbox.Add(btn5)
        
        #pngs = '/Users/dylangregersen/Desktop/Astrophysics/documents/eyeSpec/eyeSpec/icons/'
        #okButton = wx.Button(self, label = 'Average Selected+Apply')
        #self.okButton = wx.BitmapButton(self,-1,wx.Bitmap(pngs+'plus_sign.png'), style = wx.NO_BORDER)
        #self.okButton.SetSize((2,2))
                                        
        vbox.Add(pnl, proportion=1,flag=wx.ALL|wx.EXPAND, border=5)
        vbox.Add(hbox,flag=wx.ALIGN_CENTER|wx.TOP|wx.BOTTOM, border=10)
                    
        self.SetSizer(vbox)

        self.force_refresh_plot()
        
        # 18.03 is the width of a new line
        hei = 18.03*len(options)+100
        return hei

    def __len__ (self):
        return len(self.options)

    def _check_variables (self):
        all_possible = self.parent.order_ctm.keys()
        if len(all_possible) == 0: 
            print "No orders have continuum determined"
            return

        ordi = self.parent.dp.plot_data.get_current_i()
        if ordi is None: raise StandardError("Whoops, !! no order selected")

        check = np.array(map(iget(1),self._check_evts.values()))
        return (all_possible, ordi, np.any(check))

    def Average_Continuum (self):
        """
        Find the average of all the selected continua
        """
        
        # check variables
        check_var = self._check_variables()
        if check_var is None: return
        
        all_possible, ordi, check = check_var
        if not check:
            print "No orders have been clicked on"
            return None

        number_averaged  = 0
        xyfit = None
        func_type = None

        # the current order is in the options
        ordi = self.parent.dp.plot_data.get_current_i()
        # if the order is checked
        if ordi in self.options.keys():
           val = self.options[ordi]
           if val.find('DELETED') != -1 and self._check_evts[val][1]:
               # then start with that as the average
               order = self.parent.order_ctm[ordi]
               xyfit = order.get_fit_xydata()
               func_type = order.get_func_type()
               number_averaged += 1


        # find the match
        needed_interp = False
        for i in all_possible:
            val = self.options[i]
            if val.find('DELETED') != -1: continue

            # if the value is not checked skip
            if not self._check_evts[val][1]: continue

            # get the xy data from the other order
            if i == ordi: continue # it would have caught this case above
            elif ordi in self.parent._selected_order_ctm.keys():
                if i in self.parent._selected_order_ctm[ordi].keys():
                    oplot_lines = self.parent._selected_order_ctm[ordi][i]
                    oplot_xyfit = oplot_lines.get_xydata()
                else:
                    raise ValueError("Whoops, there should be a value for:"+str(i)+" because it is checked")
            else: continue
            
            func_type = None
            if xyfit == None: 
                xyfit = oplot_xyfit
                number_averaged = 1
            else:
                if not (xyfit.shape == oplot_xyfit.shape):
                    needed_interp = True
                    oplot_xyfit[1] = scipy.interp(xyfit[0],oplot_xyfit[0],oplot_xyfit[1])
                    
                xyfit[1] = (xyfit[1] + oplot_xyfit[1])/(np.ones(len(xyfit[1]))*2.0)
                number_averaged +=1
 
        if needed_interp: print "HeadsUp: applied interpolation to get average"
        return xyfit, func_type

    def Average_Continuum_Apply (self,event):
        avg_output= self.Average_Continuum()
        if avg_output is None: return
        # xyfit, func_type = avg_output
        self.Apply_Ctm(*avg_output)

#    def Interpolate_Outer (self):
#        print "!! currently not implemented"
#        if True: return
#
#        # !! ADVANCED: Create a homotopy between the points
#        all_possible,ordi,check = self._check_variables()
#        if not check:
#            print "No orders have been clicked on"
#            return
#
#        mask = np.array(map(iget(1),self._check_evts.values()))
#        all_checked = np.array(map(iget(0),self._check_evts.values()))[mask]
#
#        imin, imax = np.min(all_checked), np.max(all_checked)
#        ilen = len(np.arange(imin,imax+1))
#        
#
#        # simple:
#        # weight based on the order index so if you have the 0th ctm and the 10th 
#        # then the 7th is more like the 10th then the 0th
#        top_wei = np.arange(0,1,ilen)
#        bot_wei = np.arange(1,0,ilen)
#
#        avg_ctm = np.empty(0)
#        all_ctm = []
#        for i in np.arange(imin,imax+1):
#            if i not in self.options.keys(): continue
#            val = self.options[i]
#            if val.find("DELETED") != -1: continue
#            if not self._check_evts[val][1]: continue
#                
#            if i == ordi: continue # it would have caught this case above
#            elif ordi in self.parent._selected_order_ctm.keys():
#                if i in self.parent._selected_order_ctm[ordi].keys():
#                    ctm_xy = deepcopy(self.parent._selected_order_ctm[ordi][i].get_xydata().T)
#                else:
#                    raise ValueError("Whoops, there should be a value for:"+str(i)+" because it is checked")
#            else: continue
#
#            # !! how do I incorporate the weights to make it a 
#            if len(avg_ctm) == 0:
#                avg_ctm = deepcopy(ctm_xy)
#            else:
#                ctm_y = scipy.interp(avg_ctm[0],ctm_xy[0],ctm_xy[1])
#                avg_ctm[1] = (bot_wei[i]*avg_ctm[1] + bot_wei[i]*ctm_y)/(np.ones(len(avg_ctm[1]))*2.0)
#           
#        print "!! this needs atleast two continua, and will interpolate between them"
#        return 0,[None,None]
#
#    def Interpolate_Outer_Apply (self,event):
#        print "!! currently not implemented"
#        if True: return
#        interp_ctm,ctm_func = self.Interpolate_Outer()
#        print "!! this will apply the interp_ctm to the current order"
#        if interp_ctm is None: return
#        self.Apply_Ctm(interp_ctm,ctm_func)

    def Apply_Ctm (self,xyfit,func_type):
        if xyfit is None: return
        xyfit = np.axarray(xyfit,dtype=float)
        if xyfit.ndim != 2: raise ValueError("Continuum Dimension must be 2")

        self.parent.adopt_ctm_fit(xyfit=xyfit,
                                  func_type = func_type,
                                  method = 'modify_new',
                                  addhist = True)
        self.parent.update()
                
    def Select_All (self,event=None):
        all_possible = self.parent.order_ctm.keys()
        ordi = self.parent.dp.plot_data.get_current_i()
        for i in all_possible:
            #----------------------------------------------------#
            if i == ordi:
                # print "Treat this differently, turn on/off"
                pass
            #----------------------------------------------------#
            else:
                # add data to plot
                self.parent.add_order_continuum_to_plot(i)
            val = self.options[i]
            self._check_evts[val][1] = True
            cb = self._checkboxes[i]
            if cb.GetLabel() != val: 
                print "Whoops, check the logic"
                continue
            cb.SetValue(True)
            check = np.array(map(iget(1),self._check_evts.values()))
            if np.all(check): break
        self._are_all_selected = True
        self.force_refresh_plot()

    def Deselect_All (self,event=None):
        all_possible = self.parent.order_ctm.keys()
        ordi = self.parent.dp.plot_data.get_current_i()
        for i in all_possible:
            #----------------------------------------------------#
            if i == ordi:
                #print "Treat this differently, turn off"
                pass
            #----------------------------------------------------#
            else:
                # remove from plot
                self.parent.remove_order_continuum_to_plot(i)
            val = self.options[i]
            self._check_evts[val][1] = False

            cb = self._checkboxes[i]
            if cb.GetLabel() != val: 
                print "Whoops, check the logic"
                continue
            cb.SetValue(False)
            check = np.array(map(iget(1),self._check_evts.values()))
            if not np.any(check): break

        self._are_all_selected = False
        self.force_refresh_plot()

    def Toggle_Select_All (self,event):
        all_possible = self.parent.order_ctm.keys()
        ordi = self.parent.dp.plot_data.get_current_i()
        if self._are_all_selected: self.Deselect_All()
        else: self.Select_All()

    def OnKeyUp (self,event):    
        key = self.key_cfg.convert_wx_code_2_key(event.GetKeyCode())
        if key is None: return

        if not self.key_cfg.check_key_press_callback(key):
            return

        if key == 'q': self.Close()
        if key == 'h': self.display_help()  

        if len(self) == 0: 
            print "No Continuum To Select From"
            return

        if key == 'a': self.Toggle_Select_All(True)
        if key == ' ': self.Average_Continuum_Apply(True)
        
    def display_help (self):
        print "="*60
        editor_name = 'Continuum Calculator'
        print ("\n<<<< "+format(editor_name,'^26')+" >>>>")
        print ("-"*60)
        lines = ["This allows you to add perviously measured continuum to your",
                 "current window and then perform calcuations. 'Average' will",
                 "take the currently selected continuum and take the mean for",
                 "each data point"]
        print "\n".join(lines)
        print ("-"*60)
        self.key_cfg.display_short_info()
                  
    def _check_are_all_selected (self):
        ordi = self.parent.dp.plot_data.get_current_i()
        if ordi is None: return

        # check how many are checked
        check = np.array(map(iget(1),self._check_evts.values()))
        if np.all(check): self._are_all_selected = True
        else: self._are_all_selected = False

        # double check my the coding
        for cb in self._checkboxes:
            val = cb.GetLabel()
            if cb.GetValue() != self._check_evts[val][1]: 
                raise StandardError("Whoops, something's wrong with the logic these should parallel one another")

    def force_refresh_plot (self):
        self.p.ax.figure.canvas.draw()

    def OnClick (self,event):
        val = event.GetEventObject().GetLabel()
        ordi =  self._check_evts[val][0]

        if event.IsChecked():
            self.parent.add_order_continuum_to_plot(ordi)
            self._check_evts[val][1] = True
        else:
            self.parent.remove_order_continuum_to_plot(ordi)
            self._check_evts[val][1] = False

        self._check_are_all_selected()
        self.force_refresh_plot()
 
    def OnClose (self, event):
        self.Close()

pass
#############################################################################
    
pass   
#print "====== adding function"
#class SelectPoint:   
#    def __init__ (self,ax,**mplkwargs):
#        self.highlight_pt = 0
#        self.highlight_pt, = ax.plot(0,0,markersize=12,marker='o',linestyle='none',color='y',alpha=.8,visible=False,**mplkwargs)
#
#        self._curi = None
#        self._previ = 0
#        
#        self._kinds = {'none':[None,0]}
#        self._kind = 'none'
#     
#    def get_id (self):
#        return (self._kind,self._curi)
#     
#    def match_id (self,id):
#        return id == self.get_id()
#     
#    def add_kind (self,kind):
#        self._kinds[str(kind)] = [None,0]
#        
#    def get_kind (self):
#        return deepcopy(self._kind)
#    
#    def set_kind (self,kind):
#        # record the current values for the current kind
#        self._kinds[self._kind] = [self._curi,self._previ]
#        self.deselect_point()
#        
#        if str(kind) not in self._kinds: 
#            self.add_kind(kind)
#            return 
#        
#        # set the new kind
#        self._kind = str(kind)
#        # get the points 
#        self._curi,self._previ = self._kinds[kind]
#          
#    def get_dist (self,xpt,ypt):
#        x,y = self.get_xy()
#        return np.sqrt((x-xpt)**2 + (y-ypt)**2)    
#    def get_current_i (self):
#        return deepcopy(self._curi)
#    
#    def get_prev_i (self):
#        return deepcopy(self._previ)
#
#    def deselect_point (self):
#        if self._curi is not None: self._previ = self.get_current_i()
#        self._curi = None
#        self.highlight_pt.set_visible(False)
#
#    def is_selected (self):
#        """ return False if current index is None or the point is invisible"""
#        
#        truth = True
#        if self.get_current_i() is None or not self.highlight_pt.get_visible(): truth = False
#        return truth
# 
#    def show (self):
#        if self._curi is None: return
#        self.highlight_pt.set_visible(True)
# 
#    def hide (self):
#        if self._curi is None: return
#        self.highlight_pt.set_visible(False)
#          
#    def select_point (self,i,x,y,kind=None,always_find=False):
#        if i == None:
#            if always_find: i = self.get_prev_i()
#            else: 
#                self.deselect_point()
#                return
#        
#        if kind is not None: self.set_kind(kind) 
#        self._curi = i
#        self.highlight_pt.set_xdata(x)
#        self.highlight_pt.set_ydata(y)
#        self.highlight_pt.set_visible(True)
#
#    def get_xy (self):
#        return self.highlight_pt.get_xydata()[0]
#
#    def translate (self,xadd,yadd):
#        x,y = self.get_xy()
#        self.highlight_pt.set_xdata(x+xadd)
#        self.highlight_pt.set_ydata(y+yadd)
#        
#    def scale (self,xmult,ymult):
#        x,y = self.get_xy()
#        self.highlight_pt.set_xdata(x*xmult)
#        self.highlight_pt.set_ydata(y*ymult)    
       
class ContinuumPoints:
    def __init__ (self,kind,ax,order_index,sel_radius=0.01,**mplkwargs):        
        # for identification        
        self._ordi = order_index
        self._kind = str(kind)

        # plot points
        self.pts = 0
        self.pts, = ax.plot([],[],**mplkwargs)
        self.mplkwargs = mplkwargs
        
        # points parameters
        self._len = 0
        
        # for selecting points within this set
        self._curi = None
        self._previ = 0

        self._sel_radius = sel_radius
           
    def get_mplkwargs (self):
        return self.mplkwargs
    
    def set_mplkwargs (self,mplkwargs):
        if type(mplkwargs) != dict: raise ValueError("Didn't receive a dictionary")
        self.mplkwargs = mplkwargs
        self.pts.set(**mplkwargs)
        
    def get_id (self):
        id = (self._ordi,self._kind)   
    
    def get_order_index (self):
        return ordi
    
    def get_current_state (self):
        curstate = {'type':'points',
                    'pts':self.get_pts(),
                    'vis':self.get_visible(),
                    '_kind':self._kind,
                    '_len': self._len,
                    'mplkwargs':self.mplkwargs,
                    '_ordi':self._ordi,
                    '_curi': self._curi,
                    '_previ': self._previ,
                    '_sel_radius':self._sel_radius}
        return curstate
        
    def apply_state (self, prevstate):
        if 'type' not in prevstate or prevstate['type'] != 'points': raise ValueError("Wrong state type")
        # turn off previous points
        self.set_visible(False)
        # set new points
        self.set_pts(prevstate['pts'][0], prevstate['pts'][1])
        self.set_visible(prevstate['vis'])
        
        self._ordi = prevstate['_ordi']
        self._kind = prevstate['_kind']
        self._len = prevstate['_len']
        self.set_mplkwargs(prevstate['mplkwargs'])
        
        self._curi = prevstate['_curi']
        self._previ = prevstate['_previ']        
        
        self._sel_radius = prevstate['_sel_radius']
          
    def get_selection_radius (self):
        return deepcopy(self._sel_radius)
    
    def set_selection_radius (self,radius):
        self._sel_radius = float(radius)    
    
    def set_current_i (self,index):
        if self._curi is not None: self._previ = self.get_current_i()
        self._curi = index
    
    def get_current_i (self):
        return deepcopy(self._curi)
    
    def get_prev_i (self):
        return deepcopy(self._previ)
          
    def __str__ (self):
        return self._kind
    
    def _multi (self,arr,arrmult):
        return arr*arrmult
    
    def _add (self,arr,arradd):
        return arr+arradd
     
    def __len__ (self):
        return self._len
                      
    def get_pts (self,index=None):
        if index is None: return self.pts.get_xydata().T
        else: return self.pts.get_xydata()[index]
    
    def add_pt (self,x,y):
        X,Y = self.get_pts()
        if x in X and y in Y:
            print "Can not give two points the same value"
            return
        X = np.append(X,x)
        Y = np.append(Y,y)
        self.set_pts(X, Y)
        self._len += 1
 
    def _get_new_wl_position (self,sortit,index):
        if index is not None: index = np.where(sortit==index)[0][0]        
        return index
 
    def sort_by_wavelength (self,index=None):
        if len(self) == 0: return index
        if index not in xrange(len(self)): index = None
        
        X,Y = self.get_pts()
        sortit = np.argsort(X)
        
        self._curi = self._get_new_wl_position(sortit,self._curi)
        self._previ = self._get_new_wl_position(sortit,self._previ)
        index = self._get_new_wl_position(sortit, index)
        
        self.set_pts(X[sortit],Y[sortit])
        return index
    
    def delete_point (self,index):
        X,Y = self.get_pts()
        if len(X) == 1: 
            print "Can't delete all points"
            return
        X = np.delete(X,index)
        Y = np.delete(Y,index)
        self.set_pts(X, Y)
        
    def set_pts (self,X,Y,hist=False):
        if hist: print "!! record history"
        
        self.pts.set_xdata(X)
        self.pts.set_ydata(Y)
        
        if X.ndim == 0: self._len = 1
        else: self._len = X.shape[0]
        
    def set_visible (self,truth):
        self.pts.set_visible(truth)

    def get_visible (self):
        return self.pts.get_visible()
   
    def get_min_distance_to (self,xpt,ypt,epsx,epsy):
        scaleeps = 10.0
        epsx *= scaleeps
        epsy *= scaleeps
        outbounds = (1e20,None)
        if len(self) == 0: return outbounds
        X,Y = self.get_pts()
        
        distx = np.abs(X-xpt)
        disty = np.abs(Y-ypt)
        
        dist = np.sqrt((distx)**2+(disty)**2)
        imin = np.argmin(dist)
        if distx[imin] < epsx and disty[imin]< epsy: return dist[imin], imin
        else: return outbounds
    
    def scale (self,xmult=1.0,ymult=1.0,point=None):
        self._scale_translate(xmult, ymult, self._multi,point)

    def translate (self,xadd=0,yadd=0,point=None):
        self._scale_translate(xadd, yadd, self._add, point)    

    def _scale_translate (self,tox,toy,oper,point=None): 
        if len(self) == None: return
        X,Y = self.get_pts()
        if point is None:
            X = oper(X,tox)
            Y = oper(Y,toy)
        else:
            if point not in xrange(len(self)): return 
            X[point] = oper(X[point], tox)
            Y[point] = oper(Y[point], toy)
        self.set_pts(X, Y)    

class ContinuumPointsManager:
    def __init__ (self,ax,order_index,selection_radius=0.01):    
        self._ordi = order_index
        
        self.ax = ax
        self.default_radius = float(selection_radius)
        self._sets = {}

        self._currentlabel = None
        self.sel_pt = 0
        self.sel_pt, = self.ax.plot(0,0,linestyle='none',marker='o',visible=False,color='y',markersize=15)


    def get_id (self):
        return deepcopy(self._ordi)
    
    def get_order_index (self):
        return deepcopy(self._ordi)

    def __getitem__ (self,label):
        return self._sets[label]

    def get_current_state (self):
        pts = {}
        vis = {}
        for key in self._sets: 
            pts[key] = self._sets[key].get_current_state()
            vis[key] = self._sets[key].get_visible()
            
        curstate = {'type':'points_manager',
                    'sets':pts,
                    'sets_vis':vis,
                    'default_radius':self.default_radius,
                    'sel_xydata':self.sel_pt.get_xydata().T,
                    'sel_vis':self.sel_pt.get_visible(),
                    '_currentlabel':self._currentlabel}
        return curstate

    def apply_state (self,prevstate):
        if prevstate['type'] != 'points_manager': raise ValueError("whoops")
        self.set_visible(False)
        
        pts = prevstate['sets']
        vis = prevstate['sets_vis']
        self._sets = {}
        for key in pts: 
            if key not in self._sets: self.add_point_set(key)
            self._sets[key].apply_state(pts[key]) 
            self._sets[key].set_visible(vis[key])
            
        self.default_radius = prevstate['default_radius']
        
        xy = prevstate['sel_xydata']
        self.sel_pt.set_xdata(xy[0])
        self.sel_pt.set_ydata(xy[1])
        self.sel_pt.set_visible(prevstate['sel_vis'])
        
        self._currentlabel = prevstate['_currentlabel']
              
    def undo_redo (self,prevstate):
        curstate = self.get_current_state()
        self.apply_state(prevstate)
        return curstate

    def add_point_set (self,label,radius=None,**mplkwargs):
        label = str(label)
        mplkwargs = self._default_mplkwargs(**mplkwargs)
        
        if label == 'all': raise ValueError("Label 'all' is reserved")
        if radius is None: radius = deepcopy(self.default_radius)
        self._sets[label] = ContinuumPoints(label,self.ax,self._ordi, **mplkwargs)
        self._sets[label].set_selection_radius(radius)
        self._prevlabel = label
      
    def add_pt (self,xpt,ypt,which):
        print "huh??",xpt,ypt,which
        if which not in self.get_labels(): return
        point_set = self._sets[which]
        self._sets[which].add_pt(xpt,ypt)

    def _default_mplkwargs (self,**mplkwargs):
        mplkwargs['linestyle'] = 'none'
        if 'marker' not in mplkwargs: mplkwargs['marker']='o'
        return mplkwargs
    
    def set_plotpars (self,which,**mplkwargs):
        new = True
        if len(mplkwargs)==0: 
            new = False
            mplkwargs = self._sets[which].mplkwargs
        else: self._sets[which].mplkwargs = mplkwargs
        
        if self._currentlabel == which and not new: return
        mplkwargs = self._default_mplkwargs()
        self.sel_pt.set(**mplkwargs)

    def get_current_label (self):
        return self._currentlabel

    def get_labels (self):
        return self._sets.keys()
    
    def check_label (self,which):
        if which is 'all': return True
        return which in self._sets
    
    def get_xydata (self,which=None):
        if which is None: return self.sel_pt.get_xydata()[0]
        if which not in self.get_labels(): return
        return self._sets[which].get_pts()
    
    def get_pickradius (self,which):
        if which not in self.get_labels(): return
        return self._sets[which].get_selection_radius()

    def set_pickradius (self,which,val):
        if which not in self.get_labels(): return
        self._sets[which].set_selection_radius(val)
 
    def set_xydata (self,which,xy,hist=False):        
        if which not in self.get_labels(): return
        self._sets[which].set_pts(xy[0],xy[1],hist)
        self.deselect()
        
    def get_current_index (self,which):
        if which not in self.get_labels(): return
        return self._sets[which].get_current_i()

    def __adopt_pts (self,which,xyfit,method='modify_new',pts_set={}):
        xy = self.pts.get_xydata(which)

        xnew,ynew = xy.copy()
        if method == 'modify_new':
            if len(X) == 0: method = 'new'
            else: method = 'modify'
        
        if method == 'modify': ynew = scip.interp(xnew,xyfit[0],xyfit[1])
        elif method == 'new': print "Gah this won't work as new"
        elif method in ('set') and which in pts_set:
            xpts,ypts = pts_set[which]
            self.pts[which].set_pts(xpts,ypts)
            

    def adopt_pts (self,xyfit,method='modify_new',pts_sets={}):
        for which in self.get_labels():
            self.__adopt_pts(which, xyfit, method, pts_set)


    pass
    ##=========> on selected           
    def get_selected (self):
        if not self.is_selected(): which = None
        else: which = self._currentlabel
        
        return (which,self.get_current_index(which))
 
    def delete_selected (self):
        label, i = self.get_selected()
        if label is None: return False

        self[label].delete_point(i)
        self.__scan_through_set_of_points('-',label)
        return True
    
    def get_selected_xy (self):
        label,i = self.get_selected()
        return self._sets[label].get_pts()[:,i]
        
    pass
    ##=========> For selecting points           
    def _get_dist_index (self,which,xpt,ypt):
        xmin, xmax, ymin, ymax = self.ax.axis()
        x,y = self._sets[which].get_pts()
        if x.shape[0] == 0 or xpt is None or ypt is None or (xmax-xmin)==0 or (ymax-ymin)==0: return 99.9,0
        
        dists = ((np.abs(x-xpt))/(xmax-xmin))**2 + ((np.abs(y-ypt))/(ymax-ymin))**2
        i = np.argmin(dists)
        
        return dists[i],i
         
    def _select_xy (self,which,xpt,ypt):
        dist,i = self._get_dist_index(which, xpt, ypt)
        radius = self._sets[which].get_selection_radius()
        if dist <= radius: return (dist,i,which) 
        return None
    
    def select_by_xy (self,xpt,ypt,which='all'):        
        if not self.check_label(which): return 'error'
        
        if which == 'all': labels = self.get_labels()
        elif type(which) == str: labels = [which]
        else: labels = which

        mins = None
        for alabel in labels:
            val = self._select_xy(alabel, xpt, ypt)
            # if it wasn't close enough then just exit
            if val is None: continue
            # if no min has been found then just set it
            if mins == None: mins = val
            
            # else check if the new value is closer
            # if so then select that one
            elif val[0] < mins[0]: mins = val

        if mins == None: 
            val = 'deselect'
            if self.is_selected(): val = 'deselect changed'
            self.deselect()
            return val
        
        return self.__select_by_index(mins[2], mins[1])
    
    def select_by_index (self,which,index):
        if not self.check_label(which): return 'error'
        return self.__select_by_index(which, index)
    
    def __select_by_index (self,which,index):
        val = 'selected'
        label, curindex = self.get_selected()
      
        if label != which or index != curindex: val += ' changed'   
              
        if label != None: 
            self._prevlabel = deepcopy(self._currentlabel)
            self.deselect()
            
         
        self._currentlabel = which
        x,y = self._sets[which].get_pts(index)
        self._sets[which].set_current_i(index)
        self.sel_pt.set_xdata(x)
        self.sel_pt.set_ydata(y)
        self.set_visible(True)
        
        return val
            
    def __scan_through_set_of_points (self,direction,which=None):
        if which == None: 
            label, i = self.get_selected() 
            if label is None: return
        else: 
            label = which
            i = self._sets[which].get_current_i()
            if i is None: i = deepcopy(self._sets[which].get_prev_i())
        
        thepts = self._sets[label]
        X = thepts.get_pts()[0]
        if X.shape[0] == 0: return self.deselect()
        if direction == "+": 
            i += 1
            if i >= len(thepts): i = len(thepts)-1
        elif direction == '-':
            i-=1
            if i < 0: i = 0
        
        self.__select_by_index(label,i)
        
    def scan_through_points (self,direction,which = None):
        if which != 'all': 
            if which is not None and not self.check_label(which): raise ValueError("which is not possible value:"+str(which))
            return self.__scan_through_set_of_points(direction,which)
        
        label, i = self.get_selected()
        if label is None: 
            label = deepcopy(self._prevlabel)
            i = self._sets[label].get_current_i()
            if i is None: i = self._sets[label].get_prev_i()
            
        curx = self._sets[label].get_pts(i)
        mins = (None,)
        
        for label in self.get_labels():
            X = self._sets[label].get_pts()[0]
            if X.shape[0] == 0: continue
            if direction == '+': nextx = X[X>curx]
            else: nextx = X[X<curx]
            nexti = np.where(nextx==X)[0][0]
            if mins[0] == None or mins[0] > nextx: mins = (nextx,nexti,label)
            
        if mins[0] == None: return self.deselect()
        self.__select_by_index(mins[2],mins[1])
                    
    def set_visible (self,visible):
        self.sel_pt.set_visible(visible)
        for label in self.get_labels(): self._sets[label].set_visible(visible)
                      
    def is_selected (self):
        return self.sel_pt.get_visible()
   
    def deselect (self):
        self.sel_pt.set_visible(False)
        which = deepcopy(self._currentlabel)
        self._currentlabel = None
        if which == None: return
        
        i = self._sets[which].get_current_i()
        self._sets[which].set_current_i(None)   

    pass
    ##=========> For translating/scaling points           
    def _selected_translate_scale (self,tox,toy,val):
        xy = self.sel_pt.get_xydata()
        if val == 1:
            xy[:,0] +=tox
            xy[:,1] +=toy
        else:
            xy[:,0] *=tox
            xy[:,1] *=toy    
        self.sel_pt.set_xdata(xy[:,0])
        self.sel_pt.set_ydata(xy[:,1])        
          
    def translate (self,xadd,yadd,which=None):
        if which == None: return self.translate_selected(xadd, yadd)
        elif which == 'all': which = self.get_labels()
        else: which = [which]

        self._selected_translate_scale(xadd, yadd, 1)
        for label in which: self._sets[label].translate(xadd,yadd)
            
    def scale (self,xmult,ymult,which=None): 
        if which == None: return self.scale_selected(xadd, yadd)
        elif which == 'all': which = self.get_labels()
        else: which = [which]

        self._selected_translate_scale(xmult, ymult, 2)
        for label in which: self._sets[label].scale(xmult,ymult)
                
    def translate_selected (self,xadd,yadd):
        label, i = self.get_selected()
        if label is None: return
        self._selected_translate_scale(xadd, yadd,1)
        self._sets[label].translate(xadd,yadd,point=i)
                                       
    def scale_selected (self,xmult,ymult):
        label, i = self.get_selected()
        if label is None: return
        self._selected_translate_scale(xadd, yadd,2)
        self._sets[label].scale(xadd,yadd,point=i)
        
    def sort_by_wavelength (self,which='all'):
        if which == 'all': 
            for label in self.get_labels(): self._sets[label].sort_by_wavelength()
            return
        elif which == None:
            which,i = self.get_selected()
            if which is None: return
        
        self._sets[which].sort_by_wavelength()
                    
class CalculateContinuum:
    def __init__ (self):
        self.func_types = ['spline3',
                           'lsq spline3',
                           'piecewise poly',
                           'default',
                           None]
        self.default_func = self.func_types[2]
    
        self.no_value = [[],[]]
    
    def check_func_type (self,func_type,boolcheck=False):
        if func_type not in self.func_types:  
            print "UNKNOWN FUNCTION TYPE", str(func_type)
            if boolcheck: return False
            return self.default_func
        else: 
            if boolcheck: return True
            return func_type
    
    def get_func_types (self):
        return deepcopy(self.func_types)

    def get_next_func_type (self,func_type):
        if func_type not in self.func_types:
            print "input has no corresponding value ", func_type
            return self.default_func
        
        for i in xrange(len(self.func_types)):
            if func_type == self.func_types[i]: break
        
        i += 1
        if i >= len(self.func_types): i = 0
        return self.func_types[i]
        
    def calc_fit (self, data, order_pts, func_type):
        
        if func_type in ('default',None): func_type = 'piecewise poly'
        
        xpts,ypts = order_pts['ctm'].get_pts()
        
        if   func_type == 'spline3': X,Y = self._calc_spline_fit(xpts, ypts, data, 3, use_lsq=False)
        elif func_type == 'lsq spline3': X,Y = self._calc_spline_fit(xpts, ypts, data, 3, use_lsq=True)
        elif func_type == 'piecewise poly': X,Y = self._calc_piecewise_poly_fit(data[0], order_pts['ctm'].get_pts(), order_pts['ctl'].get_pts()[0])
                                                                                                                                         
        else: 
            print "Unknown function type "+str(func_type)
            return [],[]
        return (X, Y, func_type)
        
    def _calc_piecewise_poly_fit (self,wls,ctm_pts,ctl_pts):
        if ctm_pts.shape[1] < 3 or len(ctl_pts) == 0: return self.no_value 
        
        xpts = ctm_pts[0]
        if len(xpts) == 0: return self.no_value
        ypts = ctm_pts[1]
        
        norm = xpts[0]
        ppol = fit_piecewise_polynomial(xpts/norm,ypts,3,ctl_pts/norm)
        Y = ppol(wls/norm)
        return wls,Y
    
    def _calc_spline_fit (self,xpts,ypts,data,spline_order,use_lsq=False):
        mask = (data[0][0] < xpts)*(xpts < data[0][-1])
        
        usex = xpts[mask]
        if len(usex) <= spline_order: return [],[]
        else:
            if use_lsq: func = scipy.interpolate.LSQUnivariateSpline(data[0],data[1],uesx,data[2],k=spline_order)
            else: func = scipy.interpolate.UnivariateSpline(usex,ypts[mask],k=spline_order)
            return data[0], func(data[0])

class OrderContinuum:
    def __init__ (self,ax,order_index,zorder=0,visible=True,history=None):
        self._ordi = order_index
        self.ax = ax
        self._visible = bool(visible)
        
        if history is None: self.history = History()
        else: self.history = history 
        
        self.pts = ContinuumPointsManager(ax,order_index,selection_radius=0.001)
        self.pts.add_point_set('ctm',markersize=7,marker = 'o',linestyle='none',color='r',zorder=zorder+1,visible=visible)
        self.pts.add_point_set('ctl',markersize=12,marker = '^',linestyle='none',color='m',zorder=zorder+2,visible=visible)
                
        self.ctm_fit= 0
        self.ctm_fit, = self.ax.plot([],[],color='r',lw=2.0,zorder=zorder,visible=visible) # plot fit
        self.func_type = None

        self.calc_ctm = CalculateContinuum()

    pass
    #========= overview information ===========#   

    def get_id (self):
        ordi =  self.get_order_index()
        return ordi
    
    def get_order_index (self):
        return self._ordi
       
    def __repr__ (self):
        return "Order Continuum : "+str(self.get_func_type())

    pass
    #========= check stuff boolean ===========#   
    def is_blank (self):
        c1 = len(self.pts['ctm'])
        c2 = len(self.pts['ctl'])
        c3 = len(self.ctm_fit.get_xydata())
        return not (c1 and c2 and c3)
    
    def is_selected (self):
        return self.pts.is_selected()
    
    
    pass
    #========= deal with state information ===========#   
    def get_current_state (self):
        curstate = {'type':'order_continuum',
                    'pts':self.pts.get_current_state(),
                    'ctm_fit_xy':self.ctm_fit.get_xydata().T,
                    'ctm_fit_vis':self.ctm_fit.get_visible(),
                    'func_type':self.func_type,
                    'orderi':self._ordi,
                    '_visible':self._visible}
        return curstate

    def apply_state (self,prevstate):
        if prevstate['type'] != 'order_continuum': raise ValueError("whoops")
        
        self.set_visible(False)
        self.pts.apply_state(prevstate['pts'])
        x,y = prevstate['ctm_fit_xy']
        self.ctm_fit.set_xdata(x)
        self.ctm_fit.set_ydata(y)
        
        self.ctm_fit.set_visible(prevstate['ctm_fit_vis'])
        self.func_type = prevstate['func_type']
        
        self._ordi = prevstate['orderi']
        self._visible = prevstate['_visible']
        
    def undo_redo (self,prevstate):
        curstate = self.get_current_state()
        self.apply_state(prevstate)
        return curstate
 
    def record_history (self,info='history tag for order continuum',curstate=None):
        assert self.history != None, "Can't use this if you don't declare the history"
        if curstate is None: curstate = self.get_current_state()
        else:
            assert curstate['type'] == 'order_continuum', "Received current state of wrong class"
        self.history.add(curstate,self.undo_redo,id=self.get_id(),info=info)
        
    pass
    #========= access variables ===========#   
    def scan_func_types (self,data,direction="+"):
        func_type = self.calc_ctm.get_next_func_type(self,func_type)
        self.func_type = func_type
        self.calculate_ctm_fit(data)

    def get_func_type (self):
        return self.func_type

    def set_func_type (self,func_type):
        if self.calc_ctm.check_func_type(func_type,boolcheck=True):
            self.func_type = func_type
        else: return
        self.calculate_ctm_fit(data)
 
    def get_visible (self):
        return self._visible
    
    def set_visible (self,truth):
        self.pts.set_visible(truth)
        self.ctm_fit.set_visible(truth)
        self._visible = truth
 
    def get_selected (self):
        return self.pts.get_selected()

    pass
    #=========  ===========#   

    def adopt_ctm_fit (self,**fit_params):
    
                   
#                   xyfit=None,func_type=None,method='new',addhist=True):
        
        translate_new = [0.0,0.0]
        scale_new = [1.0,1.0]
        
        if 'xyfit' not in fit_params: raise StandardError("Sorry, but the xy data was not given here. Good luck tracking it down")
        xyfit = fit_params['xyfit']
        if xyfit.ndim != 2: raise ValueError("Continuum Dimension must be 2")
        
        
        addhist = True
        if 'addhist' in fit_params: addhist = fit_params['addhist']
        
        if addhist: self.record_history(info='adopting a continuum fit')
        
        xfit += translate_new[0]
        yfit += translate_new[1]
        
        xfit *= scale_new[0]
        yfit *= scale_new[1]
        
        self.ctm_fit.set_xdata(xfit)
        self.ctm_fit.set_ydata(yfit)
        
        self.set_visible(False)

        received_func_type = False
        if 'func_type' in fit_params:
            self.set_func_type(fit_params['func_type'])
            received_func_type = True
            
        method = 'modify_new'
        if 'method' in fit_params: method = fit_params['method']
        pts_sets = {}
        if 'pts_sets' in fit_params: pts_sets = fit_params['pts_sets']
        self.pts.adopt_pts(xyfit, method, pts_sets)

        self.set_visible(True)


    def get_ctm_xydata (self):
        return self.ctm_fit.get_xydata().T
                                         
    def get_fit_xydata (self):
        return 

    def get_xydata (self,which='selected'):
        if which in self.pts.get_labels(): return self.pts.get_xydata(which)
        elif which == 'fit': return self.ctm_fit.get_xydata().T
        elif which == 'selected': return self.pts.get_xydata()

    def select_by_xy (self,xpt,ypt,which):
        return self.pts.select_by_xy(xpt, ypt, which)

    def add_pt (self,xpt,ypt,which='ctm'):
        self.record_history('add point at '+str((xpt,ypt)))
        self.pts.add_pt(xpt, ypt, which)
        self.pts.sort_by_wavelength(which)
    
    def delete_selected (self):
        return self.pts.delete_selected()
    
    def sort_by_wavelength (self,which='all'):
        return self.pts.sort_by_wavelength(which)
    
    def scan_through_points (self,direction,which=None):
        return self.pts.scan_through_points(direction, which)
    
    pass
    #========= access variables ===========#   

    def scale_ctm (self,xmult=1.0,ymult=1.0,which=None):
        if which == 'all':
            self.pts.scale(xmult, ymult, 'all')
            print "!! scale fit"
            return
        elif which == 'pts': return self.pts.scale(xmult,ymult,'all')
        elif which == 'selected':
            self.pts.scale_selected(xmult, ymult)
            return
        self.pts.scale(xmult,ymult,which)
        
    def translate_ctm (self,xadd=0.0,yadd=0.0,which=None): 
        if which == 'all':
            self.pts.translate(xadd, yadd, 'all')
            print "!! scale fit"
            return
        elif which == 'pts': return self.pts.translate(xadd,yadd,'all')
        elif which == 'selected':
            self.pts.translate_selected(xadd, yadd)
            return
        
        self.pts.translate(xadd,yadd,which)
     
    pass
    #========= access variables ===========#    

    def __get_ypts (self,x,y,wlpts,binsize=10):
        sel_window = binsize/2.0 # A

        fpts = []
        if len(x) == 0:
            for wl in wlpts: fpts.append(1)
            return fpts
        
        for wl in wlpts:
            mask = (wl-sel_window < x)*(x < wl+sel_window)
            if not np.any(mask): ypoint = np.median(y)
            else: ypoint = np.median(y[mask])
            fpts.append(ypoint)  
        return fpts

    def __guess_ctm_pts (self,data,binsize=10):
        wlmin,wlmax = np.min(data[0]),np.max(data[0])
        wlran = wlmax-wlmin
    
        # determine the continuum points
        wlpts = np.arange(wlmin,wlmax,10)[1:-1]
        fpts = self.__get_ypts(data[0],data[1],wlpts,binsize)
        ctmpts = np.vstack((wlpts,fpts))  
        
        # determine the control points
    #            wlpts = np.linspace(wlmin,wlmax,3)
    ##            wlpts[0] += 0.05*wlran
    ##            wlpts[-1] -= 0.05*wlran
    #            fpts = get_ypts(data[0],data[1],wlpts)
    #            ctlpts = np.vstack((wlpts,fpts))
    #            
        xpts = ctmpts[0]
        num_ctl_pts = 3 
        if len(xpts) == 1: ctlpts = np.linspace(wlmin,wlmax,num_ctl_pts)
        else:
            num_ctl_pts += int(np.log(len(xpts)+1))
            if num_ctl_pts < 3: num_ctl_pts = 3
            at_start = xpts[0]
            at_end = xpts[-1]
            if at_start != at_end:
                ctlpts = np.linspace(xpts[0],xpts[-1],num_ctl_pts)
            else:
                ctlpts = np.linspace(wlmin,wlmax,num_ctl_pts)
    
        ctlpts = ctlpts[1:-1]
        ypts = self.__get_ypts(data[0],data[1],ctlpts,binsize)
        ctlpts = np.vstack((ctlpts,ypts)) 
        return ctmpts, ctlpts
      
    def guess_ctm_fit (self,data,kind='default',func_type='piecewise poly'):
        """
        Takes data = [wl,flux,inv_var] and based on the kind of guess creates a continuum
        
        Overwrites existing points
        
        
        
        """ 
        

        binsize= 10.0
        # get the function type
        func_type = self.calc_ctm.check_func_type(func_type)
            
        # get the points, both ctl and ctm for determining the fit
        if   kind == 'synthesis based': print "Use MOOG synthesis to guess ctm pts"
        elif kind == 'interactive': print "Dialog box to enter parameters"
        else: ctmpts,ctlpts = self.__guess_ctm_pts(data,binsize)
        print "here 2" 
        self.record_history(info='guess continuum')
        
        # have something called ctmpts and ctlpts and ctm_func
        self.pts.set_xydata('ctl', ctlpts, hist=False)
        self.pts.set_xydata('ctm', ctmpts, hist=False)

        self.pts.sort_by_wavelength()    
        print "should set up"
        if func_type == None: func_type = 'piecewise poly'
        self.func_type = func_type
        self.calculate_ctm_fit(data)
        print "perhaps no update?"
        return func_type
    
    def calculate_ctm_fit (self,data):
        x,y,func_type = self.calc_ctm.calc_fit(data, self.pts, self.func_type)
        self.func_type = func_type
        self.ctm_fit.set_xdata(x)
        self.ctm_fit.set_ydata(y)
        
class OrderContinuumDict (dict):
    def __init__ (self,history=None):
        if history is None: self.history = History()
        else: self.history = history
        
     
    def __setitem__ (self,key,val):
        try: key = int(key)
        except: raise ValueError("Order must be an integer")

        if repr(val).find("Order Continuum") == -1:  raise ValueError("Can only store values of OrderContinuum")
        super(OrderContinuumDict,self).__setitem__(key,val)

    def __getitem__ (self,key):
        if key not in self.keys(): raise ValueError("No Order:",key)
        else: return super(OrderContinuumDict,self).__getitem__(key)
            
    def check_key (self,key):
        return (key in self.keys())
    
    def get_ctm_list (self):
        # Get the numbers for all the orders with continuum
        cur_ctm_i = []
        
        for ordi in self.keys():
            if self[ordi].is_blank(): continue
            cur_ctm_i.append(ordi)
        
        return cur_ctm_i

    def get_current_state (self):
        ordctm = {}
        for key in self: ordctm[key] = self[key].get_current_state()
        
        curstate = {'type':'all_ctm',
                    'ordctm':ordctm}
        return curstate
    
    def apply_state (self,prevstate):
        if prevstate['type'] != 'all_ctm': raise ValueError("whoops")
        for key in self: del self[key]
        
        ordctm = prevstate['ordctm']
        rmkeys = []
        kount = 0
        for key in self:
            if key not in ordctm:
                undone.append(key)
                continue
            self[key].apply_state(ordctm[key])
            kount += 1
            
        for key in rmkeys: del self[key]
        
        if kount != len(ordctm):
            if len(self) == 0: raise ValueError("whoops, should have atleast one or rethink your logic")
            ax = self[key].ax
            for key in ordctm: 
                self[key] = OrderContinuum(ax,key,zorder=5+key)
                self[key].apply_state(ordctm[key])

    def undo_redo (self,prevstate):
        curstate = self.get_current_state()
        self.apply_state(prevstate)
        return curstate

pass
#############################################################################
                       
class InteractiveContinuumHistory:
    def __init__ (self,parent,history=None):
        self.parent = self.p = parent
        if history is None: self.history = HistoryAdvanced()
        else: self.history = history

    def get_current_state (self):
        curstate = {'type':'interative_continuum_editor',
                    'order_ctm':self.p.order_ctm.get_current_state(),
                    '_recalc_continuum':self.p._recalc_continuum,
                    '_prev_ctm_function':self.p._prev_ctm_function,
                    'data_plot':self.p.dp.get_current_state()}
                    
        return curstate
       
    def apply_state (self,prevstate):
        if prevstate['type'] != 'interactive_continuum_editor': raise ValueError("Received incorrect prev state")
        self.p.order_ctm.apply_state(prevstate['order_ctm'])
        self.p._recalc_continuum = prevstate['_recalc_continuum']
        self.p._prev_ctm_function = prevstate['_prev_ctm_function']
        self.p.dp.apply_state(prevstate['data_plot'])

    def undo_redo_data (self,prevstate):
        curstate = self.parent.dp.undo_redo(prevstate)
        
        self.parent.spec_obj = self.parent.dp.spec_obj
        
        ordi = self.p.dp.plot_data.get_current_i()
        wls = self.p.spec_obj[ordi][0]
        data = self.p.spec_obj[ordi][1]
        self.p.dp.ax.set_xlim(np.min(wls),np.max(wls))
        self.p.dp.ax.set_ylim(np.min(data),np.max(data))
        self.p.ax.figure.canvas.draw()
        
        return curstate

class InteractiveContinuumEvents:
    def __init__ (self,parent):
        self.p = parent
        self._sel_ordi = -1
        self._pressed = {'shift': False,
                         'control': False}
        self._move_pt = (-1,None)
        self.timer= Timer(0.2)
        self._diff_from_previous = [0,0]
        self._click_timer = Timer(0.2)
        
        self._prevstate = None
        
    def _check_for_order (self):
        ordi = self.p.dp.plot_data.get_current_i()
        if ordi is None: return False, None
        if ordi not in self.p.order_ctm: self.p.create_order_ctm(ordi)
        return True, ordi
    
    def _reset_vaiables (self):
        self._move_pt = [-1,None,None] # order, label, index
        self._diff_from_previous = [0,0]
        self._deselected_pt = False
  
    def check_key_press (self,key):
        return deepcopy(self._pressed[key])
  
    def key_press_ctm (self,event):
        
        # 'left' 'right' 'up' 'down' ==> 'shift' translate continuum and pts
        # 'd' ==> delete: pt, data
        # 's' ==> save stuff
        # 'o' ==> open stuff
        # 'q' ==> quit
        # 'j' ==> toggle_through_ctm_func()
        # 'm' ==> guess at continuum
        # 'n' ==> normalize the data
        #   scale, translate data
        # 'c' ==> divide by continuum
        # scan_through ==> orders, pts
        pass
  
    def key_press_move (self,event):
        ordi = self.p.dp.plot_data.get_current_i()
        if ordi is None: return
        if ordi not in self.p.order_ctm: self.p.create_order_ctm(ordi)
        order = self.p.order_ctm[ordi]

        if not order.is_selected(): return
        
        epsx,epsy = self.p.dp.get_epsilon('xy')
        
        xadd = 0
        yadd = 0
        
        scale = 1.5
        if   event.key =='right': xadd += epsx*scale
        elif event.key == 'left': xadd -= epsx*scale
        elif event.key ==   'up': yadd += epsy*scale
        elif event.key == 'down': yadd -= epsy*scale
        
        
        order.record_history('key moved a point')
        order.translate(xadd,yadd,which='selected')
        data = self.p.spec_obj[ordi]
        self.p.order_ctm[ordi].calculate_ctm_fit(data)
        
    def key_pressed (self,event,which):
        """
        This will check to see if the key is one which can have extra pressed options
        if the key is pressed then it will return True
        
        def key_press_callback (self,event):
            if key_pressed(event): return
            
        def key_release_callback (self,event):
            if key_pressed(event): return 
        
        def some_function (self,event):
            
        """
        key = event.key
        if key not in self._pressed: return False
        
        self._clicked_key = False
        if which=='press': 
            self._click_timer.reset()
            self._pressed[key] = True
            
        elif which == 'release': 
            if not self._click_timer.check(): self._clicked_key = True
            self._pressed[key] = False
            
        if self._clicked_key: return False
        return True
        
    def btn_press_ctm (self,event):
        self._reset_vaiables()
        self.p.dp.record_pts(event.xdata,event.ydata)
       
        check, ordi = self._check_for_order()
        if not check: return False
  
        order = self.p.order_ctm[ordi]
        
        out = order.select_by_xy(event.xdata,event.ydata,which='all')
        
        self._move_pt[0] = ordi           
        self._move_pt[1], self._move_pt[2] = order.get_selected()
        
        c1 = (out.find('changed')!=-1)
        c2 = (out.find('deselect')!=-1)
        
        self._prevstate = order.get_current_state()
        self._deselected_pt = (c2 and c1)
        return c1
              
    def mot_notify_ctm (self,event):
        if event.button != 1: return False
        if self._move_pt[0] == -1: return False
        if self._move_pt[1] == None: return False
        if self._move_pt[2] == None: return False
        if event.xdata is None or event.ydata is None: return False
           
        xpt,ypt = self.p.order_ctm[self._move_pt[0]].get_xydata('selected')
        
        # get the offset
        diff_x = event.xdata - xpt
        diff_y = event.ydata - ypt
        
        # determine whether the offset designates a history record
        self._diff_from_previous[0] -= diff_x
        self._diff_from_previous[1] -= diff_y
        
        # moving one point or everything?
        if self._pressed['shift']: kind = 'pts'
        else: kind = None

        self.p.order_ctm[self._move_pt[0]].translate_ctm(xadd=diff_x, yadd=diff_y, which=kind)
        
        if self.timer.check(): 
            data = self.p.spec_obj[self._move_pt[0]]
            self.p.order_ctm[self._move_pt[0]].calculate_ctm_fit(data)
        return True
        
    def btn_release_ctm (self,event):
        if round(self._diff_from_previous[0],2) == 0 and round(self._diff_from_previous[0],2)==0: return False        
    
        ordi = self._move_pt[0]
        
        self.p.order_ctm[ordi].sort_by_wavelength()
        # record history
        self.p.order_ctm[ordi].record_history(info='moved pt',curstate=self._prevstate)        
        
        data = self.p.spec_obj[ordi]
        self.p.order_ctm[ordi].calculate_ctm_fit(data)
        return True
        
    def btn_click_ctm (self,event):
        if event.inaxes is None: return False        
                        
        check, ordi = self._check_for_order()
        if not check: return False
        
        if self._deselected_pt: return False       
        if self.p.order_ctm[ordi].is_selected(): return False
    
        # add continuum point for left click
        if event.button == 1: self.p.order_ctm[ordi].add_pt(event.xdata,event.ydata,'ctm')
        
        # add control point for right click
        if event.button == 3: self.p.order_ctm[ordi].add_pt(event.xdata,event.ydata,'ctl')
        
        if event.button in (1,3): 
            data = self.p.spec_obj[ordi]
            self.p.order_ctm[ordi].calculate_ctm_fit(data)
        return True
           
    def manipulate_multi_order_ctm (self,event):
        dlg = PlotContinuumDialog(self.p,None,-1,'Add continuum to plot:')
        dlg.ShowModal()
        dlg.Destroy()    
        
class InteractiveContinuumIO:
    def __init__ (self,parent):
        self.parent = parent
        self.iochoices = {}
    
        # get save choices for spectrum object from the dataplot
    
    def save_progress (self,sprog_dir,mod,add_to=False): pass
    
    def save_order_ctm (self,filename,clobber=True): pass
    def open_order_ctm (self,filename): pass
    
    def open_snr_regions_2_ctm (self,filename): pass
    def save_selected_order_ctm (self,filename=None,clobber=True): pass
    def open_selected_order_ctm (self,filename): pass
    
    def save_current_parameters (self,filename,clobber=True): pass
    def open_current_parameters (self,filename): pass
          
class InteractiveContinuumEditor (object):
    def __init__ (self,parent_panel,ax,spec_obj,auto_scale_focus='full',history=None):
        self.ax = ax
        self.dataplot = self.dp = eyeSpecBaseDataPlot(ax, spec_obj, auto_scale_focus=auto_scale_focus, plotclass=PlotDataSingleOrder)
        self.spec_obj = self.dp.spec_obj

        self.hist = InteractiveContinuumHistory(self,history)
        self.history = self.hist.history
        self.evts = InteractiveContinuumEvents(self)        
        self.io = InteractiveContinuumIO(self)


        self._prev_ctm_function = None

        self.ppanel = parent_panel
        self.add_to_toolbar()
        
        self._delta_scale = 1.0
        self.update_deltas()
        
        self._store_ctm_state = None
        
        # determine continuum
        self.order_ctm = OrderContinuumDict(self.history) # ordi = obj, build as you go
        self._selected_order_ctm = {} # for a given ord_sel_i which others are selected

        for i in xrange(self.spec_obj.shape[1]): self.create_order_ctm(i)    

        self._start_pts =[0,0]
        self._diff_from_prev = [0,0]
        self.text = self.ax.text(0.05, 1.05, 'Order Selected: None',
                                 transform=self.ax.transAxes, va='top')

        self._recac_continuum = True

    def shifted (self):
        return self.evts._pressed['shift']

    def create_order_ctm (self,ordi):
        self.order_ctm[ordi] = OrderContinuum(self.ax,ordi,zorder=5+ordi,history=self.history)

    def add_to_toolbar (self):
        self.ppanel.toolbar.AddSeparator()

        # !! need an icon
        #ls_png = '/Library/Frameworks/Python.framework/Versions/7.3/lib/python2.7/site-packages/wx/tools/Editra/pixmaps/theme/Tango/toolbar/'
        ls_png = path_2_eyeSpec[0]+'/eS-data/images/'

        texit = self.ppanel.toolbar.AddLabelTool(-1, '', wx.Bitmap(ls_png+'find.png'))
        self.ppanel.toolbar.Realize()
        self.ppanel.Bind(wx.EVT_TOOL, self.evts.manipulate_multi_order_ctm, texit)


        choices = self._order_choice_list()
        # !! add a drop down menu for selecting which order you want
        # ctl = wx.Choice(self.parent,-1,(85,18),choices=choices)
        #self.ppanel.toolbar.Add(ctl)
        # in self.update add a thing to update the choice list

        self.ppanel.toolbar.update()

    def _order_choice_list (self):
        out = []
        for i in range(self.spec_obj.shape[1]):
            if len(self.spec_obj.get_wl(i)) == 0: continue
            out.append(i)
        return out

    def _get_offset_ctm (self,prev_ordi):
        if prev_ordi not in self.order_ctm.keys():
            print "!! Couldn't find previous order"
            return

        # take the previous continuum
        #prev_ctm = self.order_ctm[prev_ordi]
        #prev_xmid = (np.max(prev_ctm.wl) + np.min(prev_ctm.wl))/2.0
        
        # get the current plot data
        ordi = self.dp.plot_data.get_current_i()
        if ordi is None: return 0,0,1,1
        spec_data = self.spec_obj[ordi]
        prev_data = self.spec_obj[prev_ordi]
        
        xadd,yadd,xmult,ymult = self.order_ctm[prev_ordi].get_offset_estimate(xdata,ydata)

        xadd *= -1
        yadd *= -1
        xmult = 1.0/xmult
        ymult = 1.0/ymult
        return xadd,yadd,xmult,ymult

    def _plot_objects (self,ordi,visible=True):
        if ordi not in self.order_ctm: return
        self.order_ctm[ordi].set_visible(visible)

    def insert_point (self):
        ordi = self.dp.plot_data.get_current_i()
        if ordi not in self.order_ctm: self.create_order_ctm(ordi)
        order = self.order_ctm[ordi]
        
        if not order.is_selected(): return False
        x,y = order.get_selected_xy()

        label = 'ctm' # label, i = pts.get_selected() 
        xnew = x + self.delta_x
        
        order.add_pt(xnew,y,label)
            
    def guess_ctm_fit (self,kind='default'):
        ordi = self.dp.plot_data.get_current_i()
        if ordi not in self.order_ctm: self.create_order_ctm(ordi)
        
        data = self.spec_obj[ordi]
        self.order_ctm[ordi].guess_ctm_fit(data,kind,func_type=self._prev_ctm_function)
     
    def remove_order_continuum_to_plot (self,prev_ordi):
        # get the current order
        current_ordi = self.dp.plot_data.get_current_i()
        if current_ordi is None: return

        if prev_ordi == current_ordi:
            # leave always displayed
            return

        if prev_ordi not in self._selected_order_ctm[current_ordi].keys():
            return

        self._selected_order_ctm[current_ordi][prev_ordi].set_visible(False)
        del self._selected_order_ctm[current_ordi][prev_ordi]
         
    def add_order_continuum_to_plot (self,prev_ordi):  
        # get the current order
        current_ordi = self.dp.plot_data.get_current_i()
        if current_ordi is None: return

        # see if the previous order has a continuum
        if prev_ordi not in self.order_ctm.keys(): return

        # check if the previous order selected is actually the current
        if prev_ordi == current_ordi:
            # should always be visible
            return

        # see if the current order has an entry for selected order ctm
        if current_ordi not in self._selected_order_ctm.keys():
            self._selected_order_ctm[current_ordi] = {}
            
        # see if the previous order is in the selected order ctm for this 
        if prev_ordi  in self._selected_order_ctm[current_ordi].keys(): return
        
        
        xprev,yprev = self.order_ctm[prev_ordi].get_ctm_xydata()
        self.spec_obj.set_use_cropped(False)
        data = self.spec_obj[current_ordi]
        self.spec_obj.set_use_cropped("prev")

        ymid_prev = np.median(yprev)
        ymid = np.median(data[1])
        ynew = yprev/ymid_prev*ymid

        
        # if they match for pixels then just apply a scaling
        if xprev.shape == data[0].shape: xnew = data[0]
        else:
            # shift the data then resample
            xadd = data[0][0]-xprev[0]
            xnew = xprev+xadd
            
            T = resampling.get_resampling_matrix(xnew,data[0],preserve_normalization=True)
            ynew = T*ynew
            xnew = data[0]
            
        
        add_order = self.__add_selected_ctm_plot(xnew,ynew)
        self._selected_order_ctm[current_ordi][prev_ordi] = add_order
        self.update()

    def __add_selected_ctm_plot (self,x,y,visible=True):
        add_order, = self.ax.plot(x,y,color='r',lw=2.0,linestyle='--',visible=visible)
        return add_order

    def scan_through_orders (self,direction):
        starti = self.dp.plot_data.get_current_i()
        if not self.dp.plot_data.scan_through_orders(direction): return     
        if starti is not None: self._plot_objects(starti,False)
                
        ordi = self.dp.plot_data.get_current_i()
        self._plot_objects(ordi,True)
        self.update_plot_bounds()
        self.update_text()
 
    def scan_through_points (self,direction):
        ordi = self.dp.plot_data.get_current_i()
        if ordi is None: return
        order = self.order_ctm[ordi]
        order.scan_through_points(direction)
      
    def toggle_func_types (self):
        ordi = self.dp.plot_data.get_current_i()
        if ordi is None: return
        data = self.spec_obj[ordi]
        self.order_ctm[ordi].scan_func_types(data)
        self._prev_ctm_function = self.order_ctm[ordi].get_func_type()
         
    def select_order (self,index,force=False):
        ordi = self.dp.plot_data.get_current_i()
        self.dp.plot_data.select_order(index)
        
        # if there was no change
        if index == ordi and not force: return 

        if index not in self.order_ctm: self.create_order_ctm(index)
        
        self._plot_objects(ordi,False)
        self._plot_objects(index, True)
        self.update_plot_bounds()
        self.update_text()
        
    def adopt_ctm_fit (self,ordi='current',**fit_params):
        if ordi not in self.order_ctm:
            print "Whoops, given order is outside range"
            return
    
        order_ctm[ordi].adopt_ctm_fit(**fit_params)
        
    
    def delete_selected (self):
        ordi = self.dp.plot_data.get_current_i()
        # is order selected
        if ordi is None: return False
        if self.order_ctm[ordi].delete_selected(): return True
        else: print " !!Delete order"
        
    def divide_by_ctm (self):
        ordi = self.dp.plot_data.get_current_i()
        # is order selected
        if ordi is None: return False
        
        xy = self.order_ctm[ordi].get_fit_xydata()
        if len(xy[0]) == 0: return False
        
        x.set_use_cropped(False)
        wls = self.spec_obj.get_wl(ordi)
        data = self.spec_obj.get_data(ordi)
        if xy[0].shape == wls.shape and np.all(xy[0]==wls): newflux = data/xy[1]
        else:
            T = resampling.get_resampling_matrix(xy[0],wls,preserve_normalization=True)
            ctmflux = T*xy[1]
            newflux = data/ctmflux
        
        
        self.history.add(self.dp.get_current_state(),self.hist.undo_redo_data,id='ice_spec1',info='divide by continuum')
        self.spec_obj.set_data(newflux,oid=ordi)        
        x.set_use_cropped("prev")
    
        self.dp.plot_data.update_plot_data(self.spec_obj)        
        self.ax.set_ylim((np.min(newflux),np.max(newflux)))
  
    def update_plot_bounds (self):
        self.ax.axis(self.dp.plot_data.get_data_bounds())

    def update_deltas (self):
        self.delta_x = self._delta_scale*self.dp.get_epsilon('x')
        self.delta_y = self._delta_scale*self.dp.get_epsilon('y')

    def update_text (self):
        # update visual selection
        dataind = str(self.dp.plot_data.get_current_i())
        display_text = 'Order Selected: '+dataind
        if dataind != 'None': display_text += "/"+str(self.spec_obj.shape[1]-1)
        self.text.set_text(display_text)        

    def update (self): 
        c1 = self.dp.update()
        self.update_text()
        self.update_deltas()
        return c1               

pass
#############################################################################

class EditContinuumManager (eyeSpecBaseEventManager):
    def __init__ (self,ctm_panel, spec_obj):
        eyeSpecBaseEventManager.__init__(self)
        self.ppanel = ctm_panel
        self.ax = ctm_panel.ax
        
        self.ice = InteractiveContinuumEditor(self.ppanel,self.ax,spec_obj)
        self.init_connection_callbacks(self) 
        self.key_cfg = KeyboardConfiguration()       
        self.key_cfg.add_key('d','Does same as backspace')

        self.key_cfg.set_display_order(['d'])
        
        self.key_cfg.check_display()
        # self.key_cfg.add_mpl_key_convert('<','< and >')
        # self.key_cfg.add_mpl_key_convert('>','< and >')

        T = time.gmtime()
        curdate = format(T.tm_mon,'02')+format(T.tm_mday,'02')+str(T.tm_year) 
         
        self._clicked = False
        
    def delete_item (self):
        print "Delete something"
    
    def display_tips (self): 
        l = []
        l.append("="*60)
        l.append("    TIPS FOR SETTING THE CONTINUUM")
        l.append("- "*30)
        l.append("")
        l.append("o  Control points are given as purple trianges, continuum points are red circles")
        l.append("o  Using the piecewise polynomial means you're breaking the continuum up into pieces based on the control points where each piece is a polynomial. The coefficients of each polynomial piece are determined by the continuum points.")
        l.append("    -  You should have continuum points greater than the largest control point and less than the smallest control point")
        l.append("o  Using a spline3 fit means you're using the continuum points with a 3rd order spline function")
        l.append("="*60)
        
        for line in l: print line
        
    def key_press_callback (self,event):
        if self.ice.evts.key_pressed(event,'press'): return

    def key_release_callback (self,event):
        if self.ice.evts.key_pressed(event,'release'): return
        
        print "====clicked <"+event.key+">"
        
        if   event.key == ',': self.ice.scan_through_orders('-')
        elif event.key == '.': self.ice.scan_through_orders('+')
        elif event.key == 'd': self.ice.delete_selected()
        elif event.key == 'g': self.ice.dp.toggle_grid()
        elif event.key == 'c': self.ice.divide_by_ctm()
        elif event.key == 'h': print "help! "*10
        elif event.key == 'q': self.ppanel.pframe.Close()
        elif event.key == 'v': self.ice.insert_point()
        elif event.key == 't': self.display_tips()
        #elif event.key == '\t': # tab is reserved for shifting between screens
        elif event.key in ('left','right','up','down'): self.ice.evts.key_press_move(event)
        elif event.key == ' ': self.ice.scan_through_points('+')           
        elif event.key == 'b': self.ice.scan_through_points('-')
        elif event.key == 'm': self.ice.guess_ctm_fit()
        elif event.key == '[': self.ice.history.undo()
        elif event.key == ']': self.ice.history.redo()
        
        self.update()
                
    def key_held_callback (self,event): pass
    
    def button_press_callback (self,event):
        """ CTM editor """
        self._clicked = True
        if self.ice.dp.is_toolbar_button_on(): return
        c1 = self.ice.evts.btn_press_ctm(event)
        if c1: self.update()
    
    def motion_notify_callback (self,event):
        """  CTM Editor """
        self._clicked = False
        if self.ice.dp.is_toolbar_button_on(): return
        c1 = self.ice.evts.mot_notify_ctm(event)
        if c1: self.update()

    def button_release_callback (self,event):
        if self.ice.dp.is_toolbar_button_on(): return
        if self._clicked: 
            self.button_click_callback(event)
            return 
        
        self._clicked = False
        c1 = self.ice.evts.btn_release_ctm(event)
        if c1: return True
        
    def button_click_callback (self,event):
        c1 = self.ice.evts.btn_click_ctm(event)
        if c1: self.update()
        
    def update (self):
        self.ice.update()
        self.ax.figure.canvas.draw()
        
class EditContinuumPanel (eyeSpecBaseDataPanel):
    def __init__ (self,edit_data_main_panel,edit_data_frame):
        parent_panel = edit_data_main_panel
        parent_frame = edit_data_frame
        
        eyeSpecBaseDataPanel.__init__(self, parent_panel, parent_frame)
        self.spec_obj = edit_data_frame.spec_obj
        
        self.ctmManager = EditContinuumManager(self,self.spec_obj)
        self.ctmManager.ice.select_order(0,True)
        
        self.ctmManager.disconnect()
        
        del self.canvas.callbacks.callbacks['motion_notify_event'][self.statusbar_cid]
        self.canvas.mpl_connect('motion_notify_event',self.VarianceUpdateStatusBar)

    def VarianceUpdateStatusBar (self,event):
        scale_txt = " Auto Scale "
        scale_opt,info = self.ctmManager.ice.dp.get_auto_scaling_opt()
        if scale_opt == 0: scale_txt += 'X,Y'
        elif scale_opt == 1: scale_txt += 'X'
        elif scale_opt == 2: scale_txt += 'Y'
        elif scale_opt == 3: scale_txt += 'None'
                        
        current_order = self.ctmManager.ice.dp.plot_data.get_current_i()
        ord_txt = "Order: "
        if current_order is None: ord_txt += 'None'
        else: ord_txt += str(current_order)+"/" + str(self.spec_obj.shape[1] - 1)

        st = format(scale_txt,'16')+" |  "+format(ord_txt,'17')
        self.UpdateStatusBar(event,st)

    def OnStart (self,event):
        self.post_start_delete_cids()
        
        print ""
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "QUESTIONS TO USER:"
#        questions = ["o  right now none of the data editing is saved (i.e. if you remove a point) is there any you want to retain? perhaps deleting ends of orders?",
#                     "o  how much would you like to be able to undo/redo deleting SNR selections",
#                     "o  One the start are there any things you want to happen?",
#                     "o  how much would you like to be able to adjust SNR selections after they've initially been created?",
#                     "o  should the data regions be full boxes or just bounded in wavelength (including all the data points within a given wavelength range)?",
#                     "o  should the scanning forward/back be done with space/'b' or Tab/Shift_Tab"]
#        print ("\n".join(questions))
        print "="*60
        print ""

#        ive = self.varManager.ive
#        ive.select_order_by_index(0)
#        
#        xy = ive.dataplot.plot_data.plot_data[0].get_xydata().T
#        ymin,ymax = 0,np.max(xy[1])*1.2
#        xran = 10.0
#        xmin = np.min(xy[0])-0.05*xran
#        xmax = np.min(xy[0])+xran
#        self.varManager.ive.ax.axis([xmin,xmax,ymin,ymax])
#        

#
#        x = np.arange(4020,4080,20)
#        y = np.ones(len(x))*35000
#        for i in xrange(x.shape[0]):
#            self.ctmManager.ice.order_ctm[0].pts.add_pt(x[i],y[i],'ctm')
#            
        self.ctmManager.connect()
        self.ctmManager.update() 
        self.ctmManager.ice.dp.update()    
        
class EditContinuumMainPanel (eyeSpecBaseMainPanel):
    def __init__ (self,parent_frame):
        eyeSpecBaseMainPanel.__init__(self, parent_frame,'split_top_left')        
        self.datapanel = EditContinuumPanel(self.Split1,self.pframe)
        self.randpanel = RandomPanel(self.Split1,self.pframe)
        self.split_top_left(self.datapanel,self.randpanel)

class EditContinuumFrame (eyeSpecBaseFrame):
    def __init__ (self,parent_window,inputs):
        self.spec_obj = inputs
        title = "Edit Continuum: "+os.path.basename(self.spec_obj.filename)
        eyeSpecBaseFrame.__init__(self,parent_window,title)
        self.panel = EditContinuumMainPanel(self)
    
    def OnFinish (self):
        self.Backup()
        return 1
    
    def Backup (self):
        print "!! save data"
        _ctm_what_happened()
        time.sleep(.5)

############################################################################

def edit_ctm (spec_obj,clean_up=True):
    """
PURPOSE:
    This takes a eyeSpec spectrum object and allows you to normalize by a polynomial
      
CATEGORY:
   Spectral Reductions

INPUT ARGUMENTS:
    spec_obj : (eyeSpec_spec) An eyeSpec spectrum object

INPUT KEYWORD ARGUMENTS:
    clean_up : (boolean) If true then it will remove temporary files it creates, files 'TMP_*'

OUTPUTS:
    (eyeSpec_spec) Edited eyeSpec spectrum object
    
DEPENDENCIES:
   External Modules Required
   =================================================
    numpy, os, sys, time, operator, copy, pdb,
    scipy, math, matplotlib, 
    wx, resampling   
   
   External Functions and Classes Required
   =================================================
    EventConnections, Cursor, History, KeyboardConfiguration,
    eyeSpecBaseEventManager, eyeSpecBaseDataPanel, eyeSpecBaseMainPanel, eyeSpecMainFrame     
    EditContinuumFrame, EditContinuumMainPanel, EditContinuumPanel, EditContinuumManager
    InteractiveContinuumEditor, InteractiveContinuumHistory, InteractiveContinuumIO, InteractiveContinuumEvents
    OrderContinuumDict, OrderContinuum, CalculateContinuum, ContinuumPointsManager, ContinuumPoints,
    PlotContinuumDialog, SynthesisBasedGuess
       
NOTES:
   (1) 

EXAMPLE:
   >>> spec = readin("my_data.txt")
   >>> spec_editted = edit_ctm(spec)

MODIFICATION HISTORY:
    13, Jun 2013: Dylan Gregersen
    """

    if spec_obj.__class__.__name__ != 'eyeSpec_spec': raise ValueError("obj MUST BE OF CLASS eyeSpec_spec")
        
    edit_obj = spec_obj.copy()
    save_spec(spec_obj,filename='TMP_OBJ_SAVE_ORIG',clobber=True)


   # def _run_app (spec_obj,line_data):
    set_title = 'Set Continuum for: '+os.path.basename(edit_obj.filename)

    ##########################################
    # run application
    # _app_run_rv(spec_obj)
    
    app = eyeSpecBaseApp(EditContinuumFrame,edit_obj)
    sys.stdout = SysOutListener()
    try: app.MainLoop()
    finally:
        app.ExitMainLoop()
        final_spec = app.Finish()
        del app
            
    ########################################## 
    print "-"*60
    print "-"*20+format("Set Continuum Complete",'^26')+"-"*20


    output_spec_obj = None
    # load data after app.MainLoop() has exited
    if os.path.exists('TMP_OBJ_SAVE_EDIT.pkl'): 
        output_spec_obj = load_spec('TMP_OBJ_SAVE_EDIT.pkl')

    # clean up temporary files
    if clean_up:
        if os.path.exists('TMP_WHAT_JUST_HAPPENED.txt'): 
            os.system('rm TMP_WHAT_JUST_HAPPENED.txt')
        if os.path.exists('TMP_OBJ_SAVE_ORIG.pkl'): 
            os.system('rm TMP_OBJ_SAVE_ORIG.pkl')
        if os.path.exists('TMP_OBJ_SAVE_EDIT.pkl'): 
            os.system('rm TMP_OBJ_SAVE_EDIT.pkl')
    return deepcopy(output_spec_obj)
    


