# this is to be an interactive version of plotspec with more options for moving around and inspecting

if __name__ != '__main__':
    from interactive_IO import ProgressSave, InputOutput
    from interactive_classes import (EventConnections, Cursor, History, KeyboardConfiguration,
                                     eyeSpecTextRedirectPanel, SysOutListener,RandomPanel,
                                     DataSelectionBox, OrderSelectionBox, PlotData, eyeSpecBaseDataPlot,
                                     eyeSpecBaseMainPanel,eyeSpecBaseDataPanel,eyeSpecBaseLineEditor,
                                     eyeSpecBaseFrame, eyeSpecBaseApp, eyeSpecBaseEventManager)
                                             

    # now import basic dependencies from other modules
    from dependencies import (np, os, sys, time, iget, deepcopy, pdb,
                              wx, FigureCanvas, NavigationToolbar2Wx, Figure, Button, Path)



class iPlotSpec:
    
    def __init__ (self,ax,spec_obj,add_xy,xy_mplkwargs):
        self.ax = ax
        self.spec_obj = spec_obj
        self.add_xy = add_xy
        self.xy_mplkwargs = xy_mplkwargs

        
        self.dp = self.dataplot = eyeSpecBaseDataPlot(self.ax,self.spec_obj,'center')
        
        for i in xrange(len(self.add_xy)):
            xy = self.add_xy[i]
            kwargs = {"lw":.9,"alpha":.5}
            if i < len(self.xy_mplkwargs): kwargs = self.xy_mplkwargs[i]
            self.ax.plot(xy[0],xy[1],**kwargs)

        
        wlmin,fmin = self.spec_obj.get_min()
        wlmax,fmax = self.spec_obj.get_max()

        self.ax.axis([wlmin,wlmin+10,fmin,fmax])


    def set_plot_limits (self,xmin,xmax,ymin,ymax):
        self.ax.axis((xmin,xmax,ymin,ymax))
        
    def interactive_plot_limits (self):
        xmin,xmax,ymin,ymax = self.ax.axis()
        # call wx dialog to get the xmin,xmax,ymin,ymax
        # get the output
        
        # apply the values which you get  
            
    def update (self):
        pass

class EditLineFrame (eyeSpecBaseFrame, EventConnections):

    def __init__ (self, parent_window, inputs):
        self.spec_obj = inputs[0]
        self.llist = inputs[1]
        self.add_xy = inputs[2]
        self.xy_mplkwarg = inputs[3]
        
        title = 'line test'
        eyeSpecBaseFrame.__init__(self,parent_window,title)

        self.figure = Figure(figsize=(4,.4)) # default size (8,6)
        self.ax = self.figure.add_subplot(111)

        wlmin,fmin = self.spec_obj.get_min()
        wlmax,fmax = self.spec_obj.get_max()
       
        self.canvas = FigureCanvas(self, -1, self.figure)
        # set up the sizer 
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas,1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()
        self.add_toolbar()
        self.cursor = Cursor(self.ax,toolbar = self.toolbar)
        self.cursor.connect()
        
        for i in xrange(len(self.add_xy)):
            xy = self.add_xy[i]
            kwargs = {"lw":.9,"alpha":.5}
            if i < len(self.xy_mplkwarg): kwargs = self.xy_mplkwarg[i]
            self.ax.plot(xy[0],xy[1],**kwargs)





        self.dataplot = eyeSpecBaseDataPlot(self.ax,self.spec_obj)
        self.ax.axis([wlmin,wlmin+10,fmin,fmax])
        
        self.leditor = eyeSpecBaseLineEditor(self.ax,self.llist[0],self.llist[1],lock_line_positions=True) # history?
        #self.leditor = eyeSpecBaseLineEditor(self.ax,self.llist[0],lock_line_positions=True) # history? 
        
        self.init_connection_callbacks(self)
        self.connect()
        self._dragging = False

        self.canvas.mpl_connect('motion_notify_event', self.UpdateStatusBar)
        self.ax.axis([wlmin,wlmin+10,fmin,fmax])

    def UpdateStatusBar (self,event,extra_text=''):
        if event.inaxes is None: return
        xpt,ypt = event.xdata,event.ydata

        # add x,y point movement
        st = "x= "+format(xpt,'10.3f')+"   y= "+format(ypt,'10.2f')
        st = format(st,'35')
        
        # add mouse/toolbar control
        if self.toolbar.mode == '': tbmod = 'editor'
        else: tbmod = self.toolbar.mode
        st += "   |  mode: "+format(tbmod,'10')

        # add extra text
        if extra_text != '': st += "  | "+str(extra_text)
        
        self.statusBar.SetStatusText(st,0)
        
    def add_toolbar(self):
        self.toolbar = NavigationToolbar2Wx(self.canvas)
        self.toolbar.Realize()
        
        tw, th = self.toolbar.GetSizeTuple()
        fw, fh = self.canvas.GetSizeTuple()

        self.toolbar.SetSize(wx.Size(fw,th))
        self.canvas.SetSize(wx.Size(fw,fh-th))

        p1,p2 = self.canvas.GetPositionTuple()
        self.canvas.SetPosition(wx.Point(p1,p2-th))

        self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.toolbar.update()

    def display_help (self):

        lines = [['space','move forward'],
                 ['b','move back'],
                 ['enter','add line'],
                 ['d','delete selected line'],
                 ['p','pan tool'],
                 ['o','zoom tool'],
                 ['.','scan through lines >'],
                 [',','scan through lines <'],
                 [';','toggle scatter of data'],
                 ['esc','deselect all'],
                 ['n','toggle lines'],
                 ['c','select line in plot'],
                 ['s','save lines']]
                 
        for line in lines:
            print format(line[0],'>10')+" : "+line[1]

    def key_press_callback (self,event):

        if event.key == 'h': self.display_help()

        if event.key == ' ':
            self.dataplot.scan_through_walking('+')
        
        if event.key == 'b':
            self.dataplot.scan_through_walking('-')

        if event.key == 'enter':
            self.leditor.plotlines.interactive_add_line()

        if event.key == 'd':
            self.leditor.delete_selected_line()

        if event.key == 'g':
            self.dataplot.toggle_grid()

        if event.key == 'p':
            self.dataplot.toggle_pan_tool()
        if event.key == 'o':
            self.dataplot.toggle_zoom_tool()
        
        if event.key == '.': 
            self.leditor.scan_through_lines("+")
        if event.key == ',':
            self.leditor.scan_through_lines("-")

        if event.key == ';':
            self.dataplot.plot_data.toggle_data_line_scatter()

        
        if event.key == 'esc':
            self.leditor.deselect()
        if event.key == 'n':
            vis = not self.leditor.get_visible()
            self.leditor.set_visible(vis)
        if event.key == 'c':
            self.leditor.select_line_in_plot()
    
        if event.key == 's':
            self.quick_save()
    
#        epsx = self.dataplot.get_epsilon('x')
#        if event.key == 'left':
#            self.leditor.move_selected_line(-1.5*epsx)
#        if event.key == 'right':
#            self.leditor.move_selected_line(1.5*epsx)

            
        self.update()
           
    def quick_save (self):
        dlg = wx.FileDialog(None, message="Save line list file as...",
                                defaultDir=os.getcwd(), wildcard='Text (*.txt)|*.txt|*',
                                defaultFile='savelines.txt', style=wx.SAVE | wx.OVERWRITE_PROMPT)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.leditor.plotlines.save_lines(path)
        dlg.Destroy()        
              
    def button_press_callback (self,event):
        if self.dataplot.is_toolbar_button_on(): return 
        eps_x = self.dataplot.get_epsilon('x')
        self._dragging = False
        
        changed = self.leditor.btn_press_lines(event,eps_x)
        if changed: self.update()
        
    def motion_notify_callback (self,event):
        if self.dataplot.is_toolbar_button_on(): return 

        changed = self.leditor.mot_notify_lines(event)
        self._dragging = True
        if changed: self.update()
        
    def button_release_callback (self,event):
        if self.dataplot.is_toolbar_button_on():
            self.update()
            return 

        eps_x = self.dataplot.get_epsilon('x')
        if not self._dragging:
            changed = self.leditor.btn_click_lines(event,eps_x)
        else:
            changed = self.leditor.btn_release_lines(event)

        if changed: self.update()

    def update (self):
        self.dataplot.update()
        self.leditor.update()
        self.ax.figure.canvas.draw()

class iPlotSpecManager (eyeSpecBaseEventManager): 
    def __init__ (self,iplotspec_panel,inputs):
        eyeSpecBaseEventManager.__init__(self)
        parent_panel = iplotspec_panel
        
        self.ppanel = parent_panel
        self.ax = parent_panel.ax
        
        self.spec_obj = inputs[0]
        self.linelist = inputs[1]
        self.add_xy = inputs[2]
        self.xy_mplkwargs = inputs[3]


        self.ips = iPlotSpec(self.ax,self.spec_obj,self.add_xy,self.xy_mplkwargs)
        
        self.has_lines = False
        if self.linelist is None: self.leditor = None
        else: 
            self.has_lines = True
            self.leditor = eyeSpecBaseLineEditor(self.ax,self.linelist[0],self.linelist[1],lock_line_positions=True)

        self.init_connection_callbacks(self)

        
        self.key_cfg.add_key('h','Display this screen')
        self.key_cfg.add_key('g','Toggle grid on/off')
        self.key_cfg.add_key('p',"Toggle pan/zoom tool (note: editor won't work while engaged")
        self.key_cfg.add_key('z',"Toggle zoom rect tool (note: editor won't work while engaged")
        self.key_cfg.add_key('q','Close and return')
        self.key_cfg.add_key(";",'Toggle data between scatter and line plot options')
        self.key_cfg.add_key('`','Toggle auto scaling options')
        self.key_cfg.add_key(" ","Scan redward in wavelength")
        self.key_cfg.add_key("b","Scan blueward in wavelength")
        self.key_cfg.add_key('esc',"Same as 'q'")
        # self.key_cfg.add_key('[','Undo data edit')
        # self.key_cfg.add_key(']','Redo data edit')
        
        # self.key_cfg.add_key('backspace','Delete selected radial velocity record')
        # self.key_cfg.add_key('esc',"same as 'q'")
        
        #self.key_cfg.set_display_order(['backspace','d','esc','g','h','p','q','z','`',';','[',']',
        self.key_cfg.set_display_order(['esc','g','h','p','q','z','`',';',' ','b'])
        
        self.key_cfg.check_display()
   
    def key_press_callback (self,event):
        pass
    
    def key_release_callback (self,event): 
        
        if event.key == 'h': self.display_help('Interactive Spectra Plot')

        elif event.key == 'g': self.ips.dp.toggle_grid()
                
        elif event.key == 'q' and 'Close' in dir(self.ppanel.pframe): self.ppanel.pframe.Close()
            
        elif event.key == 'z': self.ips.dp.toggle_zoom_tool()  

        elif event.key == 'p': self.ips.dp.toggle_pan_tool()
          
        elif event.key == '`': self.ips.dp.set_auto_scaling() 
            
        elif event.key == ';': self.ips.dp.plot_data.toggle_data_line_scatter()
          
        elif event.key == ' ': self.ips.dp.scan_through_walking('+')
        
        elif event.key == 'b': self.ips.dp.scan_through_walking('-')

          
                
        # elif event.key == : display stats on rv records
        # if event.key == '.': scan forward
        # if event.key == ',': scan backwards
        # elif event.key == : change line list  
#        elif event.key == '[': self.vre.history.undo()   
#        elif event.key == ']': self.vre.history.redo()        
        
        self.update()   
   
    def button_press_callback (self,event):
        pass
    
    def motion_notify_callback (self,event):
        pass
    
    def button_release_callback (self,event):
        
        if self.has_lines: self.leditor.update()
        
        
        
    def update (self):
        self.ips.update()
        if self.has_lines: self.leditor.update()
            
        self.ax.figure.canvas.draw()
        
class iPlotSpecPanel (eyeSpecBaseDataPanel): 
    def __init__ (self,iplotspec_main_panel, iplotspec_frame):
        parent_panel = iplotspec_main_panel
        parent_frame = iplotspec_frame
        self.inputs = iplotspec_frame.inputs
        
        eyeSpecBaseDataPanel.__init__(self, parent_panel, parent_frame)
        
        self.Manager = iPlotSpecManager(self,self.inputs)
        self.Manager.disconnect()
        
        del self.canvas.callbacks.callbacks['motion_notify_event'][self.statusbar_cid]
        self.canvas.mpl_connect('motion_notify_event',self.iPlotSpecUpdateStatusBar)

    def iPlotSpecUpdateStatusBar (self,event):
        st = ''
        self.UpdateStatusBar(event,st)

    def OnStart (self,event):
        del self.canvas.callbacks.callbacks['key_press_event'][self._onkeystart_cid]
        del self.canvas.callbacks.callbacks['button_press_event'][self._onbutstart_cid]
        
        print ""
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "QUESTIONS TO USER:"
        questions = [""]
        print ("\n".join(questions))
        print "="*60
        print ""
        
        ips = self.Manager.ips
        
        # update initial plot parameters
        
        self.Manager.connect()
        self.Manager.update()
        
class iPlotSpecMainPanel (eyeSpecBaseMainPanel):
    def __init__ (self,parent_frame):
        eyeSpecBaseMainPanel.__init__(self, parent_frame,'split_top_left')        
        self.iplotspecpanel = iPlotSpecPanel(self.Split1,self.pframe)
        self.randpanel = RandomPanel(self.Split1,self.pframe)
        self.split_top_left(self.iplotspecpanel,self.randpanel)

class iPlotSpecFrame (eyeSpecBaseFrame):
    def __init__ (self,parent_window, inputs):
        
        self.inputs = inputs
        title = 'Inspect Data '
        eyeSpecBaseFrame.__init__(self,parent_window, title)

        self.panel = iPlotSpecMainPanel(self)

    def OnFinish (self):
        self.Backup()
        print "Return data output"
        
    def Backup (self):
        # no backup
        #time.sleep()
        pass

def iplotspec (spec,linelist=None,add_xy=[],mpl_xy_kwargs={}):
    """
PURPOSE:
    Interactive plot inspection of eyeSpec spectrum objects
   
CATEGORY:
   Spectrum Analysis

INPUT ARGUMENTS:
    spec_obj : (eyeSpec_spec) An eyeSpec spectrum object

INPUT KEYWORD ARGUMENTS:
    linelist : (array, string or None) for array this will take lines to compare against [[wl,info],[wl,info],...]
                for string this will read in a file and extract the first column as wavelength
                for None then no lines will be plotted
    add_xy  : (array) this contains a list of x,y data sets to be plotted. add_xy = [[xarray,yarray],[xarray,yarray],...]
                it will perform plt.plot(add_xy[i][0], add_xy[i][1])
    mpl_xy_kwargs : (dictionary) These are the keyword arguments for matplotlib. 
                    e.g. mpl_xy_kwargs = {'color':'r','marker':'o',linestyle='--'}
                    plt.plot(add_xy[i][0], add_xy[i][1],**mpl_xy_kwargs)

OUTPUTS:
   None

DEPENDENCIES:
   External Modules Required
   =================================================
    numpy, os, wxpython, sys, operator, copy,
   
   External Functions and Classes Required
   =================================================
    EventConnections, Cursor, History, KeyboardConfiguration,
    eyeSpecTextRedirectPanel, SysOutListener,RandomPanel,
    DataSelectionBox, OrderSelectionBox, PlotData, eyeSpecBaseDataPlot,
    eyeSpecBaseMainPanel,eyeSpecBaseDataPanel,eyeSpecBaseLineEditor,
    eyeSpecBaseFrame, eyeSpecBaseApp, eyeSpecBaseEventManager
       
NOTES:
   (1)

EXAMPLE:
   >>> spec = readin("mydata.fits"
   >>> iplotspec(spec,'mylinelist.txt')
   
   >>> lines_in = np.loadtxt('mylinelist.txt',usecols=[0,1])
   >>> iplotspec(spec,lines_in)
   

MODIFICATION HISTORY:
    13, Jun 2013: Dylan Gregersen    
    """
    
    # is linelist a file name
    if type(linelist) in (np.string_,np.str,str): 
        
        if not os.path.exists(linelist):
            print "HeadsUp: filename for linelist does not exists '"+linelist+"'"
            llist = None
        else:
            llist = []
            for line in open(linelist):
                sline = line.rstrip().split()
                if len(sline) == 0 or line.strip()[0] == '#': continue
                llist.append((sline[0]," ".join(sline[1:])))
            llist = np.asarray(llist).T
            
    elif linelist is not None: llist = np.asarray(linelist)
    else: llist = None

    app = eyeSpecBaseApp(iPlotSpecFrame,[spec,llist,add_xy,mpl_xy_kwargs],False)
    sys.stdout = SysOutListener()
    try: app.MainLoop()
    finally:
        app.ExitMainLoop()
        final_out = app.Finish()
        del app    

