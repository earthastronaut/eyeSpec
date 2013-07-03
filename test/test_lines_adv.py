execfile("main.py")
execfile("app_edit_data_v29.py")
spec = readin("../testfiles/HI.20061228.29378_1_Flux.fits")
llist = np.loadtxt("/Users/dylangregersen/Desktop/Astrophysics/data/allen/LINELISTS/lin.long.101812",unpack=True,usecols=[0,1])


smt1 = np.loadtxt("../../smart/tests/data_smt.txt/smt",unpack=True,skiprows=2)

add_xy = [smt1]


#def EditLinePanel (eyeSpecBaseDataPanel):
#    def __init__ (self,edit_data_main_panel,edit_data_frame):
#        parent_panel = edit_data_main_panel
#        parent_frame = edit_data_frame
#        
#        eyeSpecBaseDataPanel.__init__(self, parent_panel, parent_frame)
#        self.spec_obj = edit_data_frame.spec_obj
#        
#        #-------------------------------------------------#
#        # add editors
#        self.lManager = EditLineManager(self, self.spec_obj)    
#
#
#class EditLineMainPanel (eyeSpecBaseMainPanel):
#    def __init__ (self,parent_frame):
#        eyeSpecBaseMainPanel.__init__(self, parent_frame,'split_top')
#        
#        # define top panel 
#        self.datapanel = EditDataPanel(self.Split0,self.pframe)
#        self.canvas = self.datapanel.canvas
#        self.split_top(self.datapanel)

class EditLineFrame (eyeSpecBaseFrame, EventConnections):

    def __init__ (self, parent_window, inputs):
        self.spec_obj = inputs[0]
        self.llist = inputs[1]
        self.add_xy = inputs[2]
        
        title = 'line test'
        eyeSpecBaseFrame.__init__(self,parent_window,title)

        self.figure = Figure(figsize=(4,.4)) # default size (8,6)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        # set up the sizer 
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas,1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()
        self.add_toolbar()
        
        for xy in self.add_xy:
            self.ax.plot(xy[0],xy[1],lw=.9,alpha=.5)

        self.dataplot = eyeSpecBaseDataPlot(self.ax,self.spec_obj)
        self.leditor = eyeSpecBaseLineEditor(self.ax,self.llist[0],self.llist[1]) # history? 
        
        self.init_connection_callbacks(self)
        self.connect()
        self._dragging = False
      
        self.ax.axis([3990,4150,0,100000])
      
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
       
        
    def key_press_callback (self,event):
        print "=======",event.key

        if event.key == ' ':
            self.dataplot.scan_through_walking('+')
        
        if event.key == 'b':
            self.dataplot.scan_through_walking('-')

        if event.key == 'enter':
            self.leditor.plotlines.interactive_add_line()
        if event.key == 'd':
            self.leditor.delete_selected_line()
        if event.key == '.': 
            self.leditor.scan_through_lines("+")
        if event.key == ',':
            self.leditor.scan_through_lines("-")
        if event.key == 'esc':
            self.leditor.deselect()
        if event.key == 'n':
            vis = not self.leditor.get_visible()
            self.leditor.set_visible(vis)
        if event.key == 'c':
            self.leditor.select_line_in_plot()
    
        if event.key == 's':
            print "!! save"
        if event.key == 'o':
            print '!! open'
    
        epsx = self.dataplot.get_epsilon('x')
        if event.key == 'left':
            self.leditor.move_selected_line(-1.5*epsx)
        if event.key == 'right':
            self.leditor.move_selected_line(1.5*epsx)

            
        self.update()
                
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

def lines_test (spec,llist,add_xy=[]):
    app = eyeSpecBaseApp(EditLineFrame,[spec,llist,add_xy],False)
    try: app.MainLoop()
    finally: 
        print "+++ pauses here 1"
        app.ExitMainLoop()
        del app
    
