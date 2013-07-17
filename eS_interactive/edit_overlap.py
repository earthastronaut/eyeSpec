from ..  import Params 
from ._core import History
from eyeSpec.app_edit_data import InteractiveDataEditor, CanvasFrame
# now import basic dependencies from other modules
from eyeSpec.dependencies import np, os, sys, time, iget, deepcopy, pdb
from eyeSpec.dependencies import scipy, math
from eyeSpec.dependencies import plt, FormatStrFormatter, savefig
from eyeSpec.dependencies import pyfits, pickle
from eyeSpec.dependencies import wx, FigureCanvas, NavigationToolbar2Wx, Figure, Button, Path    


#!!  if I want to add some Icons here are defaults
ls_png = '/Library/Frameworks/Python.framework/Versions/7.3/lib/python2.7/site-packages/wx/tools/Editra/pixmaps/theme/Tango'

################################################################################
def _overlap_app_what_happened ():
    # General
    f = open('TMP_WHAT_JUST_HAPPENED.txt','w')
    f.write("This file is intended to help in the case that the code crashes and you're left with only TMP_* files\n")
    f.write("If this isn't helpful please email Dylan Gregersen <dylan.gregersen@utah.edu>\n")

    f.write("\n")

    # What was saved? Why?
    f.write("This file was created from the set overlap routine\n")
    f.write("The original spectrum information is in the file: TMP_OBJ_SAVE_ORIG.pkl\n")
    f.write("The edited spectrum information is in the file: TMP_OBJ_SAVE_EDIT.pkl\n")
    f.write("The overlap line information is in the file: TMP_FILE_FOR_LINE_DATA.txt\n")
    f.write("\n")
    
    # How to recover with these files.
    f.write("To restore what you were working on:\n")
    f.write(">>> from eyeSpec import load_spec")
    f.write(">>> orig_obj = load_spec('TMP_OBJ_SAVE_ORIG.pkl')\n")
    f.write(">>> edit_obj = load_spec('TMP_OBJ_SAVE_EDIT.pkl')\n")
    f.write(">>> output_lines = np.loadtxt('TMP_FILE_FOR_LINE_DATA.txt',unpack = True)\n")

    f.close()

################################################################################
class sudo_event:
    """
    this is used to pass to a function to simulate an event
    
    example:
    >>> # instead of: ax.figure.canvas.mpl_connect('key_press_event',foo.key_press_callback)
    >>> sudo_event.key = 'n'
    >>> foo.key_press_callback(sudo_event)
    >>> print "which act like a 'n' key hit"

    """
    pass


################################################################################
class AnOrder:
    def __init__ (self,obj,order_i,ax):
        self.ax = ax
        self.order_i = order_i
        self.wl = obj.get_wl(order_i)
        self.data = obj.get_data(order_i)
        self.inv_var = obj.get_inv_var(order_i)

    def set_data (self,obj,order_i):
        self.order_i = order_i
        self.wl = obj.get_wl(order_i)
        self.data = obj.get_data(order_i)
        self.inv_var = obj.get_inv_var(order_i)
           
    def plot_order (self,ax):
        pass



class LineEditor (EventConnections):
    """
    Plots lines based on input data and allows for manipulation of those lines

    INPUTS:
    line_data : (array like) the first dimension is used to plot lines
    parent : (CanvasFrame) technically of class wx.Frame so that you can call refreshing, if None it will try to use ax
    ax : plot axes, which it can adopt from the parent or be given explicitly 
    note - if ax and parent are not given then this will create it's own plot window for the data 
    !! but currently won't create it's own connections of keys to that window

    kwargs : matplotlib kwargs for the plot of the lines


    OUTPUTS:
    to get the editted data points:
    >>> lineeditor = LineEditor(line_data)
    >>> 'manipulate the data on the plot'
    >>> new_line_data = lineeditor.line_data # it won't be sorted if some lines are moved before others

    to get keyboard options
    >>> lineeditor.keyboard_cfg


    """
    def __init__(self,line_data,parent=None,ax=None,lock_lines=False,**kwargs):

        try: line_data = np.array(line_data,dtype=float)
        except:
            raise ValueError('ERROR: INPUT line_data MUST BE ARRAY LIKE FLOATING POINT') # !! error and exit
            
        if len(line_data.shape) !=1:
            print 'WARNING: ONLY TAKING THE FIRST DIMENSION OF line_data'
            line_data = line_data[0]

        # initiate where this stuff goes
        if parent is None and ax is None:
            print "HeadsUp: DIDN'T RECIEVE parent OR ax THEREFORE I AM CREATING A NEW PLOTTING WINDOW"
            self.ax = plt.figure().add_subplot(111)
            self.parent = None
        elif parent is None:
            self.ax = ax # and parent stays None
            self.parent = None
        elif ax is None and parent is not None:
            if parent.__class__.__name__.find("CanvasFrame") == -1: raise ValueError("LineEditor GIVEN PARENT WHICH IS NOT OF CLASS CanvasFrame")
            self.parent = parent
            self.ax = self.parent.ax

        # set up axes
        self.auto_scale_x = True
        self._temp_no_scale_x = False
        self.axes_xran = 50.0 # units = [\dot{A}], !! I could do something based on the number of pixels

        self.line_data = deepcopy(line_data)

        # how high/low do the lines go
        self.vymin = self.ax.axis()[2]
        self.vymax = self.ax.axis()[3]

        # the selection range for the lines
        self.epsilon_type = 'relative' # 'set' or 'relative'
        self.epsilon = 0.01 # 1% size
        self._epsilon = deepcopy(self.epsilon*(self.ax.axis()[1]-self.ax.axis()[0])) # the acutal range in xdata space

        # !! should the be relative too, do something similar to epsilon
        self.delta_x = 1.0 # how many to move left or right by

        # put on stuff
        self._visible = True

        self.lines = self.ax.vlines(self.line_data,self.vymin,self.vymax,lw=1.5,color='r',zorder=5)
        self.sel_line, = self.ax.plot([self.line_data[0],self.line_data[0]],[self.vymin,self.vymax],lw=10,color='y',alpha=.4,zorder=4,visible=True)

        self.text = self.ax.text(0.45, 1.05, 'Line Selected: None',
                                 transform=self.ax.transAxes, va='top')
        
        # set up indexing 
        self._mindex = None # for mouse
        self._kindex = None # for keyboard
        self._previndex = 0 # running baseline

        self._lock_lines = bool(lock_lines)

        # keep track of dragging
        self._dragged = False

        # record cursor click points
        self._recorded_pts = []
        self._recorded_lim = 20 # mostly only care about the -1 indexed point



        #---------------------------------------------------------#
        # these attributes are standard to my editor classes
        # give keyboard configuration
        self.keyboard_cfg = {'display':['d','i','x','alt','up','down','left','right'],
                             'info':'The line editor maniplates verticle lines interactively',
                             'd':'delete current selected point',
                             'i':'insert a point at the last clicked location or previously selected point',
                             'r':'update current plot',
                             #'t':"toggle on/off line data (other options won't work)",
                             'x':"toggle on/off resizing window with scrolling",
                             'alt':'hide all lines',
                             'up': 'select next point index to manipulate',
                             'down': 'select previous point index to manipulate',
                             'left':'adjust current selected line to the left by value = '+str(self.delta_x),
                             'right':'adjust current selected line to the right by value = '+str(self.delta_x)}

        if self._lock_lines:
            self.keyboard_cfg['display'] = ['d','x','alt','up','down']


        self.init_connection_callbacks(self)

        self._save_choices = {}
        self._save_choices["CHOICES"] = ["Line List"]
        self._save_choices["Line List"] = None
        self._save_choices["History"] = None
        self._save_choices["Current Line Editor Parameters"] = None
        self._save_choices["PROGRESS SAVE"] = None



    def _select_in_current_plot (self):
        """  use this to select the point nearest to the current plot's center"""
        xmin = self.ax.axis()[0]
        xmax = self.ax.axis()[1]
        mid = (xmin+xmax)/2.

        dist = np.abs(np.array(self.line_data - mid))
        ind = dist.argmin()

        X = self.line_data[ind]
        if xmin < X < xmax:
            if self._kindex != None: self._previndex = deepcopy(self._kindex)        
            self._kindex = ind
            self.sel_line.set_visible(True)
            xpt = self.line_data[self._kindex]
            self.sel_line.set_xdata([xpt,xpt]) 
            self.update()
        else:
            if self._kindex != None: self._previndex = deepcopy(self._kindex)        
            self._kindex = None


    def get_data (self):
        """
        This will return the edited data
        """        
        return deepcopy(self.line_data)

    def return_current (self):
        """
        This will return the current x of the data point
        """
        if self._kindex != None: return self.line_data[self._kindex]
        elif self._mindex != None: return self.line_data[self._mindex]
        else: return self.line_data[self._previndex]
        
    def key_press_callback(self, event):
        """ Line Editor """
        if event.key == 'backspace': event.key == 'd'
        
        #if event.key == '/':  print "!!debug >>",self._mindex,self._kindex,self._previndex
        # check if key is in the configuration, if not then ignore
        if event.key not in self.keyboard_cfg['display']: return


        # !! could have shift options
        # if event.key == 'shift': self._shift = not self._shift

        # toggle line visiblity
        if event.key == 'alt':
            if self._visible: self._visibility(False)
            else: self._visibility(True)
            return
        # if it's not visible then don't allow for any other keys to work
        if not self._visible: return

        # for deleting
        if event.key == 'd': self.delete_line()

        # if inserting
        elif event.key == 'i':
            if self._kindex != None:
                # if a current line is selected then insert to the right of that
                newx = deepcopy(self.line_data[self._kindex]+2.*self.delta_x)
            elif len(self._recorded_pts)!=0:
                # if no line selected but there has been a mouse click then put the new line at that click
                newx = deepcopy(self._recorded_pts[-1])
                self._kindex = 0 # index new point at beginning
            else:
                # take the previous line an insert to the right of that
                newx = deepcopy(self.line_data[self._previndex]+2.*self.delta_x)
                self._kindex = deepcopy(self._previndex)

            data = deepcopy(self.line_data)
            ymin,ymax = self.ax.axis()[2:]
            self.line_data = np.concatenate((data[:self._kindex],np.array([newx]),data[self._kindex:]))
            self.lines._paths.insert(self._kindex,(Path([[newx,self.vymin],[newx,self.vymax]])))
            self.sel_line.set_xdata([newx,newx])
            self.sel_line.set_visible(True)
            self._sort_data()

        # if redrawing
        # if event.key=='r': self.update() # this will just happen anyway at the end of this if chain


        # undo
        # !! for insert and delete could add an undo button
        # elif event.key == 'u'
        # need to keep around more record of what was there previously

        # toggle on/off rescaling x range of window
        elif event.key == 'x':
            self.auto_scale_x = not self.auto_scale_x

        # toggle text
        #elif event.key == 'shift':
        #    if self.text.get_visible(): self.text.set_visible(False)
        #    else: self.text.set_visible(True)

        # if moving a point
        elif event.key in ('right','left'):
            # if nothing selected then return
            if self._kindex is None:return
            if event.key =='left':
                self.line_data[self._kindex] -= self.delta_x
            else:
                self.line_data[self._kindex] += self.delta_x
            # !! this next section repeats in the following command, do I care about inefficent programming?
            X = self.line_data[self._kindex]
            self.lines._paths[self._kindex] = Path([[X,self.vymin],[X,self.vymax]])
            self.sel_line.set_xdata([X,X])

        # if moving point selection
        elif event.key in ['down','up']: self.scan_through_lines(event.key)
        self.update()



    def select_line (self,index = None, always_select=True):
        if not always_select:
            if self._kindex is None: return

        input_index = deepcopy(index)
        if index is not None:
            if self._kindex is not None: self._previndex = deepcopy(self._kindex)
            index = int(index)
            self._kindex = np.clip(index,0,len(self.line_data)-1)
            
        if self._kindex is None: self._kindex = deepcopy(self._previndex)

        X = self.line_data[self._kindex]
        self.lines._paths[self._kindex] = Path([[X,self.vymin],[X,self.vymax]])
        self.sel_line.set_xdata([X,X])


    def delete_line (self,index=None):
        input_index = deepcopy(index)
        if index is None:
            if self._kindex is None: return
            index = self._kindex
        else:
            index = int(index)
            index = np.clip(index,0,len(self.line_data)-1)
            
        # take out the value from the data array
        self.line_data[index]
        data = deepcopy(self.line_data)
        self.line_data = np.concatenate((data[:index],data[index+1:]))
        del self.lines._paths[index]

        if self._kindex is not None:
            if input_index is not None: self._kindex = index
            # re-index, selecting next lowest line
            self._kindex -= 1
            self._kindex = np.clip(self._kindex, 0, len(self.line_data)-1)

        if self._previndex != 0:
            self._previndex -= 1

        # select next lowest line
        self.select_line()

    def scan_through_lines (self,direction):
        if direction == 'up': direction = '+'
        elif direction == 'down':direction = '-'
        elif direction not in ['+','-']: raise TypeError("Whoops, given direction is not correct")


        if self._kindex is None:
            # if no line is selected then start at the previous index
            self._kindex = deepcopy(self._previndex)
        else:
            # else move up or down based on index
            if direction == '+': inc = 1
            if direction == '-': inc = -1
            # store previous index
            self._previndex = deepcopy(self._kindex)
            self._kindex += inc
            # max/min out at the ends
            self._kindex = np.clip(self._kindex, 0, len(self.line_data)-1)

        X = self.line_data[self._kindex]
        self.lines._paths[self._kindex] = Path([[X,self.vymin],[X,self.vymax]])
        self.sel_line.set_xdata([X,X])
        self.sel_line.set_visible(True)


    def _toolbar_button_on (self):
        if self.parent is not None and 'toolbar' in dir(self.parent) and self.parent.toolbar.mode != '': return True
        else: return False
    

    def button_press_callback(self, event):
        """ line editor """
        if event.inaxes==None: return
        if event.button != 1: return
        if not self._visible: return
        if self._toolbar_button_on(): return

        self._recorded_pts.append(event.xdata)
        if len(self._recorded_pts) > self._recorded_lim: del self._recorded_pts[0]

        self._mindex = self._get_ind_under_point(event)        
        self._dragged = False
        self._temp_no_scale_x = True # don't update until release, if appropriate

        # if clicked on line
        if self._mindex != None: # if clicked on line
            self._previndex = deepcopy(self._mindex)
            self._kindex = deepcopy(self._mindex)
            self.sel_line.set_visible(True)
            xpt = self.line_data[self._mindex]
            self.sel_line.set_xdata([xpt,xpt])
            self.update()
            
    def motion_notify_callback (self,event):
        """ line editor """
        if self._lock_lines: return
        # if point is not selected
        self._dragged=True
        if self._mindex is None: return
        if event.inaxes is None: return
        if event.button != 1: return 
        if not self._visible: return        
        if self._toolbar_button_on(): return



        # if not clicked on line then doesn't execute the following commands

        # if clicked on line and dragging
        #if self._get_ind_under_point(event) == None: return

        # adjust the X of the line below
        X = event.xdata
        self.line_data[self._mindex] = X
        self.lines._paths[self._mindex] = Path([[X,self.vymin],[X,self.vymax]])
        self.sel_line.set_xdata([X,X])
        self.update()


    def button_release_callback(self, event):
        """ line editor """
        if event.inaxes is None: return
        if event.button != 1: return
        if not self._visible: return
        if self._toolbar_button_on(): return

        self._temp_no_scale_x = True
        # clicked not on line and not dragging
        if self._mindex == None and not self._dragged:
            print "found this"
            self.sel_line.set_visible(False)
            self._kindex = None
            self._temp_no_scale_x = True
            
        # clicked on line => then self._kindex is set, line is visible
        if self._mindex != None: # and not self._dragged:
            self._temp_no_scale_x = False
            self._mindex == None

        # clicked on line and dragging
        #if self._mindex != None and self._dragged:
        #    self._mindex = None
            
        # clicked not on line and dragging
        # do nothing, don't erase the line 
        #if self._mindex == None and self._dragged:

        #if self._kindex == None: self._temp_no_scale_x = True

        self._dragged = False
        self.update()


    def _deselect_all (self):
        """ deselect all selected lines """
        self.sel_line.set_visible(False)
        self._mindex = None
        self._kindex = None


    def _visibility (self,truth):
        """ turn on or off the lines """
        if truth is not None: self._visible = truth
        else: self._visible = not self._visible # switch the truth of whether visible
        
        if self._visible:
            self.lines.set_visible(True)
            self.text.set_visible(True)
        else:
            if self._kindex != None:
                self._previndex = deepcopy(self._kindex)
                self._kindex = None
                
            self.lines.set_visible(False)
            self.sel_line.set_visible(False)
            self.text.set_visible(False)

        self.update()

    def _sort_data (self):
        'sorts the current line data'
        # need to adjust self._kindex
        current_mind = [0,0]
        current_kind = [0,0]
        if self._mindex != None: 
            current_mind[0] = self._mindex
            current_mind[1] = self.line_data[self._mindex]
        if self._kindex != None:
            current_kind[0] = self._kindex
            current_kind[1] = self.line_data[self._kindex]
            
        prev_ind = [self._previndex,self.line_data[self._previndex]]

        # sort the regular data
        self.line_data.sort()

        # where are the indexes now?
        if self._mindex != None: 
            new_ind = np.where(self.line_data==current_mind[1])[0]
            #if len(new_ind)!=1:
            #    print "WARNING: TWO OVERLAP" # !! perhaps be smarter with this
            self._mindex = new_ind[0]
        if self._kindex != None:
            new_ind = np.where(self.line_data==current_kind[1])[0]
            #if len(new_ind)!=1:
            #    print "WARNING: TWO OVERLAP" # !! perhaps be smarter with this
            self._kindex = new_ind[0]


        new_ind = np.where(self.line_data==prev_ind[1])[0]
        if len(new_ind)!=1:
            print "HeadsUp: TWO LINES OVERLAP" # !! perhaps be smarter with this
        self._prevind = new_ind[0]

        # sort the paths
        paths_dat = deepcopy(self.lines._paths)
        sort_it = []
        for i in range(len(paths_dat)):
            for j in paths_dat[i].iter_segments():
                sort_it.append([j[0][0],paths_dat[i]])
                break

        sort_it = sorted(sort_it)
        for i in range(len(sort_it)):
            self.lines._paths[i] = sort_it[i][1]


    def _get_ind_under_point(self, event):
        'find nearest point to the current mouse click'
        dist = np.abs(np.array(self.line_data - float(event.xdata)))
        ind = dist.argmin()

        if dist[ind] >= self._epsilon:
            return None
        else: return ind


    def _update_line_yrange (self):
        ymin,ymax = self.ax.axis()[2:]
        yran = ymax-ymin

        ymin -= yran
        ymax += yran

        self.vymin = ymin
        self.vymax = ymax

        self.sel_line.set_ydata([ymin,ymax])
        
        for i in range(len(self.lines._paths)):
            line_path = self.lines._paths[i]
            X = line_path.vertices[0][0]
            self.lines._paths[i] = Path([[X,ymin],[X,ymax]])

    def update(self):
        """ Update the line editor """
        #sorts the data and updates the plot
        self._sort_data()
        self.select_line(always_select=False)

        self._update_line_yrange()

        # update the selection epsilon
        if self.epsilon_type not in ['relative','set']: print "WARNING: self.epsilon_type MUST BE EITHER: 'relative' or 'set'" # !! error and exit?
        elif type(self.epsilon).__name__ not in ['float','int','int32','int64']: print "WARNING: self.epsilon MUST BE FLOAT OR INT"
        else:
            self.epsilon = float(self.epsilon)
            if self.epsilon_type == 'relative': self._epsilon = deepcopy(self.epsilon*(self.ax.axis()[1]-self.ax.axis()[0]))
            else: self._epsilon = deepcopy(self.epsilon)

        # update visual selection
        dataind = str(self._kindex)
        display_text = 'Line Selected: '+dataind
        if self._kindex != None: 
            display_text += "/"+str(len(self.line_data)-1)
            display_text += "  ("+format(self.line_data[self._kindex],'0.3f')+")"   
        self.text.set_text(display_text)


        #if self._kindex != None and self.auto_scale_x:
        if self.auto_scale_x and not self._temp_no_scale_x:
            current_x = self.return_current()
            #print current_x
            xmin,xmax = self.ax.axis()[0],self.ax.axis()[1]
            ran = (xmax - xmin)
            xran = (self.axes_xran)/2.0
#            if (current_x < xmin or current_x > xmax) or (ran > 3.*xran) or (ran< 3.*xran/2.0):
#                self.ax.set_xlim(current_x-xran,current_x+xran)
            if (current_x < xmin or current_x > xmax):
                self.ax.set_xlim(current_x-ran/2.,current_x+ran/2.)

        # temporarily don't update the x axis
        if self._temp_no_scale_x:
            self._temp_no_scale_x = False # if on turn it off
            return

        # make the actual plot update
        if self.parent is None: self.ax.figure.canvas.draw()
        else:
            if 'OnRefresh' in dir(self.parent): self.parent.OnRefresh(True)
            else: self.parent.Refresh()





class EditOverlapCanvasFrame (CanvasFrame):
    def __init__ (self,parent_window,title,spec_obj,line_data=None,**mpl_kwargs):
        if parent_window is None: id = -1
        else: id = parent_window.GetId()
        CanvasFrame.__init__(self,parent_window,id,title)

        #super(EditDataCanvasFrame,self).__init__(parent,id,title)
        #-------------------------------------------------#
        # add editors
        self._data_editor = InteractiveDataEditor(spec_obj,parent=self)
        self._data_editor.auto_scale_opt = 3
        self._data_editor.set_allow_edit(False)
        self._data_editor.set_allow_selection_box(False)
        self._data_editor.disconnect()


        if line_data is None: self.line_data = find_overlap_pts(obj)
        else: self.line_data = line_data
        self._line_editor = LineEditor(self.line_data,parent=self)#ax=self.ax)
        self._line_editor.disconnect()


    def OnStart (self,event):
        super(EditOverlapCanvasFrame,self).OnStart(event)
        print ""
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "QUESTIONS TO USER:"
        questions = ["o  should the shift in wavelength be relative or fixed which arrowing,\n  currently fixed at 1 Angstrom for every change, should it be 5% of the range?",
                     "o  I had to make some choices about how the window scales as\n  you move through, any suggestions?",
                     "o  How do you like the initial starting? any suggestions?",
                     "o  Let me know if you have any problems with freezing after\n  you close the window (dylan.gregersen@utah.edu), please note whether you exited by\n  closing the window or pressing 'q'",
                     "--"*20]
        print ("\n".join(questions))

        print ""


        first_order = self._data_editor.plot_data[0].get_xydata()

        xmin,xmax = np.min(first_order),np.max(first_order)
        ran = (xmax - xmin)
        self.ax.set_xlim(xmin+.1*ran,xmax+.1*ran) # semi-arbitrary starting point
        evt = sudo_event()
        #evt.key = 'up'
        #self._line_editor.key_press_callback(evt)
        #self._line_editor._kindex = 0

        self._data_editor.connect()
        self._line_editor.connect()
        self.canvas.draw()


    def OnDestroy (self,event):
        super(EditOverlapCanvasFrame,self).OnDestroy(event)
        self.Refresh()
        print "Closing: If this hangs up look at the file TMP_WHAT_JUST_HAPPENED.txt"
        _overlap_app_what_happened()
        save_spec(self._data_editor.edit_obj,filename='TMP_OBJ_SAVE_EDIT',clobber=True)
        np.savetxt("TMP_FILE_FOR_LINE_DATA.txt",self.line_data)
        event.Skip()






############################################################################

def edit_overlap (spec,line_data=None,clean_up=True):
    """
    This takes an eyeSpec spectrum object and optionally line data and allows you to 
    set the points where the orders overlap

    INPUTS:
    ===========  ===============================================================
    keyword      (type) Description
    ===========  ===============================================================
    spec_obj     (eyeSpec_spec) eyeSpec spectrum class object
    line_data    (array) Gives the lines where the overlap should occur
                 If None then it will guess at the overlap points
    clean_up     (bool) If true then it will remove temporary files it creates
    ===========  ===============================================================

    """

    if spec.__class__.__name__ != 'eyeSpec_spec': raise ValueError("spec MUST BE OF CLASS eyeSpec_spec")
        
    edit_spec = spec.copy()
    save_spec(spec,filename='TMP_OBJ_SAVE_ORIG',clobber=True)


   # def _run_app (edit_obj,line_data):
    set_title = 'Set Overlap for: '+os.path.basename(edit_spec.filename)


    ##########################################
    # run application
    # _app_run_rv(edit_obj)
    
    app = wx.App(redirect=False)
    app.SetExitOnFrameDelete(True)
    try:
        frame = EditOverlapCanvasFrame(None,set_title,edit_spec,line_data=line_data)
        frame.SetSize((1000,600))
        frame.Show(True)
        app.MainLoop()
    finally:
        app.ExitMainLoop()
        del app
    ##########################################
    print "-"*60
    print "-"*20+format("Set Overlap Complete",'^26')+"-"*20 

    output_lines = None
    output_spec = None

    # load data after app.MainLoop() has exited
    if os.path.exists('TMP_FILE_FOR_LINE_DATA.txt'): output_lines = np.loadtxt("TMP_FILE_FOR_LINE_DATA.txt",unpack = True)
    if os.path.exists('TMP_OBJ_SAVE_EDIT.pkl'): output_spec = load_spec('TMP_OBJ_SAVE_EDIT.pkl')

    # clean up temporary files
    if clean_up:
        if os.path.exists('TMP_WHAT_JUST_HAPPENED.txt'): os.system('rm TMP_WHAT_JUST_HAPPENED.txt')
        if os.path.exists('TMP_OBJ_SAVE_ORIG.pkl'): os.system('rm TMP_OBJ_SAVE_ORIG.pkl')
        if os.path.exists('TMP_OBJ_SAVE_EDIT.pkl'): os.system('rm TMP_OBJ_SAVE_EDIT.pkl')
        if os.path.exists('TMP_FILE_FOR_LINE_DATA.txt'): os.system('rm TMP_FILE_FOR_LINE_DATA.txt')

    #return outpfut_lines,obj
    #if output_lines is None and output_obj is None: return None
    #elif output_lines is None and output_obj is not None: return output_obj
    #elif output_lines is not None and output_obj is None: return output_lines
    #else: return output_obj,output_lines
    return deepcopy(output_spec),deepcopy(output_lines)
    


# RANDOM ERRORS I DON'T UNDERSTAND::

# sometimes when closing the window both using the X and 'q' the window freezes and I can't reaccess the interactive terminal without exiting

# sometimes after I run the setoverlap() it gives me:
# PyNoAppError: The wx.App object must be created first!
