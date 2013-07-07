# These are functions used by the wx applications to edit data

# Main todo: separate the DataPlot class to have one class which controls all the plot attributes and another which adds the plot data 

################################################################################
# Import Modules

import pdb #@UnusedImport
from eyeSpec.plotting import figure_adjust_borders, alt_order_colors, data_line_scatter #@UnresolvedImport
from eyeSpec.core import Timer #@UnresolvedImport
from eyeSpec.IO import save, readin_spec #@UnresolvedImport
from eyeSpec.dependencies import (os, sys, time, deepcopy, np, math, pickle, np_vstack_delete, #@UnresolvedImport
                                  plt, Figure, FormatStrFormatter, threading, Queue, NewEvent, wx, FigureCanvas, NavigationToolbar2Wx) #@UnresolvedImport
from eyeSpec import __path__ as path_2_eyeSpec #@UnresolvedImport
from _IO import SaveOpenChoices, InputOutput

pass
################################################################################

class State (dict):
    def __init__ (self,ID):
        self.id = ID
        
    def get_id (self):
        return deepcopy(self.id)    
        
    def check_id (self,ID):
        if not self.is_id(ID): raise ValueError("Received wrong state id :"+str(id))
        
    def is_id (self,ID):
        return (self.id == ID)   

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

class RandomPanel(wx.Panel):

    #----------------------------------------------------------------------
    def __init__(self, parent,color="#99CCB2"):
        """Constructor"""
        wx.Panel.__init__(self, parent)
        self.SetBackgroundColour(color)

class HistoryAdvanced:
    """
    This class creates a history of actions and allows for recalling
    """

    def __init__ (self):
        self._past = [] # [[prevstate, func_id, info]]
        self._future = []
        
        self._functions = {}
        
        self.max_history = 30
        self.start_deleting_at = 1

    def __repr__ (self):
        return "eyeSpecHistoryAdvanced"
        
    def add_function (self,undo_redo_func,ID='skip', unique=True):
        """
        The function is relative to a class
        it will take the previous state of the class, apply that state, and return the current state
        
        The ID is included to help you identify it later
        
        
        example:
        
        def undo_redo (self,prevstate):
            curstate = self.get_current_state()
            self.apply_previous_state(prevstate)
            return curstate
        
        """
        
        if ID in ('skip', None): 
            func_id == 'skip' #@UndefinedVariable
            undo_redo_func = None
        else:
            func_string = str(undo_redo_func)
            if func_string.find('bound method') == -1 : raise ValueError("This should have been a bound method")
            ID = str(ID)
            func_id = (func_string,str(ID))
            
            if unique:
                j = 0
                newid = deepcopy(ID)+"_"+str(j)
                crash_at = 40
                while func_id in self._functions:
                    j += 1
                    if j > crash_at: raise StandardError("Crash in loop")
                    newid = deepcopy(ID)+"_"+str(j)
                    func_id = (func_string,newid)
        
        self._functions[func_id] = undo_redo_func
        return func_id
        
    def get_function_ids (self):
        return self._functions.keys()
        
    def add (self, prevstate, undo_redo_func, ID='', info=''):
        func_id = self.add_function(undo_redo_func,ID,unique=False)
        self._past.append([prevstate, func_id, info])
        self._future = []
        
        while len(self._past) > 30:
            i = int(self.start_deleting_at)
            del self._past[i]

    def undo_n (self,n):
        for _ in xrange(int(n)): 
            out = self._undo_redo('undo')
            if out == 0: break

    def redo_n (self,n): 
        for _ in xrange(int(n)): 
            out = self._undo_redo('redo')
            if out == 0: break
        
    def _undo_redo (self, which):
        if which == 'undo': record = self._past
        elif which == 'redo': record = self._future
        else: raise ValueError("whoops, gave other than undo or redo")    
        
        if (len(record) <= 0):
            print "Can't "+which+" more"
            return
        
        prevstate, func_id, info = record[-1]
        
        if func_id is 'skip': return
        
        if func_id not in self._functions: raise StandardError("received unknown function id "+str(func_id))
        undo_redo_func = self._functions[func_id]
        
        try: current_state = undo_redo_func(prevstate)
        except:
            print "="*60
            print "Error with function for func_id: "+str(func_id)
            print "Using pdb, type 'c' to get the error from python"
            print "="*60
            pdb.set_trace()
            current_state = undo_redo_func(prevstate)
                        
        if which == 'undo':
            print "UNDO: "+str(info)+"  with "+str(func_id)
            if current_state != None:
                self._future.append([current_state,func_id,info])
                del self._past[-1] 
            return len(self._past)           
        else:
            print "REDO: "+str(info)+"  with "+str(func_id)
            if current_state != None:
                self._past.append([current_state,func_id,info])
                if len(self._future) > 0: del self._future[-1]            
            return len(self._future)

    def undo (self):
        """
        Apply a previous recorded state to be the current using the given function
        """
        self._undo_redo('undo')

    def redo (self):
        """
        Apply a state to be the current using the given function
        """
        self._undo_redo('redo')
        
    def display (self):
        print "\n======================================="
        print "="*20+format("PAST",'^9')+"="*20
        for i in xrange(len(self._past)): print i,self._past[i][2],self._past[i][1]

        print ""
        print "="*20+format("FUTURE",'^9')+"="*20
        for i in reversed(xrange(len(self._future))): print i,self._future[i][2],self._future[i][1]
   
    def save_history (self,filename,clobber=True):
        filename = str(filename)
        if os.path.exists(filename) and not clobber: raise IOError("File exists: "+filename)
        saved_obj = [self.max_history,self.start_deleting_at,self._past,self._future]
        pickle.dump(saved_obj,open(filename,'wb'))

    def open_history (self,filename):
        if not os.path.exists(filename): raise IOError("File does not exist:"+filename)
        if len(self._functions) == 0: raise IOError("functions must be set first so that they can be referenced")

        saved_obj = pickle.load(open(filename,'rb'))
        # !! check the save _ obj
        past = saved_obj[2]
        future = saved_obj[3]
        
        def _checkit_ (timetravel):
            timetravel_out = []
            for i in xrange(len(timetravel)):
                _prevstate,func_id,_info = timetravel[i]
                if func_id in self._functions: timetravel_out.append(timetravel[i])
                else: print "Error: I don't have the function to deal with:"+str(func_id)
            return timetravel

        self._past = _checkit_(past)
        self._future = _checkit_(future)
        self.max_history = saved_obj[0]
        self.start_deleting_at = saved_obj[1]
        
class History:
    """
    This class creates a history of actions and allows for recalling
    """

    def __init__ (self):
    
    
        self._past = [] # [[hist_info,hist_data,info]]
        self._future = []
        self._options = {} # add_func, undo_func, redo_func
        self._option_names = {} # !! not currently necessary
        self.max_history = 30
        self.start_deleting_at = 1

    def __repr__ (self):
        return "eyeSpecHistory"

    def get_options (self):
        return deepcopy(self._options)
    
    def set_options (self,hist_info,undo_func,redo_func=None):
        """
        for which routine define how to undo and redo 

        def undo_func (hist_data):
             # takes hist_data and makes it the current state
             # creates redo_data from current state
             return redo_data

        redo_func does the same only returning a past_data 
             
        """
        hist_info = str(hist_info)
        if redo_func is None:
            redo_func = undo_func

        self._options[hist_info] = [undo_func,redo_func]
        self._option_names[hist_info] = [str(undo_func),str(undo_func)]
        
    def add (self,hist_info,hist_data,add_to='undo',info=''):
        add_to = str(add_to)
        if add_to not in ['past','future','undo','redo']: raise ValueError("Variable add_to must be either 'past' or 'future'")


        hist_info = str(hist_info)
        if hist_info not in self._options.keys():
            print "history info not in options: "+str(hist_info)
            return

        if add_to in ['past','undo']:
            self._past.append([hist_info,hist_data,str(info)])
            self._future = []
        
            while len(self._past) > 30:
                i = int(self.start_deleting_at)
                del self._past[i]
            
        elif add_to in ['future','redo']:
            self._future.append([hist_info,hist_data,str(info)])

    def undo (self):
        if len(self._past) <= 0:
            print "Can't undo more"
            return

        hist_info, prev_state, info = self._past[-1]
        
        if hist_info not in self._options.keys():
            print "Whoops, don't know what to do with: "+hist_info
            return

        if hist_info == 'skip':
            return

        undo_func = self._options[hist_info][0]
        
        try: current_state = undo_func(prev_state)
        except:
            print "Error with undo function for hist_info: "+hist_info
            print "Using pdb, type 'cont' to get the error from python"
            pdb.set_trace()
            current_state = undo_func(prev_state)

        print "UNDO: "+hist_info+" : "+str(info)
        if current_state != None:
            self._future.append([hist_info,current_state,info]) # [hist_info,hist_data]
            del self._past[-1]

    def redo (self):
        if len(self._future) == 0: 
            print "Can't redo more"
            return
        
        hist_info, prev_state, info = self._future[-1]
        
        if hist_info not in self._options.keys():
            print "Whoops, don't know what to do with:"+hist_info
            return
     
        if hist_info == 'skip':
            print "!! skipping over history"
            return

        print "REDO: "+hist_info+" : "+str(info)
        redo_func = self._options[hist_info][1]

        try: current_state = redo_func(prev_state)
        except:
            print "Error with redo function for hist_info: "+hist_info+str(info)
            print "Using pdb, type 'cont' to get the error from python"
            pdb.set_trace()
            current_state = redo_func(prev_state)

        if current_state != None:
            self._past.append([hist_info,current_state,info])
            if len(self._future) > 0:
                del self._future[-1]
       
    def display (self):
        print ""
        print "======================================="
        print "="*20+format("PAST",'^9')+"="*20
        for i in range(len(self._past)): print i,self._past[i][0],self._past[i][2]

        print ""
        print "="*20+format("FUTURE",'^9')+"="*20
        ran = range(len(self._future))
        ran.reverse()
        for i in ran: print i,self._future[i][0],self._future[i][2]

    def save_history (self,filename,clobber=True):
        filename = str(filename)
        if os.path.exists(filename) and not clobber: raise IOError("File exists: "+filename)
        saved_obj = [self.max_history,self.start_deleting_at,self._past,self._future]
        pickle.dump(saved_obj,open(filename,'wb'))

    def open_history (self,filename):
        if not os.path.exists(filename): raise IOError("File does not exist:"+filename)
        if len(self._options) == 0: raise IOError("Options must be set first so that the functions can be input")

        saved_obj = pickle.load(open(filename,'rb'))
        # !! check the save _ obj
        self.max_history = saved_obj[0]
        self.start_deleting_at = saved_obj[1]
        past = saved_obj[2]
        future = saved_obj[3]
        
        def _checkit_ (timetravel):
            t2 = (np.ones(len(timetravel)) < 0)
            for i in range(len(timetravel)):
                hist_info = timetravel[i][0]
                if hist_info not in self._options.keys():
                    print "Error: Don't have functions to deal with :"+hist_info
                    t2[i] = True

            offset = 0
            for i in range(len(t2)):
                if t2[i]:
                    del timetravel[i-offset]
                    offset += 1

            return timetravel

        self._past = _checkit_(past)
        self._future = _checkit_(future)
        
    def post_open_history (self,class_list):
        class_list = list(class_list)
        
        class_convert = {}
        for cclass in class_list:
            class_convert[str(cclass.__class__)] = cclass
        
    def pre_pickle (self):
        tmp_options = deepcopy(self._options)
        del self._options
        return tmp_options

    def post_pickle (self,options):
        self._options = options

pass
################################################################################

class EventConnections:
    """ this bundles some connection abilities for eyeSpec editor classes """

# #     - 'button_press_event'
# #     - 'button_release_event'
# #     - 'draw_event'
# #     - 'key_press_event'
# #     - 'key_release_event'
# #     - 'motion_notify_event'
# #     - 'pick_event'
# #     - 'resize_event'
# #     - 'scroll_event'
# #     - 'figure_enter_event',
# #     - 'figure_leave_event',
# #     - 'axes_enter_event',
# #     - 'axes_leave_event'
    def __init__ (self,ax,connections):
        self.ax = ax
        #if connections is None: self.init_basic_connections()
        self.connections = connections
        self.is_connected = False
        self.store_connections = {}

    #================================================#

    def init_connection_callbacks (self,parent):
        """
        This will look for the standard events with the name replaced as callback and use that for a connection 

        """
        possible_events = ['button_press_event', # => button_press_callback
                           'button_release_event', 
                           'draw_event',
                           'key_press_event',      
                           'key_release_event',    
                           'motion_notify_event',  
                           'pick_event',           
                           'resize_event',         
                           'scroll_event',         
                           'figure_enter_event',  
                           'figure_leave_event',  
                           'axes_enter_event',    
                           'axes_leave_event']

        _func = None
        self.connections = {}
        for con in possible_events:
            con_call = con.replace('event','callback')
            # if you find the name as a method in the parent then 
            # create connection entry
            if con_call in dir(parent):
                exec('func = parent.'+con_call)
                self.connections[con] = [None,_func] # self.connections['key_press_event'] = [None,self.key_press_callback]

    def _do_connect (self,con_event):
        """ do the work of connecting a particular con_event for all the self.connection.keys()"""
        _callbacks = self.ax.figure.canvas.callbacks
        if con_event in self.connections.keys():
            cid = self.connections[con_event][0]
            if cid is None:
                fxn = self.connections[con_event][1]
                new_cid = _callbacks.connect(con_event,fxn)
                self.connections[con_event][0] = new_cid
            
            else:
                fxn = self.connections[con_event][1]
                bmproxy = _callbacks.BoundMethodProxy(fxn)
                _callbacks.callbacks[con_event][cid] = bmproxy
        else: print "Connection Event:",con_event,"not in current connections" # !!

    def connect (self,con_event=None):
        """ for connecting event callbacks to the current plot """
        self.is_connected = True
        for con_event in self.connections.keys():
            self._do_connect (con_event)

    def _do_disconnect (self,con_event):
        """ do the work of disconnecting a particular con_event for all the self.connection.keys()"""
        _callbacks = self.ax.figure.canvas.callbacks
        if con_event in self.connections.keys():
            cid = self.connections[con_event][0]
            if cid is not None:
                if con_event not in _callbacks.callbacks.keys(): print "whoops, there's no key:",con_event
                elif cid in _callbacks.callbacks[con_event].keys():
                    del _callbacks.callbacks[con_event][cid]
                else: pass #print "whoops, cid does not appear:",cid # happens if you try to disconnect twice

        else: print "Connection Event:",con_event,"not in current connections" # !!

    #================================================#

    def disconnect(self):
        """ for disconnecting event callbacks to the current plot """
        self.is_connected = False
        for con_event in self.connections.keys():
            self._do_disconnect(con_event)

    def toggle_connections (self):
        if self.is_connected: self.disconnect()
        else: self.connect()        

class EventBindings:
    
    def __init__ (self,parent,bindings):
        self.parent = parent
        self.bindings = bindings
#        self.bindings['mouse_any'] = [None,wx.EVT_MOUSE_EVENTS,self.mouse_any_callback]

    def bind (self):
        for key in self.bindings.keys():
            bind = self.bindings[key]
            store_out = self.parent.Bind(bind[1],bind[2])
            self.bindings[key][0] = store_out

class KeyboardConfiguration:
    """
    This helps define how keys are mapped to capabilities in the program

    """
    def __init__ (self):
        self._info = ''
        self._short_info = {}
        self._long_info = {}
        self._display_order = []
        
        self._convert_mpl_key_2_key = {'\t':'Tab'}
        self._init_convert_wx_code()

        # !! I COULD CREATE A REMAPPER
        # a function which takes a dictionary for remapping
        # orig_key => new_key
        self._remap = {}
        self._remapped = False
        

#    def add_key_cfg (self,key_cfg):
#        if key_cfg.__class__.__name__ != 'KeyboardConfiguration': raise StandardError("keyboard configuration is not of the correct class")
#        new_keys = key_cfg._short_info.keys()

    def _init_convert_wx_code (self):
        cwc2mk = {8:'delete',9:'\t',13:'return',32:' ',27:'escape',
                  39:"'",44:',',46:'.',
                  49:'1', 50:'2', 51:'3', 52:'4', 53:'5', 
                  54:'6', 55:'7', 56:'8', 57:'9', 58:'0',
                  59:';',
                  65:'a', 66:'b', 67:'c', 68:'d', 69:'e', 
                  70:'f', 71:'g', 72:'h', 73:'i', 74:'j', 75:'k', 
                  76:'l', 77:'m', 78:'n', 79:'o', 80:'p', 81:'q', 
                  82:'r', 83:'s', 84:'t', 85:'u', 86:'v', 87:'w',
                  88:'x', 89:'y', 90:'z',
                  # various other keys
                  91:'[',92:']',93:'\\',
                  96:'`',306:'shift',307:'alt',308:'control'}
        self._convert_wx_code_2_mpl_key = cwc2mk 

    def convert_wx_code_2_key (self,code):
        if code in self._convert_wx_code_2_mpl_key.keys():
            key = deepcopy(self._convert_wx_code_2_mpl_key[code])
            # !! ?? something to translate '\\' to '\'
            return key
        else: return None

    def check_display (self):
        for key in self._short_info.keys():
            if key not in self._display_order:
                print "Couldn't find key:",key," to display"

    def set_configure_info (self,info):
        self._info = str(info)
    
    def get_configure_info (self):
        return deepcopy(self._info)

    def add_mpl_key_convert (self,from_mpl_key,to_key):
        from_key = str(from_mpl_key)
        to_key = str(to_key)
        self._convert_mpl_key_2_key[from_key] = to_key

    def combine_key_cfg (self,key_cfg,insert_mode='protect'):
        # append keys
        for key in key_cfg:
            self.add_key(key,
                         key_cfg._short_info[key],
                         key_cfg._long_info[key],
                         insert_mode = insert_mode)

    def change_key (self,key,short_info,long_info=None,value=None,add_key_if_missing=True):
        if key not in self._short_info: 
            if add_key_if_missing: pass
            else: return
        self.add_key(key,short_info,long_info=long_info,value=value,insert_mode='replace')
        
    def add_key (self,key,short_info,long_info=None,insert_mode='append'):
        """
        INPUTS :
        insert_mode : ('protect','replace','append')

        """
        insert_mode = str(insert_mode)
        if insert_mode not in ['protect','replace','append']: insert_mode = 'append'


        if key == 'display': return
        key = str(key)
        short_info = str(short_info)
        if key in self._short_info.keys():
            if insert_mode == 'protect': return
            if insert_mode == 'append':
                print "HeadsUp: key is already used:"+key
            return

        self._short_info[key] = short_info

        if long_info is not None:
            long_info = str(long_info)
            self._long_info[key] = long_info
        else:
            self._long_info[key] = 'No long info display'

    def load_long_info (self,fname):
        f = open(str(fname),'r')
        # !! look through for specific formatting and then import
        # the longer written help info for each key
        f.close()

    def check_key_press_callback (self,key_to_check):
        key = key_to_check
        if key_to_check not in self._short_info.keys():
            if key_to_check in self._convert_mpl_key_2_key.keys():
                key = self._convert_mpl_key_2_key[key_to_check]
        return (key in self._short_info.keys())

    def keys (self):
        return deepcopy(self._short_info.keys())
    
    def set_display_order (self,key_list=None):
        if key_list is None: key_list = self._short_info.keys()
        else: key_list = list(key_list)
        self._display_order = []
        for key in key_list:
            if key not in self._short_info.keys():
                print "Can't Display Key:",key
            else:
                self._display_order.append(key)

    def get_display_order (self):
        return deepcopy(self._display_order)

    def _get_key_display (self,key,key_dict,max_chr=80):

        output = []
        info = key_dict[key]
        split_info = info.split()
        
        line = format(str(key),'>10')+'   : '
        length_2line_tab = len(line)

        for i in range(len(split_info)):
            word = split_info[i]+" "
            if len(line+word) > max_chr:
                output.append(line)
                line = ' '*length_2line_tab+word
            else:
                if i == len(split_info)-1:output.append(line+word)
                else: line += word

        return "\n".join(output)

    def display_short_info (self,display=True):
        max_number_characters = 80
        out_display = []
        out_display.append(format('key','>10')+'   : Description')
        out_display.append(' '+'-'*79)

        for key in self._display_order:
            out_display.append(self._get_key_display(key,self._short_info,max_number_characters))
            # out_display.append('  '+'-'*78)

        if display:
            for line in out_display:
                print(line)
        else: return iter("\n".join(out_display))
        
    def display_long_info (self,key,display=True):
        max_number_characters = 80
        if key not in self._long_info.keys(): return
        output = self._get_key_display(key,self._long_info,max_number_characters)
        full_output = [format('key','>10'),' : Description\n',' '*'-'*79,output]
        line_full_output = '\n'.join(full_output).strip('\n')
        if display:
            print(line_full_output)
        else: return line_full_output

class Cursor (EventConnections):
    """
    Adds horizontal and verical lines to the the cursor position, and has some event options
    INPUTS:
    ax : which axes to attach the Cursor to 
    !! I could probably be smarter with this
    """
    def __init__(self,ax,toolbar=None):
        self.ax = ax
        self.zorder = 30
        self.lx = ax.axhline(color='k',zorder=self.zorder)  # the horiz line
        self.ly = ax.axvline(color='k',zorder=self.zorder)  # the vert line

        self.timer = Timer(0.01)
        
        self._num_scrolls = 0
        self._max_num_scrolls = 3 #outside of axes
        self.percent_to_scroll = 2

        self.recorded_limit = 7000 
        self.recorded_pts = [] # !! if one instance of this is attached to different axes and clicked on each how will this record_pts behave?
        # text location in axes coords
        self.txt = ax.text( 0.7, 0.9, '', transform=ax.transAxes,visible=False,zorder=self.zorder)
        
        # configure options for linking click and keyboard
        # self.keyboard_cfg = {} 
        # would need to bind, have a something record when key_press_event happens, then does button_press_event happen, then if key_release_event happens you can execute something

        # create a better box drawing tool
        self._better_zoom = True


        # !! I could make this better: self.ax.figure.canvas.toolbar
        if type(toolbar).__name__ not in ['NavigationToolbar2Wx']: toolbar = None
        self.toolbar= toolbar

        self._dragging = False
        self._clicked = False


        if self._better_zoom and self.toolbar is not None:
            self._selection_box = self.ax.add_patch(plt.Rectangle((0.,0.),.1,.1,facecolor='none',edgecolor='k',lw=.3,alpha=1,visible=False,zorder=100))


        #------------------------------------------------#
        self.init_connection_callbacks(self)
        
    def button_press_callback (self,event):  
        if not event.inaxes: return
        if 'button' in dir(event):
            x, y = event.xdata, event.ydata
            self.recorded_pts.append([x,y])

        if len(self.recorded_pts) > self.recorded_limit: del self.recorded_pts[0]

        self._dragging = False
        self._clicked = True
        self.update()
        
    def motion_notify_callback (self, event):
        if not event.inaxes: return
        if self.timer.check(reset=False): return
        x, y = event.xdata, event.ydata

        

        if self._clicked and self.toolbar is not None and self._better_zoom and self.toolbar.mode == 'zoom rect':
            start_corner = self.recorded_pts[-1] # x,y

            xc1,yc1 = start_corner[0], start_corner[1]
            xc2, yc2 = event.xdata, event.ydata

            self._selection_limits = np.array([xc1,xc2,yc1,yc2])
            self._update_selection_box (self._selection_limits) # !! the box must at least be visible, maybe...
            self.update()

        self._dragging = True
        # update the line positions
        self.lx.set_ydata(y)
        self.ly.set_xdata(x)
        self.txt.set_text('x=%1.2f, y=%1.2f'%(x,y))
        self.update()

    def button_release_callback (self,event): #@UnusedVariable
        if self.toolbar is not None and self._better_zoom: self._update_selection_box(None)
        self._dragging=False
        self._clicked = False
        self.update()

    def hide (self):
        self.disconnect()
        self.set_visible(False)
        
    def show (self):
        self.connect()
        self.set_visible(True)

    def set_visible (self,truth):
        self.lx.set_visible(truth)
        self.ly.set_visible(truth)
            
    def update (self):
        if self.timer.check():
            self.ax.figure.canvas.draw()
            self._last_update = time.time()

    def _update_selection_box (self,selection_limits):
        """
        uses selection_limits to define the box and which data points of self._selected_data are within
        """
        # if all the limits are the same turn everything off
        if selection_limits is None:
            self._selection_box.set_visible(False)
            return 

        # create the selection box
        xc = selection_limits[0] # xmin
        yc = selection_limits[2] # ymin
        wid = selection_limits[1] - xc 
        hei = selection_limits[3] - yc

        self._selection_box.set_bounds(xc,yc,wid,hei)
        self._selection_box.set_visible(True)

    #    if self._better_box:
    #        pass

    #def _def_bbox (self,X,Y):
    #    if '_bbox' not in dir(self): self._bbox = self.ax.plot(

    def cursor_scroll (self,event):
        #!!! this currently doesn't work!!!!
        # !! currently doesn't work, search: curse_scroll
                    
        if True: return
        if event.inaxes:
            self._num_scrolls = 0
            return
        self._num_scrolls += 1
        if self._num_scrolls > self._max_num_scrolls: return
        xran = self.ax.axis()[1] - self.ax.axis()[0]
        yran = self.ax.axis()[3] - self.ax.axis()[2]

        new_xmin,new_xmax,new_ymin,new_ymax = self.ax.axis()
        del_x = self.percent_to_scroll/100.*xran
        del_y = self.percent_to_scroll/100.*yran

        if event.xdata <= self.ax.axis()[0]:
            new_xmin = self.ax.axis()[0]-del_x
            new_xmax = self.ax.axis()[1]-del_x
        if event.xdata >= self.ax.axis()[1]:
            new_xmin = self.ax.axis()[0]+del_x
            new_xmax = self.ax.axis()[1]+del_x
        if event.xdata <= self.ax.axis()[2]:
            new_ymin = self.ax.axis()[2]-del_y
            new_ymax = self.ax.axis()[3]-del_y
        if event.xdata >= self.ax.axis()[3]:
            new_ymin = self.ax.axis()[2]+del_y
            new_ymax = self.ax.axis()[3]+del_y

        self.ax.axis((new_xmin,new_xmax,new_ymin,new_ymax))
        self.ax.figure.canvas.draw()

pass
################################################################################

ID_BEGIN=100
wxStdOut, EVT_STDDOUT= NewEvent()
wxWorkerDone, EVT_WORKER_DONE= NewEvent()

class Worker(threading.Thread):
    requestID = 0
    def __init__(self, parent, requestQ, resultQ, **kwds):
        threading.Thread.__init__(self, **kwds)
        self.setDaemon(True)
        self.requestQ = requestQ
        self.resultQ = resultQ
        self.start()

    def beginTest(self, call_able, *args, **kwds):
        Worker.requestID +=1
        self.requestQ.put((Worker.requestID, call_able, args, kwds))
        return Worker.requestID

    def run(self):
        while True:
            requestID, call_able, args, kwds = self.requestQ.get()
            self.resultQ.put((requestID, call_able(*args, **kwds)))
            evt = wxWorkerDone()
            wx.PostEvent(wx.GetApp().frame, evt)

class SysOutListener:
    def write(self, string):
        sys.__stdout__.write(string)
        evt = wxStdOut(text=string)
        wx.PostEvent(wx.GetApp().frame.panel.redirect_text_panel, evt)
           
class eyeSpecTextRedirectPanel (wx.Panel):
    def __init__ (self,parent):
        wx.Panel.__init__(self, parent)
        self.requestQ = Queue.Queue() #create queues
        self.resultQ = Queue.Queue()

        self.SetBackgroundColour("#99CCB2")

        self.output_window = wx.TextCtrl(self,-1, style=wx.TE_AUTO_SCROLL|wx.TE_MULTILINE|wx.TE_READONLY)
        self.output_window.SetBackgroundColour("#99CCB2")
        self.output_window.SetForegroundColour("#312500")

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.output_window,10,wx.EXPAND)
        self.SetSizer(sizer)
        self.Fit()
        #events
        
        wx.EVT_BUTTON(self, ID_BEGIN, self.OnBeginTest)
        self.output_window.Bind(EVT_STDDOUT, self.OnUpdateOutputWindow)
        self.output_window.Bind(wx.EVT_TIMER, self.OnProcessPendingOutputWindowEvents)
        self.Bind(EVT_WORKER_DONE, self.OnWorkerDone)

        #thread
        self.worker = Worker(self, self.requestQ, self.resultQ)

    def OnUpdateOutputWindow(self, event):
        value = event.text
        self.output_window.AppendText(value)

    def OnBeginTest(self, event): #@UnusedParameter
        lines_of_output=7
        self.go.Disable()
        self.worker.beginTest(LongRunningProcess, lines_of_output) #@UndefinedVariable
        self.output_window_timer.Start(50)

    def OnWorkerDone(self, event): #@UnusedParameter
        self.output_window_timer.Stop()
        self.go.Enable()

    def OnProcessPendingOutputWindowEvents(self, event): #@UnusedParameter
        self.output_window.ProcessPendingEvents()

pass
################################################################################
    
class OrderSelectionBox:
    """ Graphical selection box around the data """
    def __init__ (self, ax, visible=False):
        self._ord_i = None
        self.ax = ax
        xc1, xc2 = [0, 1]
        X_sel = [xc1, xc1, xc2, xc2, xc1]
        Y_sel = self._get_ybounds_axis()
        
        self.order_box, = self.ax.plot(X_sel, Y_sel, lw=4, color='#66FF99', alpha=0.45, zorder= -20, linestyle='--', visible=visible)  # left,top,right,bottom

    def set_color (self, color):
        self.order_box.set_color(color)

    def set_zorder (self, zorder):
        zorder = int(zorder)
        self.order_box.set_zorder(zorder)
        
    def get_zorder (self):
        return self.order_box.get_zorder()

    def set_xbounds (self, xbounds):
        xc1, xc2 = xbounds
        X_sel = [xc1, xc1, xc2, xc2, xc1]
        self.order_box.set_xdata(X_sel)        

    def _get_ybounds_axis (self):
        ymin, ymax = self.ax.axis()[2], self.ax.axis()[3]
        if ymin == ymax or ymin is None:
            return
        yran = ymax - ymin 
        yc1, yc2 = self.ax.axis()[2] + .01 * yran, self.ax.axis()[3] - .01 * yran
        Y_sel = [yc1, yc2, yc2, yc1, yc1]
        return Y_sel
    
    def update_ybounds (self):
        Y_sel = self._get_ybounds_axis()
        self.order_box.set_ydata(Y_sel)

    def update_box (self, pltdata, index=None):
        if pltdata is None or index is None:
            self.set_visible(False)
            return
        xbounds = (np.min(pltdata.get_xdata()), np.max(pltdata.get_xdata()))
        self.set_xbounds(xbounds)
        self.update_ybounds()
        self.set_order_index(index)
        self.set_color(pltdata.get_color())
        self.set_visible(True)

    def set_order_index (self, index):
        self._ord_i = index
        
    def get_order_index (self):
        return deepcopy(self._ord_i)

    def set_visible (self, truth):
        self.order_box.set_visible(truth)
               
    def get_visible (self):
        return self.order_box.get_visible()

    def get_bounds (self):
        xdata = self.order_box.get_xdata()
        ydata = self.order_box.get_ydata()
        return [xdata[0], xdata[2], ydata[0], ydata[1]]

class DataSelectionBox:
    """ Graphical box for selecting data"""
    
    def __init__ (self, ax, zorder=100):
        self.zorder = zorder
        self.ax = ax
        self.box = self.ax.add_patch(plt.Rectangle((0., 0.), .1, .1, facecolor='r', alpha=.3, visible=False, zorder=self.zorder))

        self.selection_limits = np.zeros(4)
        
        self._allow_selection_box = True
        self.clear_box_variables()

        self._ptime = time.time()
        self._dtime = 0.08

    def __contains__ (self, xy):
        if not self.get_visible():return False
        xpt, ypt = xy
        bb = self.get_bounds()
        check = ((xpt >= bb[0]) and (xpt <= bb[1]) and (ypt >= bb[2]) and (ypt <= bb[3]))
        return check
      
    def display_box_variables (self):
        print "="*30
        print "visible:", self.get_visible()
        print 'create box:', self._create_box
        print 'create box left:', self._create_box_from_left
        print 'create box right:', self._create_box_from_right
        print 'move box:', self._move_box
        print 'start limits:', self._box_start_lims 
        
    def clear_box_variables (self):
        self._create_box = False
        self._create_box_from_left = False
        self._create_box_from_right = False
        self._move_box = False
        self._box_start_lims = np.zeros(4)

    def set_box_start (self, start_lims=None, and_update=False):
        if start_lims is None: start_lims = self.get_bounds()
        x1, x2, y1, y2 = start_lims
        xc1 = min(x1, x2)
        xc2 = max(x1, x2)
        yc1 = min(y1, y2)
        yc2 = max(y1, y2)
        self._box_start_lims = [xc1, xc2, yc1, yc2]     
        if and_update:
            self.update_box(None, False)
        
    def set_box_in_window (self, start_lims):
        sl = start_lims
        xmin, xmax, ymin, ymax = self.ax.axis()
        if (sl[0] <= xmax) and (sl[1] >= xmin) and (sl[2] <= ymax) and (sl[3] >= ymin):
            self._box_start_lims = start_lims  
        else:
            self.set_visible(False)
            self.clear_box_variables()
                    
    def get_box_start (self):
        return deepcopy(self._box_start_lims)

    def reset_box (self):
        limits = np.zeros(4)
        self.update_box(limits, True)

    def move_box (self, tocreate=None):
        if tocreate is None: 
            return deepcopy(self._move_box)
        self.clear_box_variables()
        if not self.get_visible(): return
        self._move_box = True

    def create (self, tocreate=None):
        if tocreate is None:          
            return deepcopy(self._create_box)  
        self.clear_box_variables()
        self._create_box = True

    def create_from_left (self, tocreate=None):
        if tocreate is None: 
            return deepcopy(self._create_box_from_left)
        self.clear_box_variables()
        self._create_box_from_left = True
        self._create_box_from_right = False
        
    def create_from_right (self, tocreate=None):
        if tocreate is None:
            return deepcopy(self._create_box_from_right)
        self.clear_box_variables()
        self._create_box_from_left = False
        self._create_box_from_right = True        

    def modify_bounds (self, data_bounds=None, scale=0.005):
        if data_bounds is None: return (0, 0, 0, 0)
        xmin, xmax, ymin, ymax = data_bounds
        yran = ymax - ymin      
        xran = xmax - xmin
        scale = float(scale)
        xs1 = xmin - scale * xran
        xs2 = xmax + scale * xran
        ys1 = ymin - scale * yran
        ys2 = ymax + scale * yran  
        return (xs1, xs2, ys1, ys2)                    

    def during_btn_press (self, event, epsilon_x, data_bounds=None):
        """ Selection Box
        Returns True if plot needs to be updated
        
        if data_bounds is None then you can't select from the edge of orders
        only create free box and move else data_bounds = [xmin,xmax,ymin,ymax]
        
        """
        if event.inaxes is None: return
        ib = (data_bounds is not None)
        
        if ib:
            xmin, xmax, _ymin, _ymax = data_bounds
        
        # if the selection box is visible
        if self.get_visible():
            if (event.xdata, event.ydata) in self:
                self.move_box(True)
                self.set_box_start()
                return False
        

        # selection box is not visible and order is selected
        if event.inaxes is None: return False
        if ib:
            xs1, xs2, ys1, ys2 = self.modify_bounds(data_bounds)

        # if you click within epsilon of the left edge
        if np.abs(xmin - event.xdata) < 1.5 * epsilon_x and ib:
            self.create_from_left(True)
            self.set_box_start([xs1, xmin, ys1, ys2])
        # if you click within epsilon of the right edge
        elif np.abs(xmax - event.xdata) < 1.5 * epsilon_x and ib:
            self.create_from_right(True)
            self.set_box_start([xmax, xs2, ys1, ys2])

        else:
            self.create(True)

        self.update_box(None,False)
        self._ptime = time.time()
        return True

    def during_mot_notify (self, event, start_pt, data_bounds=None, pressed_dict=None):
        """ Selection Box
        Returns True if plot needs to be updated
        
        start_pt = the event.xdata and event.ydata recorded when button press occurred
        data_bounds = [xmin,xmax,ymin,ymax] if None then you can't do the button press region select
        """
        db = True
        _xmin, _xmax, ymin, ymax = np.ones(4)
        try: _xmin, _xmax, ymin, ymax = data_bounds
        except: db = False
        can_press = (pressed_dict is not None)
        
        xs1, xs2, ys1, ys2 = self.get_box_start()
        xc1, xc2, yc1, yc2 = deepcopy((xs1, xs2, ys1, ys2))
        
        def _get_10_rounded (val):
            # this will round the number to the nearest tens place in increasing direction
            # e.g.  9960 => 10000,   -340 => -1000,  0.0003 => 0.001
            sign = 1
            if val < 0.0: sign = -1.0
            if val == 0.0: return -0.1            
            rval = round(math.log10(abs(val)))
            rval += 1.0
            return sign * 10.0 ** (rval)
        
        # pdb.set_trace()
        # self.display_box_variables()
        if self.move_box():
            xc1 = xs1 + (event.xdata - start_pt[0])
            xc2 = xs2 + (event.xdata - start_pt[0])
            yc1 = ys1 + (event.ydata - start_pt[1])
            yc2 = ys2 + (event.ydata - start_pt[1])
        elif self.create_from_left():
            xc2 = xs1 + np.abs(event.xdata - start_pt[0])
        elif self.create_from_right(): 
            xc1 = xs2 - np.abs(event.xdata - start_pt[0])
        elif self.create():
            xc1, yc1 = deepcopy(start_pt)
            xc2, yc2 = deepcopy((event.xdata, event.ydata))
            
            if (can_press and db and 'x' in pressed_dict) and pressed_dict['x']:
                yc1 = _get_10_rounded(ymin)
                yc2 = _get_10_rounded(ymax)
        else: return False
        
        self.set_visible(True)
        if (time.time()-self._ptime) > self._dtime:
            self.update_box([xc1, xc2, yc1, yc2]) 
            self._ptime      
        return True
                  
    def during_btn_click (self, event):
        """ Selection Box
        Returns True if plot needs to be updated
        """
        if not self.get_visible(): return False
        
        if self.move_box(): return True

        if (event.xdata, event.ydata) not in self:          
            self.set_visible(False)
            return True 
        return False
            
    def during_key_release (self, event, deltax=1.0):
        """ Selection Box
        Returns True if plot needs to be updated
        """   
        if not (self.create_from_left() or self.create_from_right()): return False
    
        xb1, xb2, _yb1, _yb2 = self.get_bounds()
        xc1, xc2, yc1, yc2 = self._box_start_lims
        
        if self.create_from_left():
            if event.key == 'left':  xb2 = np.max((xc2, xb2 - deltax))
            if event.key == 'right': xb2 = np.max((xc2, xb2 + deltax))
        elif self.create_from_right():
            if event.key == 'left':  xb1 = np.min((xc1, xb1 - deltax))
            if event.key == 'right': xb1 = np.min((xc1, xb1 + deltax))           
        
        self.update_box([xb1, xb2, yc1, yc2], True)        
             
    def set_zorder (self, zorder):
        zorder = int(zorder)
        self.box.set_zorder(zorder)
        
    def get_zorder (self):
        return self.box.get_zorder()
        
    def _check_in_window_visible (self):
        sl = self.get_bounds()
        xmin, xmax, ymin, ymax = self.ax.axis()
        if ((sl[0] > xmax) or (sl[1] < xmin)) and ((sl[3] < ymin) or (sl[1] > ymax)):        
            return False
        else: return True
        
    def set_visible (self, truth):
        truth = bool(truth)  
        self.box.set_visible(truth)
             
    def get_window_visible (self):         
        if not self._check_in_window_visible(): self.set_visible(False)
        return self.box.get_visible()        
        
    def get_visible (self):
        return self.box.get_visible()
        
    def get_bounds (self):
        get_bounds = [self.box.get_bbox().xmin, self.box.get_bbox().xmax,
                      self.box.get_bbox().ymin, self.box.get_bbox().ymax]
        return get_bounds
    
    def selection_box_focus (self, focus):
        self.clear_box_variables()
        if focus.find('edge') == -1: return
        
        if focus == 'edge_right':
            self._create_box = False
            self._create_box_from_left = False
            self._create_box_from_right = True
            
        elif focus == 'edge_left':
            self._create_box = False
            self._create_box_from_left = True
            self._create_box_from_right = False             

    def update_box (self, selection_limits=None, visible=True): 
        # if all the limits are the same turn everything off
        if selection_limits is None:
            selection_limits = self._box_start_lims

        # create the selection box
        xc = selection_limits[0]  # xmin
        yc = selection_limits[2]  # ymin
        wid = selection_limits[1] - xc 
        hei = selection_limits[3] - yc

        self.box.set_bounds(xc, yc, wid, hei)
        self.box.set_visible(visible)

class LineSelection:
    """ Line Selection """
    def __init__ (self,ax,wl0=0.0,color='y',alpha=0.3,visible=False):
        self.ax = ax
        self._linei = None
        self._prev_linei = 0
        ymax,ymin = self.ax.axis()[2:]
        self.highlight_line, = ax.plot([wl0,wl0], [ymax,ymin], color=color, alpha=0.65, lw=4, visible=visible, alpha=alpha, zorder=-50)

    def reset (self):
        self._linei = None
        self._prev_linei = 0
          
    def get_current_state (self):
        curstate = State('line_selection')
        curstate['_linei'] = self._linei
        curstate['_prev_linei'] = self._prev_linei
        curstate['xpt'] = self.get_xdata()
        curstate['vis'] = self.is_selected()
        
        return curstate
        
    def set_state (self,prevstate):
        prevstate.check_id('line_selection')
        
        self._linei = prevstate['_linei']
        self._prev_linei = prevstate['_prev_linei']
        self.set_xdata(prevstate['xpt'])
        self.highlight_line.set_visible(prevstate['vis'])
    
    def set_ydata (self,ybounds):
        ymin,ymax = ybounds
        self.highlight_line.set_ydata([ymin,ymax])
    
    def select_line (self,index,xpt,always_find=False):
        if index is None: 
            if always_find: self._linei = deepcopy(self._prev_linei)
            else: 
                self.deselect()
                return
        
        if index != self._linei: self._prev_linei = deepcopy(self._linei)
        
        self.set_xdata(xpt)
        self._linei = index
        self.highlight_line.set_visible(True)
        
    def deselect_line (self):
        if self._linei is not None: self._prev_linei = deepcopy(self._linei)
        self._linei = None
        self.highlight_line.set_visible(False)
        
    def is_selected (self):
        return self.highlight_line.get_visible()
    
    def get_line_index (self):
        return deepcopy(self._linei)
    
    def get_prev_line_index (self):
        return deepcopy(self._prev_linei)
    
    def get_xdata (self):
        return self.highlight_line.get_xdata()[0]

    def set_xdata (self,xpt):
        if xpt is None: raise ValueError("xpt is None")
        self.highlight_line.set_xdata([xpt,xpt])

pass
################################################################################

class PlotLines:
    def __init__ (self, ax, line_data, line_info=None, line_info_vis=False, crop_xbounds=None, lock_order=False):
        self.ax = ax
                
        def dropping_function (index):
            i = index % 3 # number of changes
            return 0.1 + 0.05*i 
        # !! or I could create a dictionary of the line drops for when you add/remove lines       
        self._dropby = dropping_function
        self._overshoot_y = 0.5 #%
        self._lock_order = lock_order
          
        self._redshift_info = 0.05
          
        self._prev_bounds = ax.axis()
        # NOTE: quite a few variables are defined in change_linelist
        self.change_linelist(line_data, line_info, crop_xbounds)
       
       
        T = time.gmtime()
        curdate = format(T.tm_mon,'02')+format(T.tm_mday,'02')+str(T.tm_year) 

        self.iochoices = SaveOpenChoices()
        self.iochoices.add('Line List',
                           '(*.dat)|*.dat'\
                           'Text File (*.txt)|*.txt'\
                            ' All Files (*)|*',
                           'Linelist_'+curdate+".txt", 
                           'Line List',
                           self.save_lines,
                           self.open_lines)
        
        self.set_info_visible(line_info_vis)
        
        # !! self.iochoices.add(parameters

    def get_current_state (self):
        curstate = State('plot_lines')
        curstate['_overshoot_y'] = self._overshoot_y
        curstate['lines'] = self.lines
        curstate['line_info'] = self.line_info
        curstate['lock_order'] = self._lock_order
        return curstate
    
    def apply_state (self,prevstate):
        prevstate.check_id('plot_lines')
        
        self._overshoot_y = prevstate['_overshoot_y']
        lines = prevstate['lines']
        info = prevstate['info']
        self.change_linelist(lines, info)
        self._lock_order = prevstate['lock_order']
        
    def __len__ (self):
        return len(self.lines)
                       
    def __contains__ (self,wl):
        return np.any(self.lines == wl)
    
    def contains (self,wl,tol=0):
        return np.argwhere(np.abs(self.lines-wl) <= tol)
     
    def get_bounds (self):
        xmin = np.min(self.lines)
        xmax = np.max(self.lines)
        return (xmin,xmax)
     
    def change_linelist (self,line_data,line_info = None,crop_xbounds=None,color='#990033',zorder=0,alpha=.8):
        self._has_info = (line_info is not None)
        if self._has_info:
            assert len(line_data) == len(line_info)  
                 
        if 'lines' in dir(self): self.hide()
                    
        # These should match exactly
        line_data = np.array(line_data,dtype=float)
        self.lines = np.sort(line_data)
        self.plines = []
        self._plines_visible = np.array([],dtype=bool)
        
        self.line_info = []
        self._info_visible = np.array([],dtype=bool)
        
        if crop_xbounds is not None:
            xc1,xc2 = crop_xbounds
            xmask = (xc1 < self.plines)*(self.plines < xc2)        
            self.plines = self.plines[xmask]          
 
        xmin,xmax,ymin,ymax = self.ax.axis()
        yran = ymax - ymin
        lymin = ymin - self._overshoot_y*yran
        lymax = ymax + self._overshoot_y*yran

        for i in np.argsort(line_data):
            wl = line_data[i]
            if not self._has_info: linfo = None
            else: linfo = line_info[i]
            self._add_pline(wl, lymin, lymax,[xmin,xmax], color, zorder, alpha)
            if self._has_info: self._add_text(i,wl,linfo,[zorder,alpha,color])         
                        
    def change_line (self,index,line_data):
        # !! checks, is index in range, is line_data a float
        wl = line_data
        self.lines[index] = wl
        self.plines[index].set_xdata([wl,wl])
         
    def change_info (self,index,line_info=None):
        wl = self.lines[index]
        linfo = line_info
        if line_info is None:
            line_info = self.line_info[index].get_text()
            linfo = line_info.split("\n")[0]
        
        if not self._has_info: return
        xpt = wl + self._redshift_info
        self.line_info[index].set_x(xpt)        
        self.line_info[index].set_text(self._get_text_to_add(linfo, wl))
    
    def _add_pline (self,wl,lymin=0,lymax=1,xbounds=None,color=None,zorder=None,alpha=None,visible=None):
        if color is None: color = self.plines[0].get_color()
        if zorder is None: zorder = self.get_zorder()
        if alpha is None: alpha = self.plines[0].get_alpha()
        if visible is None:
            if xbounds is None: xmin,xmax = self.ax.axis()[:2]
            else: xmin,xmax = xbounds
            
            visible = False
            if xmin < wl < xmax: visible = True
        
        plotline, = self.ax.plot([wl,wl],[lymin,lymax],lw=1.5,color=color,zorder=zorder,alpha=alpha,visible=visible)
        self.plines.append(plotline)

        self._plines_visible = np.append(self._plines_visible,visible)
        
    def _get_text_to_add (self,linfo,wl):
        return format(linfo,'<10')+"\n"+format(wl,'<10.3f') 
      
    def _add_text (self,index,wl,linfo,mpl_kwarg=None):
        if not self._has_info: return
        if mpl_kwarg is not None: zorder,alpha,color = mpl_kwarg
        else:
            zorder = self.plines[0].get_zorder()
            alpha = self.plines[0].get_alpha()
            color = self.plines[0].get_color()
            
        if linfo is None: return
        
        xmin,xmax,ymin,ymax = self.ax.axis()
        yran = ymax - ymin
    
        # !! could have wl only if no line info given
        linfo = self._get_text_to_add(linfo, wl)
        scl = self._dropby(index)
        ypt = ymax - scl*yran
        xpt = wl + self._redshift_info
        visible = False
        if (xmin < xpt < xmax): visible = True
        
        txt = self.ax.text(xpt,ypt,linfo,zorder=zorder,alpha=alpha,color=color,size='small',visible=visible)
        self.line_info.append(txt)
        
        self._info_visible = np.append(self._info_visible,visible)
                   
    def get_zorder (self):
        return self.plines[0].get_zorder()
    
    def set_zorder (self,zorder):
        for i in xrange(len(self.plines)):
            self.plines[i].set_zorder(zorder)
            if self._has_info: self.line_info[i].set_zorder(zorder)    
 
    def get_lock_order (self):
        return deepcopy(self._lock_order)
    
    def set_lock_order (self,lock_order):
        self._lock_order = bool(lock_order)
 
    def _get_visible_in_ranges (self,xbounds):
        # set lines visible
        xmin,xmax = xbounds        
        L = self.lines
        vis_ran = np.argwhere((xmin < L)*(L < xmax)).T[0]

        # set lines invisible            
        invis_ran = np.argwhere(self._plines_visible == True).T[0]

        invis_info_ran = None
        if self._has_info:
            invis_info_ran = np.argwhere(self._info_visible == True).T[0]

        # !! could be smarter and figure out which were just set to visible/invisible
        return vis_ran, invis_ran, invis_info_ran
 
    def _set_prev_visible (self,ran_lines,ran_info,truth=False):
        for i in ran_lines: self._set_line_visible(i,truth)
        if self._has_info:
            for i in ran_info: self._set_info_visible(i, truth)
        
    def set_window_visible (self,xbounds,bounds_update=True):
        if self._prev_bounds[:2] == xbounds: return False
        if bounds_update: self._prev_bounds[:2] = xbounds
        
        vis_ran, invis_ran, invis_info_ran = self._get_visible_in_ranges(xbounds)
        
        if len(vis_ran) == 0: return False
        
        #---------------------------------------------------------------------#
        # set invisible
        self._set_prev_visible(invis_ran, invis_info_ran)

        # set these visible:
        for i in vis_ran:
            self._set_line_visible(i,True)
            self._set_info_visible(i,True)

        return True
    
    def hide (self):
        self.set_visible(False)

    def set_visible (self,truth):
        self.set_info_visible(truth)
        self.set_lines_visible(truth)

    def set_lines_visible (self,truth):
        check = truth == np.all(self._plines_visible)
        if check: return

        ran_truth = np.argwhere(self._plines_visible != truth).T[0]
        for i in ran_truth:
            self._set_line_visible(i, truth)

    def set_info_visible (self,truth):
        check = truth == np.all(self._info_visible)
        if check: return
        
        ran_truth = np.argwhere(self._info_visible != truth).T[0]
        for i in ran_truth:
            self._set_info_visible(i, truth)

    def _set_line_visible (self,i,truth):
        lines_vis = deepcopy(truth)        
        if truth == self._plines_visible[i]: return
        
        self._plines_visible[i] = lines_vis
        self.plines[i].set_visible(lines_vis)        

    def _set_info_visible (self,i,truth):
        if not self._has_info: return
        if truth == self._info_visible[i]: return
        
        self._info_visible[i] = truth
        self.line_info[i].set_visible(truth)
                     
    def get_lines_visible (self):
        return deepcopy(self._plines_visible)
    
    def get_info_visible (self):
        if not self._has_info: return False
        return deepcopy(self._info_visible)
 
    def set_info_size (self,text_size):
        if not self._has_info: return
        if not self.get_info_visible(): return
        for i in xrange(len(self.line_info)):
            self.line_info[i].set_size(text_size)
    
    def get_index_by_wl (self,wl):
        return np.argmin(np.abs(self.lines - wl))
    
    def get_line_data (self,index):
        return self.lines[index]

    def get_line_info (self,index):
        if not self._has_info: return ''
        txt=  self.line_info[index].get_text().split("\n")
        return txt[0]
        
    def sort_by_wavelength (self):
        if self._lock_order: return
        sortit = np.argsort(self.lines)
        self.lines = self.lines[sortit]
        
        self._info_visible = self._info_visible[sortit]
        self._plines_visible = self._plines_visible[sortit]
            
        new_plines = range(len(self.plines))
        new_lineinfo = range(len(self.line_info))
        for i in xrange(len(sortit)):
            j = sortit[i]
            new_plines[i] = self.plines[j]
            if self._has_info: new_lineinfo[i] = self.line_info[j]
        
        self.plines = new_plines
        if self._has_info: self.line_info = new_lineinfo
                    
    def add_line (self,line, line_info=None, ybounds=None): 
        self.lines = np.append(self.lines,line)
        
        if ybounds is None:
            ymin,ymax = self.ax.axis()[2:]
            yran = ymax-ymin
            lymin = ymin - yran
            lymax = ymax - yran
        else:
            lymin,lymax = ybounds
        pdb.set_trace()
        self._add_pline(line, lymin, lymax)
        # This is kind of kluge, but when an update is done the code decides
        # which to update based on visible, by setting this to False it will
        # make that line visible when the update is done 
        
        index = len(self.lines)
        if self._has_info:
            if line_info == None: line_info = ' '
            self._add_text(index, line, line_info)
        
        # !! could make it smarter so that it puts the line in the right place
        if np.any(self.lines > line) and not self._lock_order:
            self.sort_by_wavelength()   
        
    def interactive_add_line (self):
        dlg = wx.TextEntryDialog(None, 'Enter line and info:',defaultValue = '6562.28 H Alpha')
        if dlg.ShowModal() == wx.ID_OK:
            add_line = dlg.Value
            sinput = add_line.split()
            try: wl = float(sinput[0])
            except: 
                print "Could not convert given string to floating point value:"+sinput[0]
                return
            info= ' '.join(sinput[1:])
        else: return
        self.add_line(wl,info)
        
    def delete_line_by_wl (self,wl,by_nearest=False):
        if by_nearest:
            index = np.argmin(self.lines-wl)
        else:
            index = np.argwhere(self.lines==wl).T[0]
            if len(index)==0: 
                print "Couldn't find wavelength:",wl
                return
            elif len(index) > 1: 
                print "Found more than one match to wavelength using the first:",len(index)
            index = index[0]
        self.delete_line_by_index(index)
        
    def delete_line_by_index (self,index):
        if len(self.lines)==1: 
            print "Can't delete all lines"
            return
                
        self.lines = np_vstack_delete(self.lines,index)
        self._plines_visible = np_vstack_delete(self._plines_visible,index)

        self.plines[index].set_visible(False)
        del self.plines[index]
        
        if self._has_info: 
            self.line_info[index].set_visible(False)
            del self.line_info[index]
            self._info_visible = np_vstack_delete(self._info_visible,index)
          
    def crop_list_by_bounds (self,xbounds):
        xmin,xmax = xbounds
        xmask = (xmin < self.lines)*(self.lines < xmax)
        
        self.lines = np_vstack_delete(self.lines,xmask)
        self._plines_visible = np_vstack_delete(self._plines_visible,xmask)
        if self._has_info: 
            self._info_visible = np_vstack_delete(self._info_visible,xmask)
        
        for i in np.argwhere(xmask==True):
            self.plines[i].set_visible(False)
            del self.plines[i]
            if not self._has_info: continue
            self.line_info[i].set_visible(False)
            del self.line_info[i]

    def save_lines (self,filename,clobber=True):
        if os.path.exists(filename) and not clobber: raise IOError("File exists:"+filename)
        output_list = zip(self.lines,self.line_info)
                
        lines = []
        lines.append("# Linelist created from eyeSpec "+time.ctime())
        for val in output_list:
            output = [format(val[0],'<15.10f'),format(val[1],"<20")]
            lines.append("  ".join(output))
            
        f = open(filename,'w')
        f.writelines("\n".join(lines))
        f.close()
     
    def open_lines (self,filename,change_list=True):
        f = open(filename,'r')
        flines = f.readlines()
        f.close()
        
        lines_data = []
        lines_info = []
        
        found_info = False
        for line in flines:
            line = line.rstrip()
            sline = line.split()
            if len(sline)==0 or line[0] == '#': continue
            wl = float(sline[0])
            info = ''
            if len(sline)>1:
                found_info = True
                info = " ".join(sline[1:])
            
            lines_data.append(wl)
            lines_info.append(info)
        
        if not found_info: lines_info = None
           
        if change_list: self.change_linelist(lines_data,lines_info) 
        else: return lines_data,lines_info
          
    def update_info (self,ybounds):
        print "!! need to update line update_lines_and_info"

        if not self._has_info: return
        if not self._info_visible: return
        xmin,xmax,ymin,ymax = ybounds
        yran = ymax - ymin
        
        # !! also turn off all line info if beyond a certain range
        # !! maybe only go through for lines which are visible?
        for i in xrange(len(self.line_info)):
            scl = self._dropby(i)
            ypt = ymax - scl*yran
            
            # change the y of the lines info
            self.line_info[i].set_y(ypt)

            wl = self.lines[i]            
            visible = False
            if xmin < wl < xmax: visible = True
            self._set_info_visible(i, visible)
            
    def update_lines (self,bounds):
        print "!! need to update line update_lines_and_info"
        
        xmin,xmax,ymin,ymax = bounds
        yran = ymax - ymin
        lymin = ymin - self._overshoot_y*yran
        lymax = ymax + self._overshoot_y*yran
        # !! maybe only go through for lines which are visible?
        for i in xrange(len(self.plines)):
            self.plines[i].set_ydata([lymin,lymax])
            wl = self.lines[i]            
            visible = False
            if xmin < wl < xmax: visible = True
            self._set_line_visible(i, visible)
            
    def update_lines_and_info (self,bounds):
        if self._prev_bounds == bounds: return False
        self._prev_bounds = bounds
        xmin,xmax,ymin,ymax = bounds

        yran = ymax - ymin
        lymin = ymin - self._overshoot_y*yran
        lymax = ymax + self._overshoot_y*yran
       
        vis_ran, invis_ran, invis_info_ran = self._get_visible_in_ranges([xmin,xmax])
        
        # set previous visible, invisible
        self._set_prev_visible(invis_ran, invis_info_ran)

        # update for lines which are visible
        for i in vis_ran:            
            self.plines[i].set_ydata([lymin,lymax])
            self._set_line_visible(i, True)

            if not self._has_info: continue
            scl = self._dropby(i)
            ypt = ymax - scl*yran       
            self.line_info[i].set_y(ypt)
            self._set_info_visible(i, True)
        return True

    def update (self):
        """
        Plot Lines
        only need to update when bounds change
        """
        changed = self.update_lines_and_info(self.ax.axis())
        return changed
 
pass
#!! Note: PlotDataSingleOrder and PlotData was a hacky way to plot only one order and have it still work as the same plotting function within the plot data
class PlotDataSingleOrder:
    def __init__ (self,ax,spec_obj,starti=0,**mpl_kwargs):
        self.spec_obj = spec_obj
        self._ordi = starti
        self._previ = 0

        X = spec_obj.get_wl(starti)
        Y = spec_obj.get_data(starti)
        
        self._bounds = (np.min(X),np.max(X),np.min(Y),np.max(Y))
        
        c,_ = alt_order_colors(0)
        self.plotdata, = ax.plot(X,Y,color=c, markersize=7, **mpl_kwargs)
            
    def get_current_state (self):
        curstate = {"type":"PlotDataSingleOrder",
                    '_ordi':self._ordi,
                    '_previ':self._previ,
                    'vis':self.get_visible(),
                    'spec_obj':self.spec_obj.copy()}
        
        return curstate

    def apply_state (self,prevstate):
        if prevstate['type'] != 'PlotDataSingleOrder': raise ValueError("Whoops, the given prev state is wrong")
        self.set_visible(False)
        ordi = prevstate['_ordi']
        previ = prevstate['_previ']

        self.update_plot_data(prevstate['spec_obj'])
        self.set_visible(prevstate['vis'])        
        
        # !! check to make sure the ordi and previ are resonable values
        if ordi is None: 
            self._ordi = None
            self._previ = 0
            self.set_visible(False)
        else:
            
            self.select_order(ordi,direction='-')
            self._previ = previ 
        
    def get_current_i (self):
        return deepcopy(self._ordi)
    
    def get_prev_i (self):
        return deepcopy(self._previ)
    
    def get_zorder_bounds (self):
        return self.plotdata.get_zorder()
    
    def get_data_bounds (self):
        x,y = self.get_all_plot_data()
        return (np.min(x),np.max(x),np.min(y),np.max(y))
           
    def scan_through_orders (self,direction,always_find = True):
        # returns True if order changed
        
        if direction not in ('-','+',None): raise StandardError("Direction needs to be +/-")
        
        inc = 0
        if direction == '-': inc = -1
        elif direction == '+': inc = 1
        
        if self._ordi is None:
            if always_find: index = self.get_prev_i()
            else: return False
        else: index = self._ordi
        
        index += inc
        if index < 0: index = 0
        if index >= self.spec_obj.shape[1]: index = self.spec_obj.shape[1]-1
        
        return self.select_order(index)
        
    def _find_non_deleted (self,index,direction):
        if direction not in ('-','+',None): return
        # find the next highest order which hasn't been deleted
        if direction == '+':
            ran = range(index,self.spec_obj.shape[1])

        # find the next lowest order which hasn't been deleted
        elif direction == '-':
            ran = range(index+1)
            ran.reverse()

        for i in ran:
            if len(self.spec_obj.get_wl(i)) == 0: continue
            return (True,deepcopy(i))

        if True: return (False,0)  
    
    def select_order (self,index,direction=None): 
        index = int(index)
        if index not in xrange(self.spec_obj.shape[1]): raise StandardError("Whoops, index out of range for the spectrum")

        #====== Find the next order which hasn't been deleted
        found_order = False
        if direction == None:
            found_order,i = self._find_non_deleted(index,'+')
            if not found_order: found_order,i = self._find_non_deleted(index,'-')
            if not found_order: raise StandardError("Whoops, looks like all orders were deleted!")
        else:
            found_order,i = self._find_non_deleted(index,direction) 
            if not found_order:
                if direction == '+': found_order,i = self._find_non_deleted(index,'-')
                else: found_order,i = self._find_non_deleted(index,'+') 
            if not found_order: raise StandardError("Whoops, looks like all orders were deleted!")    
        index = i                

        # if the order is the same make no change
        if index == self._ordi: return False

        # turn off previous if changing
        if (self._ordi is not None):
            self._previ = self.get_current_i()
            self.set_visible(False)
    
        self._ordi = index
        X = self.spec_obj.get_wl(index)
        Y = self.spec_obj.get_data(index)
        self.plotdata.set_xdata(X)
        self.plotdata.set_ydata(Y)
        self.set_visible(True)
        return True
    
    pass
    ####################### These are used by DataPlot from PlotData:

    def get_all_plot_data (self):
        return self.plotdata.get_xydata().T
    
    def update_plot_data (self,spec_obj):   
        i = self._ordi
        if i is None: self.set_visible(False)
        self.spec_obj = spec_obj
        X = spec_obj.get_wl(i)
        Y = spec_obj.get_data(i)
        
        self.plotdata.set_xdata(X)   
        self.plotdata.set_ydata(Y)   
    
    def get_visible (self):
        return self.plotdata.get_visible()
    
    def set_visible (self,truth):
        return self.plotdata.set_visible(truth)
    
    def set_visible_window (self,bounds): #@UnusedParameter
        return False
             
class PlotData:
    """
    Adds data to the plot and creates a dictionary for each ax.plot using the index
    """
    def __init__ (self, ax, spec_obj, zorder=0, **mpl_kwargs):
        self.ax = ax
        self.ax.set_xlabel("WAVELENGTH")
        self.ax.set_ylabel("INTENSITY")
        self.ax.xaxis.set_major_formatter(FormatStrFormatter('%10.1f'))

        self.zorder = int(zorder)
        self.delta_zorder = 1

        xmin = spec_obj.get_min()[0]
        xmax = spec_obj.get_max()[0]
        if xmin != None: self.ax.set_xlim(xmin, xmax)
        
        self._ptime = time.time()
        self._dtime = 0.1
                
        self.plot_data = {}
        oi = 0
        for i in range(spec_obj.shape[1]):
            wl = spec_obj.get_wl(i)
            data = spec_obj.get_data(i)

            ocolor, oi = alt_order_colors(oi)
            oi += 1
            zo = self.zorder + self.delta_zorder * i
            dat, = self.ax.plot(wl, data, color=ocolor, markersize=7, zorder=zo, **mpl_kwargs)
            self.plot_data[i] = dat        
   
    def get_current_state (self):
        curstate = {'type':'PlotData',
                    'vis':self.get_visible(),
                    'datastyle':self.get_current_data_style()}
   
        return curstate
    
    def apply_state (self,prevstate):
        if prevstate['type'] != 'PlotData': raise ValueError("Got the wrong prevstate")
        self.set_visible(False)
        
        self.apply_prev_data_style(prevstate['datastyle'])
        self.set_visible(prevstate['vis'])
                
    def __getitem__ (self, key):
        if key not in self.plot_data:
            print "Order is not in plot data: " + str(key)
            return None

        return self.plot_data[key]
      
    def __contains__ (self, val):
        return val in self.plot_data
      
    def __iter__ (self):
        keylist = self.keys()
        return iter(keylist)
      
    def keys (self):
        return deepcopy(self.plot_data.keys())
          
    def get_zorder_bounds (self):
        zo_max = self.zorder + self.delta_zorder * len(self.plot_data)
        return (deepcopy(self.zorder), zo_max + 1)
      
    def bring_order_to_top (self, index):
        self.reset_order_zorder()
        zo_max = self.zorder + self.delta_zorder * len(self.plot_data)
        if index not in self.plot_data: return
        self.plot_data[index].set_zorder(zo_max + 1)
        
    def reset_order_zorder (self):
        for i in self.plot_data:
            zo = self.zorder + self.delta_zorder * i
            self.plot_data[i].set_zorder(zo)
        
    def update_plot_data (self, spec_obj):
        for ordi in self.plot_data:
            wl = deepcopy(spec_obj.get_wl(ordi))
            dat = deepcopy(spec_obj.get_data(ordi))

            if len(wl) == 0:
                spec_obj.set_use_cropped(False)
                all_wl = spec_obj.get_wl(ordi)[0]
                all_dat = spec_obj.get_dat(ordi)[0]
                wl = [all_wl[0], all_wl[-1]]
                dat = [np.mean(all_dat), np.mean(all_dat)]
                spec_obj.set_use_cropped('previous')

            self.plot_data[ordi].set_xdata(wl)
            self.plot_data[ordi].set_ydata(dat)

    def set_visible_window (self, bounds, force_set=False):
        """
        Will set all data outside of the bounds as invisible
        """
        xmin, xmax, _ymin, _ymax = bounds

        if not force_set:
            if time.time() - self._ptime < self._dtime: return
            else: self._ptime = time.time()
        
        for ordi in self.plot_data:
            xdata = self.plot_data[ordi].get_xdata()
            visible = False
            if len(xdata) > 2 and xmin < xdata[-1] and xmax > xdata[0]: 
                visible = True
            self.plot_data[ordi].set_visible(visible)
            
    def get_visible (self):
        check = np.ones(self.plot_data.keys()) < 0
        for key in self.plot_data:
            check = self.plot_data[key].get_visible()
        # !! could add check for some being visible 

        if check[0]: return True
        
    def set_visible (self, truth):
        truth = bool(truth)
        for key in self.plot_data:
            self.plot_data[key].set_visible(truth)

    def get_current_data_style (self):
        if len(self.plot_data)==0: return None
        cur_linestyle = self.plot_data[0].get_linestyle()
        cur_marker = self.plot_data[0].get_marker()
        alpha = 1.0
        prevstyle = (cur_linestyle,cur_marker,alpha)
        return prevstyle 
    
    def apply_prev_data_style (self,prevstyle):
        if prevstyle is None: return
        new_linestyle,new_marker,alpha = prevstyle
        for i in self.plot_data:
            self.plot_data[i].set_linestyle(new_linestyle)
            self.plot_data[i].set_marker(new_marker)
            self.plot_data[i].set_alpha(alpha)        

    def toggle_data_line_scatter (self): 
        """
        Change between line, scatter, line+scatter plot of the data
        """
        if len(self.plot_data)==0: return
        cur_linestyle, cur_marker, alpha = self.get_current_data_style()
        new_marker, new_linestyle, alpha = data_line_scatter(cur_linestyle,cur_marker,alpha)
        self.apply_prev_data_style(new_linestyle,new_marker,alpha)

    def get_all_plot_data (self):
        orders = self.plot_data.keys()
        all_xdata = self.plot_data[orders[0]].get_xdata()
        all_ydata = self.plot_data[orders[0]].get_ydata()
        
        for i in range(1, len(orders)):
            ordi = orders[i]
            all_xdata = np.concatenate((all_xdata, self.plot_data[ordi].get_xdata()))
            all_ydata = np.concatenate((all_ydata, self.plot_data[ordi].get_ydata()))
        return all_xdata, all_ydata

pass
################################################################################

class EditSubplotDialog (wx.Dialog, EventConnections):
    """
    Either set or return a dictionary axes_kwargs = {} which can be given to self.ax.set(**axes_kwargs)
    
    """
    
    def __init__ (self,ax,*args,**kwargs):
        self.ax = ax # xmin,xmax,ymin,ymax current
        
        super(EditSubplotDialog,self).__init__(None,-1,*args,**kwargs)
        
        self.InitUI()
        
        self.SetSize((400,600))
        self.SetPosition((860,30))

        self.SetTitle("Set the plot Parameters:")
        
        self._invalid_float = False
        
        
        self.Bind(wx.EVT_WINDOW_DESTROY,self.OnClose)
        self.Bind(wx.EVT_KEY_UP,self.OnKeyUp)

        self.key_cfg = KeyboardConfiguration()
        self.key_cfg.add_key('a','Toggle select all then deselect all')
        self.key_cfg.add_key('Space','Average checked continuum and apply to current order')
        self.key_cfg.add_key('h','Display this help screen')
        self.key_cfg.add_key('q','Quit sub-window')

        self.key_cfg.set_display_order(['a','Space','h','q'])
        self.key_cfg.check_display()
        self.key_cfg.add_mpl_key_convert(' ','Space')

    def InitUI (self):
        self.SetPlotLimits()

    def OnClose (self,event):
        print "CLOSING"
        pass

    def OnKeyUp (self,event):
        
        pass

    def OnKeyDown (self,event):
    
        
    
        # if 'q' then quit
        
        # if 'enter' then confirm
        
    
    
    
        # is this a typing in the box:
        if str(event.GetEventObject()) in self._plot_limits: self._plot_limits_entry(event)

    def _check_invalid_float (self,text):
        
        if len(text) == 0 or text == '-' or text[-1] == 'e': 
            # special case of typing a float
            # special case when you start typing '-'
            self._invalid_float = False
            return False            
            
        try:
            float(text or 0.0)
            self._invalid_float = False
        except: self._invalid_float = True
        return self._invalid_float

    def _invalid_float_display (self):
        if self._invalid_float: print "Invalid float"
        else: print "float ok"

    def _plot_limits_entry (self,event):
        which = str(event.GetEventObject())
        if which not in self._plot_limits: return
        
        # get text dialog
        txtctl = self._plot_limits[which][2]
        text = txtctl.GetValue()
        new_text = deepcopy(text)
        cursori = txtctl.InsertionPoint
        sel = txtctl.Selection
        is_selected = not (sel[0] == sel[1])

        
        # get key values
        #arrows = {314:'left', 316:'right'}
        
        keycode = event.GetKeyCode()
                
        if keycode in xrange(256): key = chr(keycode).lower()
        else: key = None
        
        if key is None:
            if keycode in (314,316): 
                i = txtctl.GetInsertionPoint()  
                
                if is_selected:
                    if keycode == 314: i = sel[0] - 1
                    if keycode == 316: i = sel[1] + 1
                else:       
                    if keycode == 314: i -= 1
                    if keycode == 316: i += 1
                i = np.clip(i,0,len(text))                
                txtctl.SetInsertionPoint(i)
            return
            
        # pressed ctl-k then delete rest of the line
        elif event.ControlDown():
            if key == 'e': insert_point = len(new_text)
            if key == 'a': insert_point = 0
            if key == 'k': 
                new_text = text[:cursori]
                insert_point = len(new_text)                      
        # delete stuff with backspace
        elif keycode == 8: 
            if is_selected: 
                pdb.set_trace()
                new_text = text[:sel[0]]+text[sel[1]:]
                insert_point = sel[0]
            else: 
                new_text = text[:cursori-1]+text[cursori:]
                insert_point = cursori-1
            
        # else add the text
        else:
            new_text = text[:cursori]+key+text[cursori:]
            insert_point = cursori+1

        self._check_invalid_float(new_text)
        self._invalid_float_display()
        
        if not self._invalid_float: 
            txtctl.SetValue(new_text)
            txtctl.SetInsertionPoint(insert_point)
        
    def SetPlotLimits (self):
        self.ax.figure.canvas.draw()
        
        xmin, xmax, ymin, ymax = self.ax.axis()        
#            def __init__(self):
#        wx.Frame.__init__(self, None, -1, 'Text Entry Example', size=(300, 100))
#        panel = wx.Panel(self, -1) 
#        basicLabel = wx.StaticText(panel, -1, "Basic Control:")
#        basicText = wx.TextCtrl(panel, -1, "I've entered some text!", size=(175, -1))
#        basicText.SetInsertionPoint(0)
#
#        pwdLabel = wx.StaticText(panel, -1, "Password:")
#        pwdText = wx.TextCtrl(panel, -1, "password", size=(175, -1),style=wx.TE_PASSWORD)
#        sizer = wx.FlexGridSizer(cols=2, hgap=6, vgap=6)
#        sizer.AddMany([basicLabel, basicText, pwdLabel, pwdText])
#        panel.SetSizer(sizer)
        
    
        panel = wx.Panel(self, -1)

        # plotlimits_box = wx.StaticBox(panel, label = "Plot limits:")
        # sbs = wx.StaticBoxSizer(plotlimits_box, orient=wx.VERTICAL)

        plot_limits = [['Xmin',xmin],
                       ['Xmax',xmax],
                       ['Ymin',ymin],
                       ['Ymax',ymax]]

        self._plot_limits = {}
        plot_sizer = []
        
        for i in xrange(len(plot_limits)):
            arr = plot_limits[i]
            which = arr[0]
            val = arr[1]
            
            label = wx.StaticText(panel,-1,which+"=")
            txt = wx.TextCtrl(panel,-1,format(val,'<1.2f'),size=(175,-1))
            
            txt.SetInsertionPoint(0)
            # txt.CanCopy(True)
            # txt.Bind(wx.EVT_KEY_DOWN, self.OnKeyDown, id2=9000+i)
            
            self._plot_limits[repr(txt)] = [which,label,txt,val]
            
            plot_sizer += [label,txt]
#        
#       
#        plt.SetSizer(self.xmin_txt)
#        plt.SetSizer(self.xmax_txt)
#       
#        btn1 = wx.Button(self,label='Ok')
#        # btn1.Bind(wx.EVT_BUTTON, self.Average_Continuum_Apply)
#        hbox.Add(btn1)
#        
#
#        btn3 = wx.Button(self,label='Cancel')
#        # btn3.Bind(wx.EVT_BUTTON, self.
#        hbox.Add(btn3)

                    
#        self.SetSizer(panel)

        sizer = wx.FlexGridSizer(cols=2, hgap=6, vgap=12)
        sizer.AddMany(plot_sizer)
        panel.SetSizer(sizer)
        panel.SetSize((300,500))
        
pass
################################################################################

class eyeSpecBaseLineEditor:
    def __init__ (self,ax,line_data,line_info=None,history=None,crop_xbounds=None,lock_line_positions=False,lock_order=False,color='#990033',zorder=0,alpha=.8): 
        self.ax = ax
        if crop_xbounds != None: crop_xbounds = (np.min(crop_xbounds),np.max(crop_xbounds))
        
        if history is None: self.history = History()
        else:
            if repr(history) != "eyeSpecHistory": raise ValueError("History of wrong type")
            self.history = history
        
        #!! add history
        
        self.plotlines = PlotLines(self.ax,line_data,line_info,line_info_vis=False,crop_xbounds=crop_xbounds,lock_order=lock_order)
        # could then set color
        self.sel_line = LineSelection(self.ax)
        self._start_xpt = None

        # !! need epsilon?
        # self.dataplot.get_epsilon('x')
        
        self._active = True

        self._visible = True

        # display selected text
        self.text = self.ax.text(0.45, 1.05, 'Line Selected: None',
                                 transform=self.ax.transAxes, va='top')
        
        # keep lines from being moved
        self._lock_line_positions = bool(lock_line_positions)
        
        self.timer = Timer(0.2)
        self._prev_bounds = self.ax.axis()
  
    def get_current_state (self):
        curstate = State('line_editor')
        curstate['plotlines'] = self.plotlines.get_current_state()
        curstate['sel_line'] = self.sel_line.get_current_state()
        # ? curstate['data_plot'] = self.dataplot.get_current_state() 
        curstate['vis'] = self.get_visible()
        curstate['_active'] = self._active
        curstate['lock'] = self._lock_line_positions
        return curstate
  
    def apply_state (self,prevstate):
        prevstate.check_id('line_editor')
        self.plotlines.apply_state(prevstate['plotlines'])
        self.sel_line.apply_state(prevstate['sel_line'])
        self.set_visible(prevstate('vis'))
        self._active = prevstate['active']
        self._lock_line_positions = prevstate['lock']
          
    def change_linelist (self,line_data,line_info=None,crop_xbounds=None):
        self.plotlines.hide()
        self.plotlines.change_linelist(line_data, line_info, crop_xbounds)
        self.plotlines.set_visible(True)
        
        self.sel_line.deselect_line()
        self.sel_line.reset()
  
    def set_lock_line_positions (self,truth):
        self._lock_line_positions = bool(truth)
        
    def select_line_in_plot (self):
        xmin,xmax = self.ax.axis()[:2]
        xmid = (xmin+xmax)/2.0
        L = self.plotlines.lines
        if np.any(((xmin < L)*(L < xmax))):
            self.select_line_by_x(xmid)
        else:
            self.deselect()
    
    def get_index_from_x (self,xpt):
        return np.argmin(np.abs(self.plotlines.lines - xpt))
    
    def select_line_by_x (self,xpt,epsilon_x=None):       
        dists = np.abs(self.plotlines.lines-xpt)
        index = np.argmin(dists)
        min_dist = dists[index]
        
        if epsilon_x is not None and (min_dist > epsilon_x): 
            print "deselected",index,min_dist,epsilon_x
            self.sel_line.deselect_line()
            return 

        self.select_line_by_index(index)

    def select_line_by_index (self,index,always_find=False,verbose = True):        
        if index is None:
            if always_find: index = self.sel_line.get_prev_line_index()
            else: return
            
        curi = self.sel_line.get_line_index()
        if curi == index: return
        elif curi is not None: self.sel_line.deselect_line()
            
        N = len(self.plotlines.lines)
        if index not in xrange(N):
            if always_find: index = N-1
            else:
                if verbose: print "index out of range:  "+str(index)+" not in "+str((0,N-1))
                return

        # !! put a find in data range
        # index = self.line_in_range(index,bounds)

        xpt = self.plotlines.lines[index]        
        self.sel_line.select_line(index, xpt)
        return xpt
          
    def get_selected_x (self):
        if not self.sel_line.is_selected(): return None
        else: return self.sel_line.get_xdata()
             
    def scan_through_lines (self,direction,verbose=True): 
        if direction not in ['-','+',None]: raise ValueError("direction not recognized")
        linei = self.sel_line.get_line_index()
        if linei is None: 
            linei = self.sel_line.get_prev_line_index()
            direction = None
            
        inc = 0
        if direction == '-': inc = -1
        elif direction == '+': inc = +1
        
        # !! find line over data
        
        index = linei + inc
        return self.select_line_by_index(index,verbose=verbose)
        
    def delete_selected_line (self):
        linei = self.sel_line.get_line_index()
        if linei is None: return
        self.plotlines.delete_line_by_index(linei)
        self.sel_line.deselect_line()
             
    def deselect (self): 
        self.sel_line.deselect_line()
        
    def set_visible (self,truth):
        truth = bool(truth)
        self._visible = truth
        self.text.set_visible(truth)
        if truth: 
            self.plotlines.set_window_visible(self.ax.axis(),False)
            self.plotlines.update()
        else:
            self.plotlines.hide()
            
    def get_visible (self):
        return deepcopy(self._visible)
           
    def move_selected_line (self,delta_x):
        index = self.sel_line.get_line_index()
        if index is None: return
        
        delta_x = float(delta_x)
        xcur = self.sel_line.get_xdata()
        xpt = xcur+delta_x
        
        self._change_the_sel_line(xpt)
                   
    def btn_press_lines (self,event,epsilon_x=None):
        """
        Interactive Line Editor
        returns True if something changed
        """

        if event.inaxes is None: return False
        if event.button != 1: return False
        
        # get the current
        curi = self.sel_line.get_line_index()
        now_sel = self.sel_line.is_selected()
        
        # use the xdata to select line
        self.select_line_by_x(event.xdata, epsilon_x)
        
        # check to see if you selected/deselected similar
        p1 = curi == self.sel_line.get_line_index()
        p2 = now_sel == self.sel_line.is_selected()
        
        if self.sel_line.is_selected(): self._start_xpt = self.sel_line.get_xdata()
        
        if p2 and p1: return False
        else: return True         
        
    def mot_notify_lines (self,event):
        """
        Interactive Line Editor
        returns True if something changed
        """

        if event.inaxes is None: return False
        if event.button != 1: return False
        if not self.sel_line.is_selected(): return False
        if self._lock_line_positions: return False  

        # !! todo remove these if they don't crash anything
        # curx = self.sel_line.get_xdata()
        # xnew = curx + (event.xdata - self._start_xpt)
        
        if self.timer.check():
            self._change_the_sel_line(event.xdata)
            return True
        return False
    
    def btn_click_lines (self,event,epsilon_x=None):
        """
        Interactive Line Editor
        returns True if something changed
        """  
        if event.inaxes is None: return False
        if event.button != 1: return False
        
        if not self.sel_line.is_selected(): return False

        # use the xdata to see if you need to deselect line
        self.select_line_by_x(event.xdata, epsilon_x)
        if not self.sel_line.is_selected(): return True
        return False # you reselected the line
          
    def btn_release_lines (self,event):
        # ok not to be in the axis any more
        if event.button != 1: return False
        if not self.sel_line.is_selected(): return False
        if event.xdata is None: return False

        self._change_the_sel_line(event.xdata)
        if event.inaxes is None:
            self.plotlines.set_window_visible(self.ax.axis()[:2])
        return True
        # !! should I put the self.plotlines.update() here?
        # !! what about a self.update()
 
    def _change_the_sel_line (self,xpt):
        if self._lock_line_positions: return
        self.sel_line.set_xdata(xpt)
        index = self.sel_line.get_line_index()
        self.plotlines.change_line(index,xpt)  
        self.plotlines.change_info(index)
    
    def save_current_parameters (self): pass
    
    def open_current_parameters (self): pass
               
    def update_text (self):
        if not self.text.get_visible(): return False
        
        linei = self.sel_line.get_line_index()

        prev_text = self.text.get_text()
        display_text = 'Line Selected: '+str(linei)
        if linei != None: 
            lines = self.plotlines.lines
            display_text += "/"+str(len(lines)-1)
            display_text += "  ("+format(lines[linei],'0.3f')+")"   
            
        if prev_text == display_text: return False
        self.text.set_text(display_text)
        return True
    
    def update (self):
        c1 = self.update_text()
        c2 = self.plotlines.update()
        return c1 or c2
        
class eyeSpecBaseDataPlot:
    """
    This plots a spectrum object onto a canvas and adds some capabilities
    """
    def __init__ (self, ax, spec_obj, auto_scale_focus='full', plotclass=PlotData, **mpl_kwargs):
        #----------------------------------------------------------------------#
        # initiate where this stuff goes
        self.ax = ax
        
        if repr(spec_obj).find('eyeSpec_spec') != 0: raise ValueError("Input spectrum object must be eyeSpec_spec class")
        self.orig_spec = spec_obj.copy()
        self.spec_obj = spec_obj.copy()
        self.spec_obj.edit.set_protection(False)

        self.plot_data = plotclass(self.ax, self.spec_obj, **mpl_kwargs)
        self.overplot_data = []

        #----------------------------------------------------------------------#
        # define how close to data to select
        self._epsilon_type = 'relative'  # 'set' or 'relative'
        self._epsilon = 0.01  # 1% size
        # !! if self.epslion_type is 'set' i'll probably need self.epsilon = [x,y]
        self._epsilon_x = deepcopy(self._epsilon * (self.ax.axis()[1] - self.ax.axis()[0]))  # the acutal range in xdata space
        self._epsilon_y = deepcopy(self._epsilon * (self.ax.axis()[3] - self.ax.axis()[2]))  # the acutal range in ydata space
        
        #----------------------------------------------------------------------#
        # set up axes
        self.set_auto_scaling(3)
        self.set_auto_scale_focus(auto_scale_focus) 
        self._tmp_no_scale = False
        self._auto_xrange = 50.0  # units = [\dot{A}], !! I could do something based on the number of pixels
    
    
        self._record_pts = []
        self._max_record_pts = 10
    
        self._ptime = time.time()
        self._dtime = 0.02
    
        self._prev_bounds = self.ax.axis()
        
        #----------------------------------------------------------------------#
        # plot objects
        self._is_grid_on = False

#        #----------------------------------------------------------------------#       
#        # keyboard configuration
#        self.key_cfg = KeyboardConfiguration()
#        self.key_cfg.add_key('h', 'Dispay this help screen')
#        self.key_cfg.add_key('g', 'Toggle on/off visual grid')
#        self.key_cfg.add_key('o', 'Open')
#        self.key_cfg.add_key('p', 'Toggle Matplotlib pan/zoom tool')
#        self.key_cfg.add_key('q', 'Quit')
#        self.key_cfg.add_key('s', 'Save')
#        self.key_cfg.add_key('z', 'Toggle Matplotlib zoom rect tool')
#
#        self.key_cfg.add_key('`', 'Toggle auto scaling options')
#        self.key_cfg.add_key(';', 'Toggle showing the data as a line, scatter, or both')
#
#        self.key_cfg.add_key('>', 'Scan through orders in positive wavelength direction')
#        self.key_cfg.add_key('<', 'Scan through orders in negative wavelength direction')
#
#        self.key_cfg.set_display_order(['g', 'h', 'o', 'p', 'q', 's', 'z', '`', ';', '>', '<'])
#        self.key_cfg.check_display()
#        self.key_cfg.add_mpl_key_convert('.', '>')
#        self.key_cfg.add_mpl_key_convert(',', '<')
#        
        #----------------------------------------------------------------------#       
        # input/output
        T = time.gmtime()
        curdate = format(T.tm_mon,'02')+format(T.tm_mday,'02')+str(T.tm_year) 
        
        self.iochoices = SaveOpenChoices()
        self.iochoices.add('Input Spectrum',
                           '(*.spec)|*.spec',
                           'InputSpec_'+curdate+".spec",
                           'eyeSpec Spectrum File',
                           self.save_orig_obj,
                           self.open_orig_obj)
        
        self.iochoices.add('Edited Spectrum',
                           '(*.spec)|*.spec',
                           'EditSpec_'+curdate+'.spec',
                           'eyeSpec Spectrum File',
                           self.save_spec_obj,
                           self.open_orig_obj)

    def __repr__ (self):
        return 'eyeSpecBaseDataPlot'

    def record_pts (self,x,y):
        self._record_pts.append([x,y])
        if len(self._record_pts) > self._max_record_pts: del self._record_pts[0]
        
    def get_record_pts (self,index):
        return self._record_pts[index]

    def overplot_from_file (self,filename,Zorder='top',**mpl_kwargs):
        try: xy = np.loadtxt(filename,unpack=True,usecols=[0,1])
        except: raise IOError("Couldn't import X,Y pairs from the first to columns of '"+filename+"' using np.loadtxt")
        self.overplot(xy[0],xy[1],Zorder,**mpl_kwargs)
        
    def overplot (self,X,Y,Zorder='top',**mpl_kwargs):
        zmin,zmax = self.plot_data.get_zorder_bounds()[0]
        if Zorder=='top': mpl_kwargs['zorder'] = zmax+2
        elif Zorder == 'bot':mpl_kwargs['zorder'] = zmin-2

        pdat, = self.ax.plot(X,Y,**mpl_kwargs)
        self.overplot_data.append(pdat)

    def match_up_spec_obj (self):
        return self.spec_obj

    def save_orig_obj (self, filename, clobber=True):
        filename = str(filename)
        if os.path.exists(filename) and not clobber: raise IOError("File exists: " + filename)
        save(self.orig_obj, filename, clobber=True)
        
    def save_spec_obj (self, filename, clobber=True):
        filename = str(filename)
        if os.path.exists(filename) and not clobber: raise IOError("File exists: " + filename)
        save(self.spec_obj, filename, clobber=True)
                
    def open_orig_obj (self, filename):
        filename = str(filename)
        if not os.path.exists(filename): raise IOError("File does not exist:" + filename)
        self.orig_obj = readin_spec(filename)
        
    def open_spec_obj (self, filename):
        filename = str(filename)
        if not os.path.exists(filename): raise IOError("File does not exist:" + filename)
        self.spec_obj = readin_spec(filename)
                
    def save_parameters_data_plot (self, filename, clobber=True):pass
    
    def open_parameters_data_plot (self, filename):pass

    def scan_through_walking (self,direction=None):
        """
        Takes the current x range and steps it forward||backwards 
        """
        if direction is None: return
        if direction not in ['+','-']: return
        cur_xmin = self.ax.axis()[0]
        cur_xmax = self.ax.axis()[1]
        cur_xran = cur_xmax - cur_xmin
        all_xdata,_all_ydata = self.plot_data.get_all_plot_data()

        if direction == '+': # move forwards
            new_xmin = cur_xmax
            new_xmax = cur_xmax + cur_xran
            check = (all_xdata > new_xmin)
            if len(np.ones(len(check))[check]) > 5: self.ax.set_xlim(new_xmin,new_xmax)

        elif direction == '-': # move backwards
            new_xmin = cur_xmin - cur_xran
            new_xmax = cur_xmin
            check = (all_xdata < new_xmax)
            if len(np.ones(len(check))[check]) > 5: self.ax.set_xlim(new_xmin,new_xmax)

    def sprog_open (self):
        sprog = self.Input_From_File(self._open_choices)
        if sprog is not None:
            print "=== READING IN PROGRESS SAVE FOR CTM_EDITOR"
            self._open_snr_progress_save(sprog)
            print "=== PROGRESS READ IN FOR CTM_EDITOR SAVE NUMBER:" + str(sprog.get_current_save_num())

    def display_help (self):
        print "="*60
        editor_name = 'Data Editor'
        print "="*10 + " " + editor_name
        # print ("-"*60)
        # lines = ["This program allows editing data",
        # print "\n".join(lines)
        print ("-"*60)
        self.key_cfg.display_short_info() 
        print ("-"*60)
        print "Continue, press 'q' to finish and quit\n"

    def toggle_grid (self, truth=None):
        if '_one_line' not in dir(self):
            self._one_line = self.ax.axhline(y=1.0, color='k', lw=2, alpha=.4, zorder= -20, visible=False, linestyle='-')

        if truth is not None:
            if truth: self._is_grid_on = False  # will be forced on
            else: self._is_grid_on = True  # will be forced off

        if self._is_grid_on:
            self.ax.xaxis.grid(False)
            self.ax.yaxis.grid(False)
            self._one_line.set_visible(False)
            self._is_grid_on = False
        else:
            self.ax.xaxis.grid(b=True, color='k', alpha=.5, zorder= -21, linestyle=':', lw=1.5)
            self.ax.yaxis.grid(b=True, color='k', alpha=.5, zorder= -21, linestyle=':', lw=1.5)
            self._one_line.set_visible(True)
            self._is_grid_on = True

    def toggle_pan_tool (self, truth=None): 
        # !! currently can't force set this on or off
        self.ax.figure.canvas.toolbar.pan()

    def toggle_zoom_tool (self, truth=None):
        # !! currently can't force set this on or off
        self.ax.figure.canvas.toolbar.zoom()

    def get_current_state (self):
        plot_data_state = self.plot_data.get_current_state()
        if 'spec_obj' in plot_data_state: spec_obj = None
        else: spec_obj = self.spec_obj.copy()
        
        curstate = {'type':'DataPlot',
                    'spec_obj':spec_obj,
                    'plot_data':plot_data_state,
                    '_auto_scale_opt':self._auto_scale_opt,
                    '_auto_scale_focus':self._auto_scale_focus,
                    '_epsilon_type':self._epsilon_type}
        
        return curstate

    def apply_state (self, prevstate, replot=True):
        if prevstate['type'] != 'DataPlot': raise ValueError("Received wrong previous state")
        
        
        
        spec_obj = prevstate['spec_obj']
        plot_data_prev_state = prevstate['plot_data']
        if spec_obj is None:
            spec_obj = plot_data_prev_state['spec_obj']
        
        self.plot_data.apply_state(plot_data_prev_state)
        self.spec_obj = spec_obj

        self.set_auto_scaling(prevstate['_auto_scale_opt'])
        self.set_auto_scale_focus(prevstate['_auto_scale_focus'])
        self.set_epsilon_type(prevstate['_epsilon_type'])
        
        self.update_epsilon() 
        if replot and self.update(): self.ax.figure.canvas.draw()
        
        return self.spec_obj

    def undo_redo (self,prevstate):
        curstate = self.get_current_state()
        self.apply_state(prevstate)
        return curstate

    def skip_auto_scale (self):
        self._tmp_no_scale = True
        
    def continue_auto_scale (self):
        self._tmp_no_scale = False

    def get_auto_scaling_opt (self):
        if self._auto_scale_opt == 0  : return (0,"Auto scale X and Y")
        elif self._auto_scale_opt == 1: return (1,"Auto scale X axis only")
        elif self._auto_scale_opt == 2: return (2,"Auto scale Y axis only")
        elif self._auto_scale_opt == 3: return (3,"Auto scale axis off")

    def set_auto_scaling (self, auto_scale_opt=None):
        """
        Set the auto scaling options
        if auto_scale_opt is None then it will toggle through
        """
        txt = True
        if auto_scale_opt is not None:
            if int(auto_scale_opt) not in range(4):
                raise TypeError("Auto scale option must be in:" + np.array(range(4), dtype='a1'))
            self._auto_scale_opt = auto_scale_opt
            txt = False
        else:
            self._auto_scale_opt += 1
            if self._auto_scale_opt == 4: self._auto_scale_opt = 0

        if self._auto_scale_opt == 0: # auto x and y
            if txt: print "Auto scale X and Y axis"
            self._auto_scale_x = True
            self._auto_scale_y = True
        elif self._auto_scale_opt == 1: # auto x
            if txt: print "Auto scale X axis only"
            self._auto_scale_x = True
            self._auto_scale_y = False
        elif self._auto_scale_opt == 2: # auto y
            if txt: print "Auto scale Y axis only"
            self._auto_scale_x = False
            self._auto_scale_y = True
        elif self._auto_scale_opt == 3: # auto off
            if txt: print "Auto scale axis off"
            self._auto_scale_x = False
            self._auto_scale_y = False

    def set_auto_scale_focus (self, auto_scale_focus):
        if auto_scale_focus not in ['edges', 'center', 'full', 'edge_right', 'edge_left']:
            print "HeadsUp: auto_scale_focus must be either 'edges','center', or 'full'. Setting to 'full' from:" + str(auto_scale_focus)
            auto_scale_focus = 'full'
        self._auto_scale_focus = str(auto_scale_focus)

    def set_auto_scale_range (self, xran):
        self._auto_xrange = float(xran)

    def get_auto_xrange (self):
        return self._auto_xrange

    def get_auto_scale_xy (self):
        return deepcopy((self._auto_scale_x, self._auto_scale_y))
    
    def get_auto_scale_focus (self):
        return deepcopy(self._auto_scale_focus)

    def is_toolbar_button_on (self):
        if self.ax.figure.canvas.toolbar.mode != '': return True
        else: return False

    def timed_update (self,del_time=None):
        if del_time is None: del_time = self._dtime
        if time.time() - self._ptime > del_time:
            self._ptime = time.time()
            return True
        else: return False

    def get_visible (self):
        return self.plot_data.get_visible()    
        
    def set_visible (self, truth):
        truth = bool(truth)
        self.plot_data.set_visible()

    def set_epsilon_type (self, epsilon_type=None):
        if epsilon_type not in ['relative', 'absolute']: raise ValueError("epsilon_type must be either relative or absolute")
        self._epsilon_type = epsilon_type

    def set_epsilon (self, epsilon_value=None):
        if epsilon_value is not None: 
            self._epsilon = float(epsilon_value)
            
    def update_epsilon (self):
        if self._epsilon_type == 'relative':
            self._epsilon_x = self._epsilon * (self.ax.axis()[1] - self.ax.axis()[0])
            self._epsilon_y = self._epsilon * (self.ax.axis()[3] - self.ax.axis()[2])
        elif self._epsilon_type == 'absolute':
            self._epsilon_x = self._epsilon
            self._epsilon_y = self._epsilon

    def get_epsilon_type (self):
        return deepcopy(self._epsilon_type)

    def get_epsilon (self, which=None):
        epsx = deepcopy(self._epsilon_x)
        epsy = deepcopy(self._epsilon_y)
        if which is 'x': return epsx
        elif which is 'y': return epsy
        elif which is 'xy': return (epsx, epsy)
        else: return deepcopy(self._epsilon)

    def scale_x (self, ordi=None, focus='default'):
        if not self._auto_scale_x: return
        if ordi is None: return
        ordi = int(ordi)
        if ordi not in range(self.spec_obj.shape[1]): return
        if not self.bounds_changed(): return
        
        # possible to scale x to: center, full, edge_left, edge_right
        auto_focus_options = ('default', 'center', 'full', 'edge_left', 'edge_right', 'edges')
        if focus not in auto_focus_options: raise ValueError("focus must be in: " + ", ".join(auto_focus_options) + "\n instead got:" + str(focus))


        xmin = self.spec_obj.get_min(ordi)[0]
        xmax = self.spec_obj.get_max(ordi)[0]
        if xmin == xmax: return
        xmid = (xmax + xmin) / 2.
        xran = xmax - xmin

        cur_xmin = self.ax.axis()[0]
        cur_xmax = self.ax.axis()[1]

        if (xmin == cur_xmin and xmax == cur_xmax) or xmin is None:
            return

        if focus == 'default': 
            focus = deepcopy(self._auto_scale_focus)
            if focus == 'edges': focus = 'edge_right'

        auto_xran = self.get_auto_xrange()

        if focus is 'edge_left': self.ax.set_xlim(xmin - auto_xran, xmin + auto_xran)
        elif focus is 'edge_right': self.ax.set_xlim(xmax - auto_xran, xmax + auto_xran)
        elif focus is 'center': self.ax.set_xlim(xmid - auto_xran, xmid + auto_xran)
        elif focus is 'full': self.ax.set_xlim(xmin - 0.05 * xran, xmax + 0.05 * xran)
        
    def scale_y (self):
        """
        Scale the y axis to match the data
        ""        """

        if not self._auto_scale_y: return
        if not self.bounds_changed(): return
        
        over_scale = 0.01

        all_xdata, all_ydata = self.plot_data.get_all_plot_data()

        xmin, xmax = self.ax.axis()[0], self.ax.axis()[1]
        mask = (xmin < all_xdata) * (all_xdata < xmax)
        ydata = all_ydata[mask]
        if len(ydata) == 0: return
        ymin, ymax = np.min(ydata), np.max(ydata)

        cur_ymin = self.ax.axis()[2]
        cur_ymax = self.ax.axis()[3]

        # if the bounds are the same don't change
        if ymin == ymax: return
        yran = ymax - ymin
        ymin -= over_scale * yran
        ymax += over_scale * yran
        if ymin == cur_ymin and ymax == cur_ymax: return
        
        self.ax.set_ylim(ymin,ymax)

    def bounds_changed (self,allow_reset=False):
        # print "= bounds called"
        if self._prev_bounds != self.ax.axis():
            if allow_reset:self._prev_bounds = self.ax.axis()
            return True
        else: return False
        
    def update (self,ordi=None):   
        c1 = False
        c2 = False
        
        if self._tmp_no_scale: self._tmp_no_scale = False
        else:
            self.scale_x(ordi)
            self.scale_y()
            c1 = True
            
        if self.bounds_changed(): 
            self.update_epsilon()  
            self.plot_data.set_visible_window(self.ax.axis())
            c2 = True
        return (c1 or c2)
                      
class eyeSpecBaseMainPanel (wx.Panel):
    """ The main panel of the frame which is split into the important sub-panels"""
    def __init__ (self,parent_frame, which_split):
        self.pframe = parent_frame
        wx.Panel.__init__(self,parent_frame)  

        self._init_splits(which_split)
        self._which_split = which_split
        
        self.textpanel = eyeSpecTextRedirectPanel(self.Split0)
        self.redirect_text_panel = self.textpanel.output_window

    def _init_splits (self,which_split):
        split_opts = ['split_top','split_top_left']
        if which_split not in split_opts: raise ValueError("Which split not in : "+",".join(split_opts))
        
        # default 'split_top'
        self.Split0 = wx.SplitterWindow(self)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        
        if which_split == 'split_top_left':
            self.Split1 = wx.SplitterWindow(self.Split0)

    def split_top (self,panel):
        if self._which_split is not 'split_top': return        
        self.Split0.SplitHorizontally(panel,self.textpanel)
        _wid,hei = self.pframe.GetSizeTuple()
        self.Split0.SetSashPosition(1.8*hei)
        self.Split0.SetMinimumPaneSize(1)
        self.Split0.SetSashGravity(1)
                
        self.sizer.Add(self.Split0,1,wx.EXPAND)
        self.SetSizer(self.sizer)

    def split_top_left (self,panel1,panel2):
        if self._which_split != 'split_top_left': return
                                
        self.Split1.SplitVertically(panel1,panel2)
        self.Split0.SplitHorizontally(self.Split1,self.textpanel)

        wid,hei = self.pframe.GetSizeTuple()
        self.Split0.SetSashPosition(1.8*hei)
        self.Split0.SetMinimumPaneSize(1)
        self.Split0.SetSashGravity(1)
        
        self.Split1.SetSashPosition(3*wid)
        self.Split1.SetMinimumPaneSize(1)
        self.Split1.SetSashGravity(1)
                
        self.sizer.Add(self.Split0,1,wx.EXPAND)
        self.SetSizer(self.sizer)
                
class eyeSpecBaseDataPanel (wx.Panel):
    def __init__ (self,parent_panel,parent_frame,include_cursor=False):
        wx.Panel.__init__(self,parent_panel)
                
        self.ppanel = parent_panel
        self.pframe = parent_frame
                
        self.SetBackgroundColour("#FFFFFF")
 
        # initiate the figure and canvas that the data will go onto
        self.figure = Figure(figsize=(4,.4)) # default size (8,6)
        self.ax = self.figure.add_subplot(111)

        scale_fig = {'top'   : 4,#20,
                     'bottom': 4,#30,
                     'right' : 5,
                     'left'  : 25,#30,
                     'hspace': 0,#30,
                     'wspace': 0}#30}

        figure_adjust_borders(self.figure,scale_fig)
        self.canvas = FigureCanvas(self, -1, self.figure)
       
        # set up the sizer 
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas,1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()
        
        # add the toolbar 
        self.add_toolbar()

        # connect the cursor to this panel
        self._include_cursor = not include_cursor
        self.toggle_cursor()
            
        # Change cursor as it enters/exits the screen
        if wx.Platform == '__WXMAC__':
            self.canvas.mpl_connect('axes_enter_event',self.axes_enter_callback)
            self.canvas.mpl_connect('axes_exit_event',self.axes_exit_callback)
        
                
        # Update the cursor position
        self.statusbar_cid = self.canvas.mpl_connect('motion_notify_event', self.UpdateStatusBar)
        
        # to disconnect:
        # del self.canvas.callbacks.callbacks['motion_notify_event'][self.statusbar_cid]


        #-------------------------------------------------#
        # fancy show first
        self._onkeystart_cid = self.canvas.mpl_connect('key_press_event',self.OnStart)
        self._onbutstart_cid = self.canvas.mpl_connect('button_press_event',self.OnStart)
        print "------ CLICK ON PLOT TO BEGIN -------"

    def post_start_delete_cids (self):
        del self.canvas.callbacks.callbacks['key_press_event'][self._onkeystart_cid]
        del self.canvas.callbacks.callbacks['button_press_event'][self._onbutstart_cid]

    def toggle_cursor (self):
        if 'cursor' not in dir(self): self.cursor= Cursor(self.ax,toolbar = self.toolbar) 
        if self._include_cursor: self.cursor.hide()
        else: self.cursor.show()
        self._include_cursor = not self._include_cursor 

    def OnStart (self,_event):
        # this runs on first keystroke or mouse click
        self.post_start_delete_cids()
        
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
        
        self.pframe.statusBar.SetStatusText(st,0)
        
    def OnRefresh(self,_event):
        self.canvas.draw()

    def axes_enter_callback (self,_event):        
        self.canvas.SetCursor(wx.StockCursor(wx.CURSOR_CROSS))#wx.CURSOR_BLANK))
        self.canvas.Update()

    def axes_exit_callback (self,_event):
        self.canvas.SetCursor(wx.StockCursor(wx.CURSOR_ARROW))    
        self.canvas.Update()

    def add_toolbar(self):
        self.toolbar = NavigationToolbar2Wx(self.canvas)
        self.toolbar.Realize()
        
        _tw, th = self.toolbar.GetSizeTuple()
        fw, fh = self.canvas.GetSizeTuple()

        self.toolbar.SetSize(wx.Size(fw,th))
        self.canvas.SetSize(wx.Size(fw,fh-th))

        p1,p2 = self.canvas.GetPositionTuple()
        self.canvas.SetPosition(wx.Point(p1,p2-th))

        self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.toolbar.update()
   
class eyeSpecBaseFrame (wx.Frame):
    """
    basic plot initialization
    INPUTS:
    parent: the parent plot
    id : (int) the frame id
    title : (str) Put along the top bar

    """
    def __init__ (self, parent_window, title=''):     
        if parent_window is None: ID = -1
        else: ID = parent_window.GetId()
        title = 'eyeSpec: ' + title
        #-------------------------------------------------#
        # Initiate Frame
        wx.Frame.__init__(self, parent_window, ID, title=title)
        self.SetBackgroundColour(wx.NamedColour("WHITE"))

        # what to do when it's all over
        self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)

        #-------------------------------------------------#
        # add status bar
        self.statusBar = self.CreateStatusBar(1)
        self.statusBar.SetFieldsCount(1)
        self.SetStatusBar(self.statusBar)
     
    def OnFinish (self):
        # get the final data and put into output
        final_out = None
        return final_out
        
    def OnRefresh (self,_event):
        print "aaaaaaaahhhhhh!!!!, Refreshing" #!! todo remove this 
        pass

    def OnClose (self, event):
        # wind = wx.Window(None,-1,size=(200,100))
        dlg = wx.MessageDialog(None, 'Finish and close?', style=wx.OK | wx.CANCEL)
        result = dlg.ShowModal()
        dlg.Destroy()
        if result == wx.ID_OK: self.Destroy()
        else:event.Veto()

    def OnDestroy (self,_event):
        print "="*60
        wx.GetApp().PreFinish(self.OnFinish())
        self.Refresh()

class eyeSpecBaseApp (wx.App):
    """ Basic application for wx"""
    def __init__ (self,FrameClass,inputs,redirect=True):        
        self.FrameClass = FrameClass
        self.inputs = inputs
        self.redirect = redirect # send text to window
        
        redirect = False 
        filename = None
        wx.App.__init__(self, redirect, filename)

        self.Bind(wx.EVT_WINDOW_DESTROY,self._exit)

    def _exit (self,event):
        if event.GetWindow().GetId() != self.window.GetId(): return
        self.ExitMainLoop()

    def PreFinish (self,final_out):
        self.final_out = final_out
        
    def Finish (self):
        if 'final_out' not in dir(self):
            raise ValueError("Need to PreFinish before you run finish")
        return self.final_out
    
    def OnInit (self):
        self.window = wx.Window(None)
        self.window.SetId(10)
        
        self.frame = self.FrameClass(self.window,self.inputs)
        self.frame.SetSize((850, 700))
        self.frame.Show(True)
        
        if self.redirect: self.redirect_text_panel = self.frame.panel.redirect_text_panel
        
        self.SetExitOnFrameDelete(True)
        return True        
        
class eyeSpecBaseEventManager (EventConnections, InputOutput):
    """ Template of a Event Manager """
    
    def __init__ (self):
        
        self._dragging = False
        self._clicked = False
        
        self._pressing = False
        self._pressing_time = time.time()
        self._pressing_dtime = 2.5  # seconds

        self.key_cfg = KeyboardConfiguration()
        # self.init_connection_callbacks(self)

        self._enter_command = False
        self._command = ''
    
    def enter_command (self):
        return self._enter_command
        
    def get_command_entry (self,_event):
        print "Use WX python to get an entry" #!! todo add this option, maybe?
        
    def key_press_callback (self,event):
        pass
    def key_release_callback (self,event):
        pass    
    def button_press_callback (self,event):
        pass
    def motion_notify_callback (self,event):
        pass
    def button_release_callback (self,event):
        pass
    
    def display_help (self,editor_name,info=None):
        print "="*60
        print "="*10+" "+str(editor_name)
        if info is not None:
            print "-"*60
            print info
            print "-"*60
        self.key_cfg.display_short_info()
        print "-"*60
        print "Continue, press 'q' to finish and quit\n"
        
    def update (self):
        pass        

pass
################################################################################

class OverplotStandardStar:
    def __init__ (self,data_dir):
        self.data_path = path_2_eyeSpec+"/eS-data/spectrum-stds/"+data_dir
        
class StandardStars:
    def __init__ (self):
        self.stars = {'sun':OverplotStandardStar("Solar/"),
                      'arcturus':OverplotStandardStar("Arcturus/")}






