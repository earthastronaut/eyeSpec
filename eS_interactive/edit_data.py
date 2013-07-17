from ._core import ( History, SysOutListener,
                     DataSelectionBox, OrderSelectionBox, eyeSpecBaseDataPlot,
                     eyeSpecBaseMainPanel,eyeSpecBaseDataPanel,
                     eyeSpecBaseFrame, eyeSpecBaseApp, eyeSpecBaseEventManager)                                     
import pdb #@UnusedImport
from ..dependencies import (np, os, sys, time, deepcopy, math)
from ..io import save_txt

################################################################################
################################################################################

def _app_edit_data_what_happened ():
    # General
    f = open('TMP_WHAT_JUST_HAPPENED.txt', 'w')
    f.write("This file is intended to help in the case that the code crashes and you're left with only TMP_* files\n")
    f.write("If this isn't helpful please email Dylan Gregersen <dylan.gregersen@utah.edu>\n")
    f.write("\n")

    # What was saved? Why?
    f.write("This file was created from the set overlap routine\n")
    f.write("The original spectrum information is in the file: TMP_OBJ_SAVE_ORIG.pkl\n")
    f.write("The edited spectrum information is in the file: TMP_OBJ_SAVE_EDIT.pkl\n")
    f.write("\n")
    
    # How to recover with these files.
    f.write("To restore what you were working on:\n")
    f.write(">>> orig_obj = load_obj('TMP_OBJ_SAVE_ORIG.pkl')\n")
    f.write(">>> spec_obj = load_obj('TMP_OBJ_SAVE_EDIT.pkl')\n")

    f.close()

################################################################################

class EditDataPanel (eyeSpecBaseDataPanel):
    
    def __init__ (self,edit_data_main_panel,edit_data_frame):
        parent_panel = edit_data_main_panel
        parent_frame = edit_data_frame
        
        eyeSpecBaseDataPanel.__init__(self, parent_panel, parent_frame)
        self.spec_obj = edit_data_frame.spec_obj
        
        #-------------------------------------------------#
        # add editors
        self.edManager = EditDataManager(self, self.spec_obj)  # InteractiveDataEditor(spec_obj,parent=self)
        self.edManager.disconnect()
        
        del self.canvas.callbacks.callbacks['motion_notify_event'][self.statusbar_cid]
        self.canvas.mpl_connect('motion_notify_event',self.DataUpdateStatusBar)
        
    def DataUpdateStatusBar (self,event):
        scale_txt = "Auto Scale "
        scale_opt,_ = self.edManager.ide.dp.get_auto_scaling_opt()
        if scale_opt == 0: scale_txt += 'X,Y'
        elif scale_opt == 1: scale_txt += 'X'
        elif scale_opt == 2: scale_txt += 'Y'
        elif scale_opt == 3: scale_txt += 'None'
        
                
        current_order = self.edManager.ide.get_order_index()
        ord_txt = "Order: "
        if current_order is None: ord_txt += 'None'
        else: ord_txt += str(current_order)+"/" + str(self.spec_obj.shape[1] - 1)
        
        st = format(scale_txt,'16')+" |  "+format(ord_txt,'17')
        
        self.UpdateStatusBar(event,st)
        
    def OnStart (self, event):
        super(EditDataPanel,self).OnStart(event)
        
        del self.canvas.callbacks.callbacks['key_press_event'][self._onkeystart_cid]
        del self.canvas.callbacks.callbacks['button_press_event'][self._onbutstart_cid]

        
        
        print ""
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "QUESTIONS TO USER:"
        questions = ["o  should the shift in wavelength be relative or fixed which arrowing, currently fixed at 1 Angstrom for every change, should it be 5% of the range?",
                     "o  I had to make some choices about how the window scales as  you move through, any suggestions?",
                     "o  How do you like the initial starting? any suggestions?",
                     "o  Let me know if you have any problems with freezing after  you close the window (dylan.gregersen@utah.edu), please note whether you exited by closing the window or pressing 'q'",
                     "--"*20]
        print ("\n".join(questions))

        print ""
        
        first_order = self.edManager.ide.dp.plot_data[0].get_xdata()

        xmin, xmax = np.min(first_order), np.max(first_order)
        ran = (xmax - xmin)
        self.ax.set_xlim(xmin + .1 * ran, xmax + .1 * ran)  # semi-arbitrary starting point

        self.edManager.connect()
        self.edManager.update()

class EditDataMainPanel (eyeSpecBaseMainPanel):
    def __init__ (self,parent_frame):
        eyeSpecBaseMainPanel.__init__(self, parent_frame,'split_top')
        
        # define top panel 
        self.datapanel = EditDataPanel(self.Split0,self.pframe)
        self.canvas = self.datapanel.canvas
        self.split_top(self.datapanel)
         
class InteractiveDataEditor:
    """
    
    """
    def __init__ (self, ax, spec_obj, auto_scale_focus='edges'):
        # WHICH DO I WANT:
        # DataPlot.__init__(stuff)
        # self.dataplot = DataPlot(stuff)
        # self.spec_obj = self.dataplot.spec_obj
        # self.ax = self.dataplot.ax
        self.dataplot = self.dp = eyeSpecBaseDataPlot(ax, spec_obj, auto_scale_focus=auto_scale_focus)

        self.spec_obj = self.dataplot.spec_obj
        self.ax = self.dataplot.ax
        
        min_zorder, max_zorder = self.dataplot.plot_data.get_zorder_bounds()
        self._selection_pts, = self.ax.plot(self.spec_obj.get_wl(0), self.spec_obj.get_data(0), linestyle='none', marker='o', markersize=4, color='r', visible=False, zorder=max_zorder + 5)

        self._allow_edit = True

        self._previ = 0
        self._ordi = None

        self.history = History()
        self.history.set_options('edit data', self.undo_redo_edit_data)
        self.history.add('edit data', self.dataplot.get_current_state_data(), info='original data')

        self._pressed = {'x':False}


        self.order_box = OrderSelectionBox(self.ax) 
        self.selection_box = DataSelectionBox(self.ax)

        self.selection_box.set_zorder(min_zorder - 2)
        self.order_box.set_zorder(min_zorder - 3)
        self._allow_selection_box = True
        
        
        self._moving_right = False

        self._recorded_pts = []
        self._recorded_lim = 20

        # information
        self.text = self.ax.text(0.05, 1.05, 'Order Selected: None',
                                 transform=self.ax.transAxes, va='top')

    def add_to_recorded_pts (self, x, y):
        x, y = float(x), float(y)
        self._recorded_pts.append([x, y])
        rl = int(self._recorded_lim)
        diff = len(self._recorded_pts) - rl
        if diff > 0: del self._recorded_pts[:diff]

    def get_recorded_xy (self,key):
        return self._recorded_pts[key]

    def get_order_index (self):
        return deepcopy(self._ordi)

    def get_prev_order_index (self):
        return deepcopy(self._previ)

    def set_allow_selection_box (self, truth):
        truth = bool(truth)
        self._allow_selection_box = truth

    def get_allow_selection_box (self):
        return deepcopy(self._allow_selection_box)

    def set_allow_edit (self, truth):
        self._allow_edit = bool(truth)

    def get_allow_edit (self):
        return deepcopy(self._allow_edit)

    def save_params_data_editor (self, filename, clobber=True):
        filename = str(filename)
        if os.path.exists(filename) and not clobber: raise IOError("File exists: " + filename)
        f = open(filename, 'w')
        lines = self._current_params_data()
        f.write("\n".join(lines))
        f.close()

    def open_params_data_editor (self, filename):
        filename = str(filename)
        if not os.path.exists(filename): raise IOError("File does not exist:" + filename)
        f = open(filename)
        lines = f.readlines()
        self._apply_params_data_edit(lines)
        f.close()                
   
    def _current_params_data_edit (self): 
        print "!! need to update"
        return ['one line only']

    def _apply_params_data_edit (self, lines):
        print "!! need to update", lines
        self.update()


    def undo_redo_edit_data (self, prev_state):
        current_state = self.dataplot.get_current_state_data()
        self.dataplot.apply_prev_state_data(prev_state)  
        self.spec_obj = self.dataplot.spec_obj      
        return current_state

    def deselect_selection (self):
        self.selection_box.set_visible(False)
        self._selection_pts.set_visible(False)

    def deselect_all (self):
        if self._ordi is not None: self._previ = self.get_order_index()
        self._ordi = None
        self.order_box.set_visible(False)
        self.deselect_selection()
        self.dataplot.plot_data.reset_order_zorder()

    def force_select_order (self):
        if self._ordi is None: self._ordi = self.get_prev_order_index()
        self.select_order_by_index(self._ordi)

    def _select_order_by_plot_divisions (self, xpt):
        
        xmin, xmax, ymin, ymax = self.ax.axis()
     
        pd = self.dataplot.plot_data
        orders = {}
        for i in pd:
            xdata = pd[i].get_xdata() 
            ydata = pd[i].get_ydata()
            if xdata is None or ydata is None or len(xdata) == 0: continue
            xmask = (xmin < xdata) * (xdata < xmax)
            ymask = (ymin < ydata) * (ydata < ymax)
            mask = xmask * ymask
            if np.any(mask):
                orders[i] = [np.min(xdata), np.max(xdata)]
        
        if len(orders) == None: return None
        ind = np.array(orders.keys())
        xmins = np.array(orders.values()).T[0]
        if xmins.ndim != 1: raise StandardError("Whoops, check the coding:" + str(xmins))
        sorti = np.argsort(xmins)
        ind = ind[sorti]
        
        bins = np.linspace(xmin, xmax, len(ind) + 1)
        ordi = None
        for i in range(len(bins) - 1):
            if xpt >= bins[i]:
                ordi = int(ind[i])
        return ordi

    def _min_dist_to_data_xy (self, xpt, ypt):
        all_dist = {}
        pd = self.dataplot.plot_data
        eps_xy = self.dataplot.get_epsilon('xy')

        for i in pd:
            xdata = pd[i].get_xdata() 
            ydata = pd[i].get_ydata()
            
            if xdata is None or ydata is None or len(xdata) == 0: continue
            
            x_dist = np.min(np.abs(xdata - xpt))
            y_dist = np.min(np.abs(ydata - ypt))

            if x_dist < eps_xy[0] and y_dist < eps_xy[1]:
                the_dist = (math.sqrt(x_dist ** 2 + y_dist ** 2))
                all_dist[the_dist] = [i, x_dist, y_dist]
                
        if len(all_dist) == 0: return None, None


        min_dist = np.min(np.array(all_dist.keys()))
        min_i = all_dist[min_dist][0]
        min_i = int(min_i)

        return min_dist, min_i
  
    def _select_order (self, index, direction=None):
        # turn off previous
        if index is None or index not in self.dataplot.plot_data: 
            self.deselect_all()
            return
        self._ordi = index  
        self.dataplot.plot_data.bring_order_to_top(index)
        
        if direction is not None:
            self._scanning_auto_focus(direction)
        
        pltdata = self.dataplot.plot_data[index]
        self.order_box.update_box(pltdata, index)

    def select_order_by_index (self, index=None, direction=None, always_find=False):
        if index is None:
            if always_find: 
                self.force_select_order()
                index = deepcopy(self._ordi)
            else: 
                self.deselect_all()
                return

        #====== Find the next order which hasn't been deleted
        given_direction = (direction in ['-','+'])
        if not given_direction: direction = '+'
        
        # look through in one direction starting at index point
        found_non_deleted, _i = self._find_non_deleted_orders(index, direction) 
        if not found_non_deleted:
            other_direct = '-'
            if direction == '-': other_direct = '+' 
            found_non_deleted, _i = self._find_non_deleted_orders(index, other_direct)
            if not found_non_deleted: raise StandardError("Whoops, looks like all orders were deleted!")   
             
        if given_direction: self._select_order(index,direction)            
        else: self._select_order(index)
                  
    def select_order_by_xy (self, xypts, also_by_divisions=False):
        xpt, ypt = xypts
        if xpt is np.NaN or ypt is np.NaN:
            print "Whoops, received a nan:", xypts
            self.deselect_all()
            return False
        
        # choose based on distance to point
        _min_dist, index = self._min_dist_to_data_xy(xpt, ypt)
        
        if index is None:
            # if no points are close enough
            if also_by_divisions: 
                # check divisions
                index = self._select_order_by_plot_divisions(xpt, ypt)
            
        self.select_order_by_index(index)

    def select_order_in_plot (self, which_edge, near_to_xpt=None):
        xmin = self.ax.axis()[0]
        xmax = self.ax.axis()[1]

        if near_to_xpt is None: mid = (xmin + xmax) / 2.
        else: mid = float(near_to_xpt)
        
        # find the closest order to the mid depending on the direction you gave
        closest_i = None
        closest_d = 1.e6

        for i in self.plot_data:
            if len(self.plot_data[i].get_xdata()) == 0: continue
            dat_xmax = np.max(self.plot_data[i].get_xdata())
            dat_xmin = np.min(self.plot_data[i].get_xdata())
            xdat = self.plot_data[i].get_xdata()
            
            prev_d = deepcopy(closest_d)
            if which_edge == 'right' and  xmin < dat_xmax < xmax:
                closest_d = min(((dat_xmax - mid), closest_d))
            elif which_edge == 'left' and xmin < dat_xmin < xmax:
                closest_d = min(((dat_xmin - mid), closest_d))
            elif which_edge == 'center' and np.any((xdat > xmin) * (xdat < xmax)):
                closest_d = min((np.min(xdat - mid), closest_d))

            if closest_d != prev_d: closest_i = i

        # if a closest was found then select it
        if closest_i is not None:

            if self._ordi is not None: self._previ = deepcopy(self._ordi)
            self._ordi = np.clip(closest_i, 0, self.spec_obj.shape[1] - 1)                

            self.selection_box.clear_box_variables()
            if which_edge == 'right':
                self.selection_box.create_from_right(True)
            elif which_edge == 'left':
                self.selection_box.create_from_left(True)
  
            self.update_order_box()

    def _find_non_deleted_orders (self, index, direction):
        if direction not in ['-', '+']: return (False,0)
        # find the next highest order which hasn't been deleted
        index = int(index)
        if index == self.spec_obj.shape[1]: index = self.spec_obj.shape[1] - 1
        
        if direction == '+':
            ran = xrange(index, self.spec_obj.shape[1])
        # find the next lowest order which hasn't been deleted
        elif direction == '-':
            ran = reversed(xrange(index + 1))
            
        for i in ran:
            if len(self.spec_obj.get_wl(i)) == 0: continue
            return (True, deepcopy(i))
        if True: return (False, 0)
     
    def _scanning_auto_focus (self, direction): 
        if self.dataplot.get_auto_scale_focus().find('edge') == -1: return
        if self._ordi is None: return
        
        index = self._ordi
        
        xmin, ymin = self.dataplot.spec_obj.get_min(index)
        xmax, ymax = self.dataplot.spec_obj.get_max(index)
        xs1, xs2, ys1, ys2 = self.selection_box.modify_bounds([xmin, xmax, ymin, ymax])    
        if direction == '-':
            self.dataplot.set_auto_scale_focus('edge_left')
            self.selection_box.create_from_left(True)
            self.selection_box.set_box_start([xs1, xmin, ys1, ys2], True)                 
                
        elif direction == '+':
            self.dataplot.set_auto_scale_focus('edge_right')
            self.selection_box.create_from_right(True)
            self.selection_box.set_box_start([xmax, xs2, ys1, ys2], True)
            
        self.dataplot.scale_x(index)
             
    def scan_through_orders (self, direction, always_find=True):
        if direction not in ['-', '+', None]: return

        if self._ordi is None:
            if always_find: self._ordi = deepcopy(self._previ)
            else: return

        if direction == '-': inc = -1
        elif direction == '+': inc = 1

        # if you are focusing on edges then for a particular order
        # you want to be able to switch which side before switching 
        # orders
        if self.dataplot.get_auto_scale_focus().find('edge') != -1:
            if direction == '-':
                if self._moving_right:
                    self._moving_right = False
                    inc = 0
                    
            elif direction == '+':
                if not self._moving_right:           
                    self._moving_right = True
                    inc = 0       
                
        index = self._ordi + inc
        self.select_order_by_index(index, direction, True)
  
    def btn_press_selection_box (self, event):
        """ for data editor """
        # returns True if something changed
        # if no order is selected then check then select order
        if event.button != 1: return False
        if self._ordi is None:
            self.dataplot.skip_auto_scale()
            self.select_order_by_xy((event.xdata, event.ydata))
            return True

        # if order is deleted don't update        
        wl = self.spec_obj.get_wl(self._ordi)
        if len(wl) == 0: return False

        # if you're not allowed to have a selection box skip        
        if not self._allow_selection_box: 
            self.deselect_selection()
            return True           


        xmin, ymin = self.spec_obj.get_min(self._ordi)
        xmax, ymax = self.spec_obj.get_max(self._ordi)
        data_bounds = [xmin, xmax, ymin, ymax]
        should_update = self.selection_box.during_btn_press(event, self.dataplot.get_epsilon('x'), data_bounds)
        return should_update

    def mot_notify_selection_box (self, event):
        """ for data editor """        
        # returns true if edited
        if event.button != 1: return False
        if self._ordi is None: return False
        wl = self.spec_obj.get_wl(self._ordi)
        if len(wl) == 0: return False
        
        if not self.get_allow_selection_box():
            self.deselect_selection() 
            return False
        
        if len(self._recorded_pts) == 0: return False
    
        start_pt = self._recorded_pts[-1]

        xmin, ymin = self.spec_obj.get_min(self._ordi)
        xmax, ymax = self.spec_obj.get_max(self._ordi)
        data_bounds = [xmin, xmax, ymin, ymax]
        
        changed = self.selection_box.during_mot_notify(event, start_pt, data_bounds, self._pressed)
        return changed
        
    def btn_click_data_editor (self, event):
        """ for data editor """        
        # returns true if edited
        # print "\n button release"
        # self.selection_box.display_box_variables()

        if not self.get_allow_selection_box(): 
            self.deselect_selection()
            return True       

        if event.button != 1: return False
        
        changed = False
        if self.selection_box.get_visible():
            changed = self.selection_box.during_btn_click(event)    
        else:
            self.selection_box.clear_box_variables()
            self.dataplot.skip_auto_scale()
            self.select_order_by_xy((event.xdata, event.ydata))
            changed = True
            
            
        if not self.selection_box.get_visible():
            self._selection_pts.set_visible(False)
        self.selection_box.clear_box_variables()
        return changed 
    
    def key_held_callback (self, event, which):
        """ 
        controls key held events
        Should be done early in a key_press_callback and returns boolean whether to continue onto other key press events
        """

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

    def key_release_selection_box (self, event):
        eps_x = self.dataplot.get_epsilon('x')
        delta_x = 0.8 * eps_x
        self.selection_box.during_key_release(event, delta_x)

    def delete_data_selection (self):
        if not self.selection_box.get_visible(): return
        box_bounds = self.selection_box.get_bounds()
        self.delete_data_in_bounds(box_bounds)

    def delete_data_in_bounds (self, box_bounds):
        if self._ordi is None: return
        if not self._allow_edit: return

        xmin, xmax, ymin, ymax = self.ax.axis()
        # check if the box bounds are visible on the plot
        if (box_bounds[0] > xmax and (box_bounds[2] > ymax or box_bounds[3] < ymin))\
                or (box_bounds[1] < xmin and (box_bounds[2] > ymax or box_bounds[3] < ymin)):
            print "Box outside of current window. Can't delete points"
            return
        
        bb = box_bounds
        self.history.add('edit data', self.dataplot.get_current_state_data(), info='Cropped data ' + str((bb[0], bb[1], bb[2], bb[3])))
        # !! might need some prechecks so that obj.edit.crop doesn't raise errors
        # I know it'll be ok with self._order_sel_i, and box_bounds is given/ndarray/float
        # but will it find points, if it doesn't how will it complain? probably badly
        self.spec_obj.edit.crop(order=deepcopy(self._ordi), include=True, box_bounds=box_bounds)
        self.dataplot.plot_data.update_plot_data(self.spec_obj)

    def translate_data (self,to_wl=0,to_data=0):
        self.history.add('edit data', self.dataplot.get_current_state_data(), info='translate data +(delx, dely) = '+str((to_wl,to_data)))
        self.spec_obj.edit.translate(add_to_wl=to_wl,add_to_data=to_data)
        self.dataplot.plot_data.update_plot_data(self.spec_obj)
       
    def scale_data (self,to_wl=1, to_data=1):
        self.history.add('edit data', self.dataplot.get_current_state_data(), info='scale the data *(delx, dely) = '+str((to_wl,to_data)))
        self.spec_obj.edit.scale(mult_to_wl=to_wl,mult_to_data=to_data)
        self.dataplot.plot_data.update_plot_data(self.spec_obj)
        
    def radial_velocity_shift (self,rv=1):
        self.history.add('edit data', self.dataplot.get_current_state_data(), info='applying radial velocity shift to data = '+str(rv))
        self.spec_obj.edit.apply_rv(rv)
        self.dataplot.plot_data.update_plot_data(self.spec_obj)        
                        
    
    ##################################################

    def _add_2_toolbar (self):pass
     
    def update_highlight_selected_points (self, box_bounds=None):
        if self._ordi is None: return False
        
        # if True:
        #    self._selection_pts.set_visible(False)
        #    return
        
        
        if box_bounds is None:
            # check for the selection box
            if not self.selection_box.get_visible():
                self._selection_pts.set_visible(False)
                return  False
            else: box_bounds = self.selection_box.get_bounds()

        self._selection_pts.set_visible(True)
             
        wl = self.spec_obj.get_wl(self._ordi)
        data = self.spec_obj.get_data(self._ordi)

        xmask = (wl > box_bounds[0]) * (wl < box_bounds[1])
        ymask = (data > box_bounds[2]) * (data < box_bounds[3])
        mask = xmask * ymask

        max_zorder = self.dataplot.plot_data.get_zorder_bounds()[1]
        self._selection_pts.set_zorder(max_zorder)
        self._selection_pts.set_xdata(wl[mask])
        self._selection_pts.set_ydata(data[mask])
        return True

    def update_text (self): 
        # update visual selection
        prev_text = self.text.get_text()
        dataind = str(deepcopy(self._ordi))
        display_text = 'Order Selected: ' + dataind
        if self._ordi != None: 
            display_text += "/" + str(self.spec_obj.shape[1] - 1)
        
        if prev_text == display_text: return False
        self.text.set_text(display_text)
        return True
        
    def update_order_box (self):
        if self._ordi is None: pltdata = None
        else: pltdata = self.dataplot.plot_data[self._ordi]
        self.order_box.update_box(pltdata, self._ordi)

    def update (self):   
        self.spec_obj = self.dataplot.spec_obj
        
        c1 = self.update_highlight_selected_points()
        c2 = self.update_text()
        
        # update_highlight_ps
        c3 = self.dataplot.update(self._ordi)
        c4 = False
        if self.dataplot.bounds_changed(True):
            self.update_order_box()
            c4 = True
        return (c1 or c2 or c3 or c4)

class EditDataFrame (eyeSpecBaseFrame):
    def __init__ (self, parent_window, inputs):
        """
        inputs must be a single eyeSpec spectrum object
        """
        
        self.spec_obj = inputs
        
        title = 'Edit Data: '+os.path.basename(self.spec_obj.filename)
        eyeSpecBaseFrame.__init__(self, parent_window, title)
    
        self.panel = EditDataMainPanel(self)    

    def OnFinish (self):
        self.Backup()
        spec_obj = self.panel.datapanel.edManager.ide.dp.spec_obj
        return spec_obj
       
    def Backup (self):
        print "Closing: If this hangs up look at the file TMP_WHAT_JUST_HAPPENED.txt"
        _app_edit_data_what_happened()
        spec_obj = self.panel.datapanel.edManager.ide.dp.spec_obj
        save_txt(spec_obj,filename='TMP_OBJ_SAVE_EDIT.pkl',clobber=True)
        time.sleep(.5)        

class EditDataManager (eyeSpecBaseEventManager):
    def __init__ (self, edit_data_panel, spec_obj):
        eyeSpecBaseEventManager.__init__(self)
        parent_panel = edit_data_panel
        self.ppanel = parent_panel
        self.ax = parent_panel.ax
        
        self.ide = InteractiveDataEditor(self.ax, spec_obj, 'edges')
        
        self.init_connection_callbacks(self)


        self.key_cfg.add_key('d',"Same as 'backspace")
        self.key_cfg.add_key('h','Display this screen')
        self.key_cfg.add_key('g','Toggle grid on/off')
        self.key_cfg.add_key('p',"Toggle pan/zoom tool (note: editor won't work while engaged")
        self.key_cfg.add_key('z',"Toggle zoom rect tool (note: editor won't work while engaged")
        self.key_cfg.add_key('q','Close and return')
        self.key_cfg.add_key(";",'Toggle data between scatter and line plot options')
        self.key_cfg.add_key('`','Toggle auto scaling options')
        self.key_cfg.add_key('[','Undo data edit')
        self.key_cfg.add_key(']','Redo data edit')
        self.key_cfg.add_key('< and >','Scan through orders')
        self.key_cfg.add_key('left/right','Edit box from the edge of order')
        self.key_cfg.add_key('backspace','Delete Selected Data')
        self.key_cfg.add_key('esc',"same as 'q'")
        
        self.key_cfg.set_display_order(['backspace','d','esc','g','h','p','q','z','`',';','[',']',
                                        '< and >','left/right'])
        
        self.key_cfg.check_display()
        self.key_cfg.add_mpl_key_convert('left','left/right')
        self.key_cfg.add_mpl_key_convert('right','left/right')
        self.key_cfg.add_mpl_key_convert('<','< and >')
        self.key_cfg.add_mpl_key_convert('>','< and >')
        
    def key_press_callback (self, event):
        self._pressing_time = time.time()
        self._pressing = True
        
        keypress = self.ide.key_held_callback(event, 'press')
        if keypress: return           
        # self.key_held_callback(event)
    
    def key_release_callback (self, event):
        self._pressing = False
        self._key_call(event)
 
    def _key_call (self, event):
        
        if event.key == 'backspace': event.key = 'd'
        if event.key == 'esc': event.key = 'q'

        keypress = self.ide.key_held_callback(event, 'release')
        if keypress: return
        
        if event.key == 'd':
            self.ide.delete_data_selection()
        elif event.key == 'h':
            self.display_help('Data Editor')
        elif event.key == 'g':
            self.ide.dp.toggle_grid()
        elif event.key == 'p':
            self.ide.dp.toggle_pan_tool()
        elif event.key == 'q':
            if 'Close' in dir(self.ppanel.pframe): self.ppanel.pframe.Close()
        elif event.key == 'z':
            self.ide.dp.toggle_zoom_tool()        
        elif event.key == '`':
            self.ide.dp.set_auto_scaling() 
        elif event.key == ';':
            self.ide.dp.plot_data.toggle_data_line_scatter()
        elif event.key == '[':
            self.ide.history.undo()
        elif event.key == ']':
            self.ide.history.redo()
        elif event.key == ',': 
            self.ide.scan_through_orders("-")
        elif event.key == '.': 
            self.ide.scan_through_orders("+")
        elif self.ide.get_allow_edit() and self.ide.get_order_index() is not None:
            self.ide.key_release_selection_box(event)
            
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
        if self.ide.dp.is_toolbar_button_on():
            self.ide.dp.skip_auto_scale()
            self.plot_change_update()
            return
        
        if event.button == 1:
            self.ide.add_to_recorded_pts(event.xdata, event.ydata)
            
        self._dragging = False
        self._clicked = True
        
        if self.ide.btn_press_selection_box(event): self.update()
        
    def motion_notify_callback (self, event):
        """ Data Edit Manager  """
        if self.ide.dp.is_toolbar_button_on():
            self.ide.dp.skip_auto_scale()
            if self.ide.dp.timed_update(0.02):
                self.plot_change_update()
            return
        
        self._dragging = True
        self._clicked = False
                
        if self.ide.mot_notify_selection_box(event): self.update()
                  
    def button_release_callback (self, event):
        """ Data Edit Manager  """
        if self.ide.dp.is_toolbar_button_on():
            self.ide.dp.skip_auto_scale()
            self.plot_change_update()
            return

        changed = False
        # !! selection box
        if self._clicked:
            changed = self.ide.btn_click_data_editor(event)
    
        self._dragging = False  
        self._clicked = False      
        if changed: self.update()
        
    def plot_change_update (self):
        bounds = self.ide.ax.axis()
        if "_prev_bounds" not in dir(self): self._prev_bounds = bounds
        if self._prev_bounds != bounds: 
            self.update()
            self._prev_bounds = bounds       
        
    def update (self):
        self.ide.update()
        self.ide.dp.bounds_changed(True)
        self.ax.figure.canvas.draw()
    

############################################################################

def edit_data (spec, clean_up=True):
    """
    This performs the basic data editing including deleting points

    INPUTS:
    ============  ==============================================================
    keyword       (type) Description
    ============  ==============================================================
    spec_obj      (eyeSpec_spec) A eyeSpec spectrum object
    clean_up      (bool) If true then it will remove temporary files it creates
    ============  ==============================================================
    """
    # !! should also let you perform radial velocity shifts, translations, scalings, normalization(base on plot)

    if spec.__class__.__name__ != 'eyeSpec_spec': raise ValueError("spec MUST BE OF CLASS eyeSpec_spec")
        
    edit_spec = spec.copy()
    # save_spec(spec,filename='TMP_OBJ_SAVE_ORIG',clobber=True)

    ##########################################
    # run application
    app = eyeSpecBaseApp(EditDataFrame,edit_spec)
    sys.stdout = SysOutListener()
    try: app.MainLoop()
    finally:
        app.ExitMainLoop()
        final_spec = app.Finish()
        del app

    print "-"*60
    print "-"*20+format("Edit Data Complete",'^26')+"-"*20

    # load data after app.MainLoop() has exited
    if os.path.exists('TMP_OBJ_SAVE_EDIT.pkl'): pass
        #backup_spec = load_spec('TMP_OBJ_SAVE_EDIT.pkl')

    # clean up temporary files
    if bool(clean_up):
        if os.path.exists('TMP_WHAT_JUST_HAPPENED.txt'): os.system('rm TMP_WHAT_JUST_HAPPENED.txt')
        if os.path.exists('TMP_OBJ_SAVE_ORIG.pkl'): os.system('rm TMP_OBJ_SAVE_ORIG.pkl')
        if os.path.exists('TMP_OBJ_SAVE_EDIT.pkl'): os.system('rm TMP_OBJ_SAVE_EDIT.pkl')
        
    return final_spec
    
