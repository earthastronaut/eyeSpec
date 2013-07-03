execfile("main.py")
print "--- read data ---"
x = readin('../testfiles/multispectest.fits')


ord10 = x[10]
ord11 = x[11]

# mask = ord10[0]> 5531.
# part10 = np.vstack((ord10[0][mask],ord10[1][mask]))




# mask = ord11[0] < 5541.
# part11 = np.vstack((ord11[0][mask],ord11[1][mask]))






overlap_pt = 5535.2

new_i = np.abs(new_wl - overlap_pt).argmin()
cur_i = np.abs(cur_wl - overlap_pt).argmin()


# because already cropped this will be appropriate
new_a = np.abs(new_wl-np.min(cur_wl)).argmin()
cur_a = 0 

new_b = len(new_wl)-1
cur_b = np.abs(cur_wl - np.max(new_wl)).argmin()


new_wl_ai = new_wl[new_a:new_i] # not including i
cur_wl_ai = cur_wl[cur_a:new_i] # not including i
        #!! print "wl_ai",len(new_wl_ai),len(cur_wl_ai), "< don't need to be the same, only care about new_wl_ai (#1)"

new_data_ai = new_data[new_a:new_i] # not including i
cur_data_ai = cur_data[cur_a:new_i] # not including i
        
new_inv_var_ai = new_inv_var[new_a:new_i] # not including i
cur_inv_var_ai = cur_inv_var[cur_a:new_i] # not including i
        #!! print "inv_var_ai",len(new_inv_var_ai),len(cur_inv_var_ai)
        
# segment [i to b]
new_wl_ib = new_wl[new_i:-1] # including i
cur_wl_ib = cur_wl[cur_i:cur_b+1] # including b and i
        #!! print "wl_ib",len(new_wl_ib),len(cur_wl_ib), "< don't need to be the same, only care about cur_wl_ib (#2)"

new_data_ib = new_data[new_i:-1] # including i
cur_data_ib = cur_data[cur_i:cur_b+1] # including b and i

new_inv_var_ib = new_inv_var[new_i:-1] # including i
cur_inv_var_ib = cur_inv_var[cur_i:cur_b+1] # including b and i
        #!! print "inv_var_ib",len(new_inv_var_ib),len(cur_inv_var_ib)
        
