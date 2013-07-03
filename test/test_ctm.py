execfile('main.py')
execfile("../app_edit_ctm.py")

print "==== data ===="
x = readin("../../test/HI.20061228.29378_1_Flux.fits")

#x = readin('../testfiles/many_1d/4600-4670.svn')

# y = edit_ctm(x)



















# wl = x.get_wl(0)#np.arange(5000,5510,.65)
# data = x.get_data(0)#np.random.randn(len(fake_wl))/10.0+1.0
# inv_var = x.get_inv_var(0)# 1.0/(np.random.randn(len(fake_wl))/100.0)


# ax = plt.figure().add_subplot(111)
# ax.plot(wl,data,alpha=.6,zorder=0)


# continuum_length = 10.0

# guess_wl = np.arange(np.min(wl)+10.0,np.max(wl)-10.0,20)
# def get_guess_data (continuum_length):
#     guess_data = np.ones(len(guess_wl))
#     for i in range(len(guess_wl)):
#         mask = (wl > guess_wl[i]-continuum_length/2.0)*(wl < guess_wl[i]+continuum_length/2.0)
#         histo = np.histogram(data[mask],bins=30)
#         guess_data[i] = histo[1][histo[0].argmax()]

#     return guess_data


# guess_data = get_guess_data(10)
# guess_data1 = get_guess_data(50)


# dat_pts, = ax.plot(guess_wl,guess_data,marker = 'o',linestyle='none',color='r',zorder=5,alpha=.5)

# dat_pts1, = ax.plot(guess_wl,guess_data1,marker = '*',linestyle='none',color='r',zorder=5)

# add_dat = [dat_pts]




# #Y = scipy.interpolate.UnivariateSpline(wl,data,w=inv_var)
# Y = scipy.interpolate.LSQUnivariateSpline(wl,data,guess_wl,inv_var,k=3)
# X = np.arange(wl[4],wl[-4],.1)

# spline_line, = ax.plot(X,Y(X),color='g',lw = 3,zorder=1)



# def add_ctm_point (pltplot,xy):
#     # get the current data
#     current_x = pltplot.get_xdata()
#     current_y = pltplot.get_ydata()

#     # concatenate the new x and y data
#     new_x = np.concatenate((current_x,np.array([xy[0]])))
#     new_y = np.concatenate((current_y,np.array([xy[1]])))

#     # sort all the data by wavelength
#     sort_i = new_x.argsort()
#     new_x = new_x[sort_i]
#     new_y = new_y[sort_i]

#     # replace the current data with the new stuff
#     pltplot.set_xdata(new_x)
#     pltplot.set_ydata(new_y)
    
#     return pltplot


