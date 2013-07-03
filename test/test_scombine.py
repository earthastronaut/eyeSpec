execfile("main.py")
execfile("base_combine_v1.py")

#x = readin('../testfiles/multispectest.fits')
#obj = x

# create a list of files
# filelist = ['../testfiles/many_1d/4430-4500.svn',
#             '../testfiles/many_1d/4540-4610.svn',
#             '../testfiles/many_1d/4850-4930.svn',
#             '../testfiles/many_1d/4660-4720.svn',
#             '../testfiles/many_1d/4785-4860.svn',
#             '../testfiles/many_1d/4775-4800.svn',
#             '../testfiles/many_1d/4920-5000.svn',
#             '../testfiles/many_1d/4600-4670.svn',
#             '../testfiles/many_1d/4990-5070.svn']

# you would provide a file list
filelist = '../testfiles/many_1d/filelist'

# use the file list or whatever to import data 
x = readin_spectre_files(filelist)

# set by the interactive editor
line_data = np.array([ 4490.28307543,  4547.48506984,  4607.15427563,  4666.36801856,
                       4791.84204391,  4850.90353645,  4928.89171572,  4989.76207378])


#x,line_data = setoverlap(x)


# NOW COMBINE THEM



# # for the third order cut off the first
# x.set_inv_var()[2][:100] = np.zeros(100)

# x.set_inv_var()[6] = 1/(np.random.normal(0,.001,len(x.set_inv_var()[6]))**2)
# x.set_inv_var()[7] = 1/(np.random.normal(0,.001,len(x.set_inv_var()[7]))**2)
# x.set_inv_var()[8] = 1/(np.random.normal(0,.001,len(x.set_inv_var()[8]))**2)


#nwl,ndat,nvar = scombine(x,line_data)
nobj = spec_combine(x,line_data,include_weights=True,linear_interpolate=True)
#nobj_nw = spec_combine(x,line_data,scaling_weight=False)
#nobj_nv = spec_combine(x,line_data,include_varience=False)
#nobj_nwv= spec_combine(x,line_data,scaling_weight=False,include_varience=False)


ax = plt.figure(10).add_subplot(111)
plot_spec(x,fig_num=10,alpha=.4,lw=2,zorder=-2)
#'y','k'

ax.plot(nobj[0][0],nobj[0][1],c='k',linestyle='-',alpha=1,lw=.5,zorder=0)

                                      
ymin,ymax = ax.axis()[2],ax.axis()[3]

ax.vlines(line_data,ymin,ymax,color='y',zorder=-5,linestyle='-.')

 # because the one order has small error bars
#ax.plot(nobj_nw[0][0],nobj_nw[0][1],c='b',linestyle='-',alpha=.6)
#ax.plot(nobj_nv[0][0],nobj_nv[0][1],c='g',linestyle='--',alpha=.6)
#ax.plot(nobj_nwv[0][0],nobj_nwv[0][1],c='c',linestyle='-',alpha=.6)


ord = np.vstack((nobj.get_wl(0),nobj.get_data(0),nobj.get_inv_var(0)))
f = open("DELETE_TEST.txt",'w')

for i in range(len(ord[0])):
    output = []
    for j in range(3):
        output.append(format(ord[j][i],'<20.7f'))
    f.write("  ".join(output))
    f.write("\n")

f.close()


#ax = plt.figure(10).add_subplot(111)

