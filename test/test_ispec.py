execfile("main.py")
execfile("app_iplotspec_v1.py")

# spec = readin_hst("../testfiles/hst/o5f607010_x1d.fits")
spec = readin("../testfiles/HI.20061228.29378_1_Flux.fits")
# spec = readin("../testfiles/many1d/4990-5070.svn")
x,y,invar = readin("../testfiles/HI.20061228.29378_1_Flux.fits",get_data=True)
y *= .90
x = x[0]+10.0
y = y[0]


linelist = np.loadtxt("/Users/dylangregersen/Desktop/Astrophysics/data/allen/LINELISTS/lin.long.101812",unpack=True,usecols=[0,1])

smt1 = np.loadtxt("../../smart/tests/data_smt.txt/smt",unpack=True,skiprows=2)


add_xy = [(x[0],y[0])]
for i in xrange(1,len(x)):
    add_xy.append((x[i],y[i]))

xy_mplkwargs = [{'color':'r','linestyle':'-.','lw':2}]

#y = iplotspec(spec,llist,add_xy,xy_mplkwargs)
