# To inspect star data and normalize it with a planck function

import os
import eyeSpec as es
import numpy as np
import pdb
import matplotlib.pyplot as plt

#import pdb; pdb.set_trace()

###############################################
# Defenition Library
###############################################

# Extracting star's name along with its set #
def Star_Name(stardatafile):
    starname =stardatafile.replace('/media/solidtonyubuntu/PROJECT_HST/StarData/','')
    starname = starname.replace('/st.avg.fits','')
    starname = starname.replace('_x1d.fits','')
    starname = starname.replace('/','')
    return starname

#-------------------------------------------------#

# Smoothing
def windowed_smoothing(synth, smoothing):
    # synth assuming it has wavelength = synth[0]  flux = synth[1]
    window=np.hamming(smoothing)
    synth[1] =  np.convolve(synth[1], window, mode = 'same')/np.convolve(np.ones(len(synth[1])), window, mode = 'same')
    return synth

#-------------------------------------------------#

# Reading in Model and extracting xy values
def read_model(model_dir):
    model_dir = os.path.abspath(model_dir)+"/"    
    range_dirs = os.listdir(model_dir)
    sum_files = []
    for directory in range_dirs:
        datapath = model_dir+directory+"/"
        fnames = os.listdir(datapath)
        for fname in fnames:
            sum_files.append(datapath+fname)

    add_xy = []
    for fname in sum_files:
        xy, header = es.parse_synth_summary_out(fname)
        # smooth the synthesis
        xy = windowed_smoothing(xy,smoothing)
        # add to the list of data we're going to overplot
        add_xy.append(xy)
        
    return(add_xy)

#-------------------------------------------------#

# Data from star

# Readin data from star flux
def data_flux(stardata_file):
    spec = es.readin(stardata_file)
    data = spec.access_data()
    data = data/np.median(data)
    spec.set_data(np.array(data))
    return data

# Read in star flux divided by average
def data_flux_avg(stardata_file):
    spec = es.readin(stardata_file)
    data0 = spec.access_data()
    data0 = data0/np.median(data0)
    spec.set_data(np.array(data0))
    return data0

# Normalizing data by planck function
def data_flux_norm(stardata_file, y):
    spec = es.readin(stardata_file)
    data = spec.access_data()
    data = data/np.median(data)
    data_norm = data/y
    spec.set_data(np.array(data_norm)) 
    return data_norm

#-------------------------------------------------#

# Get Wavelength
def wavelength(stardata_file):
    spec = es.readin(stardata_file)
    wvl = spec.get_wl()
    return wvl

# Velocity Shift
def velocity_shift(v):
    spec = es.readin(stardata_file)
    wvl = spec.get_wl()
    shift =wvl*(v)*100000/c+wvl
    return shift

#-------------------------------------------------#

# Planck function
def planck_func(x,C,U,v):
    # Planck function constants and Temp for normalization
    return ((2*h*(c**2)/((x*1*10**(-8))**5))/(np.exp(h*c/((x*1*10**(-8))*k*T))-1)*C+U)

#-------------------------------------------------#

# Gathers all values for the star
def readstar(stardata_file, v, model_dir, y):

    data = data_flux(stardata_file)

    data_norm = data_flux_norm(stardata_file, y)

    shift = velocity_shift(v)

    add_xy = read_model(model_dir)
    
    return (data, data_norm, shift, add_xy)

#-------------------------------------------------#

# Save out constant valuse such as v, U, and C
def save_const_out(starname):
    f = open('radial_velocities/'+starname+'.txt','w')
    f.write('%s   %s    %s    %s' % (T, v, C, U) + '\n' + '\n')
    f.write("# Radial Velocity and Constants for star:"+starname+'\n'+'\n')
    f.write('T = %s   v = %s    C = %s    U = %s' % (T, v, C, U) + '\n' + '\n')
    f.write("###############################"+"\n"+'Terms'+'\n')
    f.write("# T is the temperature taken from the xls file for the stars""\n")
    f.write("# v is the velocity shift that best fits the synth"+"\n")
    f.write("# C is the constant for amount of flux"+"\n")
    f.write("# U is used to shift up or down the Curve fit"+'\n')
    f.close()

#-------------------------------------------------#

# Readin set velocity shift values
def velocity_shift_read(starname):
    with open('radial_velocities/%s.txt' % starname, 'r') as f:
        first_line = f.readline()
        T, v, C, U = first_line.split()
        T, v, C, U = float(T), float(v), float(C), float(U)
    return (T, v, C, U)

##############################################
# File locations of needed files from user
##############################################

stardata_file = '/media/solidtonyubuntu/PROJECT_HST/StarData/HD108317/Set1/st.avg.fits'

model_dir = '/media/solidtonyubuntu/PROJECT_HST/hst_atmo_models/hstoutphos/5100p250m235p200/'

linelist_file = '/media/solidtonyubuntu/PROJECT_HST/StarInspection/hst_inspection_linelistb.ln'

smoothing = 13.0

##############################################
# Constants/Global variables
##############################################

starname = Star_Name(stardata_file)
print('starname is ' + starname)

# Constants
c = 2.99792458*10**10
k = 1.38066*10**(-16)
h = 6.62608*10**(-27)

# Star information

if os.path.isfile('radial_velocities/%s.txt' % starname):
    T, v, C, U = velocity_shift_read(starname)
    Ti, vi, Ci, Ui = velocity_shift_read(starname)
    print(T, v, C, U)
else:
    Ti = T = float(raw_input('Input best guess for Teff (example: 5100): '))
    vi = v = float(raw_input('Input best guess for velocity (example: -97.8): '))
    Ci = C = float(raw_input('Input best guess for constant C (example: 4.7619047619e-14): '))
    Ui = U = float(raw_input('Input best guess for constant U (example: 0.37): '))
    print(T, v, C, U)

# Plotting information
x = wavelength(stardata_file)
y = planck_func(x,C,U,v)

################################################
# Normalization and more
################################################

data, data_norm, shift, add_xy = readstar(stardata_file, v, model_dir, y)

################################################
# Event handeling to move curve and shift v
################################################

def plot():

    fig1 = plt.figure(starname)

    def key_press_U (event):

        # Global values
        global U
        global C
        global v
        global T

        # Delta values
        delta_U = 0.01
        delta_C =  0.1*10**(-14)
        delta_v_01 = 0.1
        delta_v_1 = 1.0
        delta_v_10 = 10.0
        delta_T = 10

        # For U
        if event.key == 'down':
            U -= delta_U
        if event.key == 'up':
            U += delta_U

        # For C
        if event.key == 't':
            C -= delta_C
        if event.key == 'y':
            C += delta_C

        # For v shift by .1
        if event.key == 'u':
            v -= delta_v_01
        if event.key == 'i':
            v += delta_v_01

        # For v shift by 1
        if event.key == 'v':
            v -= delta_v_1
        if event.key == 'b':
            v += delta_v_1

        # For v shift by 10
        if event.key == 'n':
            v -= delta_v_10
        if event.key == 'm':
            v += delta_v_10

        # For T
        if event.key == 'q':
            T -= delta_T
        if event.key == 'w':
            T += delta_T

        # Updates each needed value
        y = planck_func (x,C,U,v)        
        data_norm = data_flux_norm(stardata_file, y)
        shift = velocity_shift(v)

        # Updates the Plots
        ctm_fit.set_ydata(y)        
        ctm_norm.set_ydata(data_norm[0, 0])
        ctm_data.set_xdata(shift)
        ctm_norm.set_xdata(shift)        

        fig1.canvas.draw()

        # Save info out
        save_const_out(starname)
        os.system('clear')
        print('Ti = %s   vi = %s    Ci = %s    Ui = %s' % (Ti, vi, Ci, Ui) + '\n')
        print('T = %s   v = %s    C = %s    U = %s' % (T, v, C, U) + '\n')

    # Setting Up first subplot
    plt.subplot(121)
    plt.xlabel('Wavelength A')
    plt.ylabel('Flux')
    plt.title('Curve Fitting')

    # First plot information
    plt.plot(x, data[0,0], 'r--')
    ctm_data, = plt.plot(shift, data[0,0], color='g')
    ctm_fit, = plt.plot(x, y, color='black')

    # Legend for first plot
    plt.legend(('Non-Shifted Star Flux', 'Shifted Star Flux', 'Curve Fit'), 'upper center', shadow=True, fancybox=True)  

    plt.draw()

    # Setting up second subplot
    plt.subplot(122, sharex=plt.subplot(121), sharey=plt.subplot(121))
    plt.title('Normalized Star and Synth Flux')

    # Second plot information
    ctm_norm, = plt.plot(shift, data_norm[0,0], 'ro-', color='g')
    plt.plot(add_xy[0][0], add_xy[0][1], color='r')
    plt.plot(add_xy[0][0], add_xy[1][1])
    plt.plot(add_xy[0][0], add_xy[2][1],color='m')
    plt.plot(add_xy[0][0], (add_xy[0][1]-add_xy[1][1])*10+1,color='black')

    # Plot legend for second plot
    plt.legend(('Shifted Star Flux', 'Synth +1Dex', 'Synth 0Dex', 'Synth -1Dex'), 'upper center', shadow=True, fancybox=True)

    plt.draw()

    # Show plot
    fig1.canvas.mpl_connect("key_press_event",key_press_U)
    plt.show()

plot()
