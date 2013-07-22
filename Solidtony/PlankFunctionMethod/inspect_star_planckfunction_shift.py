# for inspecting HST data 
import os
import eyeSpec as es
import numpy as np
import pdb
import matplotlib.pyplot as plt

#import pdb; pdb.set_trace()

def windowed_smoothing(synth, smoothing):
    # synth assuming it has wavelength = synth[0]  flux = synth[1]
    window=np.hamming(smoothing)
    synth[1] =  np.convolve(synth[1], window, mode = 'same')/np.convolve(np.ones(len(synth[1])), window, mode = 'same')
    return synth

def inspec_star(starname, stardata_file, model_dir, smoothing):
    # /uufs/astro.utah.edu/common/astro_data/astrougrad/u0724846/test/hstoutphos/
    # get a list of all the syntheses summary files
    model_dir = os.path.abspath(model_dir)+"/"    
    range_dirs = os.listdir(model_dir)
    sum_files = []
    for directory in range_dirs:
        datapath = model_dir+directory+"/"
        fnames = os.listdir(datapath)
        for fname in fnames:
            sum_files.append(datapath+fname)

    # Read in the data
    spec = es.readin(stardata_file)
    #spec = es.readin_hst(stardata_file)  # Currently not working for x1d files
    data0 = spec.access_data()
    data0 = data0/np.median(data0)
    spec.set_data(np.array(data0))

##############################################################################
# All values with # after the value are values that very between each star.
# Notes to self, the method of determining velocity is somewhat pathetic. Playing with it and visually fitting is somewhat difficult, so I just use Dylan's velocity shift code and plug in that value. Determining the curve constant C is also pretty difficult.

    # Plank function constants and Temp for normalization
    c = 2.99792458*10**10
    k = 1.38066*10**(-16)
    T = 5100 #
    h = 6.62608*10**(-27)

    # Velocity shifting
    wvl = spec.get_wl()
    v = -97.8 #
    shift =wvl*(v)*100000/c+wvl

    # Normalization
    C = 1/(10**(13)*2.1) #
    U = .37 #
    x = spec.get_wl()
    y = (2*h*(c**2)/((x*1*10**(-8))**5))/(np.exp(h*c/((x*1*10**(-8))*k*T))-1)*C+U
    data = spec.access_data()
    data_norm = data/y
    spec.set_data(np.array(data_norm)) 

    # Ploting

    # Setting up sub plot
    plt.subplot(121)
    plt.xlabel('Wavelength A')
    plt.ylabel('Flux')
    plt.title('Curve Fitting')

    # Plotted information
    plt.plot(wvl, data[0,0], 'r--')
    plt.plot(shift, data[0,0], color='g')
    plt.plot(x, y, color='black')

    plt.legend(('Non-Shifted Star Flux', 'Shifted Star Flux', 'Curve Fit'), 'upper center', shadow=True, fancybox=True)    
    plt.draw()

##############################################################################

    # To get the overplot synthesis data
    add_xy = []
    for fname in sum_files:
        xy, header = es.parse_synth_summary_out(fname)
        # smooth the synthesis
        xy = windowed_smoothing(xy,smoothing)
        # add to the list of data we're going to overplot
        add_xy.append(xy)

    # Plotting

    # Setting up sub plot
    plt.subplot(122, sharex=plt.subplot(121), sharey=plt.subplot(121))
    plt.title('Normalized Star and Synth Flux')

    # Plotted information
    plt.plot(shift, data_norm[0,0], 'ro-', color='g')
    plt.plot(add_xy[0][0], add_xy[0][1], color='r')
    plt.plot(add_xy[0][0], add_xy[1][1])
    plt.plot(add_xy[0][0], add_xy[2][1],color='m')
    plt.plot(add_xy[0][0], -add_xy[0][1]+add_xy[1][1]+1.09,color='black')

    plt.legend(('Shifted Star Flux', 'Synth +1Dex', 'Synth 0Dex', 'Synth -1Dex'), 'upper center', shadow=True, fancybox=True)
    plt.draw()

    plt.show()
##################################################################################################
    # Saving values for constants
    f = open('radial_velocities/'+starname+'.txt','w')
    f.write("# Radial Velocity and Constants for star:"+starname+'\n'+'\n')
    f.write('T='+'%s' % T+'   '+'v='+'%s' % v+'   '+'C='+'%s' % C+'   '+'U='+'%s' % U+'\n'+'\n')
    f.write("##################################################################"+"\n"+'Terms'+'\n')
    f.write("# T is the temperature taken from the xls file for the stars""\n")
    f.write("# v is the velocity shift that best fits the synth"+"\n")
    f.write("#C is the constant for amount of flux, determined by best observed fit"+"\n")
    f.write("#U is used to shift up or down the Curve fit"+'\n')
    f.close()
##################################################################################################
# this runs if you use execfile("inspect_star.py")
if __name__ == "__main__":
##################################################################################################
# for each observation:

    stardata_file = '/media/solidtonyubuntu/PROJECT_HST/StarData/HD108317/Set1/st.avg.fits'

    model_dir = '/media/solidtonyubuntu/PROJECT_HST/hst_atmo_models/hstoutphos/5100p250m235p200/'

##################################################################################################
    starname =stardata_file.replace('/media/solidtonyubuntu/PROJECT_HST/StarData/','')
    starname = starname.replace('/st.avg.fits','')
    starname = starname.replace('_x1d.fits','')
    starname = starname.replace('/','')

    linelist_file = '/media/solidtonyubuntu/PROJECT_HST/StarInspection/hst_inspection_linelistb.ln'

    smoothing = 13.0

    inspec_star(starname, stardata_file, model_dir, smoothing)

