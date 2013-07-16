import numpy as np
#import matplotlib.pyplot as plt
import re
import os

# Written by Brian Kimmig 09/2012
#-----------------------------------------------------------------------
#load in the data
#file = raw_input("\n Please enter the name or the path to your file and press (text only) 'Enter': \n")

#wl,flux = np.loadtxt(file, unpack = 'True', usecols = [0,1])
#-----------------------------------------------------------------------

# changed so that it finds continuum not lines
def continuum_finder (obj_xy,line_depth,num_sigma=4):

	wl = obj_xy[0]
	flux = obj_xy[1]


	# !! through in check for normalized

	#make the function
	# !! line_depth is threshold, must be normalized,
	def find_minima (in_vector, line_depth = line_depth):
	    isminima = np.zeros(len(in_vector))
	    for i in range(2, len(in_vector)-2):
	        if in_vector[i] < line_depth:
        	    if in_vector[i-1] > in_vector[i] and in_vector[i+1] > in_vector[i] and in_vector[i-2] > in_vector[i] and in_vector[i+2] > in_vector[i] :
			    isminima[i] = 1
	    return isminima

		
	
	#-----------------------------------------------------------------------
	#use the function
	mins = find_minima(flux)
	places = np.where(mins)

	waves = wl[places]
	fluxes = flux[places]

	#define the flux and half width regions (Af flux, Aw width)
	# if you're Af depth then you have Aw width
	# !! could add way to make it based on characteristic width FWHM
	# !! guess FWHM based on half the distance to the depth
	l_fwhm = []
	l_wl = []
	r_fwhm = []
	r_wl = []
	for i in range(len(places[0][1:-3])):
    #array = places[0][1:-3]
    #here = array[i]
            here = places[0][i]
	    min_fl = flux[here]
	    min_wv = wl[here]
    #plt.scatter(min_wv,min_fl,c= 'm',s = 60)
            fwhm = (1-min_fl)*0.5 + min_fl
  
	    j = 0
	    while True:
		    if (here-j) < 0: break
		    if  flux[here-j] >= fwhm: break
		    j = j + 1
	    l_fwhm.append(flux[here - j +1])
	    l_wl.append(wl[here-j +1])
        
    
	    k = 0
	    while True:
		    if (here+k) >= len(flux): break
		    if flux[here + k] >= fwhm: break
		    k = k + 1
	    r_fwhm.append(flux[here + k -1])
	    r_wl.append(wl[here+k-1])
        

	l_fwhm = np.array(l_fwhm[2:])
	l_wl = np.array(l_wl[2:])
	r_fwhm = np.array(r_fwhm[2:])
	r_wl = np.array(r_wl[2:])

	#plot the fwhm
	#plt.plot(wl,flux,'b')
	#plt.plot(l_wl,l_fwhm,'go')
	#plt.plot(r_wl,r_fwhm,'ro')
	f = []
	for i in range(len(l_wl)):
		fwhm_f = (l_fwhm[i] + r_fwhm[i])/2
		fwhm_l = r_wl[i] - l_wl[i]
		#    plt.hlines(fwhm_f,l_wl[i],r_wl[i],'m')
		f.append(fwhm_l)
	
	# 	plt.show()

	fw = []
	for i in range(len(f)):
		full_width = num_sigma*(f[i]/(2.3548))
		fw.append(full_width)



	fluxes = flux[places[0][1:-3]]
	waves = wl[places[0][1:-3]]
	regions = []
	for i in range(len(fw)):
	    region_start = waves[i]
	    region_end = waves[i+1]
	    flux_start = fluxes[i]
	    flux_end = fluxes[i+1]
	    wvb = region_start + fw[i]
	    wve = region_end - fw[i]
	    regions.append([wvb,wve])
 
    
	region = np.array(regions)
	pairs = []
	for i in range(len(region)):
		if region[i][0] < region[i][1]:
			pairs.append([region[i][0],region[i][1]])


	return pairs
# 	wvs = []
# 	fls =[]
# 	for i in range(len(pairs)):
# 		for j in range(len(wl)):
# 			if wl[j] >= pairs[i][0] and wl[j] <= pairs[i][1]:
# 				wvs.append(wl[j])
# 				fls.append(flux[j])

# 	wvs = np.array(wvs)
# 	fls = np.array(fls)


	#plot the continuum regions
	#plt.plot(wl,flux,'b')
	#plt.plot(wvs,fls,'k')
	#plt.show()

