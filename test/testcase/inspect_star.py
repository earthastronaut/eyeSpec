# for inspecting HST data 
import os
import eyeSpec as es
import pdb

def inspec_star (starname, stardata_file, model_dir, linelist_file, rvstat_file):
    # get a list of all the syntheses summary files
    model_dir = os.path.abspath(model_dir)+"/"    
    range_dirs = os.listdir(model_dir)
    sum_files = []
    for directory in range_dirs:
        datapath = model_dir+directory+"/"
        fnames = os.listdir(datapath)
        for fname in fnames:
            sum_files.append(datapath+fname)


    # to get the overplot synthesis data
    add_xy = []
    for fname in sum_files:
        xy,header = es.parse_synth_summary_out(fname)
        add_xy.append(xy)

    # read in the data
    spec = es.readin_hst(stardata_file)
    spec.edit.quick_continuum_normalize()
    
    # radial velocity shift
    spec_shift,rv,rvstat = es.edit_rv(spec)
    
    # record radial velocity
    f = open(rvstat_file,'w')
    f.write("# Radial Velocity stats for star:"+stardata_file+"\n")
    for i in xrange(len(rv)):
        wl = rv[i][0]
        rvel = rv[i][1]
        info = rv[i][2]
        f.write("  ".join((format(wl,">10.3f"),format(rvel,'>10.3f'),info))+"\n")
    f.write("#"+"-"*40+"\n")
    f.write("#"+" STATS:  "+"{0:>10.3f} +- {1:<10.3f}   STD={2:<10.3f}   N={3:<10}".format(*rvstat)+"\n")
    f.close()

    # inspect the data
    es.iplotspec(spec_shift,linelist_file,add_xy)
    


# this runs if you use execfile("inspect_star.py")
if __name__ == "__main__":
    starname = 'HD200790'
    stardata_file = '/uufs/astro.utah.edu/common/astro_data/astrougrad/u0805376/Projects/06082013/o6e615010_x1d.fits'
    model_dir = '/uufs/astro.utah.edu/common/astro_data/astrougrad/u0724846/test/hstoutphos/6115p398p002p100/'
    linelist_file = '/uufs/astro.utah.edu/common/home/u0565532/projects/hst_archive/linelists/hst_inspection_linelist.ln'

    rvstat_file = 'RVSTAT_'+starname+".txt"
    inspec_star(starname, stardata_file, model_dir, linelist_file, rvstat_file)

