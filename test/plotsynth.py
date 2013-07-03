execfile("test_rv.py")

wlrange = (3100,3150)
# 1850-1900
# 2080-3000
model = [5000,4.00,-2.0,1.0]

create_kurucz_atmo_model(*model)

linelistfile = "/Users/dylangregersen/Desktop/Astrophysics/applications/eyeSpec/for_testing/hst_linelist_1850-3200.ln"


results = moog_synth(linelistfile,'FINALMODEL',wlrange=wlrange)


linelist,synth,tol=process_moog_synth_output(results['summary_out'],results['standard_out'])



def whereis (linelist,tol,wl):
    for i in xrange(len(linelist)):
        val = linelist[i]
        wlll = val[0]
        if abs(wl-wlll) < tol:
            print i, linelist[i]
    
