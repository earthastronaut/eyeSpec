execfile("main.py")
execfile("base_IO_v2.py")

def get_hdr (filename):
    hdus = pyfits.open(filename)
    return hdus[0].header

target = "../testfiles/"
filelist = ['HI.20110212.58049_1_Flux.fits',
            'apStar-2M18561431+4431052.fits',
            'mult_band_spectest.fits',
            'longspec.fits',
            'multispectest.fits',
            'many_1d/4430-4500.svn']


# for file in filelist:
#     file = target+file
#     cc = disp_coeff_from_header(get_hdr(file))
#     pdb.set_trace()
