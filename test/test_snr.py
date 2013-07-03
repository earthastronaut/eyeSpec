execfile('main.py')

execfile("app_edit_data_v29.py")

import time
print "==== data ===="
#x = readin("../testfiles/multispectest.fits")
x = readin_spectre_files("../testfiles/many_1d/filelist")
#x = readin("../testfiles/many_1d/4430-4500.svn")
#x = readin("../testfiles/SUBARU_AOKI_2183-175n.fits")

#ll = np.loadtxt("/Users/dylangregersen/Desktop/Astrophysics/data/allen/LINELISTS/lin.long.101812",unpack=True,usecols=[0,1])

ll = np.loadtxt("../testfiles/SEGUE.ln2.dat",unpack=True,usecols=[0,1])





#orig = readin_makee("/Users/dylangregersen/Desktop/Astrophysics/data/allen/KOA/data_download_all/G167-50/HIRES_sci_32674_1")

#orig = readin_apogee("../testfiles/apStar-2M18561431+4431052.fits")

#x = orig.copy()

# linelist = []
# ll = []
# f = open('/Users/dylangregersen/Desktop/Astrophysics/documents/eyeSpec/testfiles/cdat_MB.lines')
# for line in f:
#     line=line.rstrip()
#     if line[:1] == '#': continue
#     wl = line.split()[0]
#     info = "  ".join(line.split()[1:])
#     linelist.append([wl,info])
#     ll.append(wl)
# f.close()

execfile("app_snr_editor_v11.py")


# y = load_spec("TMP_OBJ_SAVE_EDIT.pkl")

# ax = plt.figure(10).add_subplot(111)
# ax1 = plt.figure(11).add_subplot(111)
# for i in range(4):
#     ax.plot(y.get_wl(i),y.get_inv_var(i),color='r')
#     ax.plot(orig.get_wl(i),100*orig.get_inv_var(i),color='b',alpha=.7,lw=2)
#     ax.plot(orig.get_wl(i),orig.get_inv_var(i),color='c',alpha=.7,lw=2)
#     ax1.plot(y.get_wl(i),y.get_inv_var(i)/orig.get_inv_var(i))    








# class Base :
#     def __init__ (self) :
#         self.A = 35

#     def go (self) :
#         print self.A


# class Child (object,Base) :
#     def __init__ (self) :
#         super(Child, self).__init__()

#     def go (self) :
#         self.A = self.A * 2
#         super(Child, self).go()
