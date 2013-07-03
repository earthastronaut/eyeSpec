print "--- external modules ---"
execfile("../dependencies.py")
    
print "--- classes and functions ---"
execfile("../base_classes.py")
execfile("../base_functions.py")

execfile("../base_IO.py")
execfile("../extended_IO.py")

execfile("../moog_functions.py")
execfile("../MOOG.py")

    
print "--- 3rd Party ---"
execfile("../resampling.py")
execfile("../piecewise_poly.py")
execfile("../linefinder_dsg.py")

print "--- combining function ---"
execfile("../base_combine.py")

    
print "--- plotting ---"
execfile("../base_plotting.py")
execfile("../interactive_IO.py")
execfile("../interactive_classes.py")
    
Params = {'verbose':True}

path_2_eyeSpec = ["/Users/dylangregersen/Desktop/Astrophysics/applications/eyeSpec/eyeSpec"]
