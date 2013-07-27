
from ..dependencies import np, np_recfunc
from ..core import int_to_roman

import pdb #@UnusedImport

class AbundanceNormalization (object):
    
    # TODO: impliment this chaning of normalization  
    def __init__ (self):
        self._scale = "H"
        self._scale_opts = ('H','Si')

    def convert_to (self,instance,scale):
        # convert the input to that scale
        pass
    
    def __get__ (self,_instance,_owner):
        return self._scale
    
    def __set__ (self,_instance,value):
        
        if not isinstance(value,str):
            raise TypeError("given scale option must be string")
        
        if value not in self._scale_opts:
            print("no scale option "+str(value))
        # set instance._current = [abundances]
        
        self._scale = value
        pass

class PeriodicTable (object):
    # this holds the Z,el,element_name
    # converters to get elements, take Z=>el and el=>Z all with arrays
    # TODO : make this a data base to load in
    _table_data = np.array([(1,"H","Hydrogen"),
                     (2,"He","Helium"),
                     (3,"Li","Lithium"),
                     (4,"Be","Beryllium"),
                     (5,"B","Boron"),
                     (6,"C","Carbon"),
                     (7,"N","Nitorgen"),
                     (8,"O","Oxygen"),
                     (9,"F","Florine"),
                     (10,"Ne","Neon"),
                     (11,"Na","Sodium"),
                     (12,"Mg","Magnesium"),
                     (13,"Al","Aluminium"),
                     (14,"Si","Silicon"),
                     (15,"P","Phosphorus"),
                     (16,"S","Slufer"),
                     (17,"Cl","Chlorine"),
                     (18,"Ar","Argon"),
                     (19,"K","Potassium"),
                     (20,"Ca","Calcium"),
                     (21,"Sc","Scandium"),
                     (22,"Ti","Titanium"),
                     (23,"V","Vanadium"),
                     (24,"Cr","Chromium"),
                     (25,"Mn","Maganese"),
                     (26,"Fe","Iron"),
                     (27,"Co","Cobalt"),
                     (28,"Ni","Nickel"),
                     (29,"Cu","copper"),
                     (30,"Zn","Zinc"),
                     (31,"Ga","Gallium"),
                     (32,"Ge","Germanium"),
                     (33,"As","Arsnic"),
                     (34,"Se","Selenium"),
                     (35,"Br","Bormine"),
                     (36,"Kr","Krypton"),
                     (37,"Rb","Rubidium"),
                     (38,"Sr","Strontium"),
                     (39,"Y","Yttrium"),
                     (40,"Zr","Zirconium"),
                     (41,"Nb","Niobium"),
                     (42,"Mo","Molybdenum"),
                     (43,"Tc","Technetium"),
                     (44,"Ru","Ruthenium"),
                     (45,"Rh","Rhodium"),
                     (46,"Pd","Palladium"),
                     (47,"Ag","Silver"),
                     (48,"Cd","Cadmium"),
                     (49,"In","Indium"),
                     (50,"Sn","Tin"),
                     (51,"Sb","Antimony"),
                     (52,"Te","Tellurium"),
                     (53,"I","Iodine"),
                     (54,"Xe","Xenon"),
                     (55,"Cs","Caesium"),
                     (56,"Ba","Barium"),
                     (57,"La","Lanthanum"),
                     (58,"Ce","Cerium"),
                     (59,"Pr","Praseodymium"),
                     (60,"Nd","Neodymium"),
                     (61,"Pm","Promethium"),
                     (62,"Sm","Samarium"),
                     (63,"Eu","Europium"),
                     (64,"Gd","Gadolinium"),
                     (65,"Tb","Terbium"),
                     (66,"Dy","Dyspeosium"),
                     (67,"Ho","Holmium"),
                     (68,"Er","Erbium"),
                     (69,"Tm","Thulium"),
                     (70,"Yb","Ytterbium"),
                     (71,"Lu","Luletium"),
                     (72,"Hf","Hafnium"),
                     (73,"Ta","Tantalum"),
                     (74,"W","Tungsten"),
                     (75,"Re","Rhenium"),
                     (76,"Os","Osmium"),
                     (77,"Ir","iridium"),
                     (78,"Pt","Plantinum"),
                     (79,"Au","Gold"),
                     (80,"Hg","Mercury"),
                     (81,"Tl","Thallium"),
                     (82,"Pb","Lead"),
                     (83,"Bi","Bismuth"),
                     (84,"Po","Polomium"),
                     (85,"At","Astatine"),
                     (86,"Rn","Radon"),
                     (87,"Fr","Francium"),
                     (88,"Ra","Radium"),
                     (89,"Ac","Actinium"),
                     (90,"Th","Thorium"),
                     (91,"Pa","Protactinium"),
                     (92,"U","Uranium"),
                     (93,"Np","Neptunium"),
                     (94,"Pu","Plutonium"),
                     (95,"Am","Americium")],
                    dtype=[('z',float),('element','a2'),('element_long','a15')])
    
    def __init__ (self):
        # TODO: write doc string
        self._by_z = {}
        self._by_el = {}
        
        for i,row in enumerate(self._table_data):
            self._by_z[row[0]] = i
            self._by_el[row[1].lower()] = i

    @property
    def table_data (self):
        return self._table_data
    
    def _single_species_name (self,spe):
        
        if isinstance(spe,int):
            return self[spe][1]
        
        elif isinstance(spe,float):
            el = self[spe][1]
            ionz = int(round((spe - round(spe-0.5))*10))+1
            return el+" "+int_to_roman(ionz)
        
        else:
            return 'unknown'

    def species_name (self,spe):
        """
        Takes a species id and converts to a name representation
        
        species_id := proton_number + 0.1*(ionization_state)
        where the ionization_state is 0 for ground, 1 for singally, etc
        
        """
        if isinstance(spe,(float,int)):
            return self._single_species_name(spe)
        
        elif isinstance(spe,(list,tuple)):    
            return [self._single_species_name(s) for s in spe]
        elif isinstance(spe,np.ndarray):
            data = self.table_data['element']
            if isinstance(spe.dtype,int):
                return data[self.table_dat['']]
        else:
            raise TypeError("must receive either a floating value of the species or an array")

    def __iter__ (self):
        return iter(self.table_data)

    def __contains__ (self,spe):
        if isinstance(spe,str):
            return spe.lower().strip().split()[0] in np.core.defchararray.lower(self.table_data['element'])
        
        if isinstance(spe,(int,float)):
            return int(spe) in self.table_data['z']
    
    def __reversed__ (self):
        return reversed(self.table_data)
    
    def __len__ (self):
        return self.table_data.shape[0]

    def _getitem_single (self,spe):
        if isinstance(spe,str):
            return self.table_data[self._by_el[spe.lower().strip().split()[0]]]
        
        if isinstance(spe,(int,float)):
            return self.table_data[self._by_z[int(spe)]]

        return None

    def _getitem (self,spe):
        """ check input and return appropriate output """
        
        value = self._getitem_single(spe)
        if value is not None:
            return value

        if isinstance(spe,(tuple,list)):
            return np.array([self._getitem_single(s) for s in spe],dtype=self.table_data.dtype)
            
        if isinstance(spe,np.ndarray):
            if spe.shape[0] == 0: 
                return np.array([],dtype=self.table_data.dtype)
            idx = np.round(spe)-1
            return self.table_data[idx.astype(int)]
      
    def __getitem__ (self,spe):
        try: 
            return self._getitem(spe)
        except KeyError:
            raise KeyError("Please enter valid element name or Z, not "+str(spe))

periodic_table = PeriodicTable()

class AbundanceSystem (PeriodicTable):
    # For every value on the periodic_table it has a corresponding (abund,error)
    # have ways to input Z,el as value or array and get back abund,error
    
    # methods for [x/fe],[x/y] and logeps conversions
        
    # have a atomic_normalization that it has.
    abundance_norm = AbundanceNormalization()   
    def __init__ (self,citation,data,normalization='H'):
        """ 
        Stores the data for a particular system of abundances to normalize to
        
        also see help(PeriodicTable)
        
        Parameters
        ----------
        citation : string
            The citation information for where the abundance standard was taken
        data : array like
            Has data [[z,abund,stdev],[...],...] you can optionally include stdev
        normalization : 'H','Si'
            This is the scale the data is on, 'H' ==> H=10**12, 'Si'==> Si=10**6
        
        Attributes
        ----------
        citation : string
            Returns a string which has the input citation information
        array : np.ndarray
            Returns an floatting point array of [[z,abund,sigma],[...],...]
        convert_to_logeps : function
            Takes values from bracket notation and converts to logeps
        convert_to_xfe : function
            Takes values from logeps and converts to bracket notation
        solar_abundance : function
            Takes a value or array of element ids (Z,element_name) and returns
            the solar abundance for those values
            
        Raises
        ------
        KeyError : if given element identication is unknown 
    
        Notes
        -----
        __1)__ Element identification can be done by proton number or element name (e.g 26 or 'Fe')
            The id is not case sensitive (e.g. 'Fe' == 'fe') and you can but in species
            identification (e.g. ionized iron is 26.1 or 'Fe I' ==> return Fe values)

    # TODO: finish doc string                 
    Examples
    --------
    >>> batom = Batom()
    
    >>> print(batom[1])
    [1,"H","Hydrogen",12.00]
    
    >>> print batom[26.1]
    [26, 'Fe', 'Iron', 7.52]
    
    >>> print(batom['Hg'])
    [80, 'Hg', 'Mercury', 1.09]
    
    >>> print(batom['au']
    [79, 'Au', 'Gold', 0.83]
    
    >>> print(batom['ba II'])
    [56, 'Ba', 'Barium', 2.13]
    
    >>> # you can also give a list/array of Z
    >>> batom[[22.1,22,20.0,25.2]]
    array([ 4.99,  4.99,  6.36,  5.39])
    
    >>> # can use to convert
    >>> batom.convert_to_xfe(22, 3.1, -2.3)
    0.4099999999999997
    
    >>> batom.convert_to_xfe([22,22,20,25],[3.1,3.2,4.4,3.1],-2.3)
    array([ 0.41,  0.51,  0.34,  0.01])
    

            
        """
        if not isinstance(citation,str):
            raise TypeError("Citation must a string for the abundance system standard")
        
        super(AbundanceSystem,self).__init__()
        if not isinstance(data,(tuple,list,np.ndarray)):
            raise TypeError("data must be iterable object with [[z,abund,uncertainty],[...],...]")

        self._citation = citation        

        # get all the data and map it to the values in the periodic table
        by_z = {}
        for row in data:
            try: spe = self[row[0]]
            except KeyError:
                raise ValueError("Reveived unknown element name or Z : "+str(row[0]))
            
            abund = float(row[1])
            if len(row)>2:
                error = float(row[2] or np.nan)**2
            else:
                error = np.nan
            
            by_z[int(spe[0]-1)] = (abund,error)
        
        # create a numpy recarray of the values, filling in missing values  
        abunds = []
        sigma = []
        for i in xrange(len(self)):
            a,v = by_z.get(i,(np.nan,np.nan))
            abunds.append(a)
            sigma.append(v)
                
        # combine the tables
        names = ('abundance','sigma')
        dtypes = (float,float)
        self._table_data = np_recfunc.append_fields(self.table_data,names,(abunds,sigma),dtypes=dtypes,usemask=False,asrecarray=True)
        
        self.abundance_norm = normalization
        
        self._init_str_representation()

    def _init_str_representation (self):
        title_str = "  {0:>4} {1:>9} {2:>16} {3:>10} {4:>10}"
        data_str =  " [{0:>4},{1:>9},{2:>16},{3:>9.3f},{4:>9.3f}]"
        
        rep = title_str.format(*tuple(self.table_data.dtype.names))+"\n"
        rowlength = len(rep)
        rep += "-"*rowlength+"\n"
        rep += "["
        
        for i in xrange(0,len(self)-1): 
            row = self.table_data[i]
            rep += data_str.format(*row)+",\n"
        
        rep += data_str.format(*self.table_data[-1])+"]\n"
        rep += "-"*rowlength+"\n"+title_str.format(*tuple(self.table_data.dtype.names))+"\n"
        
        self._str_rep = rep
        
    def __repr__ (self):
        return "Abundance System : "+str(self.citation)
        
    def __str__ (self):
        return self._str_rep
              
    @property
    def citation (self):
        return self._citation
    
    @property
    def array (self):
        return np.dstack((self.table_data['z'].astype(float),self.table_data['abundance'],self.table_data['sigma']))[0]
    
    pass
#     def convert_to_xy (self,speX,speY,logepsX,logepsY,varX=None,varY=None):
#         """
#         Convert value from logepsX and logepsY to [X/Y] 
#         
#         [X/Y] = (logeps(X)/logeps(Y)) - (logeps(X)_sun/logeps(Y)_sun)
# 
#         ----------
#         speX,speY : float or array of floats
#             Gives the specie(s) which correspond to the logepsX, logepsY values
# 
#         logepsX,logepsY : float or array of floats
#             Gives the logeps notation of the abundance for the 
#             corresponding species
#             
#         Returns
#         -------
#         if spe is a float:
#             xy : float
#                 The abundance of the species with abundance logeps
#             var : float
#                 The sigma on that abundance
#         
#         if spe is an array of floats:
#             xy : array
#                 The abudnaces for the species with abundances logeps
#         
#             var : array
#                 The sigma on that abundance
#         """
#         logepsX_sun,_varX_sun = self.solar_abundance(speX).T
#         logepsY_sun,_varY_sun = self.solar_abundance(speY).T
#         return (logepsX/logepsY)-(logepsX_sun/logepsY_sun), np.nan
    
    def convert_to_logeps (self, spe, xfe, feh):
        """
        
        Convert value in [X/Fe] to logeps solar hydrogen scale
        
        [Fe/H] = logeps(Fe) - logeps(Fe)_sun
        logeps(X) = [X/Fe] + ([Fe/H] + logeps(X)_sun)
    
        Parameters
        ----------
        spe : float or array of floats
            Gives the specie(s) which correspond to the [X/Fe] values
        xfe : float or array of floats
            Gives the [X/Fe] bracket notation of the abundance for the 
            corresponding species
        feh : float
            the [Fe/H] value to be used for the converstion
            
        Returns
        -------
        if spe is a float:
            logeps : float
                The abundance of the species with abundance [X/Fe]
        if spe is an array of floats:
            logeps : array
                The abudnaces for the species with abundances [X/Fe]

        """
        return xfe + (feh + self.solar_abundance(spe)[:,0])
        
    def convert_to_xfe (self, spe, logeps, feh):
        """
        
        Convert value in logeps solar hydrogen scale to [X/Fe]
        
        [Fe/H] = logeps(Fe) - logeps(Fe)_sun
        [X/Fe] = logeps(X) - ([Fe/H] + logeps(X)_sun)

        Parameters
        ----------
        spe : float or array of floats
            Gives the specie(s) which correspond to the logeps values
        xfe : float or array of floats
            Gives the logeps notation of the abundance for the 
            corresponding species
        feh : float
            the [Fe/H] value to be used for the converstion
            
        Returns
        -------
        if spe is a float:
            xfe : float
                The abundance of the species with abundance logeps
        if spe is an array of floats:
            xfe : array
                The abudnaces for the species with abundances logeps
                
        """
        return logeps - (feh + self.solar_abundance(spe)[:,0])

    def solar_abundance (self,spe):
        """
        Get the solar logeps abundance for a species (or list of species)
        
        Parameters
        ----------
        spe : float or array of floats
            Give the species identification for which to return the solar abundance
                    
        Returns
        -------
        if spe is float:
            logeps : float
        if spe is array of floats:
            logeps : array
        """
        
        getitem = self.__getitem__(spe)
        
        if getitem is None:
            return None
        
        if isinstance(getitem,np.void):
            return np.array([tuple(getitem)[3:5]])
        else:
            return np.dstack((getitem['abundance'],getitem['sigma']))[0] 
                
class AbundanceSystems (dict):

    def __doc__ (self):
        # TODO: write doc string
        pass
    
    def _check_abundsys_in (self,abundance_system):
        if not isinstance(abundance_system,AbundanceSystem):
            raise TypeError("received wrong type of abundance_system")

    def __setitem__ (self,name,abundance_system):
        self._check_abundsys_in(abundance_system)
        if name in self:
            return
        super(AbundanceSystems,self).__setitem__(name,abundance_system)
          
class _CurrentSystem (object):
    # TODO : write doc string
    
    def __init__ (self):
        
        """
        
        """
        self._current_system = "None"
    
    def __get__ (self,_instance,_owner):
        return self._current_system
    
    def __set__ (self,instance,name):
        if not isinstance(instance,Abundance):
            raise TypeError("This class is meant to be used by AbundanceStandards")
        
        if name not in instance._abundance_systems:
            raise KeyError("Unknown abundance system : "+str(name))
        
        self._current_system = name
        
        # check value is appropriate
        abundance_system = instance._abundance_systems[name]
        
        instance._table_data = abundance_system._table_data
        instance._citation = abundance_system._citation
        instance.abundance_norm = abundance_system.abundance_norm
        instance._str_rep = abundance_system._str_rep
        
class Abundance (AbundanceSystem):
     
    # TODO: write doc string
        
    system = _CurrentSystem() 
    def __init__ (self,abundance_systems,system=None):
        super(Abundance,self).__init__('None',[])
        
        if not isinstance(abundance_systems,AbundanceSystems):
            raise TypeError("Abundance system must be of class AbundanceSystems")
        
        if len(abundance_systems) == 0:
            raise ValueError("must initiate with none empty abundance_systems")
        
        # internal dictionary of all abundance systems
        self._abundance_systems = abundance_systems
        
        # set which system to use
        self.system = abundance_systems.keys()[-1]
        if system is not None:
            self.system = system
     
    def __doc__ (self):
        super(Abundance,self).__doc__()
     
    @property
    def systems (self):
        return self._abundance_systems.keys()

# TODO: Changed the data to be stored in the es-data as a database
def _data_sets ():
    # TODO: Write doc string
    _batom = [(1, 12.0, np.nan),
             (2, 10.99, np.nan),
             (3, 3.31, np.nan),
             (4, 1.42, np.nan),
             (5, 2.88, np.nan),
             (6, 8.56, np.nan),
             (7, 8.05, np.nan),
             (8, 8.93, np.nan),
             (9, 4.56, np.nan),
             (10, 8.09, np.nan),
             (11, 6.33, np.nan),
             (12, 7.58, np.nan),
             (13, 6.47, np.nan),
             (14, 7.55, np.nan),
             (15, 5.45, np.nan),
             (16, 7.21, np.nan),
             (17, 5.5, np.nan),
             (18, 6.56, np.nan),
             (19, 5.12, np.nan),
             (20, 6.36, np.nan),
             (21, 3.1, np.nan),
             (22, 4.99, np.nan),
             (23, 4.0, np.nan),
             (24, 5.67, np.nan),
             (25, 5.39, np.nan),
             (26, 7.52, np.nan),
             (27, 4.92, np.nan),
             (28, 6.25, np.nan),
             (29, 4.21, np.nan),
             (30, 4.6, np.nan),
             (31, 2.88, np.nan),
             (32, 3.41, np.nan),
             (33, 2.37, np.nan),
             (34, 3.35, np.nan),
             (35, 2.63, np.nan),
             (36, 2.23, np.nan),
             (37, 2.6, np.nan),
             (38, 2.9, np.nan),
             (39, 2.24, np.nan),
             (40, 2.6, np.nan),
             (41, 1.42, np.nan),
             (42, 1.92, np.nan),
             (43, 0.0, np.nan),
             (44, 1.84, np.nan),
             (45, 1.12, np.nan),
             (46, 1.69, np.nan),
             (47, 1.24, np.nan),
             (48, 1.86, np.nan),
             (49, 0.82, np.nan),
             (50, 2.0, np.nan),
             (51, 1.04, np.nan),
             (52, 2.24, np.nan),
             (53, 1.51, np.nan),
             (54, 2.23, np.nan),
             (55, 1.12, np.nan),
             (56, 2.13, np.nan),
             (57, 1.22, np.nan),
             (58, 1.55, np.nan),
             (59, 0.71, np.nan),
             (60, 1.5, np.nan),
             (61, 0.0, np.nan),
             (62, 1.0, np.nan),
             (63, 0.51, np.nan),
             (64, 1.12, np.nan),
             (65, 0.33, np.nan),
             (66, 1.1, np.nan),
             (67, 0.5, np.nan),
             (68, 0.93, np.nan),
             (69, 0.13, np.nan),
             (70, 1.08, np.nan),
             (71, 0.12, np.nan),
             (72, 0.88, np.nan),
             (73, 0.13, np.nan),
             (74, 0.68, np.nan),
             (75, 0.27, np.nan),
             (76, 1.45, np.nan),
             (77, 1.35, np.nan),
             (78, 1.8, np.nan),
             (79, 0.83, np.nan),
             (80, 1.09, np.nan),
             (81, 0.82, np.nan),
             (82, 1.85, np.nan),
             (83, 0.71, np.nan),
             (84, 0.0, np.nan),
             (85, 0.0, np.nan),
             (86, 0.0, np.nan),
             (87, 0.0, np.nan),
             (88, 0.0, np.nan),
             (89, 0.0, np.nan),
             (90, 0.12, np.nan),
             (91, 0.0, np.nan),
             (92, 0.0, np.nan),
             (93, 0.0, np.nan),
             (94, 0.0, np.nan),
             (95, 0.0, np.nan)]
    
    
    _lodders03 =[("h", 12.0, 0.0),
                ("he", 10.899, 0.01),
                ("li", 3.28, 0.06),
                ("be", 1.41, 0.08),
                ("b", 2.78, 0.04),
                ("c", 8.39, 0.04),
                ("n", 7.83, 0.11),
                ("o", 8.69, 0.05),
                ("f", 4.46, 0.06),
                ("ne", 7.87, 0.1),
                ("na", 6.3, 0.03),
                ("mg", 7.55, 0.02),
                ("al", 6.46, 0.02),
                ("si", 7.54, 0.02),
                ("p ", 5.46, 0.04),
                ("s ", 7.19, 0.04),
                ("cl", 5.26, 0.06),
                ("ar", 6.55, 0.08),
                ("k ", 5.11, 0.05),
                ("ca", 6.34, 0.03),
                ("sc", 3.07, 0.04),
                ("ti", 4.92, 0.03),
                ("v ", 4.0, 0.03),
                ("cr", 5.65, 0.05),
                ("mn", 5.5, 0.03),
                ("fe", 7.47, 0.03),
                ("co", 4.91, 0.03),
                ("ni", 6.22, 0.03),
                ("cu", 4.26, 0.06),
                ("zn", 4.63, 0.04),
                ("ga", 3.1, 0.06),
                ("ge", 3.62, 0.05),
                ("as", 2.32, 0.05),
                ("se", 3.36, 0.04),
                ("br", 2.59, 0.09),
                ("kr", 3.28, 0.08),
                ("rb", 2.36, 0.06),
                ("sr", 2.91, 0.04),
                ("y ", 2.2, 0.03),
                ("zr", 2.6, 0.03),
                ("nb", 1.42, 0.03),
                ("mo", 1.96, 0.04),
                ("ru", 1.82, 0.08),
                ("rh", 1.11, 0.03),
                ("pb", 1.7, 0.03),
                ("ag", 1.23, 0.06),
                ("cd", 1.74, 0.03),
                ("in", 0.8, 0.03),
                ("sn", 2.11, 0.04),
                ("sb", 1.06, 0.07),
                ("te", 2.22, 0.04),
                ("i", 1.54, 0.12),
                ("xe", 2.27, 0.02),
                ("cs", 1.1, 0.03),
                ("ba", 2.18, 0.03),
                ("la", 1.18, 0.06),
                ("ce", 1.61, 0.02),
                ("pr", 0.78, 0.03),
                ("nd", 1.46, 0.03),
                ("sm", 0.95, 0.04),
                ("eu", 0.52, 0.04),
                ("gd", 1.06, 0.02),
                ("tb", 0.31, 0.03),
                ("dy", 1.13, 0.04),
                ("ho", 0.49, 0.02),
                ("er", 0.95, 0.03),
                ("tm", 0.11, 0.06),
                ("yb", 0.94, 0.03),
                ("lu", 0.09, 0.06),
                ("hf", 0.77, 0.04),
                ("ta", -0.14, 0.03),
                ("w", 0.65, 0.03),
                ("re", 0.26, 0.04),
                ("os", 1.37, 0.03),
                ("ir", 1.35, 0.03),
                ("pt", 1.67, 0.03),
                ("au", 0.83, 0.06),
                ("hg", 1.16, 0.18),
                ("tl", 0.81, 0.04),
                ("pb", 2.05, 0.04),
                ("bi", 0.68, 0.03),
                ("th", 0.09, 0.04),
                ("u", -0.49, 0.04)]
    return _batom,_lodders03

_batom,_lodders03 = _data_sets()
l03 = AbundanceSystem(('Lodders, K. 2003 Solar System Abundances and Condensation ' \
                      'Temperatures of the Element, APJ, 591:1220-1247'),_lodders03)
batom = AbundanceSystem(("the set of current solar (when available) or meteorite "       
                         'abundances, scaled to log(h) = 12.00 .  The data are from Anders '
                         'and Grevesse (1989, Geochim.Cosmichim.Acta, v53, p197) and the solar '
                         "values are adopted except for a) uncertain solar data, or b) Li, Be, "
                         "and B, for which the meteoritic values are adopted. " 
                         "I was told to use the new Fe value of 7.52 as adopted in Sneden "
                         "et al. 1992 AJ 102 2001."),_batom)
abundance_systems = AbundanceSystems({'lodder03':l03,'moog':batom})
abundance = Abundance(abundance_systems,system='moog')

