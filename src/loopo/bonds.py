import numpy as np

class BOND :
    def __init__( self ) :
        self.vmin = -1
        self.vmax = -1
        self.type = -1
        self.length = -1.0

class BONDS_CONFIGURATION :
    def NEW__init__( self ) :
        self.time = -1.0
        self.vmin = np.array([])
        self.vmax = np.array([])
        self.type = np.array([])
        self.length = np.array([])
        self.neigh_bonds = []
        self.N = 0

    def __init__( self ) :
        self.time = -1.0
        self.neigh_bonds = {}
        self.N = 0
        # OLD: to cancel in future
        self.bonds = []
        ## NEW
        self.vmin = np.array([]).astype(int)
        self.vmax = np.array([]).astype(int)
        self.type = np.array([]).astype(int)
        self.length = np.array([]).astype(float)

    def read_array( self , infile , FORMAT ) :
        KEEP_GOING = 1
        line = "initial_line"
        while KEEP_GOING and line != '' :
            line = infile.readline()
            words = line.strip().split()
            if len(words) > 0 :
                if "#" in words[0] :
                    if words[1] in [ "timestep" , "tstep" , "step" ] :
                        self.time = float(words[2])
                    KEEP_GOING = 0
        words = np.loadtxt( infile , dtype = 'str' )
        self.N = len(words)
        if self.N > 0 :
            if "#" in words[0,0] :
                self.N = 0
            else :
                self.type = words[:,3].astype(int)
                self.length = words[:,2].astype(float)
                self.vmin = np.min( words[:,0:2].astype(int) , axis=1 )
                self.vmax = np.max( words[:,0:2].astype(int) , axis=1 )
                self.N = len( self.vmin )

    def read( self , infile , FORMAT ) :
        lines = infile.readlines()
        i0 = -1
        KEEP_GOING = 1
        while KEEP_GOING and i0+1 < len(lines) :
            i0 += 1
            words = lines[i0].strip().split()
            if len(words) > 0 :
                if "#" in words[0] :
                    if words[1] in [ "timestep" , "tstep" , "step" ] :
                        self.time = float(words[2])
                else :
                    KEEP_GOING = 0
        if i0 == len(lines)-1 and "#" in words[0] :
            self.N = 0
        else :
            for i in range( i0 , len(lines) ) :
                words = lines[i].strip().split()
                bond = BOND()
                bond.type = int(words[3])
                bond.length = float(words[2])
                bond.vmin = min( int(words[0]) , int(words[1]) )
                bond.vmax = max( int(words[0]) , int(words[1]) )
                self.bonds.append( bond )
            self.bonds.sort( key = lambda x : x.vmin )
            self.N = len( self.bonds )

    def discard( self ) :
        self.time = -1.0
        self.bonds.clear()
        self.neigh_bonds.clear()
        self.N = 0
