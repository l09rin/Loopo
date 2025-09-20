import numpy as np

# definition of the class particle in D=3
class molecule :

    def __init__( self , id0=0 ) :
        self.id = id0
        self.atom_ids = []
        self.atom_idxs = []
        self.properties = {}
        self.COM = np.array( [0.0,0.0,0.0] )
        self.VCOM = np.array( [0.0,0.0,0.0] )
        self.q = 0.0
        self.PERMANENT = 0
