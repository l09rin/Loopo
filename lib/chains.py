import numpy as np

# definition of the class chain
class chain :

    def __init__( self ) :
        self.N = 0         # number of monomers
        self.Q = 0.0       # total charge
        self.cm = np.array( [0.0,0.0,0.0] )       # center of mass
        self.r = 0.0       # distance of chain's cm from the molecule's center of mass
        self.avg_r = 0.0       # average of the distances  of all of the chain's monomers from the molecule's center of mass
        self.Lee = 0.0     # end-to-end distance
        self.monomers = []         # list of the monomers in the chain
        self.cosines = []          # vector with the average value of the cosine of the angle among bonds at distances 1, 2, ..., N-2
        self.LOOP = False  # True if the chain is a loop
        self.inner_monomers_fraction = 0.0
        self.rg = 0.0

    def add_monomer( self , id_new ) :
        self.monomers.append( id_new )
        self.N += 1

    def periodic_image( self , points , ref , box_sides ) :
        periodic_images = np.sign(points-ref) * (  ( np.fabs(points-ref) // ( box_sides*0.5 ) + 1 ) // 2  )
        return periodic_images

    def gyration( self , particles ) :
        # changed
        self.rg = 0.0
        if self.LOOP :
            ppos = particles[ self.monomers[ 0 : self.N-1 ] ] - self.cm
            self.rg = np.sqrt( np.add.reduce( ppos**2 , (0,1) ) / len(ppos) )
        else :
            ppos = particles[ self.monomers ] - self.cm
            self.rg = np.sqrt( np.add.reduce( ppos**2 , (0,1) ) / len(ppos) )
        return self.rg

    def average_monomers_radii( self , particles , mol_cm ) :
        # changed
        self.avg_r = 0.0
        if self.LOOP :
            ppos = particles[ self.monomers[ 0 : self.N-1 ] ] - mol_cm
            self.avg_r = np.sqrt( np.add.reduce( ppos , 0 ) / len(ppos) )
        else :
            ppos = particles[ self.monomers ] - mol_cm
            self.avg_r = np.sqrt( np.add.reduce( ppos , 0 ) / len(ppos) )
        return self.avg_r

    def com( self , particles , box_side ) :
        # computation of the centre of mass and end2end distance
        endtoend = np.array( [ 0.0 , 0.0 , 0.0 ] )
        self.cm = np.array( [ 0.0 , 0.0 , 0.0 ] )
        self.cm += particles[ self.monomers[0] ]
        if self.LOOP :
            for j in range( self.N-1 ) :
                nearest_copy = particles[ self.monomers[j] ] - ( self.periodic_image( particles[ self.monomers[j] ] , particles[ self.monomers[j-1] ] , box_side ) * box_side )
                endtoend += ( nearest_copy - particles[ self.monomers[j-1] ]  )
                self.cm += (  particles[ self.monomers[0] ] + endtoend  )
            self.cm /= ( self.N - 1 )
            self.Lee = 0.0
        else :
            for j in range( 1 , self.N ) :
                nearest_copy = particles[ self.monomers[j] ] - ( self.periodic_image( particles[ self.monomers[j] ] , particles[ self.monomers[j-1] ] , box_side ) * box_side )
                endtoend += ( nearest_copy - particles[ self.monomers[j-1] ]  )
                self.cm += (  particles[ self.monomers[0] ] + endtoend  )
            self.cm /= self.N
            self.Lee = np.sqrt( np.add.reduce( endtoend**2 , 0 ) )
        return self.cm

    def compute_cosines( self , particles , box_side ) :
        # computation of bond angles
        if self.N > 2 :
            # if the chain is a loop cosines has to be averaged differently
            if self.monomers[0] == self.monomers[-1] :
                num_bonds = self.N - 1
                num_sites = num_bonds
                d_max = 0
                if num_bonds % 2 == 1 :
                    d_max = int( ( num_bonds + 1 ) / 2 )
                else :
                    d_max = int( num_bonds / 2 )
                for d in range( 1 , d_max ) :
                    avg_cos = 0.0
                    for k in range( num_sites ) :
                        nearest_copy = particles[ self.monomers[k+1] ] - ( self.periodic_image( particles[ self.monomers[k+1] ] , particles[ self.monomers[k] ] , box_side ) * box_side )
                        bond_1 = nearest_copy - particles[ self.monomers[k] ]
                        k_plus_d = ( k + d ) % num_sites
                        k_plus_d_plus_1 = ( k + d + 1 ) % num_sites
                        nearest_copy = particles[ self.monomers[k_plus_d_plus_1] ] - ( self.periodic_image( particles[ self.monomers[k_plus_d_plus_1] ] , particles[ self.monomers[k_plus_d] ] , box_side ) * box_side )
                        bond_2 = nearest_copy - particles[ self.monomers[k_plus_d] ]
                        avg_cos += (  np.add.reduce(bond_1*bond_2,0) / np.sqrt(np.add.reduce(bond_1**2,0)) / np.sqrt(np.add.reduce(bond_2**2,0))  )
                    avg_cos /= ( num_sites )
                    self.cosines.append( avg_cos )
                if num_bonds % 2 == 0 :
                    avg_cos = 0.0
                    for k in range( int(num_sites/2) ) :
                        nearest_copy = particles[ self.monomers[k+1] ] - ( self.periodic_image( particles[ self.monomers[k+1] ] , particles[ self.monomers[k] ] , box_side ) * box_side )
                        bond_1 = nearest_copy - particles[ self.monomers[k] ]
                        nearest_copy = particles[ self.monomers[k+d_max+1] ] - ( self.periodic_image( particles[ self.monomers[k+d_max+1] ] , particles[ self.monomers[k+d_max] ] , box_side ) * box_side )
                        bond_2 = nearest_copy - particles[ self.monomers[k+d_max] ]
                        avg_cos += (  np.add.reduce(bond_1*bond_2,0) / np.sqrt(np.add.reduce(bond_1**2,0)) / np.sqrt(np.add.reduce(bond_2**2,0))  )
                    avg_cos /= ( num_sites/2 )
                    self.cosines.append( avg_cos )
            else :
                for d in range( 1 , self.N-1 ) :
                    avg_cos = 0.0
                    for k in range( self.N-1-d ) :
                        nearest_copy = particles[ self.monomers[k+1] ] - ( self.periodic_image( particles[ self.monomers[k+1] ] , particles[ self.monomers[k] ] , box_side ) * box_side )
                        bond_1 = nearest_copy - particles[ self.monomers[k] ]
                        nearest_copy = particles[ self.monomers[k+d+1] ] - ( self.periodic_image( particles[ self.monomers[k+d+1] ] , particles[ self.monomers[k+d] ] , box_side ) * box_side )
                        bond_2 = nearest_copy - particles[ self.monomers[k+d] ]
                        avg_cos += (  np.add.reduce(bond_1*bond_2,0) / np.sqrt(np.add.reduce(bond_1**2,0)) / np.sqrt(np.add.reduce(bond_2**2,0))  )
                    avg_cos /= ( self.N-1-d )
                    self.cosines.append( avg_cos )
