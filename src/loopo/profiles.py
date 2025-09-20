import numpy as np
from . import particle_2D as p2D

# definition of the class profile_point and profile
class profile_point :

    def __init__( self , r=0.0 , y=0.0 , int_y=0.0 ) :
        self.r = r         # shell mid-point
        self.rho = y       # density value
        self.cumulative = int_y       # spherically integrated density
        self.aux = 0.0     # auxiliary variable



# definition of the class density_profile
class density_profile :

    def __init__( self , shell_thickness , box_side ) :
        self.rmin = []
        self.rmax = []
        self.counts = []
        self.dens = []
        self.Nbins = 0
        if isinstance( box_side , p2D.particle ) :
            Rmax = (box_side*0.5).min_comp()
        else :
            Rmax = np.amin( box_side*0.5 )
        ri = 0.0
        while ri + shell_thickness < Rmax :
            self.rmin.append( ri )
            self.rmax.append( ri + shell_thickness )
            self.counts.append( 0 )
            self.dens.append( 0.0 )
            ri = ri + shell_thickness

        self.rmin.append( ri )
        if isinstance( box_side , p2D.particle ) :
            volume = box_side.vol()
            self.rmax.append( ( volume / np.pi )**(1./2) )
        else :
            volume = np.multiply.reduce( box_side )
            if len(box_side) == 2 or volume == 0 :
                box_dims = np.where( box_side > 0 , box_side , 1.0 )
                volume = np.multiply.reduce( box_dims )
                self.rmax.append( ( volume / np.pi )**(1./2) )
            else :
                self.rmax.append( ( 3. * volume / 4. / np.pi )**(1./3) )
        self.counts.append( 0 )
        self.dens.append( 0.0 )
        self.Nbins = len( self.rmin )

    def total_count( self ) :
        Ntot = 0
        for i in range( self.Nbins ) :
            Ntot += self.counts[i]
        return Ntot


class cylindrical_profile :

    def __init__( self , box_side , cylinder_axis , radial_bin , height_bin ) :
        self.rmin = []
        self.rmax = []
        self.hmin = []
        self.hmax = []
        self.counts = []
        self.dens = []
        self.lastHbin = 0
        self.Nrbins = 0
        # case of hidden 2D systems
        embedded2D = False
        if box_side[2] == 0.0 :
            embedded2D = True
            height_bin = 1.0
            if cylinder_axis == "x" or cylinder_axis == "y" :
                print( "*** WARNING: 2D system, no sense for cylindrical profiles along x or y axis!" )
        if cylinder_axis == "x" :
            Hmax = box_side[0]*0.5
            Rmax = min( box_side[1]*0.5 , box_side[2]*0.5 )
        elif cylinder_axis == "y" :
            Hmax = box_side[1]*0.5
            Rmax = min( box_side[0]*0.5 , box_side[2]*0.5 )
        elif cylinder_axis == "z" :
            if embedded2D :
                Hmax = 0.5
            else :
                Hmax = box_side[2]*0.5
            Rmax = min( box_side[0]*0.5 , box_side[1]*0.5 )
        else :
            print( "*** For cylindrical profiles the cylinder axis has to be x, y or z !" )
            exit(-1)

        self.generate_height( height_bin , Hmax )
        self.generate_radii( radial_bin , Rmax , box_side , cylinder_axis )
        self.generate_profile_array()

    def generate_height( self , height_bin , Hmax ) :
        hi = -0.5 * height_bin
        self.hmin.append( hi )
        self.hmax.append( -hi )
        hi = hi + height_bin
        while hi + height_bin < Hmax :
            self.hmin.append( hi )
            self.hmax.append( hi + height_bin )
            hi = hi + height_bin
        if self.hmax[-1] > Hmax :
            self.hmax[-1] = Hmax
        self.lastHbin = len( self.hmax ) - 1
        i = len( self.hmax ) - 1
        while i > 0 :
            self.hmin.append( -self.hmax[i] )
            self.hmax.append( -self.hmin[i] )
            i -= 1

    def generate_radii( self , radial_bin , Rmax , box_side , cylinder_axis ) :
        ri = 0.0
        while ri + radial_bin < Rmax :
            self.rmin.append( ri )
            self.rmax.append( ri + radial_bin )
            ri = ri + radial_bin
        if cylinder_axis == "x" :
            area = box_side[1] * box_side[2]
        elif cylinder_axis == "y" :
            area = box_side[0] * box_side[2]
        elif cylinder_axis == "z" :
            area = box_side[1] * box_side[0]
        self.rmax[-1] = ( area / np.pi )**(1./2)
        self.Nrbins = len( self.rmin )

    def generate_profile_array( self ) :
        for i in range( len(self.hmin) ) :
            self.counts.append( [] )
            self.dens.append( [] )
            for j in range( len(self.rmin) ) :
                self.counts[i].append( 0 )
                self.dens[i].append( 0.0 )

    def total_count( self ) :
        Ntot = 0
        for i in range( -self.lastHbin , self.lastHbin+1 ) :
            for j in range( self.Nrbins ) :
                Ntot += self.counts[i][j]
        return Ntot


class linear_profile :

    def __init__( self , box_side , axis , bin_width ) :
        self.Nbins = 0
        self.hmin = []
        self.hmax = []
        self.counts = []
        self.dens = []
        if axis == "x" :
            Hmax = box_side[0]
        elif axis == "y" :
            Hmax = box_side[1]
        elif axis == "z" :
            Hmax = box_side[2]
        else :
            print( "*** For linear profiles the axis has to be x, y or z !" )
            exit(-1)

        hi = 0.0
        while hi + bin_width < Hmax :
            self.hmin.append( hi )
            self.hmax.append( hi + bin_width )
            self.counts.append( 0 )
            self.dens.append( 0.0 )
            hi = hi + bin_width
        self.hmin.append( hi )
        self.hmax.append( Hmax )
        self.counts.append( 0 )
        self.dens.append( 0.0 )
        self.Nbins = len( self.hmin )

    def total_count( self ) :
        Ntot = 0
        for i in range( self.Nbins ) :
            Ntot += self.counts[i]
        return Ntot
