from action import *

import scipy.linalg as la
import scipy.integrate as integration

class CONVEX_HULL( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "convex_hull"
            self.mode = argv[i+1].strip()
            self.value = -1
            j = 2
            if self.mode != "all" :
                self.value = int( argv[i+j].strip() )
                j += 1
            if len( argv ) > i+j :
                if argv[i+j].strip() == "file" :
                    j += 1
                    self.outfile = open( argv[i+j].strip() , "w" )
                else :
                    self.outfile = open( "convex_hulls.dat" , "w" )
            else :
                self.outfile = open( "convex_hulls.dat" , "w" )
            self.outfile.write( "# timestep (semiaxis of the intertia tensor) Rg volume Rh numerical_error asphericity_1=(1.5*(a^2+b^2+c^2)/(a+b+c)^2-0.5) asphericity_2=(0.5*((l2-l1)^2+(l3-l1)^2+(l3-l2)^2)/(l1+l2+l3)^2),l1=a^2/3,l2=b^2/3,l3=c^2/3\n" )
            self.outfile.flush()

    def execute( self , equil_config , col ) :
        if self.mode == "all" :
            hull, facets_com, inertia_tensor = equil_config.convexHull( )
        elif self.mode == "mol" :
            if col.mol != -1 :
                hull, facets_com, inertia_tensor = equil_config.convexHull( self.mode , self.value )
            else :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about molecules ID !" )
                self.lock.release()
                exit()
        elif self.mode == "type" :
            if col.type != -1 :
                hull, facets_com, inertia_tensor = equil_config.convexHull( self.mode , self.value )
            else :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about type !" )
                self.lock.release()
                exit()
        else :
            self.lock.acquire()
            print( " *** mode not recognized !" )
            self.lock.release()
            exit()
        eigvals, eigvecs = la.eig( inertia_tensor )
        eigvals = np.sort( eigvals , kind='mergesort' )
        eigvals = eigvals[::-1]
        Rg = sqrt( np.real(eigvals[0]) + np.real(eigvals[1]) + np.real(eigvals[2]) )
        a = sqrt( 3.0 * np.real(eigvals[0]))
        b = sqrt( 3.0 * np.real(eigvals[1]))
        c = sqrt( 3.0 * np.real(eigvals[2]))
        f = lambda theta : 1.0 / sqrt( (a**2+theta) * (b**2+theta) * (c**2+theta) )
        ( integral , error_int ) = integration.quad( f , 0 , np.inf )
        l1 = a**2/3.
        l2 = b**2/3.
        l3 = c**2/3.
        asphericity_1 = 1.5 * (a**2+b**2+c**2)/(a+b+c)**2 - 0.5
        asphericity_2 = 0.5 * ( (l2-l1)**2 + (l3-l1)**2 + (l3-l2)**2 ) / (l1+l2+l3)**2
        self.lock.acquire()
        self.outfile.write( str( equil_config.time ) + " " + repr(a) + " " + repr(b) + " " + repr(c) + " " + repr(Rg) + " " + repr(hull.volume) + " " + repr(2.0/integral) + " " + repr(2.0/integral**2*error_int) + " " + repr(asphericity_1) + " " + repr(asphericity_2) + "\n" )
        self.outfile.flush()
        self.lock.release()

    def terminate( self ) :
        self.outfile.close()

    def return_values( self ) :
        self.lock.acquire()
        self.outfile.flush()
        self.lock.release()
        return {}

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-convex_hull\033[0m <mol|type|all> <mol_id|atom_type> { file <output_fname> }          [output_fname=convex_hulls.dat] } " )
