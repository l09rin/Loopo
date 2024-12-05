from action import *
import numpy as np

class ADD_PARTICLES( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "add_particles"
            self.PARTICLES = argv[i+1].strip()
            if self.PARTICLES not in [ "coms" ] :
                print( " *** WARNING : It is not possible to add " + self.PARTICLES + " particles to the system !!!" )
                self.PARTICLES = "none"

    def execute( self , equil_config , col ) :
        if self.PARTICLES == "coms" :
            mols = []
            max_type = np.amax( equil_config.type )
            max_id = np.amax( equil_config.id )
            mols = np.unique( equil_config.mol )
            mols = mols[ np.not_equal( mols , 0 ) ]
            max_type += 1
            max_id += 1
            coms = []
            vcoms = []
            for mol in mols :
                coms.append( equil_config.com( "mol" , mol ) )
                if len( equil_config.vel ) > 0 :
                    vcoms.append( equil_config.vcom( "mol" , mol ) )
            npos = np.array( coms )
            nvel = np.array( vcoms )
            nid = np.arange( max_id , max_id + len(mols) , 1 )
            ntype = np.array( [max_type]*len(mols) )
            nq = np.array( [0.0]*len(mols) )
            nmol = np.array( [0]*len(mols) )
            equil_config.insert( npos , nid , ntype , nmol , nq , nvel )

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-add_particles\033[0m coms } " )
