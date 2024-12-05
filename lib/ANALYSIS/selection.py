from action import *

class SELECT( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "select"
            self.TYPE = argv[i+1].strip()
            if self.TYPE in [ "mol" , "type" ] :
                self.VALUE = int( argv[i+2].strip() )
            elif self.TYPE == "q" :
                self.VALUE = float( argv[i+2].strip() )
            elif self.TYPE == "xlinkers" :
                self.bondslist_file = argv[i+2].strip()
                self.bonded = self.generate_bondslist()
                self.VALUE = self.bonded
            elif self.TYPE == "IDfile" :
                idfile = open( argv[i+2].strip() , "r" )
                self.VALUE = np.loadtxt( idfile , dtype = 'str' ).astype('int32')
                idfile.close()
            else :
                print( "No valid choice for selection type !" )
                exit( EXIT_FAILURE )

    def execute( self , equil_config , col ) :
        if self.TYPE == "mol" and col.mol == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about molecules ID !" )
                self.lock.release()
                self.lock.release()
                exit()
        elif self.TYPE == "type" and col.type == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about atom types !" )
                self.lock.release()
                exit()
        elif self.TYPE == "IDfile" and col.id == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about atom ids !" )
                self.lock.release()
                exit()
        elif self.TYPE == "q" and col.q == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about charge !" )
                self.lock.release()
                exit()
        elif self.TYPE == "xlinkers" and self.bondslist_file == "" :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about bonds !" )
                self.lock.release()
                exit()
        equil_config.select( self.TYPE , self.VALUE )

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-sel\033[0m <mol|type|q> <selecting_value> } { -sel xlinkers <lmp_bonds_list_file> }" )
        print( "\t\t\t\t{ -sel IDfile <filename> [file of 1 col with ids to select] }" )

class REMOVE( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "remove"
            self.TYPE = argv[i+1].strip()
            self.VALUE = ""
            self.max2remove = 0
            self.MODE = "none"
            self.RADIUS = 0.0
            self.bonded = {}
            if self.TYPE not in [ "mol" , "type" , "q" , "xlinkers" , "sphere" , "by_dist" ] :
                print( "No valid choice for remove type !" )
                exit( EXIT_FAILURE )
            if self.TYPE == "xlinkers" :
                self.bondslist_file = argv[i+2].strip()
                self.bonded = self.generate_bondslist()
                self.VALUE = self.bonded
            elif self.TYPE == "sphere" :
                self.MODE = argv[i+2].strip()
                self.VALUE = int( argv[i+3].strip() )
                self.RADIUS = float( argv[i+4].strip() )
            elif self.TYPE == "by_dist" :
                self.MODE = argv[i+2].strip()
                self.CUTOFF_DIST = float( argv[i+3].strip() )
                self.RM_MODE = argv[i+4].strip()
                self.VALUE = argv[i+5].strip()
                self.REFERENCE_MODE = argv[i+6].strip()
                self.REFERENCE_VALUE = argv[i+7].strip()
            elif self.TYPE in [ "mol" , "type" ] :
                self.VALUE = int( argv[i+2].strip() )
            elif self.TYPE == "q" :
                self.VALUE = float( argv[i+2].strip() )
            if self.TYPE == "type" :
                if i+3 < len(argv) :
                    if argv[i+3].strip() == "maxN" :
                        self.max2remove = int( argv[i+4].strip() )

    def execute( self , equil_config , col ) :
        if self.TYPE == "mol" and col.mol == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about molecules ID !" )
                self.lock.release()
                exit()
        elif self.TYPE == "type" and col.type == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about atom types !" )
                self.lock.release()
                exit()
        elif self.TYPE == "q" and col.q == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about charge !" )
                self.lock.release()
                exit()
        elif self.TYPE == "sphere" :
            if self.VALUE != -1 and col.type == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about atom types !" )
                self.lock.release()
                exit()
        elif self.TYPE == "xlinkers" and self.bondslist_file == "" :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about bonds !" )
                self.lock.release()
                exit()
        if self.TYPE == "by_dist" :
            if ( self.RM_MODE == "type" or self.REFERENCE_MODE == "type" ) and col.type == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about atom types !" )
                self.lock.release()
                exit()
            else :
                equil_config.remove_by_dist( self.CUTOFF_DIST , self.MODE , self.RM_MODE , self.VALUE , self.REFERENCE_MODE , self.REFERENCE_VALUE )
        else :
            equil_config.remove( self.TYPE , self.VALUE , self.max2remove , self.RADIUS , self.MODE )

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-rm\033[0m <mol|type|q> <removing_value> } { -rm xlinkers <lmp_bonds_list_file> }" )
        print( "\t\t\t\tto remove particles by particle type, an optional key-val argument can be used:" )
        print( "\t\t\t\t   > -rm type <removing_value> {maxN <N_parts_2be_removed>}  ," )
        print( "\t\t\t\tthat allows to specify the maximum number of particles to remove" )
        print( "\t\t\t{ \033[1m-rm\033[0m sphere in|out <type> <r> }                (type=-1 for all types)" )
        print( "\t\t\t{ \033[1m-rm\033[0m by_dist grtr|less <cutoff_distance> type <removing_value> type <reference_type> }" )
