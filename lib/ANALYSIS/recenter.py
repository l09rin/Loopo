from action import *

class RECENTER( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "recenter"
            self.RECENTER_mode = argv[i+1].strip()
            self.RECENTER_value = -1
            if self.RECENTER_mode != "all" :
                self.RECENTER_value = int( argv[i+2].strip() )

    def execute( self , equil_config , col ) :
        if self.RECENTER_mode not in [ "all" , "mol" , "type" ] :
            self.lock.acquire()
            print( " *** Recenter mode not recognized !" )
            self.lock.release()
            exit()
        elif self.RECENTER_mode == "mol" and col.mol == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about molecules ID !" )
                self.lock.release()
                exit()
        elif self.RECENTER_mode == "type" and col.type == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about type !" )
                self.lock.release()
                exit()
        equil_config.displace(  -equil_config.com( self.RECENTER_mode , self.RECENTER_value )  )

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-recenter\033[0m <mol|type|all> <mol_id|atom_type> } " )
