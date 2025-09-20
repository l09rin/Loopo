from .action import *

class SHIFT_TIME( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "shift_time"
            self.T0 = float( argv[i+1].strip() )

    def execute( self , equil_config , col ) :
        equil_config.time += self.T0

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-shift_time\033[0m <t0> } " )
