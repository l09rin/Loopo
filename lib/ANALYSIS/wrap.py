from action import *

class REWRAP( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "rewrap"

    def execute( self , equil_config , col ) :
        equil_config.wrap()

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-rewrap\033[m } " )


class UNWRAP( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "unwrap"
            self.bondslist_file = argv[i+1].strip()
            self.bonded = self.generate_bondslist()

    def execute( self , equil_config , col ) :
        equil_config.unwrap_v2( self.bonded )

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-unwrap\033[0m <lmp_bonds_list_file> }" )
