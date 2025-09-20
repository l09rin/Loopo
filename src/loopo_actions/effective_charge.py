from .action import *

class Q_EFF( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "effective_charge"
            self.MOL_ID = int(argv[i+1].strip())
            self.METHOD = argv[i+2].strip()
            self.CUTOFF = float(argv[i+3].strip())
            if self.METHOD not in [ "first" , "second" , "third" ] :
                print( "No valid method chosen for effective charge calculation !" )
                exit( EXIT_FAILURE )
            OUT_FNAME="Qeff.dat"
            j = 4
            if len( argv ) > i+j :
                if argv[i+j].strip() == "file" :
                    j += 1
                    OUT_FNAME = argv[i+j].strip()
            self.qeff_file = open( OUT_FNAME , "w" )
            self.qeff_file.write( "# timestep Qeff\n" )
            self.qeff_file.flush()

    def execute( self , equil_config , col ) :
        mol_eff_charge = equil_config.effective_charge( self.MOL_ID , self.METHOD , self.CUTOFF )
        self.lock.acquire()
        self.qeff_file.write( str( equil_config.time ) + " " + str( mol_eff_charge ) + "\n" )
        self.lock.release()

    def terminate( self ) :
        self.qeff_file.close()

    def return_values( self ) :
        self.lock.acquire()
        self.qeff_file.flush()
        self.lock.release()
        return {}

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-Qeff\033[0m <mol_id> <method:first|second|third> <cutoff>  { file <output_fname> }          [output_fname=Qeff.dat]  }" )
