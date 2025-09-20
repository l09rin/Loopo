from .action import *

class GYRATION( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "gyration"
            self.GYRATION_mode = argv[i+1].strip()
            self.GYRATION_value = -1
            self.RG = 0.0
            self.N = 0
            j = 1
            if self.GYRATION_mode != "all" :
                j = 2
                self.GYRATION_value = int( argv[i+j].strip() )
            j += 1
            if i+j < len(argv) :
                if argv[i+j].strip() == "files" :
                    self.RGperconf_file = open( argv[i+j+1].strip() , "w" )
                    self.RGavg_file = open( argv[i+j+2].strip() , "w" )
                    j += 2
                else :
                    self.RGperconf_file = open( "gyration_radii.dat" , "w" )
                    self.RGavg_file = open( "gyration.dat" , "w" )
                j += 1
            else :
                self.RGperconf_file = open( "gyration_radii.dat" , "w" )
                self.RGavg_file = open( "gyration.dat" , "w" )

    def execute( self , equil_config , col ) :
        if self.GYRATION_mode not in [ "all" , "mol" , "type" ] :
            self.lock.acquire()
            print( " *** Gyration mode not recognized !" )
            self.lock.release()
            exit()
        elif self.GYRATION_mode == "mol" and col.mol == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about molecules ID !" )
                self.lock.release()
                exit()
        elif self.GYRATION_mode == "type" and col.type == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about type !" )
                self.lock.release()
                exit()
        rg = equil_config.gyration( equil_config.com( self.GYRATION_mode , self.GYRATION_value ) , self.GYRATION_mode , self.GYRATION_value  )
        self.lock.acquire()
        self.RG += rg
        self.N += 1
        self.RGperconf_file.write(  str( equil_config.time ) + " " + str( rg ) + "\n"  )
        self.lock.release()

    def terminate( self ) :
        self.RGperconf_file.close()
        self.RGavg_file.write( "# Time-averaged data for gyration radius\n" )
        self.RGavg_file.write( "# number_of_confs average_gyration_radius\n" )
        self.RGavg_file.write( str(self.N) + '\t' + str(self.RG/self.N) + '\n' )
        self.RGavg_file.close()

    def return_values( self ) :
        self.lock.acquire()
        self.RGperconf_file.flush()
        self.lock.release()
        return { "RG" : self.RG , "N" : self.N }

    def merge_return_values( self , values_list ) :
        for partial in values_list :
            self.RG += partial["RG"]
            self.N += partial["N"]

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-gyration\033[0m <mol|type|all> <mol_id|atom_type> { files <fname_rg_perconf> <fname_rg_average> } } " )
