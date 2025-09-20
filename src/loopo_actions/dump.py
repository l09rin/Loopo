from .action import *

class DUMP( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "print-confs"
            self.FILENAME = argv[i+1].strip()
            self.FORMAT = argv[i+2].strip()
            if self.FORMAT not in [ "lmp" , "xyz" , "sph" ] :
                print("*** ERROR : -print option needs a valid format type !")
                exit
            if self.FORMAT == "zeno_input" :
                self.OUT_FILE = open( self.FILENAME+".traj" , "w" )
                bod_file = open( self.FILENAME+".bod" , "w" )
                bod_file.write( "TRAJECTORY " + self.FILENAME + ".traj " + self.FILENAME + ".conv\n" )
                bod_file.close()
                if i+3 >= len(argv) :
                    j = 2
                else :
                    j = 3
                if argv[i+j].strip() == "radius" :
                    conv_file = open( self.FILENAME+".conv" , "w" )
                    while argv[i+j].strip() == "radius" :
                        conv_file.write( argv[i+j+1].strip() + " " + argv[i+j+2].strip() + "\n" )
                        if i+j+3 >= len(argv) :
                            j += 2
                        else :
                            j += 3
                    conv_file.close()
                self.ATTRIBUTES = [ "type" , "pos" ]
                if argv[i+j].strip() in [ "all" , "type" , "mol" , "q" ] :
                    self.SELECTION = argv[i+j].strip()
                    if self.SELECTION == "all" :
                        self.SEL_VALUE = "-1"
                    else :
                        self.SEL_VALUE = argv[i+j+1].strip()
                else :
                    self.SELECTION = "all"
                    self.SEL_VALUE = "-1"
            else :
                self.OUT_FILE = open( self.FILENAME , "w" )
                self.ATTRIBUTES = argv[i+3].strip().split(":")
                for att in self.ATTRIBUTES :
                    if att not in [ "pos" , "vel" , "mol" , "type" , "q" , "id" ] :
                        print( "*** Attribute " + att + " not recognized!" )
                        exit( EXIT_FAILURE )
                self.SELECTION = argv[i+4].strip()
                if self.SELECTION == "all" :
                    self.SEL_VALUE = "-1"
                else :
                    self.SEL_VALUE = argv[i+5].strip()

    def execute( self , equil_config , col ) :
        print_col = COLUMNS()
        if "pos" in self.ATTRIBUTES :
            print_col.set_x( col.x )
            print_col.set_y( col.y )
            print_col.set_z( col.z )
        else :
            print_col.set_x( -1 )
            print_col.set_y( -1 )
            print_col.set_z( -1 )
        if "vel" in self.ATTRIBUTES :
            print_col.set_vx( col.vx )
            print_col.set_vy( col.vy )
            print_col.set_vz( col.vz )
        else :
            print_col.set_vx( -1 )
            print_col.set_vy( -1 )
            print_col.set_vz( -1 )
        if "id" in self.ATTRIBUTES :
            print_col.set_id( col.id )
        else :
            print_col.set_id( -1 )
        if "type" in self.ATTRIBUTES :
            print_col.set_type( col.type )
        else :
            print_col.set_type( -1 )
        if "mol" in self.ATTRIBUTES :
            print_col.set_mol( col.mol )
        else :
            print_col.set_mol( -1 )
        if "q" in self.ATTRIBUTES :
            print_col.set_q( col.q )
        else :
            print_col.set_q( -1 )
        self.lock.acquire()
        equil_config.print_selection( self.OUT_FILE , self.SELECTION , self.SEL_VALUE , self.FORMAT , print_col.id , print_col.type , print_col.mol , print_col.q , print_col.x , print_col.y , print_col.z , print_col.vx , print_col.vy , print_col.vz )
        self.OUT_FILE.flush()
        self.lock.release()

    def terminate( self ) :
        self.OUT_FILE.close()

    def return_values( self ) :
        self.lock.acquire()
        self.OUT_FILE.flush()
        self.lock.release()
        return {}

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-print\033[0m <output_file_name> <format:lmp|xyz|sph> <attributes> <all|mol|type|q> <sel_value> } " )
        print( "\t\t\t\tattributes is a list of strings separated by character \":\" among pos, vel, mol, type, q, id " )
        print( "\t\t\t{ \033[1m-print\033[0m <output_file_name_w/out_extension> zeno_input { radius <type> <r> . . . } { <all(default)|mol|type|q> <sel_value> } } " )
        print( "\t\t\t\tif radii are provided a file .conv is created, and all the type radii have to be specified! " )
