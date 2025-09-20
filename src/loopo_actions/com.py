from .action import *

class COM( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "com"
            self.COM_mode = argv[i+1].strip()
            self.COM_value = -1
            if self.COM_mode != "all" :
                self.COM_value = int( argv[i+2].strip() )
            if argv[i+2].strip() == "file" :
                self.COM_file = open( argv[i+3].strip() , "w" )
            elif argv[i+3].strip() == "file" :
                self.COM_file = open( argv[i+4].strip() , "w" )
            else :
                self.COM_file = open( "com.dat" , "w" )

    def execute( self , equil_config , col ) :
        if self.COM_mode not in [ "all" , "mol" , "type" ] :
            self.lock.acquire()
            print( " *** Com mode not recognized !" )
            self.lock.release()
            exit()
        elif self.COM_mode == "mol" and col.mol == -1 :
            self.lock.acquire()
            print( " *** ERROR: I do not have information about molecules ID !" )
            self.lock.release()
            exit()
        elif self.COM_mode == "type" and col.type == -1 :
            self.lock.acquire()
            print( " *** ERROR: I do not have information about type !" )
            self.lock.release()
            exit()
        COM = equil_config.com( self.COM_mode , self.COM_value )
        self.lock.acquire()
        self.COM_file.write(  str( equil_config.time ) + " " + str( COM[0] ) + " " + str( COM[1] ) + " " + str( COM[2] ) + "\n"  )
        self.lock.release()

    def terminate( self ) :
        self.COM_file.close()

    def return_values( self ) :
        self.lock.acquire()
        self.COM_file.flush()
        self.lock.release()
        return {}

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-com\033[0m <mol|type|all> <mol_id|atom_type> { file <fname_com_perconf> } } " )




class VCOM( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "vcom"
            self.VCOM_mode = argv[i+1].strip()
            self.VCOM_value = -1
            self.VCOM_file = ""
            if self.VCOM_mode != "all" :
                self.VCOM_value = int( argv[i+2].strip() )
            if i+2 < len(argv) :
                if argv[i+2].strip() == "file" :
                    self.VCOM_file = open( argv[i+3].strip() , "w" )
                elif i+3 < len(argv) :
                    if argv[i+3].strip() == "file" :
                        self.VCOM_file = open( argv[i+4].strip() , "w" )
            if self.VCOM_file == "" :
                self.VCOM_file = open( "com_vel.dat" , "w" )

    def execute( self , equil_config , col ) :
        if self.VCOM_mode not in [ "all" , "mol" , "type" ] :
            self.lock.acquire()
            print( " *** Vcom mode not recognized !" )
            self.lock.release()
            exit()
        elif self.VCOM_mode == "mol" and col.mol == -1 :
            self.lock.acquire()
            print( " *** ERROR: I do not have information about molecules ID !" )
            self.lock.release()
            exit()
        elif self.VCOM_mode == "type" and col.type == -1 :
            self.lock.acquire()
            print( " *** ERROR: I do not have information about type !" )
            self.lock.release()
            exit()
        VCOM = equil_config.vcom( self.VCOM_mode , self.VCOM_value )
        self.lock.acquire()
        self.VCOM_file.write(  str( equil_config.time ) + " " + str( VCOM[0] ) + " " + str( VCOM[1] ) + " " + str( VCOM[2] ) + "\n"  )
        self.lock.release()

    def terminate( self ) :
        self.VCOM_file.close()

    def return_values( self ) :
        self.lock.acquire()
        self.VCOM_file.flush()
        self.lock.release()
        return {}

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-vcom\033[0m <mol|type|all> <mol_id|atom_type> { file <fname_com_perconf> } } " )
