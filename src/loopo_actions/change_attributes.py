from .action import *

class CHANGE_BOX( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "change_box"
            self.BOX_FILENAME = argv[i+1].strip()

    def execute( self , equil_config , col ) :
        equil_config.change_box( self.BOX_FILENAME )

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-change_box\033[0m <lmp-format_data_file-name> }" )


class CHANGE( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "change_attribute"
            self.SELECTION_ATTR = argv[i+1].strip()
            self.NEW_ATTR = argv[i+3].strip()
            if self.SELECTION_ATTR in [ "all" , "mol" , "type" , "q" , "id" , "x" , "y" , "z" ] :
                self.SELECTION_VAL = argv[i+2].strip()
            elif self.SELECTION_ATTR == "xlinkers" :
                self.bondslist_file = argv[i+2].strip()
                self.bonded = self.generate_bondslist()
                self.generate_XL_list()
            elif self.SELECTION_ATTR == "IDfile" :
                IDlist_file = open( argv[i+2].strip() )
                lines = IDlist_file.readlines()
                IDlist_file.close()
                self.SELECTION_VAL = {}
                if self.NEW_ATTR in [ "mol" , "type" , "id" ] :
                    for line in lines :
                        words = line.strip().split()
                        self.SELECTION_VAL[ int( words[0] ) ] = int( words[1] )
                elif self.NEW_ATTR in [ "q" , "m" , "r" ] :
                    for line in lines :
                        words = line.strip().split()
                        self.SELECTION_VAL[ int( words[0] ) ] = float( words[1] )
            else :
                print( "No valid parameter for change attribute !" )
                exit( EXIT_FAILURE )
            if self.SELECTION_ATTR == "IDfile" :
                self.NEW_VAL = 0
            elif self.NEW_ATTR in [ "mol" , "type" , "id" ] :
                self.NEW_VAL = int( argv[i+4].strip() )
            elif self.NEW_ATTR in [ "q" , "m" , "r" ] :
                self.NEW_VAL = float( argv[i+4].strip() )
            else :
                print( "No valid parameter for change attribute !" )
                exit( EXIT_FAILURE )

    def execute( self , equil_config , col ) :
        if self.SELECTION_ATTR == "mol" and col.mol < 0 :
            self.lock.acquire()
            print( " *** ERROR: I do not have information about molecules ID !" )
            self.lock.release()
            exit()
        if self.SELECTION_ATTR == "type" and col.type < 0 :
            self.lock.acquire()
            print( " *** ERROR: I do not have information about atom types !" )
            self.lock.release()
            exit()
        if ( self.SELECTION_ATTR == "id" or self.SELECTION_ATTR == "IDfile" ) and col.id < 0 :
            self.lock.acquire()
            print( " *** ERROR: I do not have information about atom ids !" )
            self.lock.release()
            exit()
        if self.SELECTION_ATTR == "q" and col.q < 0 :
            self.lock.acquire()
            print( " *** ERROR: I do not have information about charge !" )
            self.lock.release()
            exit()
        if ( self.SELECTION_ATTR == "x" and col.x < 0 ) or ( self.SELECTION_ATTR == "y" and col.y < 0 ) or ( self.SELECTION_ATTR == "z" and col.z < 0 ) :
            self.lock.acquire()
            print( " *** ERROR: I do not have information about this position coordinate !" )
            self.lock.release()
            exit()
        if self.NEW_ATTR == "mol" and col.mol < 0 :
            col.set_mol( max( [ col.id , col.type , col.mol , col.q , col.x , col.y , col.z , col.vx , col.vy , col.vz ] ) + 1 )
        if self.NEW_ATTR == "type" and col.type < 0 :
            col.set_type( max( [ col.id , col.type , col.mol , col.q , col.x , col.y , col.z , col.vx , col.vy , col.vz ] ) + 1 )
        if self.NEW_ATTR == "id" and col.id < 0 :
            col.set_id( max( [ col.id , col.type , col.mol , col.q , col.x , col.y , col.z , col.vx , col.vy , col.vz ] ) + 1 )
        if self.NEW_ATTR == "q" and col.q < 0 :
            col.set_q( max( [ col.id , col.type , col.mol , col.q , col.x , col.y , col.z , col.vx , col.vy , col.vz ] ) + 1 )
        if self.SELECTION_ATTR in [ "all" , "mol" , "type" , "id" , "q" , "x" , "y" , "z" ] :
            equil_config.change_attribute( self.SELECTION_ATTR , self.SELECTION_VAL , self.NEW_ATTR , self.NEW_VAL )
        elif self.SELECTION_ATTR == "xlinkers" and len(self.xlinkers_ids) != 0 :
            equil_config.change_attribute( "atoms_list" , self.xlinkers_ids , self.NEW_ATTR , self.NEW_VAL )
        elif self.SELECTION_ATTR == "IDfile" and len(self.SELECTION_VAL) != 0 :
            equil_config.change_attribute( "IDlist" , self.SELECTION_VAL , self.NEW_ATTR , self.NEW_VAL )

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-change\033[0m <sel_attr=id|x|y|z|mol|type|q|xlinkers|IDfile|all> <sel_val(min:max,molID,type,q)>|<lmp_bonds_list_file(xlinkers)|attribute_filename(IDfile)|0(all)> <new_attr=id|type|mol|q|m|r> <new_val|none(IDfile)>    [interval selection convention: min<=val<=max]}" )
