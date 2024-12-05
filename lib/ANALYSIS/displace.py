from action import *
import numpy as np

class DISPLACE( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "displace_attribute"
            self.SELECTION_ATTR = argv[i+1].strip()
            self.ATTR2DISP = argv[i+3].strip()
            if self.SELECTION_ATTR in [ "mol" , "type" ] :
                self.SELECTION_VAL = int( argv[i+2].strip() )
            elif self.SELECTION_ATTR == "q" :
                self.SELECTION_VAL = float( argv[i+2].strip() )
            elif self.SELECTION_ATTR == "id" :
                self.SELECTION_VAL = argv[i+2].strip()
            elif self.SELECTION_ATTR == "xlinkers" :
                self.bondslist_file = argv[i+2].strip()
                self.bonded = self.generate_bondslist()
                self.generate_XL_list()
            elif self.SELECTION_ATTR == "all" :
                self.SELECTION_VAL = 0
            else :
                print( "No valid parameter for change attribute !" )
                exit( EXIT_FAILURE )
            if self.ATTR2DISP in [ "mol" , "type" , "id" ] :
                self.NEW_VAL = int( argv[i+4].strip() )
            elif self.ATTR2DISP == "pos" :
                point = argv[i+4].strip().strip("()").split(',')
                self.NEW_VAL = np.array( [ float(point[0]) , float(point[1]) , float(point[2]) ] )
            else :
                print( "No valid parameter for change attribute !" )
                exit( EXIT_FAILURE )

    def execute( self , equil_config , col ) :
        if self.SELECTION_ATTR in [ "mol" , "type" , "id" , "q" , "all" ] :
            if self.SELECTION_ATTR == "mol" and col.mol == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about molecules ID !" )
                self.lock.release()
                exit()
            if self.SELECTION_ATTR == "type" and col.type == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about atom types !" )
                self.lock.release()
                exit()
            if self.SELECTION_ATTR == "id" and col.id == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about atom ids !" )
                self.lock.release()
                exit()
            if self.SELECTION_ATTR == "q" and col.q == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about charge !" )
                self.lock.release()
                exit()
            if self.ATTR2DISP == "mol" and col.mol == -1 :
                self.lock.acquire()
                print( "*** ERROR: There are no molecular ids in the configuration !" )
                self.lock.release()
                exit()
            if self.ATTR2DISP == "type" and col.type == -1 :
                self.lock.acquire()
                print( "*** ERROR: There are no types in the configuration !" )
                self.lock.release()
                exit()
            if self.ATTR2DISP == "id" and col.id == -1 :
                self.lock.acquire()
                print( "*** ERROR: There are no ids in the configuration !" )
                self.lock.release()
                exit()
            if self.ATTR2DISP == "pos" and ( col.x == -1 or col.y == -1 ) :
                self.lock.acquire()
                print( "*** ERROR: There are no x,y positions in the configuration !" )
                self.lock.release()
                exit()
            elif col.z == -1 :
                self.lock.acquire()
                print( "*** ERROR: There are no z coordinates in the configuration !" )
                self.lock.release()
                exit()
            equil_config.displace_attribute( self.SELECTION_ATTR , self.SELECTION_VAL , self.ATTR2DISP , self.NEW_VAL )
        elif self.SELECTION_ATTR == "xlinkers" and len(self.xlinkers_ids) != 0 :
            equil_config.displace_attribute( "atoms_list" , self.xlinkers_ids , self.ATTR2DISP , self.NEW_VAL )

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-displace\033[0m <sel_attr=id|mol|type|q|xlinkers|all> <sel_val(IDmin:IDmax,mol,type,q,0_for_all)>|<lmp_bonds_list_file(xlinkers)>" )
        print("\t\t\t\t<attr2displace=id|type|mol|pos> <offset|dx,dy,dz> }" )
