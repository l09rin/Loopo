from action import *

# FREUD # for order parameters calculation
try :
    from freud.order import Hexatic
    from freud.locality import NeighborList
    from freud.box import Box
except :
    print("***WARNING: module freud not found! This is needed to calculate order parameters [install package: freud-analysis].")
    pass

class HEXATIC( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "hexatic"
            orders = argv[i+1].strip().split(",")
            self.ORDERS = []
            for w in orders :
                w2 = w.split(":")
                if len(w2) == 1 :
                    if 1 < int(w) < 13 :
                        self.ORDERS.append(int(w))
                elif len(w2) == 2 :
                    if 1 < int(w2[0]) < 13 and 1 < int(w2[1]) < 14 and int(w2[0]) < int(w2[1]) :
                        self.ORDERS += list(range( int(w2[0]) , int(w2[1]) ))
            if len(self.ORDERS) == 0 :
                print( "*** ERROR: Invalid order value(s)!" )
                exit(EXIT_FAILURE)
            self.CUTOFFS = np.array([])
            j = 2
            clines = []
            if i+j < len(argv) :
                while argv[i+j].strip() == "c" :
                    clines.append( argv[i+j+1].strip().split(":") )
                    if i+j+2 < len(argv) :
                        j += 2
                    else :
                        j += 1
            if len(clines) == 0 :
                print( "*** ERROR: Need to specify at least one cutoff value!" )
                exit(30)
            clines = np.array( clines )
            typemax = np.max( clines[:,0:2].astype(int) )
            # Store for each pair of types the square of the cutoff distance.
            self.CUTOFFS = np.zeros((typemax + 1, typemax + 1))
            for cline in clines :
                self.CUTOFFS[int(cline[0]), int(cline[1])] = float(cline[2])
                self.CUTOFFS[int(cline[1]), int(cline[0])] = float(cline[2])
            self.CUTOFFS = self.CUTOFFS**2
            self.MODE = "all"
            self.VALUE = -1
            if i+j < len(argv) :
                if argv[i+j].strip() in [ "mol" , "type" ] :
                    self.MODE = argv[i+j].strip()
                    self.VALUE = int( argv[i+j+1].strip() )
                    j += 2
            filename = "order_parameters.dat"
            if i+j < len(argv) :
                if argv[i+j].strip() == "file" :
                    filename = argv[i+j+1].strip()
                    j += 2
            self.FILE = open( filename , "w" )
            self.AVGD_FILE = open( filename.rsplit(".",1)[0]+"_avgd."+filename.rsplit(".",1)[1] , "w" )
            self.FILE.write( "# t" )
            self.AVGD_FILE.write( "# t" )
            for k in self.ORDERS :
                self.FILE.write( " "+str(k) )
                self.AVGD_FILE.write( " "+str(k) )
            self.FILE.write("\n")
            self.FILE.flush()
            self.AVGD_FILE.write("\n")
            self.AVGD_FILE.flush()
            self.SAVE_PARTICLE_VALUES = False
            if i+j < len(argv) :
                if argv[i+j].strip() == "save_confs" :
                    self.SAVE_PARTICLE_VALUES = True
                    self.CONF_FNAME = filename.rsplit(".",1)[0]

    def execute( self , equil_config , col ) :
        if self.MODE == "mol" and col.mol == -1 :
            self.lock.acquire()
            print( " *** ERROR: I do not have information about molecules ID !" )
            self.lock.release()
            exit()
        elif self.MODE == "type" and col.type == -1 :
            self.lock.acquire()
            print( " *** ERROR: I do not have information about type !" )
            self.lock.release()
            exit()
        box_side = equil_config.box_sup - equil_config.box_inf
        box = Box(box_side[0],box_side[1],box_side[2])
        equil_config.create_contacts( self.CUTOFFS )
        neighlist = equil_config.list_of_contacts( symmetric=True )
        uniques = np.unique( neighlist[:,0] )
        nearest_img = equil_config.pos[ neighlist[:,1] ] - ( equil_config.periodic_image( equil_config.pos[ neighlist[:,1] ] , equil_config.pos[ neighlist[:,0] ] , box_side ) * box_side )
        # neighlist = NeighborList.from_arrays( len(uniques), len(uniques), neighlist[:,0], neighlist[:,1], nearest_img - equil_config.pos[ neighlist[:,0] ] )
        neighlist = NeighborList.from_arrays( equil_config.N, equil_config.N, neighlist[:,0], neighlist[:,1], nearest_img - equil_config.pos[ neighlist[:,0] ] )

        op_values = []
        avgd_op_values = []
        for k in self.ORDERS :
            hex_order = Hexatic(k=k)
            hex_order.compute( system=(box_side[0:2], equil_config.pos) , neighbors=neighlist )
            op_values.append( hex_order.particle_order )
            avgd_particle_order = np.copy( hex_order.particle_order )
            for i in range(equil_config.N) :
                if len(equil_config.contacts[i]) > 0 :
                    avgd_particle_order[i] = (avgd_particle_order[i] + np.sum(avgd_particle_order[equil_config.contacts[i]])) / (len(equil_config.contacts[i])+1) 
            avgd_op_values.append( avgd_particle_order )
            # hex_order.plot()
        if self.MODE == "type" :
            for i in range(len(op_values)) :
                op_values[i] = op_values[i][ equil_config.type == self.VALUE ]
                avgd_op_values[i] = avgd_op_values[i][ equil_config.type == self.VALUE ]
        elif self.MODE == "mol" :
            for i in range(len(op_values)) :
                op_values[i] = op_values[i][ equil_config.mol == self.VALUE ]
                avgd_op_values[i] = avgd_op_values[i][ equil_config.mol == self.VALUE ]
        print_values = ""
        for val in np.mean( np.abs(op_values), axis=1 ) :
            print_values += (" "+repr(val))
        print_values_avgd = ""
        for val in np.mean( np.abs(avgd_op_values), axis=1 ) :
            print_values_avgd += (" "+repr(val))
        self.lock.acquire()
        self.FILE.write(  str( equil_config.time ) + print_values + "\n"  )
        self.AVGD_FILE.write(  str( equil_config.time ) + print_values_avgd + "\n"  )
        self.lock.release()

        if self.SAVE_PARTICLE_VALUES :
            np.savetxt( self.CONF_FNAME+"-"+str(equil_config.time)+".dat" , np.abs(np.transpose(op_values)) , fmt='%.16g', delimiter=' ', newline='\n', header=str(list(np.array(self.ORDERS))).strip("[]").replace(',','')+" ; timestep "+str(equil_config.time) )
            np.savetxt( self.CONF_FNAME+"_avgd-"+str(equil_config.time)+".dat" , np.abs(np.transpose(avgd_op_values)) , fmt='%.16g', delimiter=' ', newline='\n', header=str(list(np.array(self.ORDERS))).strip("[]").replace(',','')+" ; timestep "+str(equil_config.time) )

    def terminate( self ) :
        self.FILE.close()
        self.AVGD_FILE.close()

    def return_values( self ) :
        self.lock.acquire()
        self.FILE.flush()
        self.AVGD_FILE.flush()
        self.lock.release()
        return {}

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-hexatic\033[0m <order(s)_e.g. 4,8:12>"
               " { c <type1>:<type2>:<cutoff> ...}"
               " { <mol|type> <mol_id|atom_type> }"
               " { file <fname_com_perconf> }  { save_confs } } " )
