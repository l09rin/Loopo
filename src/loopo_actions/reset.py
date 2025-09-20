from .action import *

class RESET( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.chains = [] # contains atom ids
            self.operation = "reset"
            self.ATTRIBUTE = argv[i+1].strip()
            if self.ATTRIBUTE == "type" :
                self.new_value = argv[i+2].strip()
                if argv[i+3].strip() == "chain_id_type" :
                    self.ATTRIBUTE = "type_w_chainID"
                    self.bondslist_file = argv[i+4].strip()
                    self.bonded = self.generate_bondslist()
            elif self.ATTRIBUTE in [ "x" , "y" , "z" ] :
                self.new_value = float(argv[i+2].strip())

    def execute( self , equil_config , col ) :
        if self.ATTRIBUTE == "id" :
            equil_config.reset( self.ATTRIBUTE )
        elif self.ATTRIBUTE == "type" :
            equil_config.reset( self.ATTRIBUTE , self.new_value )
            col.set_type( 10 )
        elif self.ATTRIBUTE in [ "x" , "y" , "z" ] :
            equil_config.reset( self.ATTRIBUTE , self.new_value )
            if self.ATTRIBUTE == "x" and col.x < 0 :
                col.set_x( 10 )
            elif self.ATTRIBUTE == "y" and col.y < 0 :
                col.set_y( 10 )
            elif self.ATTRIBUTE == "z" and col.z < 0 :
                col.set_z( 10 )
        elif self.ATTRIBUTE == "type_w_chainID" :
            indexes = {}
            for i in range( equil_config.N ) :
                indexes[ equil_config.id[i] ] = i
            self.lock.acquire()
            if len(self.chains) == 0 and len(self.bonded) != 0 :
                blist = [ [] for i in range( equil_config.N ) ]
                for ip in self.bonded.keys() :
                    for jp in self.bonded[ip] :
                        blist[ indexes[ip] ].append( indexes[jp] )
                equil_config.build_chains( blist )
                blist.clear()
                for chain in equil_config.chain :
                    mon_ids = []
                    for mon_idx in chain.monomers :
                        mon_ids.append( equil_config.id[mon_idx] )
                    self.chains.append( mon_ids )
            self.lock.release()
            for i in range( len(self.chains) ) :
                for mon_id in self.chains[i] :
                    if len(self.bonded[mon_id]) < 3 :
                        equil_config.type[ indexes[mon_id] ] = i + int(self.new_value)
            indexes.clear()
            col.set_type( 10 )

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-reset\033[0m <id|type|x|y|z> {<new_type|new_x,y,z> {chain_id_type <lmp_bonds_list_file>} } }" )
        print( "\t\t\t\tfor     -reset type <new_type> if the additional key-val argument chain_id_type <lmp_bonds_list_file> is present" )
        print( "\t\t\t\tall the monomers are given the type chain_id + new_type ." )
