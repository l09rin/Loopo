from action import *

class CLUSTERS_CHECK( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "clusters"
            self.bondslist_file = argv[i+1].strip()
            self.bonded = self.generate_bondslist()
            self.DONE = 0

    def execute( self , equil_config , col ) :
        if self.DONE == 0 :
            bonded_IDXs = equil_config.translate_ID_2_IDX_bondslist( self.bonded , "NOtype" )
            clusters = equil_config.clusters( bonded_IDXs )
            molsXcluster = []
            for i in range( len(clusters) ) :
                molsXcluster.append( set() )
                for idx in clusters[i] :
                    molsXcluster[-1].add( equil_config.mol[idx] )
            self.lock.acquire()
            print( " Molecule IDs present in clusters :" )
            self.lock.release()
            for i in range( len(clusters) ) :
                self.lock.acquire()
                print( "    -> " + str(i+1) + " : " + str(molsXcluster[i]) )
                self.lock.release()
            self.DONE = 1

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-clusters\033[0m <lmp_bonds_list_file>  <<provvisorio>> } " )
