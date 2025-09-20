from .action import *
import numpy as np

class COM_MOL_DISTANCE( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "com_mol_dist"
            self.FIRST_CONF = True
            self.avg_pos = np.zeros( (1,3) , dtype=float )
            self.avg_distances = np.zeros( (1) , dtype=float )
            self.sorted_ids = np.zeros( (1) , dtype=int )
            self.Nconfs_analysed = 0
            fname_root = "com_mol_distances.dat"
            self.AVG_file = open( fname_root , "w" )
            self.AVG_file.write( "# id xyz_avg_pos avg_dist\n" )
            self.AVG_file.flush()

    def execute( self , equil_config , col ) :
        if len( equil_config.mol ) == 0 :
            self.lock.acquire()
            print( " *** I need molecule IDs !" )
            self.lock.release()
            exit()
        if self.FIRST_CONF :
            self.sorted_ids = np.sort( equil_config.id )
            self.avg_pos = np.zeros( (equil_config.N,3) , dtype=float )
            self.avg_distances = np.zeros( (equil_config.N) , dtype=float )
            self.FIRST_CONF = False
        mols, mols_idx = np.unique( equil_config.mol , return_inverse=True )
        mol_coms = np.zeros( (len(mols),3) )
        for i in range(len(mols)) :
            mol_coms[i,:] = equil_config.com( "mol" , mols[i] )
        pos = equil_config.unwrapped_coordinates()
        com_mol_distances = pos - mol_coms[ mols_idx ]
        self.avg_pos = self.avg_pos + com_mol_distances[ np.argsort( equil_config.id ) ]
        self.avg_distances = self.avg_distances + np.sqrt( np.add.reduce( com_mol_distances[ np.argsort( equil_config.id ) ]**2 , 1 ) )
        self.Nconfs_analysed += 1

    def terminate( self ) :
        # avg_norm = np.sqrt( np.add.reduce( self.avg_pos**2 , 1 ) )
        for i in range( len(self.avg_distances) ) :
            self.AVG_file.write( str(self.sorted_ids[i]) + " " + str(self.avg_pos[i,0]) + " " + str(self.avg_pos[i,1]) + " " + str(self.avg_pos[i,2]) + " " + str(self.avg_distances[i]) + "\n" )
        self.AVG_file.close()

    def return_values( self ) :
        return { "avg_distances" : self.avg_distances ,
                 "avg_pos" : self.avg_pos ,
                 "sorted_ids" : self.sorted_ids ,
                 "Nconfs_analysed" : self.Nconfs_analysed }

    def merge_return_values( self , values_list ) :
        self.avg_distances = np.copy( values_list[0]["avg_distances"] )
        self.avg_pos = np.copy( values_list[0]["avg_pos"] )
        self.sorted_ids = np.copy( values_list[0]["sorted_ids"] )
        self.Nconfs_analysed = values_list[0]["Nconfs_analysed"]
        for i in range( 1 , len(values_list) ) :
            self.avg_distances = self.avg_distances + values_list[i]["avg_distances"]
            self.avg_pos = self.avg_pos + values_list[i]["avg_pos"]
            self.Nconfs_analysed += values_list[i]["Nconfs_analysed"]
        self.avg_distances = self.avg_distances / self.Nconfs_analysed
        self.avg_pos = self.avg_pos / self.Nconfs_analysed

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-com_mol_dist\033[0m } " )
