#!/usr/bin/env python3

# OS INTERACTION TOOLS
import os
import sys
import glob  # It allows to expand matching pattern like in bash shell
import time

LOOPO_PATH = os.path.dirname(os.path.abspath(__file__))
sys.path.append( os.path.expanduser( LOOPO_PATH + "/lib" ) )
if 'DISPLAY' in os.environ.keys() : del os.environ['DISPLAY']
if 'XAUTHORITY' in os.environ.keys() : del os.environ['XAUTHORITY']

# MULTI THREADING
import multiprocessing as mp
from multiprocessing.managers import BaseManager

# CUSTOM LIBRARY
import configuration as cnf

# ANALYSIS TOOLS
sys.path.append( os.path.expanduser( LOOPO_PATH + "/lib/ANALYSIS" ) )
from dump import DUMP
from density_profiles import PROFILE
from com_mol_distance import COM_MOL_DISTANCE
from recenter import RECENTER
from clusters import CLUSTERS_CHECK
from add_particles import ADD_PARTICLES
from displace import DISPLACE
from rescale import RESCALE
from change_attributes import CHANGE , CHANGE_BOX
from selection import SELECT , REMOVE
from shape_overlap import SHAPE_OVERLAP
from convex_hull import CONVEX_HULL
from surface_mesh import SURFACE_MESH
from com import COM , VCOM
from gyration import GYRATION
from hexatic import HEXATIC
from effective_charge import *
from reset import *
from wrap import *
from timeshift import *
EXIT_FAILURE = 1


# custom manager to support custom classes
class CustomManager(BaseManager):
    # nothing
    pass

def configurations_analysis( MAINlock , datafiles , ACTIONS , RETURN_ACTIONS ) :
    col = COLUMNS()
    equil_config = cnf.CONFIGURATION()
    subprocessNAME = mp.current_process().name

    MAINlock.acquire()
    print( " >> Thread " + str(subprocessNAME) + " started ." )
    MAINlock.release()
    while datafiles.isnotempty() :
        MAINlock.acquire()
        datafiles.update()

        if datafiles.isnotempty() :
            t0 = time.time()
            if datafiles.get_format() in [ "xyz" , "sph" , "patch" , "ptc" ] :
                col.id , col.type , col.x , col.y , col.z , col.q , col.mol , col.vx , col.vy , col.vz = equil_config.smart_auto_read( datafiles , datafiles.get_format() , "NO" )
            else :
                col.id , col.type , col.x , col.y , col.z , col.q , col.mol , col.vx , col.vy , col.vz = equil_config.smart_auto_read( datafiles , datafiles.get_format() , "YES" )
            if not datafiles.file_isclosed() and equil_config.N > 0 :
                datafiles.increment()
                print( " >> Thread " + str(subprocessNAME) + " read configuration number " + str(datafiles.get_read_confs()) + ", at timestep " + str(equil_config.time) + " in " + str(time.time()-t0) + "s ." )

        MAINlock.release()
        if not datafiles.file_isclosed() and equil_config.N > 0 :
            # Handling of the configuration
            t0 = time.time()
            for ACTION in ACTIONS :
                ACTION.execute( equil_config , col )
            equil_config.clear_configuration()
            tf = time.time()
            MAINlock.acquire()
            print( " >>>> Thread " + str(subprocessNAME) + ", computation time : " + str(tf-t0) + "s ." , flush=True )
            MAINlock.release()
    MAINlock.acquire()
    print( " >> Thread " + str(subprocessNAME) + " finished ." )
    MAINlock.release()
    RETURN_VALUES = []
    for i in range(len( ACTIONS )) :
        RETURN_VALUES.append( ACTIONS[i].return_values() )
    RETURN_ACTIONS.append( RETURN_VALUES )





if __name__ == '__main__':
    print( "=============================================================================================================================================" )
    print( "This program select particles over configurations dumped by LAMMPS." )
    print( "    The option -recenter applies only to the selected atoms." )
    print( "    The option -attributes allows to choose the atom features to keep in the output file." )
    print( "" )
    print( "    "+sys.argv[0].strip().split("/")[-1]+" \033[1m-confs\033[0m <lammps_dumped_configurations:singlefile|\"pattern_string\">|<list <file_with_files-list>> { <lmp|xyz|sph|patch|ptc> }" )
    print( "\t\t\t{ \033[1m-parallel\033[0m <number_of_threads:1> } " )
    print( "\t\t\t{ \033[1m--help\033[0m | \033[1m-h\033[0m | \033[1m-help\033[0m } " )
    print( "\t\t\t{ <action> <options>}" )
    print( "" )
    print( "" )
    print( "                        [...] = default values" )
    print( "" )
    print( "" )

    number_of_arguments = len(sys.argv)
    Nthreads = 1

    CustomManager.register('DATA_FILES', DATA_FILES)
    manager = CustomManager()
    manager.start()
    datafiles = manager.DATA_FILES()
    ACTIONS = []
    RETURN_ACTIONS = mp.Manager().list()

    for i in range( number_of_arguments ) :

        if sys.argv[i] == "-confs" :
            j = 1
            if sys.argv[i+j].strip() == "list" :
                fileslist = open( sys.argv[i+j+1].strip() , "r" )
                lof = fileslist.readlines()
                fileslist.close()
                datafiles.set_files( [ line.strip().split()[0] for line in lof ] )
                j = 3
            else :
                single_files = []
                for fnames in sys.argv[i+j].strip().split() :
                    single_files = single_files + glob.glob( fnames )
                datafiles.set_files( single_files )
                j = 2
            if sys.argv[i+j].strip() == "lmp" :
                datafiles.set_format( "lmp" )
            elif sys.argv[i+j].strip() == "xyz" :
                datafiles.set_format( "xyz" )
            elif sys.argv[i+j].strip() == "sph" :
                datafiles.set_format( "sph" )
            elif sys.argv[i+j].strip() == "patch" :
                datafiles.set_format( "patch" )
            elif sys.argv[i+j].strip() == "ptc" :
                datafiles.set_format( "ptc" )
            else :
                datafiles.set_format( "lmp" )

        elif sys.argv[i] == "-parallel" :
            Nthreads = int( sys.argv[i+1] )
            print("Maximum number of processors: ", mp.cpu_count())
            print("Number of threads: ", Nthreads)
            if Nthreads < 1 :
                print( "The number of threads must be a positive integer greater than 1 !!!" )
                exit(EXIT_FAILURE)

        elif sys.argv[i] == "-print" :
            ACTIONS.append( DUMP( sys.argv , i ) )

        elif sys.argv[i] == "-recenter" :
            ACTIONS.append( RECENTER( sys.argv , i ) )

        elif sys.argv[i] == "-shift_time" :
            ACTIONS.append( SHIFT_TIME( sys.argv , i ) )

        elif sys.argv[i] == "-rewrap" :
            ACTIONS.append( REWRAP( sys.argv , i ) )

        elif sys.argv[i] == "-unwrap" :
            ACTIONS.append( UNWRAP( sys.argv , i ) )

        elif sys.argv[i] == "-reset" :
            ACTIONS.append( RESET( sys.argv , i ) )

        elif sys.argv[i] == "-Qeff" :
            ACTIONS.append( Q_EFF( sys.argv , i ) )

        elif sys.argv[i] == "-gyration" :
            ACTIONS.append( GYRATION( sys.argv , i ) )

        elif sys.argv[i] == "-profile" :
            ACTIONS.append( PROFILE( sys.argv , i ) )

        elif sys.argv[i] == "-com_mol_dist" :
            ACTIONS.append( COM_MOL_DISTANCE( sys.argv , i ) )

        elif sys.argv[i] == "-com" :
            ACTIONS.append( COM( sys.argv , i ) )

        elif sys.argv[i] == "-vcom" :
            ACTIONS.append( VCOM( sys.argv , i ) )

        elif sys.argv[i] == "-convex_hull" :
            ACTIONS.append( CONVEX_HULL( sys.argv , i ) )

        elif sys.argv[i] == "-surface_mesh" :
            ACTIONS.append( SURFACE_MESH( sys.argv , i ) )

        elif sys.argv[i] == "-shape_overlap" :
            ACTIONS.append( SHAPE_OVERLAP( sys.argv , i ) )

        elif sys.argv[i] == "-sel" :
            ACTIONS.append( SELECT( sys.argv , i ) )

        elif sys.argv[i] == "-rm" :
            ACTIONS.append( REMOVE( sys.argv , i ) )

        elif sys.argv[i] == "-change_box" :
            ACTIONS.append( CHANGE_BOX( sys.argv , i ) )

        elif sys.argv[i] == "-change" :
            ACTIONS.append( CHANGE( sys.argv , i ) )

        elif sys.argv[i] == "-displace" :
            ACTIONS.append( DISPLACE( sys.argv , i ) )

        elif sys.argv[i] == "-rescale" :
            ACTIONS.append( RESCALE( sys.argv , i ) )

        elif sys.argv[i] == "-add_particles" :
            ACTIONS.append( ADD_PARTICLES( sys.argv , i ) )

        elif sys.argv[i] == "-clusters" :
            ACTIONS.append( CLUSTERS_CHECK( sys.argv , i ) )

        elif sys.argv[i] == "-hexatic" :
            ACTIONS.append( HEXATIC( sys.argv , i ) )

        elif sys.argv[i] in [ "--help" , "-h" , "-help" ] :
            # here it prints help messages
            print( "" )
            DUMP( [] , -1 )
            print( "" )
            SELECT( [] , -1 )
            print( "" )
            REMOVE( [] , -1 )
            print( "" )
            ADD_PARTICLES( [] , -1 )
            print( "" )
            REWRAP( [] , -1 )
            UNWRAP( [] , -1 )
            print( "" )
            CHANGE_BOX( [] , -1 )
            CHANGE( [] , -1 )
            print( "" )
            CLUSTERS_CHECK( [] , -1 )
            print( "" )
            COM( [] , -1 )
            VCOM( [] , -1 )
            print( "" )
            CONVEX_HULL( [] , -1 )
            print( "" )
            SURFACE_MESH( [] , -1 )
            print( "" )
            PROFILE( [] , -1 )
            print( "" )
            COM_MOL_DISTANCE( [] , -1 )
            print( "" )
            RECENTER( [] , -1 )
            print( "" )
            GYRATION( [] , -1 )
            print( "" )
            DISPLACE( [] , -1 )
            print( "" )
            RESCALE( [] , -1 )
            print( "" )
            SHIFT_TIME( [] , -1 )
            print( "" )
            Q_EFF( [] , -1 )
            print( "" )
            RESET( [] , -1 )
            print( "" )
            SHAPE_OVERLAP( [] , -1 )
            print( "" )
            HEXATIC( [] , -1 )
            print( "" )
            exit()

    if number_of_arguments < 3 :
        print( "The number of arguments is too low" )
        exit(EXIT_FAILURE)

    MAINlock = mp.Lock()
    for i in range( len(ACTIONS) ) :
        ACTIONS[i].MAINlock = MAINlock
    threads = []
    for i in range( Nthreads ) :
        threads.append(  mp.Process( name=str(i) , target=configurations_analysis , args=(MAINlock , datafiles , ACTIONS , RETURN_ACTIONS) )  )
        threads[i].start()

    for i in range( Nthreads ) :
        threads[i].join()

    for i in range(len( ACTIONS )) :
        ACTIONS[i].merge_return_values( [ RETURN_ACTIONS[j][i] for j in range(len(RETURN_ACTIONS)) ] )
        ACTIONS[i].terminate()


    print( "Done .\n" )
