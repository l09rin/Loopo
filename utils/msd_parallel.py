#!/usr/bin/env python3

from math import *
import cmath
import random
import numpy as np
# interacting with the OS
import os
import sys
import glob

import loopo.configuration as cnf

import multiprocessing
from multiprocessing import Manager
from multiprocessing.managers import BaseManager, BaseProxy

EXIT_FAILURE = 1

class SharedObjectsManager(BaseManager):
    pass

# definition of the class MSD
class MSD :
    def __init__( self ) :
        self.t = 0
        self.msd = 0.0
        self.normalization = 0
        self.lock = None

    def update( self , msd_val , msd_norm ) :
        with self.lock :
            self.msd += msd_val
            self.normalization += msd_norm

SharedObjectsManager.register('MSD', MSD)

def make_action( config , action , mode , val ) :
    if action == "select" :
        config.select( mode , val )
    if action == "remove" :
        config.remove( mode , val )
    if action == "recenter" :
        config.displace(  -config.com( mode , val )  )


def calculate_ballistic_point( task_queue , msd , t , normalization , locks , conf_files , IN_FORMAT , N_confs , cycle_nconfs , ACTIONS ) :
    while True:
        try:
            dt = task_queue.get_nowait()
        except multiprocessing.queues.Empty:
            break

        config_t0 = cnf.CONFIGURATION()
        config_t0dt = cnf.CONFIGURATION()
        t0 = 0
        while t0 < N_confs :
            # loading of the configurations t0
            data_file = open( conf_files[ t0 ][1] , "r" )
            config_t0.smart_auto_read( data_file , IN_FORMAT )
            data_file.close()
            config_t0.unwrapped_coordinates( True )
            for ACTION in ACTIONS :
                make_action( config_t0 , ACTION["operation"] , ACTION["mode"] , ACTION["value"] )
            if t0+dt < N_confs :
                # loading of the configuration t0+dt
                data_file = open( conf_files[ t0+dt ][1] , "r" )
                config_t0dt.smart_auto_read( data_file , IN_FORMAT )
                data_file.close()
                config_t0dt.unwrapped_coordinates( True )
                for ACTION in ACTIONS :
                    make_action( config_t0dt , ACTION["operation"] , ACTION["mode"] , ACTION["value"] )
                if fabs( config_t0dt.time - config_t0.time - t[dt-1] ) > 0.001*t[0] :
                    print( "  ERROR:  MSD is being averaged for different time intervals!!" )
                with locks[dt-1] :
                    msd[dt-1] += np.add.reduce( ( config_t0dt.pos - config_t0.pos )**2 , (0,1) )
                    normalization[dt-1] += config_t0.N
                config_t0dt.discard_particles()
            config_t0.discard_particles()
            t0 += cycle_nconfs
        print( " Point " + str(dt) + " / " + str(cycle_nconfs-1) + " done." )
        task_queue.task_done()


def calculate_diffusive_point( task_queue , msd , t , normalization , locks , conf_files , IN_FORMAT , N_confs , cycle_nconfs , N_cycles , ACTIONS ) :
    while True:
        try:
            dt = task_queue.get_nowait()
        except multiprocessing.queues.Empty:
            break

        config_t0 = cnf.CONFIGURATION()
        config_t0dt = cnf.CONFIGURATION()
        for l in range( N_cycles-dt ) :
            # loading of the reference configuration
            data_file = open( conf_files[ l*cycle_nconfs ][1] , "r" )
            config_t0.smart_auto_read( data_file , IN_FORMAT )
            data_file.close()
            config_t0.unwrapped_coordinates( True )
            for ACTION in ACTIONS :
                make_action( config_t0 , ACTION["operation"] , ACTION["mode"] , ACTION["value"] )
            data_file = open( conf_files[ (l+dt)*cycle_nconfs ][1] , "r" )
            config_t0dt.smart_auto_read( data_file , IN_FORMAT )
            data_file.close()
            config_t0dt.unwrapped_coordinates( True )
            for ACTION in ACTIONS :
                make_action( config_t0dt , ACTION["operation"] , ACTION["mode"] , ACTION["value"] )
            if fabs( config_t0dt.time - config_t0.time - t[ cycle_nconfs-1 + dt -1 ] ) > 0.001*t[0] :
                print( "  ERROR:  MSD is being averaged for different time intervals (diffusive part) !!" )
            with locks[ cycle_nconfs-1 + dt -1 ] :
                msd[ cycle_nconfs-1 + dt -1 ] += np.add.reduce( ( config_t0dt.pos - config_t0.pos )**2 , (0,1) )
                normalization[ cycle_nconfs-1 + dt -1 ] += config_t0.N
            config_t0dt.discard_particles()
            config_t0.discard_particles()

        print( " Point " + str(dt) + " / " + str(N_cycles-1) + " done." )
        task_queue.task_done()


def main() :
    print( "=============================================================================================================================================" )
    print( "This program computes the mean squared displacement of a system of particles." )
    print( "The input file can be in the LAMMPS format (lmp), sph format, or in the following one (xyz, default) :" )
    print( "     > first line : # N <number_of_particles>" )
    print( "     > second line : # step <time_step>" )
    print( "     > N lines containing the unwrapped coordinates of atoms." )
    print( "     > blank line(s)" )
    print( "     > repetition of the last three entries . . ." )
    print( "" )
    print( "   Usage :" )
    print( "    "+sys.argv[0].strip().split("/")[-1]+"  -confs { list <files_list.dat> | single \"<files_with_single_confs>\" | \"<file_name/s>\" }" )
    print("                    { -format <lmp|xyz|sph> -chunks <spherical_thickness> <number_of_chunks> -out <output_filename> }" )
    print("                    { -select <mol|type|q|xlinkers> <val|bonds.file> -remove <mol|type|q|xlinkers> <val|bonds.file> }" )
    print("                    { -recenter mol <val> }" )
    print("                    { -parallel <n_procs> }" )
    print( "" )
    print( "" )

    number_of_arguments = len(sys.argv)
    SINGLE_CONF_FILES = False
    IN_FORMAT="xyz"
    CHUNKS_YES=0
    dr = -1.0
    N_chunks = 0
    out_filename = ""
    ACTIONS = []
    PARALLEL = 1
    if number_of_arguments < 3 :
        print( "The number of arguments is too low" )
        exit( EXIT_FAILURE )
    for i in range( number_of_arguments ) :

        if sys.argv[i] == "-confs" :
            j = 1
            if sys.argv[i+j].strip() == "list" :
                SINGLE_CONF_FILES = True
                fileslist = open( sys.argv[i+j+1].strip() , "r" )
                lof = fileslist.readlines()
                fileslist.close()
                file_names = [ line.strip().split()[0] for line in lof ]
                j = 3
            elif sys.argv[i+j].strip() == "single" :
                SINGLE_CONF_FILES = True
                file_names = []
                for fnames in sys.argv[i+j+1].strip().split() :
                    file_names = file_names + glob.glob( fnames )
                j = 3
            else :
                file_names = []
                for fnames in sys.argv[i+j].strip().split() :
                    file_names = file_names + glob.glob( fnames )
                j = 2

        elif sys.argv[i] == "-out" :
            out_filename = sys.argv[i+1].strip()

        elif sys.argv[i] == "-parallel" :
            PARALLEL = int(sys.argv[i+1].strip())

        elif sys.argv[i] == "-format" :
            IN_FORMAT = sys.argv[i+1].strip()
            if IN_FORMAT not in [ "xyz" , "lmp" , "sph" ] :
                print( " *** Unrecognized format !" )
                exit( EXIT_FAILURE )

        elif sys.argv[i] == "-chunks" :
            CHUNKS_YES = 1
            dr = float( sys.argv[i+1].strip() )
            N_chunks = int( sys.argv[i+2].strip() )

        elif sys.argv[i] == "-select" :
            SELECT = True
            SEL_MODE = sys.argv[i+1].strip()
            if SEL_MODE in [ "mol" , "type" ] :
                SEL_VAL = int( sys.argv[i+2].strip() )
            elif SEL_MODE == "q" :
                SEL_VAL = float( sys.argv[i+2].strip() )
            elif SEL_MODE == "xlinkers" :
                SEL_VAL = sys.argv[i+2].strip()
            else :
                SELECT = False
            if SELECT :
                ACTIONS.append( { "operation" : "select" , "mode" : SEL_MODE , "value" : SEL_VAL } )

        elif sys.argv[i] == "-remove" :
            REMOVE = True
            RM_MODE = sys.argv[i+1].strip()
            if RM_MODE in [ "mol" , "type" ] :
                RM_VAL = int( sys.argv[i+2].strip() )
            elif RM_MODE == "q" :
                RM_VAL = float( sys.argv[i+2].strip() )
            elif RM_MODE == "xlinkers" :
                RM_VAL = sys.argv[i+2].strip()
            else :
                REMOVE = False
            if REMOVE :
                ACTIONS.append( { "operation" : "remove" , "mode" : RM_MODE , "value" : RM_VAL } )

        elif sys.argv[i] == "-recenter" :
            RECENTER = True
            RM_MODE = sys.argv[i+1].strip()
            if RM_MODE != "mol" :
                RECENTER = False
                print( "*** No molecule recentering will be applied ." )
            else :
                REC_MOL = int( sys.argv[i+2].strip() )
            if RECENTER :
                ACTIONS.append( { "operation" : "recenter" , "mode" : "mol" , "value" : REC_MOL } )

    if out_filename == "" :
        if CHUNKS_YES :
            out_filename = "msd_chunk.dat"
        else :
            out_filename = "msd.dat"

    # loading the configurations
    conf_files = []
    tmp_conf = cnf.CONFIGURATION()
    if not SINGLE_CONF_FILES :
        for fname in file_names :
            data_file = open( fname , "r" )
            time_steps , part_numbers = tmp_conf.split_confs( data_file , ".msd_tmp" , IN_FORMAT )
            data_file.close()
            for tstep in time_steps :
                print( "I found configuration at timestep " + str(tstep) + " ." )
                conf_files.append( [ tstep , ".msd_tmp/cnf-"+str(tstep)+".dat" ] )
    else :
        for fname in file_names :
            tstep = tmp_conf.fast_timestep_read( fname , IN_FORMAT )
            print( "I found configuration at timestep " + str(tstep) + " ." )
            conf_files.append( [ tstep , fname ] )
    conf_files = sorted( conf_files )

    # calculation of the cycle's length
    N_confs = len( conf_files )
    flag = 0
    cycle_nconfs = 0
    first_dt = conf_files[1][0] - conf_files[0][0]
    while flag == 0 and cycle_nconfs < int(floor(0.5*N_confs)) :
        cycle_nconfs += 1
        flag = 1
        j = 2
        while j*cycle_nconfs < N_confs :
            if fabs( conf_files[j*cycle_nconfs][0] - conf_files[(j-1)*cycle_nconfs][0] - conf_files[cycle_nconfs][0] + conf_files[0][0] ) > 0.001*first_dt :
                flag = 0
            j += 1
    if flag == 0 :
        print( " *** ERROR: the mean squared displacement cannot be computed :" )
        print( " ***        no logarithmic time-cycles can be individuated, neither linear ones !" )
        exit( EXIT_FAILURE )
    else :
        flag = 1
        for i in range( 1 , cycle_nconfs ) :
            j = 2
            while i+j*cycle_nconfs < N_confs :
                if fabs( conf_files[i+j*cycle_nconfs][0] - conf_files[i+(j-1)*cycle_nconfs][0] - conf_files[i+cycle_nconfs][0] + conf_files[i][0] ) > 0.001*first_dt :
                    flag = 0
                j += 1
        if flag == 0 :
            print( " *** ERROR: the mean squared displacement cannot be computed :" )
            print( " ***        no logarithmic time-cycles can be individuated, neither linear ones !" )
            exit( EXIT_FAILURE )
    print( "Number of configurations:  " + str(N_confs) )
    print( "Number of configurations per cycle:  " + str( cycle_nconfs ) )
    print( "Cycle duration:  " + str( conf_files[ cycle_nconfs ][0] - conf_files[ 0 ][0] ) )

    # Create a manager
    manager = Manager()
    # manager = SharedObjectsManager()
    # manager.start()
    # Create a managed list
    # mean_sq_disp = multiprocessing.Manager().list()
    msd = manager.list()
    normalization = manager.list()
    t = manager.list()
    locks = manager.list()

    ####################################################################################################
    ####################################### CHUNKS : NO ################################################
    ####################################################################################################

    if CHUNKS_YES == 0 :
        # computation of the "ballistic" part
        if cycle_nconfs > 1 :
            print( "Computation of the points in the ballistic region . . ." )
            for dt in range( 1 , cycle_nconfs ) :
                msd.append(0.0)
                t.append( conf_files[ dt ][0] - conf_files[ 0 ][0] )
                normalization.append( 0 )
                locks.append( manager.Lock() )
            # Create a task queue
            task_queue = multiprocessing.JoinableQueue()
            for dt in range( 1 , cycle_nconfs ) :
                task_queue.put( dt )
            # creation of the processes
            processes = []
            for _ in range(PARALLEL):
                p = multiprocessing.Process(target=calculate_ballistic_point, args=(task_queue , msd , t , normalization , locks , conf_files , IN_FORMAT , N_confs , cycle_nconfs , ACTIONS))
                p.start()
                processes.append(p)

            # Wait for all tasks to be processed
            task_queue.join()
            # Ensure all worker processes are terminated
            for p in processes:
                p.join()

            for dt in range( 1 , cycle_nconfs ) :
                if normalization[dt-1] > 0 :
                    msd[dt-1] /= float( normalization[dt-1] )
        else :
            print( " *** The configurations are linearly time-spaced ." )

        # computation of the "diffusive" part
        N_cycles = int( floor(N_confs/cycle_nconfs) )
        print( "Number of cycles:  " + str( N_cycles ) )
        print( "Computation of the points in the diffusive region . . ." )
        for dt in range( 1 , N_cycles ) :
            msd.append(0.0)
            t.append( conf_files[ dt*cycle_nconfs ][0] - conf_files[ 0 ][0] )
            normalization.append( 0 )
            locks.append( manager.Lock() )

        # Create a task queue
        task_queue = multiprocessing.JoinableQueue()
        for dt in range( 1 , N_cycles ) :
            task_queue.put( dt )
        # creation of the processes
        processes = []
        for _ in range(PARALLEL):
            p = multiprocessing.Process(target=calculate_diffusive_point, args=(task_queue , msd , t , normalization , locks , conf_files , IN_FORMAT , N_confs , cycle_nconfs , N_cycles , ACTIONS))
            p.start()
            processes.append(p)

        # Wait for all tasks to be processed
        task_queue.join()
        # Ensure all worker processes are terminated
        for p in processes:
            p.join()

        for dt in range( 1 , N_cycles ) :
            if normalization[ cycle_nconfs-1 + dt -1 ] > 0 :
                msd[ cycle_nconfs-1 + dt -1 ] /= float( normalization[ cycle_nconfs-2 + dt ] )

        # creation of the output file
        print( "I'm saving data . . . ." )
        out_file = open( out_filename , "w" )
        out_file.write( "# time_step msd\n" )
        N_timesteps = len( msd )
        for i in range( N_timesteps ) :
            out_file.write( str( t[i] ) + ' ' + str( msd[i] ) + "\n" )
        out_file.close()

    ####################################################################################################
    ####################################### CHUNKS : YES ###############################################
    ####################################################################################################

    elif CHUNKS_YES == 1 :
        print( "*** ERROR: Chunk function not yet coded !" )


    if not SINGLE_CONF_FILES :
        os.system( "rm -r .msd_tmp" )








if __name__ == "__main__":
    main()
