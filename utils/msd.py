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

EXIT_FAILURE = 1

def make_action( config , action , mode , val ) :
    if action == "select" :
        config.select( mode , val )
    if action == "remove" :
        config.remove( mode , val )
    if action == "recenter" :
        config.displace(  -config.com( mode , val )  )

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
            if sys.argv[i+j].strip() == "single" :
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

    list_of_confs = []
    for conf in conf_files :
        equil_config = cnf.CONFIGURATION()
        equil_config.time = conf[0]
        equil_config.description = conf[1]
        list_of_confs.append( equil_config )
    data_file = open( list_of_confs[0].description , "r" )
    list_of_confs[0].smart_auto_read( data_file , IN_FORMAT )
    data_file.close()
    box_side = np.copy( list_of_confs[0].box_sup - list_of_confs[0].box_inf )
    list_of_confs[0].discard_particles()
    N_confs = len( list_of_confs )

    # calculation of the cycle's length
    flag = 0
    cycle_nconfs = 0
    first_dt = list_of_confs[1].time - list_of_confs[0].time
    while flag == 0 and cycle_nconfs < int(floor(0.5*N_confs)) :
        cycle_nconfs += 1
#        cycle_length = list_of_confs[cycle_nconfs].time - list_of_confs[0].time
        flag = 1
        j = 2
        while j*cycle_nconfs < N_confs :
            if fabs( list_of_confs[j*cycle_nconfs].time - list_of_confs[(j-1)*cycle_nconfs].time - list_of_confs[cycle_nconfs].time + list_of_confs[0].time ) > 0.001*first_dt :
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
                if fabs( list_of_confs[i+j*cycle_nconfs].time - list_of_confs[i+(j-1)*cycle_nconfs].time - list_of_confs[i+cycle_nconfs].time + list_of_confs[i].time ) > 0.001*first_dt :
                    flag = 0
                j += 1
        if flag == 0 :
            print( " *** ERROR: the mean squared displacement cannot be computed :" )
            print( " ***        no logarithmic time-cycles can be individuated, neither linear ones !" )
            exit( EXIT_FAILURE )
    print( "Number of configurations:  " + str(N_confs) )
    print( "Number of configurations per cycle:  " + str( cycle_nconfs ) )
    print( "Cycle duration:  " + str( list_of_confs[ cycle_nconfs ].time - list_of_confs[ 0 ].time ) )

    ####################################################################################################
    ####################################### CHUNKS : NO ################################################
    ####################################################################################################

    if CHUNKS_YES == 0 :
        mean_sq_disp = []
        # computation of the "ballistic" part
        if cycle_nconfs > 1 :
            print( "Computation of the points in the ballistic region . . ." )
            for dt in range( 1 , cycle_nconfs ) :
                msd_dt = MSD()
                msd_dt.t = list_of_confs[ dt ].time - list_of_confs[ 0 ].time
                msd_dt.normalization = 0
                mean_sq_disp.append( msd_dt )
            t0 = 0
            while t0 < N_confs :
                # loading of the configurations t0
                data_file = open( list_of_confs[ t0 ].description , "r" )
                list_of_confs[ t0 ].smart_auto_read( data_file , IN_FORMAT )
                data_file.close()
                list_of_confs[ t0 ].unwrapped_coordinates( True )
                for ACTION in ACTIONS :
                    make_action( list_of_confs[ t0 ] , ACTION["operation"] , ACTION["mode"] , ACTION["value"] )
                for dt in range( 1 , cycle_nconfs ) :
                    if t0+dt < N_confs :
                        # loading of the configuration t0+dt
                        data_file = open( list_of_confs[ t0+dt ].description , "r" )
                        list_of_confs[ t0+dt ].smart_auto_read( data_file , IN_FORMAT )
                        data_file.close()
                        list_of_confs[ t0+dt ].unwrapped_coordinates( True )
                        for ACTION in ACTIONS :
                            make_action( list_of_confs[ t0+dt ] , ACTION["operation"] , ACTION["mode"] , ACTION["value"] )
                        if fabs( list_of_confs[ t0+dt ].time - list_of_confs[ t0 ].time - mean_sq_disp[dt-1].t ) > 0.001*first_dt :
                            print( "  ERROR:  MSD is being averaged for different time intervals!!" )
                        # changed
                        mean_sq_disp[dt-1].msd += np.add.reduce( ( list_of_confs[t0+dt].pos - list_of_confs[t0].pos )**2 , (0,1) )
                        mean_sq_disp[dt-1].normalization += list_of_confs[ t0 ].N
                        list_of_confs[ t0+dt ].discard_particles()
                list_of_confs[ t0 ].discard_particles()
                t0 += cycle_nconfs
            for dt in range( 1 , cycle_nconfs ) :
                if mean_sq_disp[dt-1].normalization > 0 :
                    mean_sq_disp[dt-1].msd /= float( mean_sq_disp[dt-1].normalization )
        else :
            print( " *** The configurations are linearly time-spaced ." )

        # computation of the "diffusive" part
        N_cycles = int( floor(N_confs/cycle_nconfs) )
        print( "Number of cycles:  " + str( N_cycles ) )
        print( "Computation of the points in the diffusive region . . ." )
        for dt in range( 1 , N_cycles ) :
            msd_dt = MSD()
            msd_dt.t = list_of_confs[ dt*cycle_nconfs ].time - list_of_confs[ 0 ].time
            msd_dt.normalization = 0
            mean_sq_disp.append( msd_dt )
        for l in range( N_cycles-1 ) :
            print( " Cycle number : " + str(l+1) )
            # loading of the reference configuration
            data_file = open( list_of_confs[ l*cycle_nconfs ].description , "r" )
            list_of_confs[ l*cycle_nconfs ].smart_auto_read( data_file , IN_FORMAT )
            data_file.close()
            list_of_confs[ l*cycle_nconfs ].unwrapped_coordinates( True )
            for ACTION in ACTIONS :
                make_action( list_of_confs[ l*cycle_nconfs ] , ACTION["operation"] , ACTION["mode"] , ACTION["value"] )
            for m in range( l+1 , N_cycles ) :
                data_file = open( list_of_confs[ m*cycle_nconfs ].description , "r" )
                list_of_confs[ m*cycle_nconfs ].smart_auto_read( data_file , IN_FORMAT )
                data_file.close()
                list_of_confs[ m*cycle_nconfs ].unwrapped_coordinates( True )
                for ACTION in ACTIONS :
                    make_action( list_of_confs[ m*cycle_nconfs ] , ACTION["operation"] , ACTION["mode"] , ACTION["value"] )
                if fabs( list_of_confs[ m*cycle_nconfs ].time - list_of_confs[ l*cycle_nconfs ].time - mean_sq_disp[ cycle_nconfs-1 + m-l -1 ].t ) > 0.001*first_dt :
                    print( "  ERROR:  MSD is being averaged for different time intervals (diffusive part) !!" )
                # changed
                mean_sq_disp[ cycle_nconfs-1 + m-l -1 ].msd += np.add.reduce( ( list_of_confs[ m*cycle_nconfs ].pos - list_of_confs[ l*cycle_nconfs ].pos )**2 , (0,1) )
                mean_sq_disp[ cycle_nconfs-1 + m-l -1 ].normalization += list_of_confs[ l*cycle_nconfs ].N
                list_of_confs[ m*cycle_nconfs ].discard_particles()
            list_of_confs[ l*cycle_nconfs ].discard_particles()
        for dt in range( 1 , N_cycles ) :
            if mean_sq_disp[ cycle_nconfs-1 + dt -1 ].normalization > 0 :
                mean_sq_disp[ cycle_nconfs-1 + dt -1 ].msd /= float( mean_sq_disp[ cycle_nconfs-2 + dt ].normalization )

        # creation of the output file
        print( "I'm saving data . . . ." )
        out_file = open( out_filename , "w" )
        out_file.write( "# time_step msd\n" )
        N_timesteps = len( mean_sq_disp )
        for i in range( N_timesteps ) :
            out_file.write( str( mean_sq_disp[i].t ) + ' ' + str( mean_sq_disp[i].msd ) + "\n" )
        out_file.close()

    ####################################################################################################
    ####################################### CHUNKS : YES ###############################################
    ####################################################################################################

    elif CHUNKS_YES == 1 :
        # definition of chunks
        if dr < 0 :
            print( "ERROR: You must instert a value for the radial thickness of chunks!" )
            exit( EXIT_FAILURE )
        if N_chunks == 0 :
            shorter_side = box_side[0] if ( box_side[0]<box_side[1] and box_side[0]<box_side[2] ) else ( box_side[1] if ( box_side[1]<box_side[2] ) else box_side[2] )
            N_chunks = int( floor( 0.5*shorter_side/dr ) ) + 1           # the last chunk is not a sphere-shell, but the remaining part of the box
                                                                         #     outside the greater shell

        # computation of the "ballistic" part
        msd_chunk_t = []
        if cycle_nconfs > 1 :
            for nchk in range( N_chunks ) :
                print( "Computation of the points in the ballistic region for the chunk " + str( nchk+1 ) + " ." )
                mean_sq_disp = []
                for dt in range( 1 , cycle_nconfs ) :
                    msd_dt = MSD()
                    msd_dt.t = list_of_confs[ dt ].time - list_of_confs[ 0 ].time
                    msd_dt.normalization = 0
                    mean_sq_disp.append( msd_dt )
                msd_chunk_t.append( mean_sq_disp )
            t0 = 0
            while t0 < N_confs :
                # loading of the configurations t0
                data_file = open( list_of_confs[ t0 ].description , "r" )
                list_of_confs[ t0 ].smart_auto_read( data_file , IN_FORMAT )
                data_file.close()
                list_of_confs[ t0 ].unwrapped_coordinates( True )
                for ACTION in ACTIONS :
                    make_action( list_of_confs[ t0 ] , ACTION["operation"] , ACTION["mode"] , ACTION["value"] )
                list_of_confs[ t0 ].create_chunks( "sphere" , dr , N_chunks , box_side )
                for dt in range( 1 , cycle_nconfs ) :
                    if t0+dt < N_confs :
                        # loading of the configuration t0+dt
                        data_file = open( list_of_confs[ t0+dt ].description , "r" )
                        list_of_confs[ t0+dt ].smart_auto_read( data_file , IN_FORMAT )
                        data_file.close()
                        list_of_confs[ t0+dt ].unwrapped_coordinates( True )
                        for ACTION in ACTIONS :
                            make_action( list_of_confs[ t0+dt ] , ACTION["operation"] , ACTION["mode"] , ACTION["value"] )
                        if fabs( list_of_confs[ t0+dt ].time - list_of_confs[ t0 ].time - msd_chunk_t[nchk][dt-1].t ) > 0.001*first_dt :
                            print( "  ERROR:  MSD is being averaged for different time intervals!!" )
                        ## efficientare con numpy, chunks array !
                        for nchk in range( N_chunks ) :
                            msd_chunk_t[nchk][dt-1].msd += np.add.reduce( ( list_of_confs[t0+dt].pos[ list_of_confs[ t0 ].chunk[nchk] ] - list_of_confs[t0].pos[ list_of_confs[ t0 ].chunk[nchk] ] )**2 , (0,1) )
                            msd_chunk_t[nchk][dt-1].normalization += len( list_of_confs[ t0 ].chunk[nchk] )
                        list_of_confs[ t0+dt ].discard_particles()
                list_of_confs[ t0 ].discard_particles()
                t0 += cycle_nconfs
            for dt in range( 1 , cycle_nconfs ) :
                for nchk in range( N_chunks ) :
                    if msd_chunk_t[nchk][dt-1].normalization > 0 :
                        msd_chunk_t[nchk][dt-1].msd /= float(msd_chunk_t[nchk][dt-1].normalization)
        else :
            print( " *** The configurations are linearly time-spaced ." )
            for nchk in range( N_chunks ) :
                mean_sq_disp = []
                msd_chunk_t.append( mean_sq_disp )

        # computation of the "diffusive" part
        N_cycles = int( floor(N_confs/cycle_nconfs) )
        print( "Number of cycles:  " + str( N_cycles ) )
        for nchk in range( N_chunks ) :
            print( "Computation of the points in the diffusive region for the chunk " + str( nchk+1 ) + " ." )
            for dt in range( 1 , N_cycles ) :
                msd_dt = MSD()
                msd_dt.t = list_of_confs[ dt*cycle_nconfs ].time - list_of_confs[ 0 ].time
                msd_dt.normalization = 0
                msd_chunk_t[nchk].append( msd_dt )
        for l in range( N_cycles-1 ) :
            # loading of the reference configuration
            data_file = open( list_of_confs[ l*cycle_nconfs ].description , "r" )
            list_of_confs[ l*cycle_nconfs ].smart_auto_read( data_file , IN_FORMAT )
            data_file.close()
            list_of_confs[ l*cycle_nconfs ].unwrapped_coordinates( True )
            for ACTION in ACTIONS :
                make_action( list_of_confs[ l*cycle_nconfs ] , ACTION["operation"] , ACTION["mode"] , ACTION["value"] )
            list_of_confs[ l*cycle_nconfs ].create_chunks( "sphere" , dr , N_chunks , box_side )
            for m in range( l+1 , N_cycles ) :
                data_file = open( list_of_confs[ m*cycle_nconfs ].description , "r" )
                list_of_confs[ m*cycle_nconfs ].smart_auto_read( data_file , IN_FORMAT )
                data_file.close()
                list_of_confs[ m*cycle_nconfs ].unwrapped_coordinates( True )
                for ACTION in ACTIONS :
                    make_action( list_of_confs[ m*cycle_nconfs ] , ACTION["operation"] , ACTION["mode"] , ACTION["value"] )
                if fabs( list_of_confs[ m*cycle_nconfs ].time - list_of_confs[ l*cycle_nconfs ].time - msd_chunk_t[nchk][ cycle_nconfs-1 + m-l -1 ].t ) > 0.001*first_dt :
                    print( "  ERROR:  MSD is being averaged for different time intervals (diffusive part) !!" )
                for nchk in range( N_chunks ) :
                    msd_chunk_t[ nchk ][ cycle_nconfs-1 + m-l -1 ].msd += np.add.reduce( ( list_of_confs[ m*cycle_nconfs ].pos[ list_of_confs[ l*cycle_nconfs ].chunk[nchk] ] - list_of_confs[ l*cycle_nconfs ].pos[ list_of_confs[ l*cycle_nconfs ].chunk[nchk] ] )**2 , (0,1) )
                    msd_chunk_t[ nchk ][ cycle_nconfs-1 + m-l -1 ].normalization += len( list_of_confs[ l*cycle_nconfs ].chunk[nchk] )
                list_of_confs[ m*cycle_nconfs ].discard_particles()
            list_of_confs[ l*cycle_nconfs ].discard_particles()
        for dt in range( 1 , N_cycles ) :
            for nchk in range( N_chunks ) :
                if msd_chunk_t[ nchk ][ cycle_nconfs-1 + dt -1 ].normalization > 0 :
                    msd_chunk_t[ nchk ][ cycle_nconfs-1 + dt -1 ].msd /= float( msd_chunk_t[ nchk ][ cycle_nconfs-2 + dt ].normalization )

        # creation of the output file
        print( "I'm saving data . . . ." )
        out_file = open( out_filename , "w" )
        out_file.write( "# BOX_EDGES " + str(box_side[0]) + ' ' + str(box_side[1]) + ' ' + str(box_side[2]) + "\n" )
        out_file.write( "# NUMBER_OF_CHUNKS " + str(N_chunks) + "\n" )
        out_file.write( "# THICKNESS_OF_CHUNKS " + str(dr) + "\n" )
        out_file.write( "# time_step msd_of_chunks . . .\n" )
        N_timesteps = len( msd_chunk_t[0] )
        for i in range( N_timesteps ) :
            out_file.write( str( msd_chunk_t[0][i].t ) )
            for ck in range( N_chunks ) :
                out_file.write( ' ' + str( msd_chunk_t[ck][i].msd ) )
            out_file.write( "\n" )
        out_file.close()


    if not SINGLE_CONF_FILES :
        os.system( "rm -r .msd_tmp" )







# definition of the class MSD
class MSD :

    def __init__( self ) :
        self.t = 0
        self.msd = 0.0
        self.normalization = 0









main()
