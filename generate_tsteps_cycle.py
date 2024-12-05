#!/usr/bin/env python3

import cmath
import random
from math import *
# interacting with the OS
import os
import sys
import argparse

LOOPO_PATH = os.path.dirname(os.path.abspath(__file__))
sys.path.append( os.path.expanduser( LOOPO_PATH + "/lib" ) )

EXIT_FAILURE = 1

parser = argparse.ArgumentParser(
    prog='create_msd_timesteps_cycles',
    description="This script outputs a series of INTEGER time steps approximately logarithmically spaced."
    "  It needs the logarithmic base (progression factor), the cycle length and the number of cycles.")

parser.add_argument( '-start', default=0, type=float,
                     help='Gives just the starting step of the progression.')

parser.add_argument( '-dt0', default=1.0, type=float,
                     help='Determines the smallest time interval at the beginning of the cycle.')

parser.add_argument( '-cyL', default=-1, type=float,
                     help='Number of timesteps for each logarithmic cycle.')

parser.add_argument( '-cyN', default=1, type=int,
                     help='Number of logarithmic cycles.')

parser.add_argument( '-base', default=-1, type=float,
                     help='logarithmic base (progression factor for subsequent timesteps).')

parser.add_argument( '-cyNsteps', default=-1, type=int,
                     help='Number of steps within each logarithmic cycle (alternative to -base option).')

parser.add_argument( '-int_steps', default=False, action='store_true',
                     help='To output a series of integer time steps.')

args = parser.parse_args()

def main() :
    if args.start < 0 :
        args.start = 0
        print( "*** WARNING: starting step should be a positive integer !" )
    if args.base <= 1 and args.cyNsteps <= 1 :
        print( "ERROR: Or base either cycle's number of steps have to be greater than 1 !!" )
    elif args.base > 1 :
        args.cyNsteps = int(round( log( args.cyL/args.dt0 ) / log( args.base ) ))
        args.cyL = args.dt0 * (args.base ** args.cyNsteps)
    else :
        args.base = ( args.cyL/args.dt0 ) ** (1.0/args.cyNsteps)
    if args.int_steps :
        args.start = int(args.start)
        args.cyL = int(args.cyL)
    print( "Cycle starting point : " + str(args.start) )
    print( "Logarithmic base : " + str(args.base) )
    print( "Cycle max length : " + str(args.cyL) )
    print( "Total number of cycles : " + str(args.cyN) )

    tsteps_cycle = []
    tsteps_cycle.append( 0 )
    if args.int_steps :
        for i in range(args.cyNsteps) :
            step = int(round( args.dt0 * args.base**i ))
            if step > tsteps_cycle[ len(tsteps_cycle)-1 ] :
                tsteps_cycle.append( step )
    else :
        for i in range(args.cyNsteps) :
            step = args.dt0 * args.base**i
            tsteps_cycle.append( step )

    # creation of the output file
    print( "I'm saving data . . . ." )
    out_file = open( "tsteps_log_prog.dat" , "w" )
    step = args.start
    for i in range( 0 , args.cyN ) :
        for j in range( len(tsteps_cycle) ) :
            out_file.write( str( step + tsteps_cycle[j] ) + "\n" )
        step += args.cyL
    out_file.write( str( step ) + "\n" )
    out_file.close()


main()
