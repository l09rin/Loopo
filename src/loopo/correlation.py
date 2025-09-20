# Library defining correlation points and a function that computes the time correlation function over a series of configurations

from math import *
import cmath
import random
import numpy as np



# definition of the class CORRELATION
class CORRELATION :

    def __init__( self ) :
        self.t = 0.0
        self.correlation = 0.0
        self.correlation_v3D = np.array( [ 0.0 , 0.0 , 0.0 ] )
        self.normalization = 0

    def normalize( self ) :
        if self.normalization > 0 :
            self.correlation /= float( self.normalization )
            self.correlation_v3D /= float( self.normalization )
            self.normalization = 0
        else :
            self.correlation = 0.0
            self.correlation_v3D[:] = 0.0


def build_correlation_function( N_cycles , cycle_nconfs , timesteps ) :
    correlation_function = []
    for dt in range( cycle_nconfs ) :
        correlation_dt = CORRELATION()
        correlation_dt.t = timesteps[ dt ] - timesteps[ 0 ]
        correlation_dt.normalization = 0
        correlation_function.append( correlation_dt )
    for dt in range( 1 , N_cycles+1 ) :
        correlation_dt = CORRELATION()
        correlation_dt.t = timesteps[ dt*cycle_nconfs ] - timesteps[ 0 ]
        correlation_dt.normalization = 0
        correlation_function.append( correlation_dt )
    return correlation_function



def calculate_cycles( timesteps ) :
    N_confs = len( timesteps )
    # calculation of the cycle's length
    flag = 0
    starting_step = timesteps[0]
    first_step = timesteps[1] - timesteps[0]
    previous_step = starting_step
    cycle_nconfs = 0
    for i in range( 1 , len(timesteps) ) :
        if flag == 0 :
            if timesteps[i] - previous_step > first_step :
                flag = 1
        elif flag == 1 :
            if timesteps[i] - previous_step == first_step :
                cycle_nconfs = i-1
                flag = -1
        previous_step = timesteps[i]
    if flag == 0 :
        for i in range( 1 , N_confs ) :
            if timesteps[i] - timesteps[i-1] != timesteps[1] - timesteps[0] :
                flag = 10
    if flag == 10 :
        print( " *** ERROR: the averaged velocity correlation function cannot be computed :" )
        print( " ***        no logarithmic time-cycles can be individuated, neither linear ones !" )
        exit( EXIT_FAILURE )
    if cycle_nconfs == 0 :
        cycle_nconfs = 1
    print( "Number of configurations:  " + str( len(timesteps) ) )
    print( "Number of configurations per cycle:  " + str( cycle_nconfs ) )
    print( "Cycle duration:  " + str( timesteps[ cycle_nconfs ] - timesteps[ 0 ] ) )
    N_cycles = 0
    while (N_cycles+1)*cycle_nconfs+1 <= N_confs :
        N_cycles += 1
    print( "Number of cycles:  " + str( N_cycles ) )
    if cycle_nconfs == 1 :
        print( " *** The configurations are linearly time-spaced ." )
    return cycle_nconfs , N_cycles



def correlate_cycles( conf_files , correlator , CONFIG_CLASS , FORMAT_PARAMETERS = [] , ARGS = [] ) :
    conf_files = sorted( conf_files )
    N_confs = len( conf_files )
    timesteps = [ conf_files[i][0] for i in range(len( conf_files )) ]
    cycle_nconfs , N_cycles = calculate_cycles( timesteps )
    correlation_function = build_correlation_function( N_cycles , cycle_nconfs , timesteps )

    bra_config = CONFIG_CLASS()
    ket_config = CONFIG_CLASS()
    delta0 = conf_files[1][0] - conf_files[0][0]
    # computation of the "ballistic" part
    print( "Computation of the points in the ballistic region . . ." )
    t0 = 0
    while t0 < N_confs :
        # loading of configuration t0
        bra_file = open( conf_files[ t0 ][1] , "r" )
        bra_config.read( bra_file , FORMAT_PARAMETERS )
        bra_file.close()
        for dt in range( cycle_nconfs ) :
            if t0+dt < N_confs :
                # loading of the configuration t0+dt
                ket_file = open( conf_files[ t0 + dt ][1] , "r" )
                ket_config.read( ket_file , FORMAT_PARAMETERS )
                ket_file.close()
                if fabs( ket_config.time - bra_config.time - correlation_function[dt].t ) > 10**(-3) * delta0 :
                    print( "  ERROR:  CORRELATION is being averaged for different time intervals!!" )

                correlation_function[dt].correlation += correlator( bra_config , ket_config , ARGS )
                correlation_function[dt].normalization += 1

                ket_config.discard()
        bra_config.discard()
        t0 += cycle_nconfs
    for dt in range( cycle_nconfs ) :
        correlation_function[dt].normalize()

    # computation of the "diffusive" part
    print( "Computation of the points in the diffusive region . . ." )
    for l in range( N_cycles ) :
        print( " Cycle number : " + str(l+1) )
        # loading of the reference configuration
        bra_file = open( conf_files[ l*cycle_nconfs ][1] , "r" )
        bra_config.read( bra_file , FORMAT_PARAMETERS )
        bra_file.close()
        for m in range( l+1 , N_cycles+1 ) :
            ket_file = open( conf_files[ m*cycle_nconfs ][1] , "r" )
            ket_config.read( ket_file , FORMAT_PARAMETERS )
            ket_file.close()
            if fabs( ket_config.time - bra_config.time - correlation_function[ cycle_nconfs + m-l -1 ].t ) > 10**(-3) * delta0 :
                print( "  ERROR:  CORRELATION is being averaged for different time intervals (diffusive part) !!" )

            correlation_function[ cycle_nconfs + m-l -1 ].correlation += correlator( bra_config , ket_config , ARGS )
            correlation_function[ cycle_nconfs + m-l -1 ].normalization += 1

            ket_config.discard()
        bra_config.discard()
    for dt in range( 1 , N_cycles+1 ) :
        correlation_function[ cycle_nconfs + dt -1 ].normalize()
    return correlation_function
