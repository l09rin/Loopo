from .action import *
import numpy as np
from loopo import profiles as pf

class PROFILE( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "profile"
            self.PROFILE_TYPE = argv[i+1].strip()
            self.PROFILES = []
            self.last_step = -1
            self.first_step = -1
            self.last_box = np.array( [0.0,0.0,0.0] )
            self.CYLINDER_AXIS = "z"
            self.VARIANCE = 0
            self.PROFILE_MODE = "all"
            self.PROFILE_VALUE = -1
            self.COM_MODE = "all"
            self.COM_VALUE = -1
            self.FIXED_DIAMETER = 0.0
            if self.PROFILE_TYPE == "sphere" :
                self.RADIAL_BIN = float( argv[i+2].strip() )
                j=3
            elif self.PROFILE_TYPE == "cylinder" :
                if argv[i+2].strip() not in [ "x" , "y" , "z" ] :
                    print( "*** For cylinder profiles axis has to be x, y or z !" )
                    exit( EXIT_FAILURE )
                else :
                    self.CYLINDER_AXIS = argv[i+2].strip()
                self.HEIGHT_BIN = float( argv[i+3].strip() )
                self.RADIAL_BIN = float( argv[i+4].strip() )
                j=5
            elif self.PROFILE_TYPE == "linear" :
                if argv[i+2].strip() not in [ "x" , "y" , "z" ] :
                    print( "*** For cylinder profiles axis has to be x, y or z !" )
                    exit( EXIT_FAILURE )
                else :
                    self.CYLINDER_AXIS = argv[i+2].strip()
                self.HEIGHT_BIN = float( argv[i+3].strip() )
                j=4
            else :
                print( "*** ERROR: The profile type has not been recognized !" )
                exit()
            fname_root = "atoms_profile.dat"
            if i+j < len(argv) :
                while argv[i+j].strip() in [ "specie" , "center" , "file" , "fixed_diameter" , "variance" ] :
                    if argv[i+j].strip() == "specie" :
                        j += 1
                        self.PROFILE_MODE = argv[i+j].strip()
                        j += 1
                        self.PROFILE_VALUE = int( argv[i+j].strip() )
                    elif argv[i+j].strip() == "center" :
                        j += 1
                        self.COM_MODE = argv[i+j].strip()
                        j += 1
                        self.COM_VALUE = argv[i+j].strip()
                    elif argv[i+j].strip() == "file" :
                        j += 1
                        fname_root = argv[i+j].strip()
                    elif argv[i+j].strip() == "fixed_diameter" :
                        j += 1
                        self.FIXED_DIAMETER = float( argv[i+j].strip() )
                    elif argv[i+j].strip() == "variance" :
                        j += 1
                        if argv[i+j].strip() == "yes" :
                            self.VARIANCE = 1
                            self.VARIANCE_file = open( "variance_" + fname_root , "w" )
                    j += 1
                    if i+j == len(argv) :
                        j -= 1
            self.AVG_file = open( fname_root , "w" )
            self.PERCONF_file = open( "singleconfs_" + fname_root , "w" )
            if self.PROFILE_TYPE == "sphere" :
                self.PERCONF_file.write( "# Single-configuration density profiles\n" )
                self.PERCONF_file.write( "# Timestep Number-of-chunks Total-count\n" )
                self.PERCONF_file.write( "# Chunk midpoint Ncount density/number\n" )
            elif self.PROFILE_TYPE == "cylinder" :
                self.PERCONF_file.write( "# Single-configuration density profiles\n" )
                self.PERCONF_file.write( "# Timestep Number-of-chunks Total-count\n" )
                self.PERCONF_file.write( "# Chunk axial-midpoint radial-midpoint Ncount density/number\n" )
            elif self.PROFILE_TYPE == "linear" :
                self.PERCONF_file.write( "# Single-configuration density profiles\n" )
                self.PERCONF_file.write( "# Timestep Number-of-chunks Total-count\n" )
                self.PERCONF_file.write( "# Chunk midpoint Ncount density/number\n" )
            self.PERCONF_file.flush()

    def execute( self , equil_config , col ) :
        if ( self.PROFILE_MODE not in [ "mol" , "type" , "all" ] ) or ( self.COM_MODE not in [ "mol" , "type" , "all" , "fixed" ] ) :
            self.lock.acquire()
            print( " *** Profile modes not recognized !" )
            self.lock.release()
            exit()
        if ( self.PROFILE_MODE == "mol" or self.COM_MODE == "mol" ) and col.mol == -1 :
            self.lock.acquire()
            print( " *** ERROR: I do not have information about molecules ID !" )
            self.lock.release()
            exit()
        if ( self.PROFILE_MODE == "type" or self.COM_MODE == "type" ) and col.type == -1 :
            self.lock.acquire()
            print( " *** ERROR: I do not have information about type !" )
            self.lock.release()
            exit()
        if self.COM_MODE in [ "all" , "mol" , "type" ] :
            center_of_mass = equil_config.com( self.COM_MODE , int( self.COM_VALUE ) )
        elif self.COM_MODE == "fixed" :
            point = self.COM_VALUE.strip("()").split(',')
            center_of_mass = np.array( [ float(point[0]) , float(point[1]) , float(point[2]) ] )
        last_box = np.array( [0.0,0.0,0.0] )
        if self.FIXED_DIAMETER == 0.0 :
            last_box = equil_config.box_sup - equil_config.box_inf
        else :
            last_box = np.array( [ self.FIXED_DIAMETER , self.FIXED_DIAMETER , self.FIXED_DIAMETER ] )
        if self.PROFILE_TYPE == "sphere" :
            actual_profile = equil_config.radial_profile( center_of_mass , self.RADIAL_BIN , last_box , self.PROFILE_MODE , self.PROFILE_VALUE )
        elif self.PROFILE_TYPE == "cylinder" :
             actual_profile = equil_config.cylinder_profile( center_of_mass , last_box , self.CYLINDER_AXIS , self.RADIAL_BIN , self.HEIGHT_BIN , self.PROFILE_MODE , self.PROFILE_VALUE )
        elif self.PROFILE_TYPE == "linear" :
             actual_profile = equil_config.linear_profile( last_box , self.CYLINDER_AXIS , self.HEIGHT_BIN , self.PROFILE_MODE , self.PROFILE_VALUE )
        self.lock.acquire()
        self.last_step = max( equil_config.time , self.last_step )
        if self.first_step >= 0 :
            self.first_step = min( equil_config.time , self.first_step )
        else :
            self.first_step = equil_config.time
        self.last_box = np.copy( last_box )
        if self.PROFILE_TYPE == "sphere" :
            self.PROFILES.append( actual_profile )
            self.PERCONF_file.write( str( equil_config.time ) + '\t' + str( self.PROFILES[-1].Nbins ) + '\t' + str( self.PROFILES[-1].total_count() ) + '\n' )
            for i in range( self.PROFILES[-1].Nbins ) :
                self.PERCONF_file.write( str( i+1 ) + '\t' + repr( 0.5*(self.PROFILES[-1].rmin[i]+self.PROFILES[-1].rmax[i]) ) + '\t' + repr( self.PROFILES[-1].counts[i] ) + '\t' + str( self.PROFILES[-1].dens[i] ) + '\n' )
        elif self.PROFILE_TYPE == "cylinder" :
             self.PROFILES.append( actual_profile )
             self.PERCONF_file.write( str( equil_config.time ) + '\t' + str( self.PROFILES[-1].Nrbins*(2*self.PROFILES[-1].lastHbin+1) ) + '\t' + str( self.PROFILES[-1].total_count() ) + '\n' )
             for j in range( self.PROFILES[-1].Nrbins ) :
                 for i in range( -self.PROFILES[-1].lastHbin , self.PROFILES[-1].lastHbin+1 ) :
                    self.PERCONF_file.write( str( (i+self.PROFILES[-1].lastHbin)+1+j*(2*self.PROFILES[-1].lastHbin+1) ) + '\t' + repr( 0.5*(self.PROFILES[-1].hmin[i]+self.PROFILES[-1].hmax[i]) ) + '\t' + repr( 0.5*(self.PROFILES[-1].rmin[j]+self.PROFILES[-1].rmax[j]) ) + '\t' + repr( self.PROFILES[-1].counts[i][j] ) + '\t' + str( self.PROFILES[-1].dens[i][j] ) + '\n' )
        elif self.PROFILE_TYPE == "linear" :
            self.PROFILES.append( actual_profile )
            self.PERCONF_file.write( str( equil_config.time ) + '\t' + str( self.PROFILES[-1].Nbins ) + '\t' + str( self.PROFILES[-1].total_count() ) + '\n' )
            for i in range( self.PROFILES[-1].Nbins ) :
                self.PERCONF_file.write( str( i+1 ) + '\t' + repr( 0.5*(self.PROFILES[-1].hmin[i]+self.PROFILES[-1].hmax[i]) ) + '\t' + repr( self.PROFILES[-1].counts[i] ) + '\t' + str( self.PROFILES[-1].dens[i] ) + '\n' )
        self.PERCONF_file.flush()
        self.lock.release()

    def terminate( self ) :
        if self.PROFILE_TYPE == "sphere" :
            self.PERCONF_file.close()
            self.AVG_file.write( "# Time-averaged density profile\n" )
            self.AVG_file.write( "# Timesteps Number-of-chunks Total-count\n" )
            self.AVG_file.write( "# Chunk midpoint Ncount density/number\n" )
            FLAG_DIFF_BOXES = 0
            rmax0 = self.PROFILES[0].rmax[-1]
            for profile in self.PROFILES :
                if rmax0 != profile.rmax[-1] :
                    FLAG_DIFF_BOXES = 1
            if FLAG_DIFF_BOXES :
                # if you want to take into account a series of configurations with variable box sides !!
                min_boxside = np.amin( self.last_box )
                max_boxside = np.multiply.reduce( self.last_box , 0 )**(1./3)
                for profile in self.PROFILES :
                    max_boxside = max( max_boxside , profile.rmax[-1]*(4./3*np.pi)**(1./3) )
                    min_boxside = min( min_boxside , profile.rmin[-1]*(4./3*np.pi)**(1./3) )
                avg_prof = pf.density_profile( self.RADIAL_BIN , np.array( [ min_boxside , min_boxside , min_boxside ] ) )
                avg_prof.rmax[-1] = max_boxside/(4./3*np.pi)**(1./3)
                midpoint = []
                for i in range( avg_prof.Nbins ) :
                    midpoint.append( 0.5 * ( avg_prof.rmin[i] + avg_prof.rmax[i] ) )
                N = len( self.PROFILES )
                for profile in self.PROFILES :
                    for i in range( avg_prof.Nbins-1 ) :
                        avg_prof.counts[i] += profile.counts[i]
                        avg_prof.dens[i] += profile.dens[i] / N
                    for i in range( avg_prof.Nbins-1 , profile.Nbins ) :
                        avg_prof.counts[-1] += profile.counts[i]
                        avg_prof.dens[-1] += profile.dens[i] / N / (profile.Nbins-avg_prof.Nbins+1)
            else :
                avg_prof = pf.density_profile( self.RADIAL_BIN , self.last_box )
                midpoint = []
                for i in range( avg_prof.Nbins ) :
                    midpoint.append( 0.5 * ( avg_prof.rmin[i] + avg_prof.rmax[i] ) )
                N = len( self.PROFILES )
                for profile in self.PROFILES :
                    for i in range( profile.Nbins ) :
                        avg_prof.counts[i] += profile.counts[i]
                        avg_prof.dens[i] += profile.dens[i] / N
            self.AVG_file.write( str( self.last_step-self.first_step ) + '\t' + str( avg_prof.Nbins ) + '\t' + str( avg_prof.total_count() ) + '\n' )
            for i in range( avg_prof.Nbins ) :
                self.AVG_file.write( str( i+1 ) + '\t' + repr( midpoint[i] ) + '\t' + repr( avg_prof.counts[i] ) + '\t' + repr( avg_prof.dens[i] ) + '\n' )
            self.AVG_file.close()
            if self.VARIANCE :
                profile_variance = pf.density_profile( self.RADIAL_BIN , self.last_box )
                self.VARIANCE_file.write( "# Chunk midpoint density/number dev_stand\n" )
                for profile in self.PROFILES :
                    for i in range( profile.Nbins ) :
                        profile_variance.dens[i] += ( profile.dens[i] - avg_prof.dens[i] )**2 / N
                for i in range( avg_prof.Nbins ) :
                    self.VARIANCE_file.write( str( i+1 ) + '\t' + repr( midpoint[i] ) + '\t' + repr( avg_prof.dens[i] ) + '\t' + repr( np.sqrt(profile_variance.dens[i]) ) + '\n' )
                self.VARIANCE_file.close()

        elif self.PROFILE_TYPE == "cylinder" :
            self.PERCONF_file.close()
            self.AVG_file.write( "# Time-averaged density profile\n" )
            self.AVG_file.write( "# Timesteps Number-of-chunks Total-count\n" )
            self.AVG_file.write( "# Chunk h_midpoint r_midpoint Ncount density/number\n" )
            avg_prof = pf.cylindrical_profile( self.last_box , self.CYLINDER_AXIS , self.RADIAL_BIN , self.HEIGHT_BIN )
            rmidpoint = []
            hmidpoint = []
            for i in range( avg_prof.Nrbins ) :
                rmidpoint.append( 0.5 * ( avg_prof.rmin[i] + avg_prof.rmax[i] ) )
            hmidpoint.append( 0.0 )
            for i in range( 1 , avg_prof.lastHbin+1 ) :
                hmidpoint.append( 0.5 * ( avg_prof.hmin[i] + avg_prof.hmax[i] ) )
            for i in range( -avg_prof.lastHbin , 0 ) :
                hmidpoint.append( 0.5 * ( avg_prof.hmin[i] + avg_prof.hmax[i] ) )
            N = len( self.PROFILES )
            for profile in self.PROFILES :
                for i in range( len(hmidpoint) ) :
                    for j in range( len(rmidpoint) ) :
                        avg_prof.counts[i][j] += profile.counts[i][j]
                        avg_prof.dens[i][j] += ( profile.dens[i][j] / N )
            self.AVG_file.write( str( self.last_step-self.first_step ) + '\t' + str( len(hmidpoint)*len(rmidpoint) ) + '\t' + str( avg_prof.total_count() ) + '\n' )
            for j in range( len(rmidpoint) ) :
                for i in range( len(hmidpoint) ) :
                    ii = i - ( ( len(hmidpoint) - 1 ) / 2 )
                    self.AVG_file.write( str( i+1+j*len(hmidpoint) ) + '\t' + repr( hmidpoint[ii] ) + '\t' + repr( rmidpoint[j] ) + '\t' + repr( avg_prof.counts[ii][j] ) + '\t' + repr( avg_prof.dens[ii][j] ) + '\n' )
            self.AVG_file.close()
            if self.VARIANCE :
                profile_variance = pf.cylindrical_profile( self.last_box , self.CYLINDER_AXIS , self.RADIAL_BIN , self.HEIGHT_BIN )
                self.VARIANCE_file.write( "# Chunk hmidpoint rmidpoint density/number dev_stand\n" )
                for profile in self.PROFILES :
                    for i in range( len(hmidpoint) ) :
                        for j in range( len(rmidpoint) ) :
                            profile_variance.dens[i][j] += ( profile.dens[i][j] - avg_prof.dens[i][j] )**2 / N
                for j in range( len(rmidpoint) ) :
                    for i in range( len(hmidpoint) ) :
                        ii = i - ( ( len(hmidpoint) - 1 ) / 2 )
                        self.VARIANCE_file.write( str( i+1+j*len(hmidpoint) ) + '\t' + repr( hmidpoint[ii] ) + '\t' + repr( rmidpoint[j] ) + '\t' + repr( avg_prof.dens[ii][j] ) + '\t' + repr( np.sqrt(profile_variance.dens[ii][j]) ) + '\n' )
                self.VARIANCE_file.close()

        elif self.PROFILE_TYPE == "linear" :
            self.PERCONF_file.close()
            self.AVG_file.write( "# Time-averaged density profile\n" )
            self.AVG_file.write( "# Timesteps Number-of-chunks Total-count\n" )
            self.AVG_file.write( "# Chunk midpoint Ncount density/number\n" )
            FLAG_DIFF_BOXES = 0
            hmax0 = self.PROFILES[0].hmax[-1]
            for profile in self.PROFILES :
                if hmax0 != profile.hmax[-1] :
                    FLAG_DIFF_BOXES = 1
            if FLAG_DIFF_BOXES :
                # if you want to take into account a series of configurations with variable box sides !!
                print( "*** ERROR: This part of the code should still be written! Within density_profile module, linear density profile" )
                exit(3)
            else :
                avg_prof = pf.linear_profile( self.last_box , self.CYLINDER_AXIS , self.HEIGHT_BIN )
                midpoint = []
                for i in range( avg_prof.Nbins ) :
                    midpoint.append( 0.5 * ( avg_prof.hmin[i] + avg_prof.hmax[i] ) )
                N = len( self.PROFILES )
                for profile in self.PROFILES :
                    for i in range( profile.Nbins ) :
                        avg_prof.counts[i] += profile.counts[i]
                        avg_prof.dens[i] += profile.dens[i] / N
            self.AVG_file.write( str( self.last_step-self.first_step ) + '\t' + str( avg_prof.Nbins ) + '\t' + str( avg_prof.total_count() ) + '\n' )
            for i in range( avg_prof.Nbins ) :
                self.AVG_file.write( str( i+1 ) + '\t' + repr( midpoint[i] ) + '\t' + repr( avg_prof.counts[i] ) + '\t' + repr( avg_prof.dens[i] ) + '\n' )
            self.AVG_file.close()
            if self.VARIANCE :
                profile_variance = pf.linear_profile( self.last_box , self.CYLINDER_AXIS , self.HEIGHT_BIN )
                self.VARIANCE_file.write( "# Chunk midpoint density/number dev_stand\n" )
                for profile in self.PROFILES :
                    for i in range( profile.Nbins ) :
                        profile_variance.dens[i] += ( profile.dens[i] - avg_prof.dens[i] )**2 / N
                for i in range( avg_prof.Nbins ) :
                    self.VARIANCE_file.write( str( i+1 ) + '\t' + repr( midpoint[i] ) + '\t' + repr( avg_prof.dens[i] ) + '\t' + repr( np.sqrt(profile_variance.dens[i]) ) + '\n' )
                self.VARIANCE_file.close()

    def return_values( self ) :
        self.lock.acquire()
        self.PERCONF_file.flush()
        self.lock.release()
        return { "PROFILES" : self.PROFILES ,
                 "last_step" : self.last_step ,
                 "first_step" : self.first_step ,
                 "last_box" : [ self.last_box[0] , self.last_box[1] , self.last_box[2] ] }

    def merge_return_values( self , values_list ) :
        for partial in values_list :
            self.PROFILES = self.PROFILES + partial["PROFILES"]
            self.last_step = partial["last_step"]
            self.first_step = partial["first_step"]
            self.last_box[0] = partial["last_box"][0]
            self.last_box[1] = partial["last_box"][1]
            self.last_box[2] = partial["last_box"][2]

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-profile\033[0m [ sphere <bin_width> | cylinder <axis:x|y|z> <height_width> <radial_thickness> | linear <axis:x|y|z> <bin_width> ] " )
        print( "\t\t\t\t{ specie <mol|type> <mol_id|atom_type> [default:all] }" )
        print( "\t\t\t\t{ center <mol|type|fixed> <to_compute_com:mol_id|atom_type|x,y,z> [default:all] }" )
        print( "\t\t\t\t{ file <fname_prof_perconf:atoms_profile.dat> } { fixed_diameter <box_side_value> } { variance <yes|no> } } " )
