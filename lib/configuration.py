import os
import math
import random
import molecule as molobj
import profiles as pf
import chains as ch
from tiles import Tile
import numpy as np
from scipy.spatial import ConvexHull
from itertools import islice
EXIT_FAILURE = -1
X = 0
Y = 1
Z = 2

# OVITO # for surface meshes
try :
    from ovito.io import *
    from ovito.modifiers import *
    from ovito.data import *
    from ovito.pipeline import *
except :
    print("***WARNING: module ovito not found! This is needed to use configuration.surface_mesh() method.")
    pass

class BOND :
    def __init__( self , bid = 0 , btype = 0 , p1 = 0 , p2 = 0 ) :
        self.id = bid
        self.type = btype
        self.p1 = p1
        self.p2 = p2

class ANGLE :
    def __init__( self , aid = 0 , atype = 0 , p1 = 0 , p2 = 0 , p3 = 0 ) :
        self.id = aid
        self.type = atype
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.K = 0
        self.theta0 = 180.0  ## degree units !!!!!!!

# definition of the class configuration
class CONFIGURATION :

    def __init__( self ) :
        self.time = 0.0
        self.N = 0
        self.COM = np.float64([ 0.0 , 0.0 , 0.0 ])
        self.VCOM = np.float64([ 0.0 , 0.0 , 0.0 ])
        self.box_inf = np.float64([ 0.0 , 0.0 , 0.0 ])
        self.box_sup = np.float64([ 0.0 , 0.0 , 0.0 ])
        self.pos = np.float64( [] )
        self.vel = np.float64( [] )
        self.box_img = np.int32( [] )
        self.q = np.float64( [] )
        self.mass = np.float64( [] )
        self.mol = np.int32( [] )
        self.type = np.int32( [] )
        self.id = np.int32( [] )
        self.radius = np.float64( [] )
        self.patch_params = np.float64( [] ) # For patchy particles, Nx2 with cos(halfopening),lambda/sigma
        self.rot_matrix = np.float64( [] )   # 3D rotation matrix [[cos sin ...]...]
        self.molecules = {}
        self.Nchunks = 0
        self.chunk = []
        self.chain = []
        self.bond = []
        self.angle = []
        self.description = ""
        self.cells = []
        self.contacts = []
        self.Ncontacts = 0
        self.lofcontacts = np.float64( [] )

    def insert( self , npos = np.float64([]) , nid = np.int32([]) , ntype = np.int32([]) , nmol = np.int32([]) , nq = np.float64([]) , nvel = np.float64([]) ) :
        if ( len(self.pos) > 0 or self.N == 0 ) and len(npos) > 0 :
            if len(self.box_img) > 0 :
                ref = ( self.box_inf + self.box_sup ) * 0.5
                box_sides = self.box_sup - self.box_inf
                # in case of purely 2D systems :
                D3T2D = False
                if len(box_sides) == 3 :
                    if box_sides[Z] == 0.0 :
                        D3T2D = True
                if D3T2D :
                    box_sides[Z] = 1.0
                    nbox_img = ( np.sign(npos-ref) * (  ( np.fabs(npos-ref) // ( box_sides*0.5 ) + 1 ) // 2  ) ).astype('int32')
                    box_sides[Z] = 0.0
                    nbox_img[:,Z] = 0
                else :
                    nbox_img = ( np.sign(npos-ref) * (  ( np.fabs(npos-ref) // ( box_sides*0.5 ) + 1 ) // 2  ) ).astype('int32')
                npos = npos - ( nbox_img * box_sides )
                self.box_img = np.append( self.box_img , nbox_img , axis=0 )
            self.pos = np.append( self.pos , npos , axis=0 )
        if ( len(self.id) > 0 or self.N == 0 ) and len(nid) > 0 :
            self.id = np.append( self.id , nid )
        if ( len(self.type) > 0 or self.N == 0 ) and len(ntype) > 0 :
            self.type = np.append( self.type , ntype )
        if ( len(self.mol) > 0 or self.N == 0 ) and len(nmol) > 0 :
            self.mol = np.append( self.mol , nmol )
        if ( len(self.q) > 0 or self.N == 0 ) and len(nq) > 0 :
            self.q = np.append( self.q , nq )
        if ( len(self.vel) > 0 or self.N == 0 ) and len(nvel) > 0 :
            self.vel = np.append( self.vel , nvel , axis=0 )
        self.N = max( [ len(self.pos) , len(self.vel) , len(self.id) , len(self.type) , len(self.q) , len(self.mol) ] )

    def insert_bond( self , bid , btype , p1 , p2 ) :
        self.bond.append( BOND( bid , btype , p1 , p2 ) )

    def insert_angle( self , aid , atype , p1 , p2 , p3 ) :
        self.angle.append( ANGLE( aid , atype , p1 , p2 , p3 ) )

    def compute_periodic_images( self ) :
        if len( self.box_img ) == 0 :
            ref = ( self.box_inf + self.box_sup ) * 0.5
            box_sides = self.box_sup - self.box_inf
            # in case of purely 2D systems :
            D3T2D = False
            if len(box_sides) == 3 :
                if box_sides[Z] == 0.0 :
                    D3T2D = True
            if D3T2D :
                box_sides[Z] = 1.0
                self.box_img = ( np.sign(self.pos-ref) * (  ( np.fabs(self.pos-ref) // ( box_sides*0.5 ) + 1 ) // 2  ) ).astype('int32')
                box_sides[Z] = 0.0
                self.box_img[:,Z] = 0
            else :
                self.box_img = ( np.sign(self.pos-ref) * (  ( np.fabs(self.pos-ref) // ( box_sides*0.5 ) + 1 ) // 2  ) ).astype('int32')
            self.pos = self.pos - ( self.box_img.astype('float64') * box_sides )
        else :
            print( "WARNING : attempting to compute periodic images, but they are already present !" )

    def unwrapped_coordinates( self , MODIFY = False ) :
        if len( self.box_img ) == self.N :
            if MODIFY :
                self.pos = self.pos + ( self.box_img.astype('float64') * ( self.box_sup - self.box_inf ) )
                self.box_img = np.int32( [] )
                return self.pos
            else :
                ppos = self.pos + ( self.box_img.astype('float64') * ( self.box_sup - self.box_inf ) )
                return ppos
        else :
            return self.pos

    def periodic_image( self , points , ref=np.float64([0.0,0.0,0.0]) , box_sides=np.float64([]) ) :
        if len( box_sides ) == 0 :
            box_sides = self.box_sup - self.box_inf
        # in case of purely 2D systems :
        D3T2D = False
        if len(box_sides) == 3 :
            if box_sides[Z] == 0.0 :
                D3T2D = True
        if D3T2D :
            box_sides[Z] = 1.0
            periodic_images = np.sign(points-ref) * (  ( np.fabs(points-ref) // ( box_sides*0.5 ) + 1 ) // 2  )
            box_sides[Z] = 0.0
            periodic_images[:,Z] = 0
        else :
            periodic_images = np.sign(points-ref) * (  ( np.fabs(points-ref) // ( box_sides*0.5 ) + 1 ) // 2  )
        return periodic_images

    def create_chunks( self , type , dr , Nchunks , box_sides ) :
        if type == "sphere" :
            self.Nchunks = Nchunks
            for i in range( Nchunks ) :
                list_of_parts = []
                self.chunk.append( list_of_parts )
            periodic_images = self.periodic_image( self.pos , box_sides=box_sides )
            dist = np.sqrt( np.add.reduce( ( self.pos-(periodic_images*box_sides) )**2 , 1 ) )
            chunk_id = np.floor( dist/dr )
            chunk_id = chunk_id.astype( 'int32' )
            chunk_id = np.where( chunk_id < Nchunks , chunk_id , [Nchunks-1]*self.N )
            for i in range( self.N ) :
                self.chunk[chunk_id[i]].append( i )
        else :
            print( "ERROR: The type of chunks has not been recognised" )
            exit( EXIT_FAILURE )

    def split_confs( self , data_file , directory = "single_confs" , fmt = "lmp" , block = -1 ) :
        time_steps = []
        part_numbers = []
        line = data_file.readline()
        if block > 0 :
            out_file_tless = open( "cnf-tless-"+str(block)+".dat" , "w" )
            out_file_tgtr = open( "cnf-tgtr-"+str(block)+".dat" , "w" )
        if line == "" :
            data_file.close()
            return time_steps , part_numbers
        else :
            if block > 0 :
                directory = "."
            if not os.path.exists( directory ) :
                os.system( "mkdir " + directory )
            if fmt == "lmp" :
                if block <= 0 :
                    out_file = open( directory + "/.null.dat" , "w" )
                while ( not line.isspace() ) and line != "" :
                    line.strip()
                    words = line.split()
                    if len(words) > 1 :
                        if words[1] == "TIMESTEP" :
                            line2 = data_file.readline()
                            line2.strip()
                            words2 = line2.split()
                            timestep = float(words2[0])
                            time_steps.append( timestep )
                            if block > 0 :
                                if timestep < block :
                                    out_file = out_file_tless
                                else :
                                    out_file = out_file_tgtr
                            else :
                                out_file.close()
                                out_file = open( directory+"/cnf-"+str(timestep)+".dat" , "w" )
                            out_file.write( line )
                            out_file.write( line2 )
                        elif words[1] == "NUMBER" :
                            line2 = data_file.readline()
                            line2.strip()
                            words2 = line2.split()
                            part_N = int(words2[0])
                            part_numbers.append( part_N )
                            out_file.write( line )
                            out_file.write( line2 )
                        elif words[1] == "ATOMS" :
                            out_file.write( line )
                            particles_lines = islice( data_file , part_N )
                            for pline in particles_lines :
                                out_file.write( pline )
                        else :
                            out_file.write( line )
                    else :
                        out_file.write( line )
                    line = data_file.readline()
                if block <= 0 :
                    out_file.close()
                    os.system( "rm " + directory + "/.null.dat" )

            elif fmt == "xyz" :
                BOX_YES = 0
                while line != "" :
                    line.strip()
                    words = line.split()
                    if len(words) > 1 :
                        if words[0] == "#" :
                            if words[1] == "N" :
                                N_particles = int(words[2])
                            elif words[1] in [ "step" , "timestep" ] :
                                timestep = float(words[2])
                            elif words[1] == "box" :
                                BOX_YES = 1
                                box_sides = line
                        else :
                            time_steps.append( timestep )
                            part_numbers.append( N_particles )
                            if block > 0 :
                                if timestep < block :
                                    out_file = out_file_tless
                                else :
                                    out_file = out_file_tgtr
                            else :
                                out_file = open( directory + "/cnf-" + str(timestep) + ".dat" , "w" )
                            out_file.write( "# N " + str(N_particles) + "\n" )
                            out_file.write( "# timestep " + str(timestep) + "\n" )
                            if BOX_YES :
                                out_file.write( box_sides )
                            BOX_YES = 0
                            out_file.write( line )
                            for pline in islice( data_file , N_particles-1 ) :
                                out_file.write( pline )
                            if block <= 0 :
                                out_file.close()
                    line = data_file.readline()

            elif fmt in [ "sph" , "patch" , "ptc" ] :
                timestep = 0
                while line != "" :
                    line.strip()
                    words = line.split()
                    if len(words) > 0 and words[0][0] != "#" :
                        if len(words) in [ 1, 2 ] and words[0].strip('&').isdigit() :
                            N_particles = int(words[0].strip('&'))
                            if len(words) > 1 :
                                timestep = float(words[1])
                            else :
                                timestep += 1
                            time_steps.append( timestep )
                            box_sides = data_file.readline()
                            part_numbers.append( N_particles )
                            if block > 0 :
                                if timestep < block :
                                    out_file = out_file_tless
                                else :
                                    out_file = out_file_tgtr
                            else :
                                out_file = open( directory + "/cnf-" + str(timestep) + ".dat" , "w" )
                            if fmt == "sph" :
                                out_file.write( str(N_particles) + " " + str(timestep) + "\n" )
                            else :
                                out_file.write( "&" + str(N_particles) + " " + str(timestep) + "\n" )
                            out_file.write( box_sides )
                            particles_lines = islice( data_file , N_particles )
                            for pline in particles_lines :
                                out_file.write( pline )
                            if block <= 0 :
                                out_file.close()
                    line = data_file.readline()

            else :
                print( " *** ERROR: input format not recognized in function configuration.split_confs()" )
                exit
            if block > 0 :
                out_file_tless.close()
                out_file_tgtr.close()
            return time_steps , part_numbers




    def change_box( self , data_file_name ) :
        # reads in the box boundaries stored in a file in LAMMPS format
        success = False
        data_file = open( data_file_name , "r" )
        lines = data_file.readlines()
        data_file.close()
        for i in range( len(lines) ) :
            words = lines[i].strip().split()
            if len( words ) > 1 :
                if words[1] == "BOX" :
                    words = lines[i+1].strip().split()
                    self.box_inf[X] = float( words[0] )
                    self.box_sup[X] = float( words[1] )
                    words = lines[i+2].strip().split()
                    self.box_inf[Y] = float( words[0] )
                    self.box_sup[Y] = float( words[1] )
                    words = lines[i+3].strip().split()
                    self.box_inf[Z] = float( words[0] )
                    self.box_sup[Z] = float( words[1] )
                    success = True
        if success == False :
            print( "*** ERROR in reading box boundaries from file " + data_file_name + "!" )
            exit( EXIT_FAILURE )
        else :
            box_side = self.box_sup - self.box_inf
            return box_side




    def fast_timestep_read( self , fname , fmt = "lmp" ) :
        TIMESTEP = -1
        data_file = open( fname , "r" )
        line = data_file.readline()
        if fmt == "lmp" :
            TSTEP_READ = False
            while TSTEP_READ == False and line != '' :
                line.strip()
                words = line.split()
                if len(words) > 1 :
                    if words[1] == "TIMESTEP" :
                        line = data_file.readline()
                        line.strip()
                        words = line.split()
                        TIMESTEP = float(words[0])
                        TSTEP_READ = True
                    line = data_file.readline()

        elif fmt == "xyz" :
            line.strip()
            words = line.split()
            while len(words) == 0 and line != "" :
                line = data_file.readline()
                line.strip()
                words = line.split()
            TSTEP_READ = False
            while TSTEP_READ == False and line != '' :
                line.strip()
                words = line.split()
                if len(words) > 1 :
                    if words[1] in [ "step" , "timestep" ] :
                        TIMESTEP = float(words[2])
                        TSTEP_READ = True
                line = data_file.readline()

        elif fmt in [ "sph" , "patch" , "ptc" ] :
            line.strip()
            words = line.split()
            if len(words) == 1 :
                TIMESTEP = -1.0
            else :
                TIMESTEP = float(words[1])

        data_file.close()
        return TIMESTEP

    def smart_auto_read( self , data_file , fmt = "lmp" , VELOCITIES = "NO" ) :
        line = data_file.readline()
        if line == "" :
            data_file.close()
            return -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1
        else :
            if fmt == "lmp" :

                id_col = -1
                type_col = -1
                mol_col = -1
                q_col = -1
                x_col = -1
                y_col = -1
                z_col = -1
                ix_col = -1
                iy_col = -1
                iz_col = -1
                vx_col = -1
                vy_col = -1
                vz_col = -1
                box_side = np.float64([ 0.0 , 0.0 , 0.0 ])
                if line != "" :
                    ATOMS_READ = False
                    while ATOMS_READ == False and line != '' :
                        line.strip()
                        words = line.split()
                        if len(words) > 1 :

                            if words[1] == "NUMBER" :
                                line = data_file.readline()
                                line.strip()
                                words = line.split()
                                self.N = int(words[0])

                            elif words[1] == "BOX" :
                                line = data_file.readline()
                                line.strip()
                                words = line.split()
                                self.box_inf[X] = float( words[0] )
                                self.box_sup[X] = float( words[1] )
                                line = data_file.readline()
                                line.strip()
                                words = line.split()
                                self.box_inf[Y] = float( words[0] )
                                self.box_sup[Y] = float( words[1] )
                                line = data_file.readline()
                                line.strip()
                                words = line.split()
                                self.box_inf[Z] = float( words[0] )
                                self.box_sup[Z] = float( words[1] )
                                box_side = self.box_sup - self.box_inf

                            elif words[1] == "TIMESTEP" :
                                line = data_file.readline()
                                line.strip()
                                words = line.split()
                                self.time = float(words[0])

                            elif words[1] == "ATOMS" :
                                ATOMS_READ = True
                                for i in range( 2 , len(words) ) :
                                    if words[i] == "type" :
                                        type_col = i-2
                                    elif words[i] == "mol" :
                                        mol_col = i-2
                                    elif words[i] == "x" or words[i] == "xu" :
                                        x_col = i-2
                                    elif words[i] == "y" or words[i] == "yu" :
                                        y_col = i-2
                                    elif words[i] == "z" or words[i] == "zu" :
                                        z_col = i-2
                                    elif words[i] == "ix" :
                                        ix_col = i-2
                                    elif words[i] == "iy" :
                                        iy_col = i-2
                                    elif words[i] == "iz" :
                                        iz_col = i-2
                                    elif words[i] == "vx" :
                                        vx_col = i-2
                                    elif words[i] == "vy" :
                                        vy_col = i-2
                                    elif words[i] == "vz" :
                                        vz_col = i-2
                                    elif words[i] == "id" :
                                        id_col = i-2
                                    elif words[i] == "q" :
                                        q_col = i-2
                                # words = np.array( [ data_file.readline().strip().split() for i in range( self.N ) ] )
                                if hasattr( data_file , 'loadtxt' ) and callable( data_file.loadtxt ) :
                                    words = data_file.loadtxt( self.N )  # way faster
                                else :
                                    words = np.loadtxt( data_file , dtype = 'str' , max_rows = self.N )
                                if id_col != -1 :
                                    self.id = words[:,id_col]
                                    self.id = self.id.astype('int32')
                                if type_col != -1 :
                                    self.type = words[:,type_col]
                                    self.type = self.type.astype('int32')
                                if q_col != -1 :
                                    self.q = words[:,q_col]
                                    self.q = self.q.astype('float64')
                                if mol_col != -1 :
                                    self.mol = words[:,mol_col]
                                    self.mol = self.mol.astype('int32')
                                xarr = np.float64( [] )
                                yarr = np.float64( [] )
                                zarr = np.float64( [] )
                                if x_col != -1 :
                                    xarr = words[:,x_col]
                                    xarr = xarr.astype('float64')
                                if y_col != -1 :
                                    yarr = words[:,y_col]
                                    yarr = yarr.astype('float64')
                                if z_col != -1 :
                                    zarr = words[:,z_col]
                                    zarr = zarr.astype('float64')
                                if x_col != -1 and z_col == -1 :
                                    self.pos = np.column_stack(( xarr , yarr ))
                                else :
                                    self.pos = np.column_stack(( xarr , yarr , zarr ))
                                xarr = np.int32( [] )
                                yarr = np.int32( [] )
                                zarr = np.int32( [] )
                                if ix_col != -1 :
                                    xarr = words[:,ix_col]
                                    xarr = xarr.astype('int32')
                                if iy_col != -1 :
                                    yarr = words[:,iy_col]
                                    yarr = yarr.astype('int32')
                                if iz_col != -1 :
                                    zarr = words[:,iz_col]
                                    zarr = zarr.astype('int32')
                                if ix_col != -1 and iz_col == -1 :
                                    self.box_img = np.column_stack(( xarr , yarr ))
                                else :
                                    self.box_img = np.column_stack(( xarr , yarr , zarr ))
                                ## change : DA ORA SELF.POS NON È PIÙ AUTOMATICAMENTE UNWRAPPED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                if VELOCITIES == "YES" :
                                    xarr = np.float64( [] )
                                    yarr = np.float64( [] )
                                    zarr = np.float64( [] )
                                    if vx_col != -1 :
                                        xarr = words[:,vx_col]
                                        xarr = xarr.astype('float64')
                                    if vy_col != -1 :
                                        yarr = words[:,vy_col]
                                        yarr = yarr.astype('float64')
                                    if vz_col != -1 :
                                        zarr = words[:,vz_col]
                                        zarr = zarr.astype('float64')
                                    if x_col != -1 and z_col == -1 :
                                        self.vel = np.column_stack(( xarr , yarr ))
                                    else :
                                        self.vel = np.column_stack(( xarr , yarr , zarr ))
                        if ATOMS_READ == False :
                            line = data_file.readline()
                return id_col , type_col , x_col , y_col , z_col , q_col , mol_col , vx_col , vy_col , vz_col


            elif fmt == "xyz" :

                line.strip()
                words = line.split()
                INITIAL_INFO = 1
                DIMENSION = 3
                while INITIAL_INFO :
                    if line == "" :
                        INITIAL_INFO = 0
                    elif len(words) == 0 :
                        line = data_file.readline()
                        line.strip()
                        words = line.split()
                    elif words[0] == "#" :
                        if words[1] == "N" :
                            self.N = int(words[2])
                        elif words[1] in [ "step" , "timestep" ] :
                            self.time = float(words[2])
                        elif words[1] == "box" :
                            self.box_inf[X] = 0.0
                            self.box_sup[X] = float(words[2])
                            self.box_inf[Y] = 0.0
                            self.box_sup[Y] = float(words[3])
                            self.box_inf[Z] = 0.0
                            if len(words) > 4 :
                                self.box_sup[Z] = float(words[4])
                                DIMENSION = 3
                            else :
                                self.box_sup[Z] = 0
                                DIMENSION = 2
                        line = data_file.readline()
                        line.strip()
                        words = line.split()
                    else :
                        INITIAL_INFO = 0
                if len(words) > 1 and line != "" :
                    if hasattr( data_file , 'loadtxt' ) and callable( data_file.loadtxt ) :
                        words = np.concatenate( ( np.array([words]) , data_file.loadtxt( self.N-1 ) ) , axis=0 )  # way faster
                    else :
                        words = np.concatenate( ( np.array([words]) , np.loadtxt( data_file , dtype = 'str' , max_rows = self.N-1 ) ) , axis=0 )
                    # words = np.array( [ data_file.readline().strip().split() for i in range( self.N ) ] )
                    if VELOCITIES == "YES" :
                        self.vel = words[:,0:3].astype('float64')
                    else :
                        if DIMENSION == 3 :
                            self.pos = words[:,0:3].astype('float64')
                            if len(words[0]) > 5 :
                                self.box_img = words[:,3:6].astype('int32')
                        if DIMENSION == 2 :
                            self.pos = np.zeros((len(words),3)).astype('float64')
                            self.pos[:,0:2] = words[:,0:2].astype('float64')
                            if len(words[0]) > 3 :
                                self.box_img = np.zeros((len(words),3)).astype('int32')
                                self.box_img[:,0:2] = words[:,2:4].astype('int32')
                if VELOCITIES == "YES" :
                    return -1 , -1 , -1 , -1 , -1 , -1 , -1 , 1 , 2 , 3
                else :
                    return -1 , -1 , 1 , 2 , 3 , -1 , -1 , -1 , -1 , -1

            elif fmt in [ "sph" , "patch" , "ptc" ] :

                line.strip()
                words = line.split()
                while ( len(words) == 0 or words[0][0] == "#" ) and line != "" :
                    line = data_file.readline()
                    line.strip()
                    words = line.split()
                if len(words) > 0 and line != "" :
                    if len(words) in [ 1, 2 ] and words[0].strip('&').isdigit() :
                        self.N = int(words[0].strip('&'))
                        if len(words) > 1 :
                            self.time = float(words[1])
                        else :
                            self.time = -1
                        line = data_file.readline()
                        words = line.strip().split()
                        self.box_inf[X] = 0.0
                        self.box_inf[Y] = 0.0
                        self.box_sup[X] = float( words[0] )
                        self.box_sup[Y] = float( words[1] )
                        self.box_inf[Z] = 0.0
                        if len(words) > 2 :
                            self.box_sup[Z] = float( words[2] )
                        else :
                            self.box_sup[Z] = 0.0
                        box_side = self.box_sup - self.box_inf
                        # words = np.array( [ data_file.readline().strip().split() for i in range( self.N ) ] )
                        if hasattr( data_file , 'loadtxt' ) and callable( data_file.loadtxt ) :
                            words = data_file.loadtxt( self.N )  # way faster
                        else :
                            words = np.loadtxt( data_file , dtype = 'str' , max_rows = self.N )
                        self.pos = words[:,1:4].astype('float64')
                        self.radius = words[:,4].astype('float64')
                        types = words[:,0]
                        self.type = np.array( [ ord(types[i])-ord('a')+1 for i in range( self.N ) ] ).astype('int32')
                        self.id = np.arange( 1 , self.N+1 , dtype=int )
                        if fmt in [ "patch" , "ptc" ] :
                            self.patch_params = words[:,5:7].astype('float64')
                            self.rot_matrix = words[:,7:16].astype('float64')
                    else :
                        print( " *** ERROR unknown . " )
                        exit
                if VELOCITIES == "YES" :
                    print( " *** ERROR: the VELOCITIES option is not available for format sph" )
                    return -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1
                else :
                    return -1 , 0 , 1 , 2 , 3 , -1 , -1 , -1 , -1 , -1

            else :
                print( " *** ERROR: input format not recognized in function configuration.smart_auto_read()" )
                exit


    def smart_read( self , data_file , id_col , type_col , mol_col , q_col , x_col , y_col , z_col ) :
        # if columns id are negative they are ignored
        words = np.array( [ data_file.readline().strip().split() for i in range( self.N ) ] )
        if id_col != -1 :
            self.id = words[:,id_col]
            self.id = self.id.astype('int32')
        if type_col != -1 :
            self.type = words[:,type_col]
            self.type = self.type.astype('int32')
        if q_col != -1 :
            self.q = words[:,q_col]
            self.q = self.q.astype('float64')
        if mol_col != -1 :
            self.mol = words[:,mol_col]
            self.mol = self.mol.astype('int32')
        xarr = np.float64( [] )
        yarr = np.float64( [] )
        zarr = np.float64( [] )
        if x_col != -1 :
            xarr = words[:,x_col]
            xarr = xarr.astype('float64')
        if y_col != -1 :
            yarr = words[:,y_col]
            yarr = yarr.astype('float64')
        if z_col != -1 :
            zarr = words[:,z_col]
            zarr = zarr.astype('float64')
        if x_col != -1 and z_col == -1 :
            self.pos = np.column_stack(( xarr , yarr ))
        else :
            self.pos = np.column_stack(( xarr , yarr , zarr ))

    def print_selection( self , out_file , sel_mode , sel_val , out_format , id_col , type_col , mol_col , q_col , x_col , y_col , z_col , vx_col=-1 , vy_col=-1 , vz_col=-1 ) :
        # out_file has to be already opened
        # if columns id are negative they are ignored
        # selection of the atoms to be printed out
        if out_format != "lmp" :
            ppos = self.unwrapped_coordinates()
        else :
            ppos = self.pos
        N = self.N
        if sel_mode == "all" :
            pvel = self.vel
            pid = self.id
            ptype = self.type
            pq = self.q
            pmol = self.mol
            pimg = self.box_img
            pradius = self.radius
            ppatch_params = self.patch_params
            prot_matrix = self.rot_matrix
        elif sel_mode == "mol" :
            select_vector = np.equal( self.mol , int(sel_val) )
        elif sel_mode == "type" :
            select_vector = np.equal( self.type , int(sel_val) )
        elif sel_mode == "q" :
            select_vector = np.equal( self.q , float(sel_val) )
        elif sel_mode == "custom" :
            select_vector = np.array( sel_val ).astype(bool)
            if len(select_vector) != self.N :
                print( "*** ERROR: Invalid custom selection in print_selection() method !" )
        else :
            print( "*** ERROR: Selection option is needed! Possible values: all, mol, type, q, custom ." )
        if sel_mode != "all" :
            if ( len(self.pos) > 0 ) : ppos = ppos[ select_vector , : ]
            if ( len(self.vel) > 0 ) : pvel = self.vel[ select_vector , : ]
            if ( len(self.id) > 0 ) : pid = self.id[ select_vector ]
            if ( len(self.type) > 0 ) : ptype = self.type[ select_vector ]
            if ( len(self.q) > 0 ) : pq = self.q[ select_vector ]
            if ( len(self.mol) > 0 ) : pmol = self.mol[ select_vector ]
            if ( len(self.box_img) > 0 ) : pimg = self.box_img[ select_vector , : ]
            if ( len(self.radius) > 0 ) : pradius = self.radius[ select_vector ]
            if ( len(self.patch_params) > 0 ) : ppatch_params = self.patch_params[ select_vector , : ]
            if ( len(self.rot_matrix) > 0 ) : prot_matrix = self.rot_matrix[ select_vector , : ]
            N = np.count_nonzero( select_vector )
        # global information
        if out_format == "lmp" :
            out_file.write( "ITEM: TIMESTEP\n" )
            out_file.write( str( self.time ) + '\n' )
            out_file.write( "ITEM: NUMBER OF ATOMS\n" )
            out_file.write( str( N ) + '\n' )
            out_file.write( "ITEM: BOX BOUNDS pp pp pp\n" )
            out_file.write( str( self.box_inf[X] ) + ' ' + str( self.box_sup[X] ) + '\n' )
            out_file.write( str( self.box_inf[Y] ) + ' ' + str( self.box_sup[Y] ) + '\n' )
            out_file.write( str( self.box_inf[Z] ) + ' ' + str( self.box_sup[Z] ) + '\n' )
            out_file.write( "ITEM: ATOMS" )
            if id_col != -1 :
                out_file.write( " id" )
            if type_col != -1 :
                out_file.write( " type" )
            if x_col != -1 :
                out_file.write( " x" )
            if y_col != -1 :
                out_file.write( " y" )
            if z_col != -1 :
                out_file.write( " z" )
            if q_col != -1 :
                out_file.write( " q" )
            if mol_col != -1 :
                out_file.write( " mol" )
            if vx_col != -1 :
                out_file.write( " vx" )
            if vy_col != -1 :
                out_file.write( " vy" )
            if vz_col != -1 :
                out_file.write( " vz" )
            if len(self.box_img) > 0 :
                out_file.write( " ix iy" )
                if z_col != -1 :
                    out_file.write( " iz" )
            out_file.write( "\n" )
        elif out_format == "xyz" :
            out_file.write( "# N " + str( N ) + "\n" )
            out_file.write( "# timestep " + str( self.time ) + "\n" )
            box_sides = self.box_sup - self.box_inf
            out_file.write( "# box " + str(box_sides[X]) + ' ' + str(box_sides[Y]) + ' ' + str(box_sides[Z]) + "\n" )
        elif out_format in [ "sph" , "patch" , "ptc" ] :
            if out_format == "sph" :
                out_file.write( str( N ) + " " + str( self.time ) + "\n" )
            else :
                out_file.write( "&" + str( N ) + " " + str( self.time ) + "\n" )
            box_sides = self.box_sup - self.box_inf
            if box_sides[Z] == 0.0 :
                out_file.write( str(box_sides[X]) + ' ' + str(box_sides[Y]) + '\n' )
            else :
                out_file.write( str(box_sides[X]) + ' ' + str(box_sides[Y]) + ' ' + str(box_sides[Z]) + '\n' )
            if len(pradius) == 0 :
                pradius = np.array( [0.5]*N )
            if len(ptype) == 0 :
                ptype = np.array( [0]*N )
            if out_format in [ "patch" , "ptc" ] :
                if len(ppatch_params) == 0 :
                    ppatch_params = np.array( [[0.9925,1.12]]*N )
                if len(prot_matrix) == 0 :
                    prot_matrix = np.array( [[1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0]]*N )
        elif out_format == "zeno_input" :
            out_file.write( str( N ) + "\n" )
            out_file.write( "# timestep " + str( self.time ) + "\n" )
        else :
            print( "*** ERROR: output format not recognized!" )
            exit( EXIT_FAILURE )
        # printing out particles' positions and attributes
        if out_format == "lmp" :
            for i in range( N ) :
                if id_col != -1 :
                    out_file.write( str( pid[i] ) + " " )
                if type_col != -1 :
                    out_file.write( str( ptype[i] ) + " " )
                if x_col != -1 :
                    out_file.write( str( ppos[i,X] ) + " " )
                if y_col != -1 :
                    out_file.write( str( ppos[i,Y] ) + " " )
                if z_col != -1 :
                    out_file.write( str( ppos[i,Z] ) + " " )
                if q_col != -1 :
                    out_file.write( str( pq[i] ) + " " )
                if mol_col != -1 :
                    out_file.write( str( pmol[i] ) + " " )
                if vx_col != -1 :
                    out_file.write( str( pvel[i,X] ) + " " )
                if vy_col != -1 :
                    out_file.write( str( pvel[i,Y] ) + " " )
                if vz_col != -1 :
                    out_file.write( str( pvel[i,Z] ) + " " )
                if len(self.box_img) > 0 :
                    if z_col != -1 :
                        out_file.write( str( pimg[i,X] ) + " " + str( pimg[i,Y] ) + " " + str( pimg[i,Z] ) + " " )
                    else :
                        out_file.write( str( pimg[i,X] ) + " " + str( pimg[i,Y] ) + " " + " " )
                out_file.write( "\n" )
        elif out_format == "xyz" :
            for i in range( N ) :
                if x_col != -1 :
                    out_file.write( str( ppos[i,X] ) + " " )
                if y_col != -1 :
                    out_file.write( str( ppos[i,Y] ) + " " )
                if z_col != -1 :
                    out_file.write( str( ppos[i,Z] ) + " " )
                if vx_col != -1 :
                    out_file.write( str( pvel[i,X] ) + " " )
                if vy_col != -1 :
                    out_file.write( str( pvel[i,Y] ) + " " )
                if vz_col != -1 :
                    out_file.write( str( pvel[i,Z] ) + " " )
                out_file.write( "\n" )
            out_file.write( "\n" )
        elif out_format == "sph" :
            for i in range( N ) :
                out_file.write( chr(ptype[i]-1+ord('a')) + " " + str( ppos[i,X] ) + " " + str( ppos[i,Y] ) + " " + str( ppos[i,Z] ) + " " + str( pradius[i] ) + "\n" )
        elif out_format in [ "patch" , "ptc" ] :
            for i in range( N ) :
                out_file.write( chr(ptype[i]-1+ord('a')) + " " + str( ppos[i,X] ) + " " + str( ppos[i,Y] ) + " " + str( ppos[i,Z] ) + " " + str( pradius[i] ) + " " + str( ppatch_params[i,0] ) + " " + str( ppatch_params[i,1] ) + " " + str( prot_matrix[i,0] ) + " " + str( prot_matrix[i,1] ) + " " + str( prot_matrix[i,2] ) + " " + str( prot_matrix[i,3] ) + " " + str( prot_matrix[i,4] ) + " " + str( prot_matrix[i,5] ) + " " + str( prot_matrix[i,6] ) + " " + str( prot_matrix[i,7] ) + " " + str( prot_matrix[i,8] ) + "\n" )
        elif out_format == "zeno_input" :
            for i in range( N ) :
                if type_col != -1 :
                    out_file.write( str( ptype[i] ) )
                else :
                    out_file.write( "0" )
                if x_col != -1 :
                    out_file.write( ' ' + str( ppos[i,X] ) )
                if y_col != -1 :
                    out_file.write( ' ' + str( ppos[i,Y] ) )
                if z_col != -1 :
                    out_file.write( ' ' + str( ppos[i,Z] ) )
                out_file.write( "\n" )

    def get_valence( self , bonds ) :
        # here bonds is a list of ( bond_type , part1_index , part2_index ) , with bonds counted once
        valence = np.int32( [0] * self.N )
        for ( btype , part1 , part2 ) in bonds :
            valence[ part1 ] += 1
            valence[ part2 ] += 1
        return valence

    def get_IDX_bondslist( self , bonds ) :
        # here bonds is a list of tuples ( bond_type , part1_index , part2_index ) , with bonds counted once
        # bonds_idx is a list of indexes of bonded particles (bonds[idx1]=[idx2,idx3,...])
        bonds_idx = []
        for i in range(self.N) :
            bonds_idx.append( [] )
        for ( btype , part1 , part2 ) in bonds :
            bonds_idx[ part1 ].append( part2 )
            bonds_idx[ part2 ].append( part1 )
        return bonds_idx

    def translate_IDX_2_ID_bondslist( self , bonds_idx , kind="NOtype" ) :
        # here bonds_idx is a list of indexes of bonded particles (bonds_idx[idx1]=[idx2,idx3,...])
        # bonds_id is a dictionary containing bonded particles IDs (bonds_id[ID1]=[ID2,ID3,...])
        # bonds_id and bonds_idx are full lists !
        # if tuples with bond_types are involved, bonds are counted once
        if kind == "NOtype" :
            bonds_id = {}
            for i in range( len(bonds_idx) ) :
                bonds_id[ self.id[ i ] ] = []
                for bonded in bonds_idx[i] :
                    bonds_id[ self.id[ i ] ].append( self.id[ bonded ] )
        elif kind == "YEStype" :
            bonds_id = []
            for ( btype , part1 , part2 ) in bonds_idx :
                bonds_id.append(( btype , self.id[ part1 ] , self.id[ part2 ] ))
        else :
            print( "*** Bonds list not recognized !!" )
            exit( EXIT_FAILURE )
        return bonds_id

    def get_indexes( self ) :
        IDXs = {}
        for i in range(self.N) :
            IDXs[ self.id[i] ] = i
        return IDXs

    def translate_ID_2_IDX_bondslist( self , bonds_id , kind="NOtype" ) :
        # here bonds_id is a dictionary containing bonded particles IDs (bonds_id[ID1]=[ID2,ID3,...])
        # bonds_idx is a list of indexes of bonded particles (bonds_idx[idx1]=[idx2,idx3,...])
        # bonds_id and bonds_idx are full lists !
        # if tuples with bond_types are involved, bonds are counted once
        bonds_idx = []
        IDXs = self.get_indexes()
        if kind == "NOtype" :
            for i in range(self.N) :
                bonds_idx.append( [] )
            for part1_id in IDXs.keys() :
                if part1_id in bonds_id :
                    for part2_id in bonds_id[ part1_id ] :
                        bonds_idx[ IDXs[part1_id] ].append( IDXs[part2_id] )
        elif kind == "YEStype" :
            for ( btype , part1_id , part2_id ) in bonds_id :
                bonds_idx.append(( btype , IDXs[ part1_id ] , IDXs[ part2_id ] ))
        else :
            print( "*** Bonds list not recognized !!" )
            exit( EXIT_FAILURE )
        return bonds_idx

    def read_nico_init( self , filename ) :
        bonds = []
        try :
            data_file = open( filename , "r" )
        except :
            print( "ERROR: You have to tell the name of an existent init-file" )
            exit( EXIT_FAILURE )
        line0 = data_file.readline().strip()
        lines = data_file.readlines()
        data_file.close()
        self.N = int( lines[0].strip().split()[0] )
        valence = np.int32( [0] * self.N )
        words = lines[1].strip().split()
        self.box_inf = np.float64([ 0.0 , 0.0 , 0.0 ])
        self.box_sup = np.float64([ float( words[0] ) , float( words[1] ) , float( words[2] ) ])
        words = np.array( [ lines[i].strip().split() for i in range( 2 , self.N + 2 ) ] )
        self.id = np.arange( 1 , self.N+1 , dtype=int )
        self.mol = np.int32( [1]*self.N )
        self.q = np.float64( [0.0]*self.N )
        self.type = np.int32( [1]*self.N )
        self.pos = words[:,0:3].astype('float64')
        if len( words[0] ) > 3 :
            self.box_img = words[:,3:6].astype('int32')
        i = self.N + 2
        while i < len( lines ) :
            words = lines[i].strip().split()
            part1 = int(words[0]) - 1
            valence[ part1 ] = int(words[1])
            if len( words ) > 2 :
                self.q[part1] = float( words[2] )
                if self.q[part1] == 0 :
                    self.type[part1] = 1
                else :
                    self.type[part1] = 3
            else :
                self.type[part1] = 1
                self.q[part1] = 0.0
            i += 1
            words = lines[i].strip().split()
            for word in words :
                part2 = int(word) - 1
                if part2 > part1 :
                    bonds.append(( 1 , part1 , part2 ))
            i += 1
        return valence , bonds , line0
    # bonds contains the particles indexes

    def print_nico_init( self , out_file_name , bonds , line0 = "0 0 0" ) :
        ppos = self.unwrapped_coordinates()
        out_file = open( out_file_name , "w" )
        out_file.write( line0 + "\n" )
        out_file.write( str( self.N ) + "\n" )
        out_file.write( str( self.box_sup[X] - self.box_inf[X] ) + ' ' + str( self.box_sup[Y] - self.box_inf[Y] ) + ' ' + str( self.box_sup[Z] - self.box_inf[Z] ) + "\n" )
        ppos = ppos - self.box_inf
        for i in range( self.N ) :
            out_file.write( str( ppos[i,X] ) + '\t' + str( ppos[i,Y] ) + '\t' + str( ppos[i,Z] ) + '\n' )
        bonds_list = [[] for x in range(self.N)]
        for bond in bonds :
            bonds_list[ int(bond[1]) ].append( int(bond[2]) )
            bonds_list[ int(bond[2]) ].append( int(bond[1]) )
        for i in range( self.N ) :
            out_file.write( str( i+1 ) + '\t' + str( len(bonds_list[i]) ) + '\t' + str( self.q[i] ) + '\n' )
            if len(bonds_list[i]) > 0 :
                for bonded in bonds_list[i] :
                    out_file.write( str( bonded+1 ) + '\t' )
                out_file.write( '\n' )
        out_file.close()

    def read_lammps_init( self , filename , mode = "charge" ) :
        atom_types = 0
        bond_types = 0
        angle_types = 0
        bonds = []
        valence = np.int32( [] )
        N_bonds = 0
        N_angles = 0
        try :
            data_file = open( filename , "r" )
        except :
            print( "ERROR: You have to tell the name of an existent init-file" )
            exit( EXIT_FAILURE )
        line0 = data_file.readline().strip()
        lines = data_file.readlines()
        data_file.close()
        IDXs = {}
        i = 0
        while i < len( lines ) :
            words = lines[i].strip().split()
            if len(words) > 1 :
                if words[1] == "atoms" :
                    self.N = int(words[0])
                    valence = np.int32( [0] * self.N )
                elif words[1] == "bonds" :
                    N_bonds = int(words[0])
                elif words[1] == "angles" :
                    N_angles = int(words[0])
            if len(words) > 2 :
                if words[2] == "xlo" :
                    self.box_inf[X] = float( words[0] )
                    self.box_sup[X] = float( words[1] )
                if words[2] == "ylo" :
                    self.box_inf[Y] = float( words[0] )
                    self.box_sup[Y] = float( words[1] )
                if words[2] == "zlo" :
                    self.box_inf[Z] = float( words[0] )
                    self.box_sup[Z] = float( words[1] )
                if words[1] == "atom" and words[2] == "types" :
                    atom_types = int( words[0] )
                if words[1] == "bond" and words[2] == "types" :
                    bond_types = int( words[0] )
                if words[1] == "angle" and words[2] == "types" :
                    angle_types = int( words[0] )
            if len(words) > 0 :
                if words[0] == "Atoms" :
                    box_sides = self.box_sup - self.box_inf
                    i += 1
                    if mode == "charge" :
                        words = np.array( [ lines[j].strip().split() for j in range( i+1 , i + self.N + 1 ) ] )
                        i += self.N
                        self.id = words[:,0].astype('int32')
                        self.type = words[:,1].astype('int32')
                        self.pos = words[:,2:5].astype('float64')
                        self.q = words[:,5].astype('float64')
                        self.mol = words[:,6].astype('int32')
                        if len(words[0]) == 10 :
                            self.box_img = words[:,7:10].astype('int32')
                    elif mode == "neutral" :
                        words = np.array( [ lines[j].strip().split() for j in range( i+1 , i + self.N + 1 ) ] )
                        i += self.N
                        self.id = words[:,0].astype('int32')
                        self.mol = words[:,1].astype('int32')
                        self.type = words[:,2].astype('int32')
                        self.pos = words[:,3:6].astype('float64')
                        if len(words[0]) == 9 :
                            self.box_img = words[:,6:9].astype('int32')
                    elif mode == "atomic" :
                        words = np.array( [ lines[j].strip().split() for j in range( i+1 , i + self.N + 1 ) ] )
                        i += self.N
                        self.id = words[:,0].astype('int32')
                        self.type = words[:,1].astype('int32')
                        self.pos = words[:,2:5].astype('float64')
                        if len(words[0]) == 8 :
                            self.box_img = words[:,5:8].astype('int32')
                    IDXs = self.get_indexes()
                elif words[0] == "Velocities" :
                    i += 1
                    words = np.array( [ lines[j].strip().split() for j in range( i+1 , i + self.N + 1 ) ] )
                    i += self.N
                    ids = words[:,0].astype('int32')
                    orderedIDXs = np.int32( [ IDXs[ids[j]] for j in range( self.N ) ] )
                    self.vel = words[:,1:4].astype('float64')
                    self.vel = self.vel[ orderedIDXs ]
                elif words[0] == "Bonds" :
                    i += 1
                    for j in range( N_bonds ) :
                        i += 1
                        words = lines[i].strip().split()
                        part1 = IDXs[ int(words[2]) ]
                        valence[ part1 ] += 1
                        part2 = IDXs[ int(words[3]) ]
                        valence[ part2 ] += 1
                        bonds.append(( int(words[1]) , part1 , part2 ))
                elif words[0] == "Angles" :
                    i += 1
                    for j in range( N_angles ) :
                        i += 1
                        words = lines[i].strip().split()
                        part1 = IDXs[ int(words[2]) ]
                        part2 = IDXs[ int(words[3]) ]
                        part3 = IDXs[ int(words[4]) ]
                        self.angle.append( ANGLE( int(words[0]) , int(words[1]) , part1 , part2 , part3 ) )
            i += 1
        return valence , bonds , atom_types , bond_types , line0

    def print_lammps_init( self , bonds , atom_types , bond_types , line0 , out_file_name , mode = "charge" ) :
        # Here bonds is an array containing ( bond_type , part1_index , part2_index )
        out_file = open( out_file_name , "w" )
        out_file.write( line0 + "\n" )
        out_file.write( str( self.N ) + " atoms\n" )
        if len(bonds) > 0 :
            out_file.write( str( len(bonds) ) + " bonds\n" )
        if len( self.angle ) > 0 :
            out_file.write( str( len(self.angle) ) + " angles\n" )
        out_file.write( "\n" )
        out_file.write( str( atom_types ) + " atom types\n" )
        if bond_types > 0 :
            out_file.write( str( bond_types ) + " bond types\n" )
        if len( self.angle ) > 0 :
            angle_types = max( [ x.type for x in self.angle ] )
            out_file.write( str( angle_types ) + " angle types\n" )
        out_file.write( "\n" )
        out_file.write( str( self.box_inf[X] ) + ' ' + str( self.box_sup[X] ) + ' ' + " xlo xhi\n" )
        out_file.write( str( self.box_inf[Y] ) + ' ' + str( self.box_sup[Y] ) + ' ' + " ylo yhi\n" )
        out_file.write( str( self.box_inf[Z] ) + ' ' + str( self.box_sup[Z] ) + ' ' + " zlo zhi\n" )
        out_file.write( "\n" )
        out_file.write( "Masses\n" )
        out_file.write( "\n" )
        for i in range( 1 , atom_types+1 ) :
            out_file.write( str(i) + " 1\n" )
        out_file.write( "\n" )
        out_file.write( "Atoms\n" )
        out_file.write( "\n" )
        self.box_img = self.box_img.astype('int32')
        if mode == "charge" :
            if len( self.box_img ) > 0 :
                for i in range( self.N ) :
                    out_file.write( str( self.id[i] ) + '\t' + str( self.type[i] ) + '\t' + str( self.pos[i,X] ) + '\t' + str( self.pos[i,Y] ) + '\t' + str( self.pos[i,Z] ) + '\t' + str( self.q[i] ) + '\t' + str( self.mol[i] ) + '\t' + str( self.box_img[i,X] ) + '\t' + str( self.box_img[i,Y] ) + '\t' + str( self.box_img[i,Z] ) + '\n' )
            else :
                for i in range( self.N ) :
                    out_file.write( str( self.id[i] ) + '\t' + str( self.type[i] ) + '\t' + str( self.pos[i,X] ) + '\t' + str( self.pos[i,Y] ) + '\t' + str( self.pos[i,Z] ) + '\t' + str( self.q[i] ) + '\t' + str( self.mol[i] ) + '\n' )
        elif mode == "neutral" :
            if len( self.box_img ) > 0 :
                for i in range( self.N ) :
                    out_file.write( str( self.id[i] ) + '\t' + str( self.mol[i] ) + '\t' + str( self.type[i] ) + '\t' + str( self.pos[i,X] ) + '\t' + str( self.pos[i,Y] ) + '\t' + str( self.pos[i,Z] ) + '\t' + str( self.box_img[i,X] ) + '\t' + str( self.box_img[i,Y] ) + '\t' + str( self.box_img[i,Z] ) + '\n' )
            else :
                for i in range( self.N ) :
                    out_file.write( str( self.id[i] ) + '\t' + str( self.mol[i] ) + '\t' + str( self.type[i] ) + '\t' + str( self.pos[i,X] ) + '\t' + str( self.pos[i,Y] ) + '\t' + str( self.pos[i,Z] ) + '\n' )
        elif mode == "atomic" :
            if len( self.box_img ) > 0 :
                for i in range( self.N ) :
                    out_file.write( str( self.id[i] ) + '\t' + str( self.type[i] ) + '\t' + str( self.pos[i,X] ) + '\t' + str( self.pos[i,Y] ) + '\t' + str( self.pos[i,Z] ) + '\t' + str( self.box_img[i,X] ) + '\t' + str( self.box_img[i,Y] ) + '\t' + str( self.box_img[i,Z] ) + '\n' )
            else :
                for i in range( self.N ) :
                    out_file.write( str( self.id[i] ) + '\t' + str( self.type[i] ) + '\t' + str( self.pos[i,X] ) + '\t' + str( self.pos[i,Y] ) + '\t' + str( self.pos[i,Z] ) + '\n' )
        else :
            print( " *** ERROR: You have to specify a mode among neutral or charge or atomic !" )
            exit(1)
        out_file.write( "\n" )
        if len(self.vel) == self.N :
            out_file.write( "Velocities\n" )
            out_file.write( "\n" )
            for i in range( self.N ) :
                out_file.write( str( self.id[i] ) + '\t' + str( self.vel[i,X] ) + '\t' + str( self.vel[i,Y] ) + '\t' + str( self.vel[i,Z] ) + '\n' )
            out_file.write( "\n" )
        if len(bonds) > 0 :
            out_file.write( "Bonds\n" )
            out_file.write( "\n" )
            for i in range(len(bonds)) :
                out_file.write( str( i+1 ) + '\t' + str( bonds[i][0] ) + '\t' + str( self.id[ bonds[i][1] ] ) + '\t' + str( self.id[ bonds[i][2] ] ) + '\n' )
        if len( self.angle ) > 0 :
            out_file.write( "\nAngles\n" )
            out_file.write( "\n" )
            for i in range(len(self.angle)) :
                out_file.write( str( self.angle[i].id ) + '\t' + str( self.angle[i].type ) + '\t' + str( self.id[ self.angle[i].p1 ] ) + '\t' + str( self.id[ self.angle[i].p2 ] ) + '\t' + str( self.id[ self.angle[i].p3 ] ) + '\n' )
        out_file.close()

    def print_molsim_init( self , bonds , out_file_root ) :
        # here bonds is a list of ( bond_type , part1_index , part2_index ) , with bonds counted once
        ppos = self.unwrapped_coordinates()
        pos_file = open( out_file_root+"_positions.dat" , "w" )
        atom_props_file = open( out_file_root+"_attributes.dat" , "w" )
        valence = self.get_valence( bonds )
        for i in range( self.N ) :
            pos_file.write( str( ppos[i,X] ) + '\t' + str( ppos[i,Y] ) + '\t' + str( ppos[i,Z] ) + '\n' )
            atom_props_file.write( str( self.id[i] ) + '\t' + str( self.type[i] ) + '\t' + str( self.q[i] ) + '\t' + str( valence[i] ) + '\t' + str( self.mol[i] ) + '\n' )
        pos_file.close()
        atom_props_file.close()
        bondsNN_file = open( out_file_root+"_bondNN.dat" , "w" )
        bondsCL_file = open( out_file_root+"_bondCL.dat" , "w" )
        IDXbondslist = self.get_IDX_bondslist( bonds )
        for i in range( self.N ) :
            NNbonds = []
            CLbonds = []
            if valence[i] < 3 :
                for j in IDXbondslist[i] :
                    if valence[j] < 3 :
                        NNbonds.append( j+1 )
                    else :
                        CLbonds.append( j+1 )
            else :
                for j in IDXbondslist[i] :
                    CLbonds.append( j+1 )
            for j in range( len(NNbonds) , 2 ) :
                NNbonds.append(0)
            for j in range( len(CLbonds) , 4 ) :
                CLbonds.append(0)
            bondsNN_file.write( str( NNbonds[0] ) + '\t' + str( NNbonds[1] ) + '\n' )
            bondsCL_file.write( str( CLbonds[0] ) + '\t' + str( CLbonds[1] ) + '\t' + str( CLbonds[2] ) + '\t' + str( CLbonds[3] ) + '\n' )
            NNbonds.clear()
            CLbonds.clear()
        bondsNN_file.close()
        bondsCL_file.close()
        box_file = open( out_file_root+"_box.dat" , "w" )
        box_file.write( str( self.box_inf[X] ) + ' ' + str( self.box_sup[X] ) + "\n" )
        box_file.write( str( self.box_inf[Y] ) + ' ' + str( self.box_sup[Y] ) + "\n" )
        box_file.write( str( self.box_inf[Z] ) + ' ' + str( self.box_sup[Z] ) + "\n" )
        box_file.close()

    def discard_particles( self ) :
        self.id = np.int32( [] )
        self.type = np.int32( [] )
        self.q = np.float64( [] )
        self.mass = np.float64( [] )
        self.mol = np.int32( [] )
        self.pos = np.float64( [] )
        self.box_img = np.int32( [] )
        self.radius = np.float64( [] )
        self.N = 0
        self.vel = np.float64( [] )
        for i in range( self.Nchunks ) :
            del self.chunk[i][:]
        del self.chunk[:]
        self.Nchunks = 0
        self.molecules.clear
        del self.cells[:]
        for i in range( self.Ncontacts ) :
            del self.contacts[i][:]
        del self.contacts[:]
        self.Ncontacts = 0
        self.lofcontacts = np.float64( [] )

    def clear_configuration( self ) :
        self.time = 0.0
        self.id = np.int32( [] )
        self.type = np.int32( [] )
        self.q = np.float64( [] )
        self.mass = np.float64( [] )
        self.mol = np.int32( [] )
        self.pos = np.float64( [] )
        self.box_img = np.int32( [] )
        self.radius = np.float64( [] )
        self.N = 0
        self.vel = np.float64( [] )
        for i in range( self.Nchunks ) :
            del self.chunk[i][:]
        del self.chunk[:]
        self.Nchunks = 0
        self.COM = np.float64([ 0.0 , 0.0 , 0.0 ])
        self.VCOM = np.float64([ 0.0 , 0.0 , 0.0 ])
        self.box_inf = np.float64([ 0.0 , 0.0 , 0.0 ])
        self.box_sup = np.float64([ 0.0 , 0.0 , 0.0 ])
        del self.chain[:]
        del self.bond[:]
        del self.angle[:]
        self.molecules.clear
        del self.cells[:]
        del self.contacts[:]
        self.Ncontacts = 0
        self.lofcontacts = np.float64( [] )

    def select( self , sel_mode , sel_val ) :
        # atoms selection
        if sel_mode == "mol" :
            select_vector = np.equal( self.mol , sel_val )
        elif sel_mode == "type" :
            select_vector = np.equal( self.type , sel_val )
        elif sel_mode == "q" :
            select_vector = np.equal( self.q , sel_val )
        elif sel_mode == "xlinkers" :
            bonds_list = sel_val
            select_vector = np.array( [False]*self.N )
            for i in range( self.N ) :
                if self.id[i] in bonds_list.keys() :
                    if len( bonds_list[ self.id[i] ] ) > 2 :
                        select_vector[i] = True
        elif sel_mode == "IDfile" :
            IDS_list = sel_val
            select_vector = np.array( [False]*self.N )
            for i in range( self.N ) :
                if self.id[i] in IDS_list :
                    select_vector[i] = True
        else :
            print( "*** ERROR: Selection option is needed! Possible values: xlinkers, mol, type, q ." )
        if len( select_vector ) != self.N :
            print( "*** ERROR: Selection failed ." )
        else :
            if len( self.pos ) == self.N :
                self.pos = self.pos[ select_vector , : ]
            if len( self.vel ) == self.N :
                self.vel = self.vel[ select_vector , : ]
            if len( self.id ) == self.N :
                self.id = self.id[ select_vector ]
            if len( self.type ) == self.N :
                self.type = self.type[ select_vector ]
            if len( self.q ) == self.N :
                self.q = self.q[ select_vector ]
            if len( self.mol ) == self.N :
                self.mol = self.mol[ select_vector ]
            if len( self.box_img ) == self.N :
                self.box_img = self.box_img[ select_vector , : ]
            if len( self.radius ) == self.N :
                self.radius = self.radius[ select_vector ]
            self.N = np.count_nonzero( select_vector )

    def remove( self , rm_mode , rm_val , max_removal = 0 , radius = 0.0 , mode = "none" ) :
        ## change : CONTROL IF THE POSITIONS ARE USED AS WRAPPED OR UNWRAPPED
        # atoms selection
        if rm_mode == "mol" :
            select_vector = np.not_equal( self.mol , rm_val )
        elif rm_mode == "type" :
            select_vector = np.not_equal( self.type , rm_val )
        elif rm_mode == "q" :
            select_vector = np.not_equal( self.q , rm_val )
        elif rm_mode == "xlinkers" :
            bonds_list = rm_val
            select_vector = np.array( [True]*self.N )
            for i in range( self.N ) :
                if self.id[i] in bonds_list.keys() :
                    if len( bonds_list[ self.id[i] ] ) > 2 :
                        select_vector[i] = False
        elif rm_mode == "sphere" :
        # here unwrapped positions are used !
            if mode not in [ "in" , "out" ] :
                print( "In removing sphere mode can be in or out !" )
                exit()
            if len( self.box_img ) > 0 :
                dist = self.pos + ( self.box_img.astype('float64') * (self.box_sup-self.box_inf) )
            else :
                dist = self.pos
            dist = np.sqrt( np.add.reduce( dist**2 , 1 ) )
            if mode == "in" :
                selectBdist = np.greater( dist , radius )
            elif mode == "out" :
                selectBdist = np.less_equal( dist , radius )
            if rm_val >= 0 :
                select_vector = np.where( self.type == rm_val , selectBdist , [True]*self.N )
            else :
                select_vector = selectBdist
        else :
            print( "*** ERROR: Removing option is needed! Possible values: xlinkers, mol, type, q ." )
        if len( select_vector ) != self.N :
            print( "*** ERROR: Removing failed ." )
        else :
            if max_removal > 0 :
                i = 0
                Nrm = 0
                while i < self.N and Nrm < max_removal :
                    if select_vector[i] == False :
                        Nrm += 1
                    i += 1
                select_vector[i:self.N] = True
            if len( self.pos ) == self.N :
                self.pos = self.pos[ select_vector , : ]
            if len( self.vel ) == self.N :
                self.vel = self.vel[ select_vector , : ]
            if len( self.id ) == self.N :
                self.id = self.id[ select_vector ]
            if len( self.type ) == self.N :
                self.type = self.type[ select_vector ]
            if len( self.q ) == self.N :
                self.q = self.q[ select_vector ]
            if len( self.mol ) == self.N :
                self.mol = self.mol[ select_vector ]
            if len( self.box_img ) == self.N :
                self.box_img = self.box_img[ select_vector , : ]
            if len( self.radius ) == self.N :
                self.radius = self.radius[ select_vector ]
            self.N = np.count_nonzero( select_vector )

    def remove_by_dist( self , cutoff_dist , mode , remove_by , rm_val , reference_sel , reference_val ) :
        if mode not in [ "grtr" , "less" ] :
            print( "In removing by distance, mode can be grtr or less !" )
            exit()
        # atoms selection
        reference_atoms_idx = []
        atoms2remove_idx = []
        if reference_sel == "type" :
            ref_val = int( reference_val )
            for i in range(self.N) :
                if self.type[i] == ref_val :
                    reference_atoms_idx.append( i )
        else :
            print( "In function remove_by_dist() : reference atoms can be selected by: type" )
            exit()
        box_side = self.box_sup - self.box_inf
        if remove_by == "type" :
            type_no = int( rm_val )
            if mode == "grtr" :
                for i in range(self.N) :
                    if self.type[i] == type_no :
                        j = 0
                        flag = 0
                        while ( j < len(reference_atoms_idx) and flag == 0 ) :
                            nearest_img_ref = self.pos[ reference_atoms_idx[j] ] - ( self.periodic_image( self.pos[ reference_atoms_idx[j] ] , self.pos[i] , box_side ) * box_side )
                            r = np.sqrt( np.add.reduce( ( self.pos[i] - nearest_img_ref )**2 , 0 ) )
                            if r < cutoff_dist :
                                flag = 1
                            j += 1
                        if flag == 0 :
                            atoms2remove_idx.append( i )
            elif mode == "less" :
                for i in range(self.N) :
                    if self.type[i] == type_no :
                        j = 0
                        flag = 0
                        while j < len(reference_atoms_idx) and flag == 0 :
                            nearest_img_ref = self.pos[ reference_atoms_idx[j] ] - ( self.periodic_image( self.pos[ reference_atoms_idx[j] ] , self.pos[i] , box_side ) * box_side )
                            r = np.sqrt( np.add.reduce( ( self.pos[i] - nearest_img_ref )**2 , 0 ) )
                            if r > cutoff_dist :
                                flag = 1
                            j += 1
                        if flag == 0 :
                            atoms2remove_idx.append( i )
            print( str(len(reference_atoms_idx)) + ' ' + str(len(atoms2remove_idx)) )
            select_vector = np.array( [True]*self.N )
            select_vector[ atoms2remove_idx ] = False
            if len( self.pos ) == self.N :
                self.pos = self.pos[ select_vector , : ]
            if len( self.vel ) == self.N :
                self.vel = self.vel[ select_vector , : ]
            if len( self.id ) == self.N :
                self.id = self.id[ select_vector ]
            if len( self.type ) == self.N :
                self.type = self.type[ select_vector ]
            if len( self.q ) == self.N :
                self.q = self.q[ select_vector ]
            if len( self.mol ) == self.N :
                self.mol = self.mol[ select_vector ]
            if len( self.box_img ) == self.N :
                self.box_img = self.box_img[ select_vector , : ]
            if len( self.radius ) == self.N :
                self.radius = self.radius[ select_vector ]
            self.N = np.count_nonzero( select_vector )
        else :
            print( "In function remove_by_dist() : atoms to remove can be selected by: type" )
            exit()



    def unwrap( self , blist , box_side ) :
    # here blist contains the bonded particles INDEXES
        self.unwrapped_coordinates( True )
        check = [0] * self.N
        part2relocate = []
        for i in range( self.N ) :
            if len(blist[i]) == 0 :
                check[i] = 1
            else :
                part2relocate.append( i )
        walkers = []
        while len(part2relocate) > 0 :
            start = part2relocate[0]
            walkers.append( start )
            check[start] = 1
            del part2relocate[0]
            while len( walkers ) != 0 :
                for i in range( len(walkers) ) :
                    ip = walkers[0]
                    for j in range( len(blist[ip]) ) :
                        if check[ blist[ip][j] ] == 0 :
                            self.pos[ blist[ip][j] ] -= ( self.periodic_image( self.pos[ blist[ip][j] ] , self.pos[ip] , box_side ) * box_side )
                            walkers.append( blist[ip][j] )
                            check[ blist[ip][j] ] = 1
                            part2relocate.remove( blist[ip][j] )
                    del walkers[0]

    def unwrap_v2( self , IDblist ) :
    # here blist contains the bonded particles IDs
        self.unwrapped_coordinates( True )
        box_side = self.box_sup - self.box_inf
        part2relocate = {}
        for i in range( self.N ) :
            if self.id[i] in IDblist :
                part2relocate[ self.id[i] ] = i
        walkers = []
        while len(part2relocate) > 0 :
            walkers.append( part2relocate.popitem() )
            while len( walkers ) != 0 :
                for i in range( len(walkers) ) :
                    ip_id = walkers[0][0]
                    ip_index = walkers[0][1]
                    for j in range( len(IDblist[ip_id]) ) :
                        if IDblist[ip_id][j] in part2relocate :
                            self.pos[ part2relocate[ IDblist[ip_id][j] ] ] -= ( self.periodic_image( self.pos[ part2relocate[ IDblist[ip_id][j] ] ] , self.pos[ip_index] , box_side ) * box_side )
                            walkers.append( ( IDblist[ip_id][j] , part2relocate.pop( IDblist[ip_id][j] ) ) )
                    del walkers[0]

    def clusters( self , blist ) :
    # here blist contains the bonded particles INDEXES
        check = [0] * self.N
        part2relocate = []
        for i in range( self.N ) :
            if len(blist[i]) == 0 :
                check[i] = 1
            else :
                part2relocate.append( i )
        partsINclusters = len(part2relocate)
        walkers = []
        clusters = []
        while len(part2relocate) > 0 :
            start = part2relocate.pop(0)
            walkers.append( start )
            check[start] = 1
            clusters.append( [] )
            clusters[-1].append( start )
            while len( walkers ) != 0 :
                for i in range( len(walkers) ) :
                    ip = walkers.pop(0)
                    for j in blist[ip] :
                        if check[ j ] == 0 :
                            clusters[-1].append( j )
                            walkers.append( j )
                            check[ j ] = 1
                            part2relocate.remove( j )
        # final check
        print( " Clusters (" + str(len(clusters)) + " ; " + str(partsINclusters) + ") :" )
        for i in range(len(clusters)) :
            print( "    -> " + str(i+1) + " : " + str(len(clusters[i])) )
        # the function returns an array containing the lists of particles' index for each cluster
        return clusters

    def build_chains( self , blist ) :
    # here blist contains the bonded particles INDEXES
        loops_number = 0
        self.chain.clear()
        for i in range( self.N ) :
            if len( blist[i] ) == 1 or len( blist[i] ) > 2 :
                for j in blist[i] :
                    current_chain = ch.chain()
                    current_monomer = i
                    current_chain.add_monomer( current_monomer )
                    next_monomer = j
                    while len( blist[next_monomer] ) == 2 :
                        current_chain.add_monomer( next_monomer )
                        nextnext_monomer = -1
                        if blist[ next_monomer ][1] == current_monomer :
                            nextnext_monomer = blist[ next_monomer ][0]
                        elif blist[ next_monomer ][0] == current_monomer :
                            nextnext_monomer = blist[ next_monomer ][1]
                        else :
                            print( "ERROR occurring when building chains' structure" )
                        current_monomer = next_monomer
                        next_monomer = nextnext_monomer
                    current_chain.add_monomer( next_monomer )
                    if current_chain.monomers[ 0 ] < current_chain.monomers[ -1 ] :
                        self.chain.append( current_chain )
                    if current_chain.monomers[ 0 ] == current_chain.monomers[ -1 ] :
                        if current_chain.monomers[ 1 ] < current_chain.monomers[ -2 ] :
                            loops_number += 1
                            current_chain.LOOP = True
                            self.chain.append( current_chain )
                        if current_chain.N < 4 :
                            print( "ERROR: I detected a 2-monomers loop   " + repr(current_chain.monomers[ 0 ]) + " " + repr(current_chain.monomers[ -2 ]) + " " + repr(current_chain.monomers[ -1 ]) )
        print( "   Number of loops: " + str(loops_number) )


    def effective_charge( self , mol , method , cutoff = 2.0 ) :
        ## change : efficientare con numpy
        if len(self.mol) == 0 or len(self.q) == 0 :
            print( "*** ERROR : effective charge cannot be computed without molecule index and charge information !" )
            exit()
        box_side = self.box_sup - self.box_inf
        effective_charge = 0.0
        mol_ions = np.array( [] )
        cions = np.array( [] )
        if method == "first" or method == "third" :
            charged_beads = np.argwhere( self.q != 0.0 )[:,0]
            select_vector = np.equal( self.mol[ charged_beads ] , mol )
            mol_ions = charged_beads[ select_vector ]
            select_vector = np.logical_not( select_vector )
            cions = charged_beads[ select_vector ]
            effective_charge = np.add.reduce( self.q[ mol_ions ] , 0 )
        elif method == "second" :
            charged_beads = np.argwhere( self.q != 0.0 )[:,0]
            select_vector = np.not_equal( self.mol[ charged_beads ] , mol )
            cions = charged_beads[ select_vector ]
            mol_ions = np.argwhere( self.mol == mol )[:,0]
            effective_charge = np.add.reduce( self.q[ mol_ions ] , 0 )
        else :
            print( " *** ERROR in method CONFIGURATION().effective_charge() : computation method not recognized !" )
            exit( 1 )
        if method == "first" or method == "second" :
            for cion_id in cions :
                flag = -1
                ion_id = 0
                while ion_id < len(mol_ions) and flag == -1 :
                    ## changed
                    nearest_img_ref = self.pos[ cion_id ] - ( self.periodic_image( self.pos[ cion_id ] , self.pos[ mol_ions[ion_id] ] , box_side ) * box_side )
                    if np.sqrt( np.add.reduce( ( self.pos[ mol_ions[ion_id] ] - nearest_img_ref )**2 , 0 ) ) <= cutoff :
                        effective_charge += self.q[cion_id]
                        flag = 0
                    ion_id += 1
        elif method == "third" :
            self.com( "mol" , mol )
            for cion_id in cions :
                ## changed
                nearest_img_ref = self.pos[ cion_id ] - ( self.periodic_image( self.pos[ cion_id ] , self.COM , box_side ) * box_side )
                if np.sqrt( np.add.reduce( ( self.COM - nearest_img_ref )**2 , 0 ) ) <= cutoff :
                    effective_charge += self.q[cion_id]
        return effective_charge

    def com( self , sel_mode="all" , sel_val=0 ) :
        if sel_mode == "mol" :
            select_vector = np.equal( self.mol , int(sel_val) )
        elif sel_mode == "type" :
            select_vector = np.equal( self.type , int(sel_val) )
        elif sel_mode == "q" :
            select_vector = np.equal( self.q , float(sel_val) )
        elif sel_mode != "all" :
            print( "*** ERROR: Centre of mass can be calculated with the following modes : all, mol, type, q ." )
            exit( EXIT_FAILURE )
        ppos = self.unwrapped_coordinates()
        pmass = self.mass
        if sel_mode != "all" :
            ppos = ppos[ select_vector , : ]
            if len( self.mass ) != 0 :
                pmass = pmass[ select_vector ]
        if len(ppos) == 0 :
            print( "*** ERROR: Centre of mass cannot be calculated for a number of particles equal to 0 ." )
            exit( EXIT_FAILURE )
        if len( self.mass ) == 0 :
            self.COM = np.add.reduce( ppos , 0 ) / len(ppos)
        else :
            self.COM = np.add.reduce( ppos * np.column_stack((pmass , pmass , pmass)) , 0 ) / np.add.reduce( pmass , 0 )
        return self.COM

    def vcom( self , sel_mode = "all", sel_val = 0 ) :
        ## changed
        if sel_mode == "mol" :
            select_vector = np.equal( self.mol , int(sel_val) )
        elif sel_mode == "type" :
            select_vector = np.equal( self.type , int(sel_val) )
        elif sel_mode == "q" :
            select_vector = np.equal( self.q , float(sel_val) )
        elif sel_mode != "all" :
            print( "*** ERROR: The velocity of the centre of mass can be calculated with the following modes : all, mol, type, q ." )
            exit( EXIT_FAILURE )
        pvel = self.vel
        if sel_mode != "all" :
            pvel = pvel[ select_vector , : ]
        if len(pvel) == 0 :
            print( "*** ERROR: The velocity of the centre of mass cannot be calculated for a number of particles equal to 0 ." )
            exit( EXIT_FAILURE )
        self.VCOM = np.add.reduce( pvel , 0 ) / len(pvel)
        if len( self.vel ) != self.N :
            print( " *** WARNING: velocities are not properly loaded!" )
        return self.VCOM

    def displace( self , disp ) :
        if len( self.box_img ) == 0 :
            self.pos += disp
        else :
            self.unwrapped_coordinates( True )
            self.pos += disp
            self.compute_periodic_images()

    def wrap( self ) :
        self.box_img = np.int32( [] )
        self.compute_periodic_images()
        self.box_img = np.int32( [] )

    def reset( self , attribute , val=0 ) :
        if attribute == "id" :
            self.id[:] = np.arange(1,self.N+1)
        elif attribute == "type" :
            self.type[:] = val
        elif attribute == "x" :
            self.pos[:,X] = val
        elif attribute == "y" :
            self.pos[:,Y] = val
        elif attribute == "z" :
            self.pos[:,Z] = val
        else :
            print( " *** ERROR: Insert a valid attribute to reset!" )
            exit( EXIT_FAILURE )

    def change_attribute( self , sel_att , sel_val , change_att , change_val ) :
        if sel_att == "IDlist" :
            indexes = []
            change_val = []
            for i in range( self.N ) :
                if self.id[i] in sel_val.keys() :
                    indexes.append( i )
                    change_val.append( sel_val[ self.id[i] ] )
        elif sel_att == "mol" :
            if len( str(sel_val).split(":") ) == 2 :
                Vmin = int( sel_val.split(":")[0] )
                Vmax = int( sel_val.split(":")[1] )
                select_vector = np.greater_equal( self.mol , Vmin )
                select_vector = np.where( self.mol <= Vmax , select_vector , [False]*self.N )
            else :
                select_vector = np.equal( self.mol , int(sel_val) )
        elif sel_att == "type" :
            if len( str(sel_val).split(":") ) == 2 :
                Vmin = int( sel_val.split(":")[0] )
                Vmax = int( sel_val.split(":")[1] )
                select_vector = np.greater_equal( self.type , Vmin )
                select_vector = np.where( self.type <= Vmax , select_vector , [False]*self.N )
            else :
                select_vector = np.equal( self.type , int(sel_val) )
        elif sel_att == "q" :
            if len( str(sel_val).split(":") ) == 2 :
                Vmin = float( sel_val.split(":")[0] )
                Vmax = float( sel_val.split(":")[1] )
                select_vector = np.greater_equal( self.q , Vmin )
                select_vector = np.where( self.q <= Vmax , select_vector , [False]*self.N )
            else :
                select_vector = np.equal( self.q , float(sel_val) )
        elif sel_att == "id" :
            IDmin = int( sel_val.split(":")[0] )
            IDmax = int( sel_val.split(":")[1] )
            select_vector = np.greater_equal( self.id , IDmin )
            select_vector = np.where( self.id <= IDmax , select_vector , [False]*self.N )
        elif sel_att in [ "x" , "y" , "z" ] :
            Vmin = float( sel_val.split(":")[0] )
            Vmax = float( sel_val.split(":")[1] )
            if sel_att == "x" :
                select_vector = np.greater_equal( self.pos[:,X] , Vmin )
                select_vector = np.where( self.pos[:,X] <= Vmax , select_vector , [False]*self.N )
            elif sel_att == "y" :
                select_vector = np.greater_equal( self.pos[:,Y] , Vmin )
                select_vector = np.where( self.pos[:,Y] <= Vmax , select_vector , [False]*self.N )
            elif sel_att == "z" :
                select_vector = np.greater_equal( self.pos[:,Z] , Vmin )
                select_vector = np.where( self.pos[:,Z] <= Vmax , select_vector , [False]*self.N )
        elif sel_att == "atoms_list" :
            select_vector = np.array( [False]*self.N )
            for i in range( self.N ) :
                if self.id[i] in sel_val :
                    select_vector[i] = True
        elif sel_att == "all" :
            select_vector = [True]*self.N
        else :
            print( "*** ERROR : No valid selection mode to change attributes. Valid choices: mol, q, type, id, IDlist, atoms_list" )
        if sel_att == "IDlist" :
            if change_att == "mol" :
                if len(self.mol) == 0 :
                    self.mol = np.int32( [0]*self.N )
                self.mol[ indexes ] = change_val
            elif change_att == "q" :
                if len(self.q) == 0 :
                    self.q = np.float64( [0]*self.N )
                self.q[ indexes ] = change_val
            elif change_att == "m" :
                if len(self.mass) == 0 :
                    self.mass = np.float64( [1.0]*self.N )
                self.mass[ indexes ] = change_val
            elif change_att == "type" :
                if len(self.type) == 0 :
                    self.type = np.int32( [0]*self.N )
                self.type[ indexes ] = change_val
            elif change_att == "id" :
                if len(self.id) == 0 :
                    self.id = np.int32( [0]*self.N )
                self.id[ indexes ] = change_val
            elif change_att == "r" :
                if len(self.radius) == 0 :
                    self.radius = np.float64( [0.5]*self.N )
                self.radius[ indexes ] = change_val
            else :
                print( "*** WARNING: attribute to change not recognized!" )
        else :
            if change_att == "type" :
                if len(self.type) == 0 :
                    self.type = np.int32( [0]*self.N )
                self.type = np.where( select_vector , [change_val]*self.N , self.type )
            elif change_att == "q" :
                if len(self.q) == 0 :
                    self.q = np.float64( [0.0]*self.N )
                self.q = np.where( select_vector , [change_val]*self.N , self.q )
            elif change_att == "m" :
                if len(self.mass) == 0 :
                    self.mass = np.float64( [1.0]*self.N )
                self.mass = np.where( select_vector , [change_val]*self.N , self.mass )
            elif change_att == "mol" :
                if len(self.mol) == 0 :
                    self.mol = np.int32( [0]*self.N )
                self.mol = np.where( select_vector , [change_val]*self.N , self.mol )
            elif change_att == "r" :
                if len(self.radius) == 0 :
                    self.radius = np.float64( [0.0]*self.N )
                self.radius = np.where( select_vector , [change_val]*self.N , self.radius )
            else :
                print( "*** WARNING: attribute to change not recognized!" )

    def displace_attribute( self , sel_att , sel_val , change_att , change_val ) :
        if sel_att == "mol" :
            select_vector = np.equal( self.mol , int(sel_val) )
        elif sel_att == "type" :
            select_vector = np.equal( self.type , int(sel_val) )
        elif sel_att == "q" :
            select_vector = np.equal( self.q , float(sel_val) )
        elif sel_att == "id" :
            IDmin = int( sel_val.split(":")[0] )
            IDmax = int( sel_val.split(":")[1] )
            select_vector = np.greater_equal( self.id , IDmin )
            select_vector = np.where( self.id <= IDmax , select_vector , [False]*self.N )
        elif sel_att == "atoms_list" :
            select_vector = np.array( [False]*self.N )
            for i in range( self.N ) :
                if self.id[i] in sel_val :
                    select_vector[i] = True
        elif sel_att != "all" :
            print( "*** ERROR : No valid selection mode to change attributes. Valid choices: mol, q, type, id, IDlist, atoms_list" )
            exit( EXIT_FAILURE )
        if sel_att == "all" :
            if change_att == "mol" :
                self.mol = self.mol + change_val
            elif change_att == "q" :
                self.q = self.q + change_val
            elif change_att == "type" :
                self.type = self.type + change_val
            elif change_att == "id" :
                self.id = self.id + change_val
            elif change_att == "pos" :
                if len( self.box_img ) == 0 :
                    self.pos = self.pos + change_val
                else :
                    self.unwrapped_coordinates( True )
                    self.pos = self.pos + change_val
                    self.compute_periodic_images()
            elif change_att == "vel" :
                self.vel = self.vel + change_val
        else :
            if change_att == "mol" :
                self.mol = np.where( select_vector , self.mol+change_val , self.mol )
            elif change_att == "q" :
                self.q = np.where( select_vector , self.q+change_val , self.q )
            elif change_att == "type" :
                self.type = np.where( select_vector , self.type+change_val , self.type )
            elif change_att == "id" :
                self.id = np.where( select_vector , self.id+change_val , self.id )
            elif change_att == "pos" :
                if len( self.box_img ) == 0 :
                    self.pos = np.where( select_vector , self.pos+change_val , self.pos )
                else :
                    self.unwrapped_coordinates( True )
                    self.pos = np.where( select_vector , self.pos+change_val , self.pos )
                    self.compute_periodic_images()
            elif change_att == "vel" :
                self.vel = np.where( select_vector , self.vel+change_val , self.vel )

    def rescale_attribute( self , sel_att , sel_val , change_att , change_val ) :
        if sel_att == "mol" :
            select_vector = np.equal( self.mol , int(sel_val) )
        elif sel_att == "type" :
            select_vector = np.equal( self.type , int(sel_val) )
        elif sel_att == "q" :
            select_vector = np.equal( self.q , float(sel_val) )
        elif sel_att == "id" :
            IDmin = int( sel_val.split(":")[0] )
            IDmax = int( sel_val.split(":")[1] )
            select_vector = np.greater_equal( self.id , IDmin )
            select_vector = np.where( self.id <= IDmax , select_vector , [False]*self.N )
        elif sel_att == "atoms_list" :
            select_vector = np.array( [False]*self.N )
            for i in range( self.N ) :
                if self.id[i] in sel_val :
                    select_vector[i] = True
        elif sel_att != "all" :
            print( "*** ERROR : No valid selection mode to change attributes. Valid choices: mol, q, type, id, IDlist, atoms_list" )
            exit( EXIT_FAILURE )
        if sel_att == "all" :
            if change_att == "mol" :
                self.mol = self.mol * change_val
            elif change_att == "q" :
                self.q = self.q * change_val
            elif change_att == "type" :
                self.type = self.type * change_val
            elif change_att == "id" :
                self.id = self.id * change_val
            elif change_att == "pos" :
                if len( self.box_img ) == 0 :
                    self.pos = self.pos * change_val
                else :
                    self.unwrapped_coordinates( True )
                    self.pos = self.pos * change_val
                    self.compute_periodic_images()
            elif change_att == "vel" :
                self.vel = self.vel * change_val
        else :
            if change_att == "mol" :
                self.mol = np.where( select_vector , self.mol*change_val , self.mol )
            elif change_att == "q" :
                self.q = np.where( select_vector , self.q*change_val , self.q )
            elif change_att == "type" :
                self.type = np.where( select_vector , self.type*change_val , self.type )
            elif change_att == "id" :
                self.id = np.where( select_vector , self.id*change_val , self.id )
            elif change_att == "pos" :
                if len( self.box_img ) == 0 :
                    self.pos = np.where( select_vector , self.pos*change_val , self.pos )
                else :
                    self.unwrapped_coordinates( True )
                    self.pos = np.where( select_vector , self.pos*change_val , self.pos )
                    self.compute_periodic_images()
            elif change_att == "vel" :
                self.vel = np.where( select_vector , self.vel*change_val , self.vel )

    def gyration( self , COM , mode = "all" , sel_val = 1 ) :
        ppos = self.unwrapped_coordinates()
        Rg = 0.0
        if mode == "all" :
            ppos = ppos - COM
        elif mode == "mol" :
            ppos = ppos[ np.equal( self.mol , int(sel_val) ) , : ] - COM
        elif mode == "type" :
            ppos = ppos[ np.equal( self.type , int(sel_val) ) , : ] - COM
        Rg = np.sqrt( np.add.reduce( ppos**2 , (0,1) ) / len(ppos) )
        return Rg

    def convexHull( self , mode = "all" , sel_val = 1 ) :
        # This function computes the convex Hull of the configuration or part of it and returns
        #       the hull object, the center of masses of facets and the inertia tensor of them.
        ppos = self.unwrapped_coordinates()
        if mode == "all" :
            np_conf = ppos
        elif mode == "mol" :
            np_conf = ppos[ np.equal( self.mol , int(sel_val) ) , : ]
        elif mode == "type" :
            np_conf = ppos[ np.equal( self.type , int(sel_val) ) , : ]
        elif mode == "indexes" :
            np_conf = ppos[ sel_val ]
        # convex hull
        hull = ConvexHull( np_conf )
        vertices = hull.points
        triangles = hull.simplices
        n_triangles , nvx = triangles.shape
        # Computation of triangles center of mass
        facets_com = ( vertices[triangles[:,0]] + vertices[triangles[:,1]] + vertices[triangles[:,2]] ) / 3.0
        # center of mass
        com = np.add.reduce( facets_com , 0 ) / len(facets_com)
        # computation of the inertia tensor
        I = np.zeros( (3,3) )
        points = facets_com - com
        psquared = points**2
        I[0,0] = np.add.reduce( psquared[:,0] , 0 )
        I[1,1] = np.add.reduce( psquared[:,1] , 0 )
        I[2,2] = np.add.reduce( psquared[:,2] , 0 )
        I[0,1] = np.add.reduce( points[:,0]*points[:,1] , 0 )
        I[1,0] = I[0,1]
        I[0,2] = np.add.reduce( points[:,0]*points[:,2] , 0 )
        I[2,0] = I[0,2]
        I[1,2] = np.add.reduce( points[:,1]*points[:,2] , 0 )
        I[2,1] = I[1,2]
        I = I / len(points)
        return hull, facets_com, I

    def surface_mesh( self , mode = "all" , sel_val = -1 , probe_radius = 5.0 ) :
        # Insert a new SimulationCell object into a data collection:
        DATA = DataCollection()
#        cell = SimulationCell(pbc = (True, True, True))
        cell = SimulationCell(pbc = (False, False, False))
        cell[:,0] = (self.box_sup[X]-self.box_inf[X],0,0)
        cell[:,1] = (0,self.box_sup[Y]-self.box_inf[Y],0)
        cell[:,2] = (0,0,self.box_sup[Z]-self.box_inf[Z])
        DATA.objects.append(cell)
        # Create a Particles object containing two particles:
        particles = Particles()
        ppos = self.unwrapped_coordinates()
        if mode == "all" :
            particles.create_property( 'Position' , data=ppos )
        elif mode == "mol" :
            particles.create_property( 'Position' , data=ppos[ np.equal( self.mol , int(sel_val) ) , : ] )
        elif mode == "type" :
            particles.create_property( 'Position' , data=ppos[ np.equal( self.type , int(sel_val) ) , : ] )
        elif mode == "indexes" :
            particles.create_property( 'Position' , data=ppos[ sel_val ] )
        DATA.objects.append(particles)
        # Create a new Pipeline with a StaticSource as data source:
        pipeline = Pipeline(source = StaticSource(data = DATA))
        # Apply a modifier:
        pipeline.modifiers.append( ConstructSurfaceModifier(
            method = ConstructSurfaceModifier.Method.AlphaShape ,
            radius = probe_radius ,
            identify_regions = True ) )
        # computation of surface mesh
        analysed_data = pipeline.compute()
        surface_mesh = analysed_data.surfaces['surface']
        vertices = surface_mesh.get_vertices()
        triangles = surface_mesh.get_faces()
        # Computation of triangles center of mass
        facets_com = ( vertices[triangles[:,0]] + vertices[triangles[:,1]] + vertices[triangles[:,2]] ) / 3.0
        # center of mass
        com = np.add.reduce( facets_com , 0 ) / len(facets_com)
        # computation of the inertia tensor
        I = np.zeros( (3,3) )
        points = facets_com - com
        psquared = points**2
        I[0,0] = np.add.reduce( psquared[:,0] , 0 )
        I[1,1] = np.add.reduce( psquared[:,1] , 0 )
        I[2,2] = np.add.reduce( psquared[:,2] , 0 )
        I[0,1] = np.add.reduce( points[:,0]*points[:,1] , 0 )
        I[1,0] = I[0,1]
        I[0,2] = np.add.reduce( points[:,0]*points[:,2] , 0 )
        I[2,0] = I[0,2]
        I[1,2] = np.add.reduce( points[:,1]*points[:,2] , 0 )
        I[2,1] = I[1,2]
        I = I / len(points)
        return analysed_data , facets_com , I

    def create_molecules( self ) :
        # This is to create all the molecules present in the configuration :
        # mol_ids , indexes = np.unique( self.mol , return_inverse=True )
        # for i in range( mol_ids ) :
        #     if mol_ids[i] != 0 :
        #         self.molecules[ mol_ids[i] ] = molobj.molecule( mol_ids[i] )
        #         self.molecules[ mol_ids[i] ].atom_idxs = np.where( self.mol == mol_ids[i] )
        #         self.molecules[ mol_ids[i] ].atom_ids = self.id[  self.molecules[ mol_ids[i] ].atom_idxs  ]
        for mol_id in np.unique( self.mol ) :
            if mol_id != 0 :
                self.molecules[ mol_id ] = molobj.molecule( mol_id )
                self.molecules[ mol_id ].atom_idxs = np.argwhere( self.mol == mol_id )[:,0]
                self.molecules[ mol_id ].atom_ids = self.id[  self.molecules[ mol_id ].atom_idxs  ]

    def radial_profile( self , origin , shell_thickness , box_side , mode="all" , sel_val=-1 ) :
        spherical_prof = pf.density_profile( shell_thickness , box_side )
        ## changed
        if mode == "all" :
            ppos = self.pos
        elif mode == "mol" :
            ppos = self.pos[ np.equal( self.mol , int(sel_val) ) , : ]
        elif mode == "type" :
            ppos = self.pos[ np.equal( self.type , int(sel_val) ) , : ]
        else :
            print( "*** The mode to compute the density profile has not been recognized !!" )
            exit()
        nearest_img_ref = ppos - ( self.periodic_image( ppos , origin , box_side ) * box_side )
        dist = np.sqrt( np.add.reduce( ( nearest_img_ref - origin )**2 , 1 ) )
        shell_idx = np.floor( dist / shell_thickness ).astype('int32')
        uniques , counts = np.unique( shell_idx , return_counts=True )
        for i in range( len(uniques) ) :
            if uniques[i] < spherical_prof.Nbins :
                spherical_prof.counts[ uniques[i] ] += counts[i]
            else :
                spherical_prof.counts[ -1 ] += counts[i]
        volume = np.multiply.reduce( box_side )
        if volume == 0 :
            DIMENSION = 2
        else :
            DIMENSION = 3
        for nbin in range( spherical_prof.Nbins ) :
            if DIMENSION == 3 :
                shell_volume = 4. / 3. * math.pi * ( spherical_prof.rmax[ nbin ]**3 - spherical_prof.rmin[ nbin ]**3 )
            else :
                shell_volume = math.pi * ( spherical_prof.rmax[ nbin ]**2 - spherical_prof.rmin[ nbin ]**2 )
            spherical_prof.dens[ nbin ] = spherical_prof.counts[ nbin ] / shell_volume
        return spherical_prof

    def cylinder_profile( self , origin , box_side , cylinder_axis , radial_bin , height_bin , mode="all" , sel_val=-1 ) :
        if cylinder_axis == "x" and height_bin > box_side[X] :
            height_bin = box_side[X]
        elif cylinder_axis == "y" and height_bin > box_side[Y] :
            height_bin = box_side[Y]
        elif cylinder_axis == "z" and height_bin > box_side[Z] :
            height_bin = box_side[Z]
        if height_bin == 0.0 :
            # case of embedded 2D systems
            height_bin = 1.0
        cylindrical_prof = pf.cylindrical_profile( box_side , cylinder_axis , radial_bin , height_bin )
        if mode == "all" :
            ppos = self.pos
        elif mode == "mol" :
            ppos = self.pos[ np.equal( self.mol , int(sel_val) ) , : ]
        elif mode == "type" :
            ppos = self.pos[ np.equal( self.type , int(sel_val) ) , : ]
        else :
            print( "*** The mode to compute the density profile has not been recognized !!" )
            exit()
        nearest_img_ref = ppos - ( self.periodic_image( ppos , origin , box_side ) * box_side ) - origin
        if cylinder_axis == "x" :
            hdist = nearest_img_ref[ : , X ]
            rdist = np.sqrt( nearest_img_ref[ : , Y ]**2 + nearest_img_ref[ : , Z ]**2 )
        elif cylinder_axis == "y" :
            hdist = nearest_img_ref[ : , Y ]
            rdist = np.sqrt( nearest_img_ref[ : , Z ]**2 + nearest_img_ref[ : , X ]**2 )
        elif cylinder_axis == "z" :
            hdist = nearest_img_ref[ : , Z ]
            rdist = np.sqrt( nearest_img_ref[ : , Y ]**2 + nearest_img_ref[ : , X ]**2 )
        ## changed
        H_idx = np.floor( hdist / height_bin + 0.5 ).astype('int32')
        H_idx = np.where( np.abs(H_idx) > cylindrical_prof.lastHbin , np.sign(H_idx)*cylindrical_prof.lastHbin , H_idx )
        R_idx = np.floor( rdist / radial_bin ).astype('int32')
        R_idx = np.where( R_idx >= cylindrical_prof.Nrbins , np.array( [cylindrical_prof.Nrbins-1]*len(R_idx) ) , R_idx )
        indexes = np.column_stack(( H_idx , R_idx ))
        uniques , counts = np.unique( indexes , axis=0 , return_counts=True )
        for i in range(len(uniques)) :
            cylindrical_prof.counts[ uniques[i,0] ][ uniques[i,1] ] = counts[i]
        for hbin in range( -cylindrical_prof.lastHbin , cylindrical_prof.lastHbin+1 ) :
            for rbin in range( cylindrical_prof.Nrbins ) :
                ring_volume = math.pi * ( cylindrical_prof.rmax[ rbin ]**2 - cylindrical_prof.rmin[ rbin ]**2 ) * ( cylindrical_prof.hmax[ hbin ] - cylindrical_prof.hmin[ hbin ] )
                cylindrical_prof.dens[hbin][rbin] = cylindrical_prof.counts[hbin][rbin] / ring_volume
        return cylindrical_prof

    def linear_profile( self , box_side , profile_axis , bin_width , mode="all" , sel_val=-1 ) :
        linear_prof = pf.linear_profile( box_side , profile_axis , bin_width )
        if mode == "all" :
            ppos = self.pos
        elif mode == "mol" :
            ppos = self.pos[ np.equal( self.mol , int(sel_val) ) , : ]
        elif mode == "type" :
            ppos = self.pos[ np.equal( self.type , int(sel_val) ) , : ]
        else :
            print( "*** The mode to compute the density profile has not been recognized !!" )
            exit()
        origin = box_side * 0.5
        nearest_img_ref = ppos - ( self.periodic_image( ppos , origin , box_side ) * box_side )
        slice_area = 0.0
        box_dims = np.where( box_side > 0 , box_side , 1.0 )
        if profile_axis == "x" :
            hdist = nearest_img_ref[ : , X ]
            slice_area = box_dims[Y] * box_dims[Z]
        elif profile_axis == "y" :
            hdist = nearest_img_ref[ : , Y ]
            slice_area = box_dims[X] * box_dims[Z]
        elif profile_axis == "z" :
            hdist = nearest_img_ref[ : , Z ]
            slice_area = box_dims[X] * box_dims[Y]
        H_idx = np.floor( hdist / bin_width ).astype('int32')
        uniques , counts = np.unique( H_idx , return_counts=True )
        for i in range( len(uniques) ) :
            if uniques[i] < linear_prof.Nbins :
                linear_prof.counts[ uniques[i] ] += counts[i]
            else :
                linear_prof.counts[ -1 ] += counts[i]
        for nbin in range( linear_prof.Nbins ) :
            slice_volume = slice_area * ( linear_prof.hmax[ nbin ] - linear_prof.hmin[ nbin ] )
            linear_prof.dens[ nbin ] = linear_prof.counts[ nbin ] / slice_volume
        return linear_prof

    def print_lammps_init_oldv( self , file_name , line0 , atom_types , bond_types , bonds , mode = "charge" ) :
        out_file = open( file_name , "w" )
        out_file.write( line0 + "\n" )
        out_file.write( str(self.N) + " atoms\n" )
        out_file.write( str(len(bonds)) + " bonds\n" )
        out_file.write( "\n" )
        out_file.write( str(atom_types) + " atom types\n" )
        out_file.write( str(bond_types) + " bond types\n" )
        out_file.write( "\n" )
        out_file.write( repr(self.box_inf[X]) + ' ' + repr(self.box_sup[X]) + " xlo xhi\n" )
        out_file.write( repr(self.box_inf[Y]) + ' ' + repr(self.box_sup[Y]) + " ylo yhi\n" )
        out_file.write( repr(self.box_inf[Z]) + ' ' + repr(self.box_sup[Z]) + " zlo zhi\n" )
        out_file.write( "\n" )
        out_file.write( "Masses\n" )
        out_file.write( "\n" )
        for i in range( atom_types ) :
            out_file.write( str(i+1) + " 1\n" )
        out_file.write( "\n" )
        out_file.write( "Atoms\n" )
        out_file.write( "\n" )
        if mode == "charge" :
            # atom_ID atom_type x y z charge molecule_ID
            for i in range( self.N ) :
                out_file.write( str(self.id[i]) + '\t' + str(self.type[i]) + '\t' + repr(self.pos[i,X]) + '\t' + repr(self.pos[i,Y]) + '\t' + repr(self.pos[i,Z]) + '\t' + repr(self.q[i]) + '\t' + str(self.mol[i]) + "\n" )
        elif mode == "neutral" :
            # atom_ID molecule_ID atom_type x y z
            for i in range( self.N ) :
                out_file.write( str(self.id[i]) + '\t' + str(self.mol[i]) + '\t' + str(self.type[i]) + '\t' + repr(self.pos[i,X]) + '\t' + repr(self.pos[i,Y]) + '\t' + repr(self.pos[i,Z]) + "\n" )
        out_file.write( "\n" )
        out_file.write( "Bonds\n" )
        out_file.write( "\n" )
        for i in range( len(bonds) ) :
            # bond_ID bond_type atom1 atom2
            btype , atom1 , atom2 = bonds[i]
            out_file.write( str(i+1) + '\t' + str(btype) + "\t" + str(atom1+1) + '\t' + str(atom2+1) + '\n' )
        out_file.close()



    def select_neutral_cluster( self , seg , icheck , max_length , IDXbondslist ) :
        if self.q[ icheck ] == 0 :
            seg.append( icheck )
            if len( seg ) < max_length :
                for j in IDXbondslist[ icheck ] :
                    if j not in seg and len(seg) < max_length :
                        seg = self.select_neutral_cluster( seg , j , max_length , IDXbondslist )
        return seg

    def generate_charge_distribution( self , dist = "none" , SEL_MODE = "all" , SEL_VAL = 0 , charges_fraction = 0.0 , charges_number = 0 , charge = 0.0 , charge_type = 0 , NO_XLINKERS = True , filename = "" , surface_fraction = 0.0 , chargesXchain = 0 , radial_power = 0 , bonds = [] ) :
        available_charge_modes = [ "radial" , "pnipam_kps" , "pnipam_kps_noNN" , "mixed_rndsrf" , "chainends_rnd" , "chainends_surf_rnd" , "surf_chains" , "surf_chains_blocks" , "random" ]
        if dist not in available_charge_modes :
            print( "*** ERROR in method configuration.generate_charge_distribution() :" )
            print( "  > Available dist modes :" )
            print( available_charge_modes )
            exit()
        SELECTION = { "MODE" : SEL_MODE , "VALUE" : SEL_VAL , "IDX_LIST" : [] }
        if SELECTION["MODE"] not in [ "all" , "type" , "mol" ] :
            print( "*** ERROR in method configuration.generate_charge_distribution() :" )
            print( "  > Available selection modes :" )
            print( [ "all" , "type" , "mol" ] )
            exit()
        if dist in [ "surf_chains" , "surf_chains_blocks" , "chainends_rnd" , "chainends_surf_rnd" ] :
            if filename == "" :
                print( "*** ERROR in method configuration.generate_charge_distribution() :" )
                print( "  > Option filename needed !" )
                exit()
            if dist in [ "surf_chains" , "surf_chains_blocks" ] :
                if chargesXchain == 0 :
                    print( "*** ERROR in method configuration.generate_charge_distribution() :" )
                    print( "  > Option chargesXchain > 0 needed with dist = surf_chains !" )
                    exit()
        if dist == "mixed_rndsrf" :
            if surface_fraction == 0.0 :
                print( "*** ERROR in method configuration.generate_charge_distribution() :" )
                print( "  > Option surface_fraction > 0.0 needed with dist = mixed_rndsrf !" )
                exit()
        if dist == "radial" :
            if radial_power == 0 :
                print( "*** ERROR in method configuration.generate_charge_distribution() :" )
                print( "  > radial_power option needs to be defined as a non-zero positive number !" )
                exit()

        # selection of atoms
        if SELECTION["MODE"] == "all" :
            SELECTION["IDX_LIST"] = np.argwhere( [True]*self.N )[:,0]
        elif SELECTION["MODE"] == "mol" :
            SELECTION["IDX_LIST"] = np.argwhere( self.mol == SELECTION["VALUE"] )[:,0]
        elif SELECTION["MODE"] == "type" :
            SELECTION["IDX_LIST"] = np.argwhere( self.type == SELECTION["VALUE"] )[:,0]

        # computation of the charges' number
        if charges_number == 0 :
            charges_number = int( charges_fraction * len(SELECTION["IDX_LIST"]) )
        xlinkers = 0
        valence = self.get_valence( bonds )
        for i in SELECTION["IDX_LIST"] :
            if valence[i] > 2 :
                xlinkers += 1
        if charges_number > ( len(SELECTION["IDX_LIST"]) - xlinkers ) :
            charges_number = len(SELECTION["IDX_LIST"]) - xlinkers
        print( "   Number of charges to be generated: " , charges_number )

        # reading of the bonds
        if dist == "pnipam_kps_noNN" or dist == "mixed_rndsrf" :
            IDXbondslist = self.get_IDX_bondslist( bonds )

        # calculation of center of mass of Microgel
        self.unwrapped_coordinates( True )
        self.com( SELECTION["MODE"] , SELECTION["VALUE"] )
        print( "Center of mass of the Microgel : " , self.COM , '\r' )
        rg = self.gyration( self.COM , SELECTION["MODE"] , SELECTION["VALUE"] )
        print( "Gyration radius of the Microgel : " , rg , '\r' )
        # sys.stdout.flush()

        if len(self.q) != self.N :
            self.q = np.array([0.0]*self.N)
        if dist == "random" :
            i = 0
            while i < charges_number :
                randn = SELECTION["IDX_LIST"][ int( random.random() * len(SELECTION["IDX_LIST"]) ) ]
                if ( self.q[ randn ] == 0 ) and ( valence[ randn ] < 3 ) :
                    self.q[ randn ] = charge
                    self.type[ randn ] = charge_type
                    i += 1

        elif dist == "chainends_rnd" or dist == "chainends_surf_rnd" :
            IDXs = self.get_indexes()
            chends_file = open( filename , "r" )
            lines = chends_file.readlines()
            chends_file.close()
            chends_id = []
            for i in range( 1 , len(lines) ) :
                chends_id.append( int( lines[i].strip().split()[0] ) )
            if dist == "chainends_surf_rnd" :
                # computation of the distance among chain ends and mgel's c.o.m.
                chends_dist = {}
                for i in chends_id :
                    ## changed
                    chends_dist[i] = np.sqrt( np.add.reduce( (self.pos[ IDXs[i] ] - self.COM)**2 , 0 ) )
                chends_id_sorted = [ k for k, v in reversed( sorted(chends_dist.items(), key=lambda item: item[1]) ) ]
                chends_id.clear()
                chends_id = chends_id_sorted
            # assignment of charges
            assigned = 0
            for i in chends_id :
                if assigned < charges_number :
                    assigned += 1
                    self.q[ IDXs[i] ] = charge
                    self.type[ IDXs[i] ] = charge_type
            i = 0
            while i < charges_number-len(chends_id) :
                randn = int( random.random()*self.N )
                if ( self.q[ randn ] == 0 ) and ( valence[ randn ] < 3 ) :
                    self.q[ randn ] = charge
                    self.type[ randn ] = charge_type
                    i += 1

        elif dist == "pnipam_kps" or dist == "pnipam_kps_noNN" :
            cutoff_radius = rg
            external_monomers = 0
            # only external neutral monomers are counted
            external_part = []
            for i in SELECTION["IDX_LIST"] :
                ## changed
                distance = np.sqrt( np.add.reduce(( self.pos[i] - self.COM )**2 , 0 ) )
                if ( distance > cutoff_radius ) and ( valence[ i ] < 3 ) and ( self.q[ i ] == 0 ) :
                    external_part.append( i )
                    external_monomers += 1
            if external_monomers < charges_number :
                print( "*** ERROR : you are trying to assign more charges than possible" )
                exit( EXIT_FAILURE )
            if charges_number < external_monomers*0.5 :
                i = 0
                if dist == "pnipam_kps" :
                    while i < charges_number :
                        randn = int( random.random()*external_monomers )
                        part_idx = external_part[randn]
                        if self.q[ part_idx ] == 0 :
                            self.q[ part_idx ] = charge
                            self.type[ part_idx ] = charge_type
                            i += 1
                elif dist == "pnipam_kps_noNN" :
                    while i < charges_number :
                        randn = int( random.random()*external_monomers )
                        part_idx = external_part[randn]
                        if self.q[ part_idx ] == 0 :
                            flag = 0
                            for mon in IDXbondslist[ part_idx ] :
                                if self.q[ mon ] == charge :
                                    flag = 1
                            if flag == 0 :
                                self.q[ part_idx ] = charge
                                self.type[ part_idx ] = charge_type
                                i += 1
            else :
                if dist == "pnipam_kps_noNN" :
                    print( " *** WARNING: charged nearest neighbours cannot be avoided !" )
                no_touch = []
                for i in external_part :
                    if self.q[ i ] != 0.0 :
                        no_touch.append( i )
                if external_monomers-len(no_touch) < charges_number :
                    print( "*** ERROR in method configuration.generate_charge_distribution() :" )
                    print( "  > No sufficient uncharged monomers in the periphery of the molecule ! :" )
                    exit()
                while external_monomers-len(no_touch) != charges_number :
                    randn = int( random.random()*external_monomers )
                    part_idx = external_part[randn]
                    if self.q[ part_idx ] == 0.0 :
                        no_touch.append( part_idx )
                for i in external_part :
                    if i not in no_touch :
                        self.q[ i ] = charge
                        self.type[ i ] = charge_type

        elif dist == "mixed_rndsrf" :
            cutoff_radius = rg
            external_monomers = 0
            # only external neutral monomers are counted
            external_part = []
            for i in SELECTION["IDX_LIST"] :
                ## changed but change
                distance = np.sqrt( np.add.reduce( (self.pos[i] - self.COM)**2 , 1 ) )
                if ( distance > cutoff_radius ) and ( valence[ i ] < 3 ) and ( self.q[i] == 0 ) :
                    external_part.append( i )
                    external_monomers += 1
            if external_monomers < charges_number :
                print( "ERROR : you are trying to assign more charges than possible" )
                exit(EXIT_FAILURE)
            rnd_chrg = int( (1.0-surface_fraction) * charges_number )
            srf_chrg = charges_number - rnd_chrg
            i = 0
            while i < rnd_chrg :
                randn = SELECTION["IDX_LIST"][ int( random.random() * len(SELECTION["IDX_LIST"]) ) ]
                if ( self.q[ randn ] == 0 ) and ( valence[ randn ] < 3 ) :
                    self.q[ randn ] = charge
                    self.type[ randn ] = charge_type
                    i += 1
            i = 0
            while i < srf_chrg :
                randn = int( random.random()*external_monomers )
                part_idx = external_part[randn]
                ## changed but change
                distance = np.sqrt( np.add.reduce( (self.pos[part_idx] - self.COM)**2 , 1 ) )
                if self.q[ part_idx ] == 0 :
                    flag = 0
                    for mon in IDXbondslist[ part_idx ] :
                        if self.q[ mon ] == charge :
                            flag = 1
                    if flag == 0 :
                        self.q[ part_idx ] = charge
                        self.type[ part_idx ] = charge_type
                        i += 1

        elif dist == "surf_chains" :
            chainfile = open( filename , "r" )
            lines = chainfile.readlines()
            chainfile.close()
            chains = []
            chains_com = {}
            for i in range( 2 , len(lines) ) :
                chain = []
                words = lines[i].strip().split()
                for j in range( 2 , len(words) ) :
                    chain.append( int( words[j] ) )
                chains.append( chain )
            if charges_number > len(chains) :
                print( " *** ERROR: The number of chains is less than the number of charges to be redistributed !" )
                exit( EXIT_FAILURE )
            # computation of the distance among the c.o.m. of chains and mgel's c.o.m.
            for i in range( len(chains) ) :
                if valence[ chains[i][0] ] > 2 :
                    chains[i].pop( 0 )
                if valence[ chains[i][-1] ] > 2 :
                    chains[i].pop( -1 )
                ## changed
                chains_com[i] = ( np.add.reduce( self.pos[ chains[i] ] , 0 ) / len(chains[i]) - self.COM )
            chains_sorted = [ k for k, v in reversed( sorted(chains_com.items(), key=lambda item: item[1]) ) ]
            # assignment of charges
            charges2assign = charges_number
            ich = 0
            while charges2assign > 0 :
                chid = chains_sorted[ich]
                stopv = 0
                if valence[ chains[ chid ][0] ] == 1 :
                    self.q[ chains[ chid ][0] ] = charge
                    self.type[ chains[ chid ][0] ] = charge_type
                    stopv += 1
                    charges2assign -= 1
                elif valence[ chains[ chid ][-1] ] == 1 :
                    self.q[ chains[ chid ][-1] ] = charge
                    self.type[ chains[ chid ][-1] ] = charge_type
                    stopv += 1
                    charges2assign -= 1
                minAssignableN = int(floor( len(chains[ chid ]) / 3.0 ))
                while stopv < chargesXchain and stopv < minAssignableN :
                    randn = int( random.random() * len( chains[ chid ] ) )
                    if self.q[ chains[ chid ][randn] ] == 0 :
                        flag = 0
                        for mon in IDXbondslist[ chains[ chid ][randn] ] :
                            if self.q[ mon ] != 0 :
                                flag = 1
                        if flag == 0 :
                            self.q[ chains[ chid ][randn] ] = charge
                            self.type[ chains[ chid ][randn] ] = charge_type
                            charges2assign -= 1
                            stopv += 1
                ich += 1

        elif dist == "surf_chains_blocks" :
            chainfile = open( filename , "r" )
            lines = chainfile.readlines()
            chainfile.close()
            chains = []
            chains_com = {}
            for i in range( 2 , len(lines) ) :
                chain = []
                words = lines[i].strip().split()
                for j in range( 2 , len(words) ) :
                    chain.append( int( words[j] ) )
                chains.append( chain )
            if charges_number > len(chains)*chargesXchain :
                print( " *** ERROR: The number of chains is less than the number of charges to be redistributed !" )
                exit( EXIT_FAILURE )
            # computation of the distance among the c.o.m. of chains and mgel's c.o.m.
            for i in range( len(chains) ) :
                if valence[ chains[i][0] ] > 2 :
                    chains[i].pop( 0 )
                if valence[ chains[i][-1] ] > 2 :
                    chains[i].pop( -1 )
                ## changed
                chains_com[i] = ( np.add.reduce( self.pos[ chains[i] ] , 0 ) / len(chains[i]) - self.COM )
            chains_sorted = [ k for k, v in reversed( sorted(chains_com.items(), key=lambda item: item[1]) ) ]
            # assignment of charges
            charges2assign = charges_number
            ich = 0
            while charges2assign > 0 :
                chid = chains_sorted[ich]
                if len(chains[ chid ]) < chargesXchain :
                    minAssignableN = len(chains[ chid ])
                else :
                    minAssignableN = chargesXchain
                if valence[ chains[ chid ][0] ] == 1 :
                    self.q[ chains[ chid ][0] ] = charge
                    self.type[ chains[ chid ][0] ] = charge_type
                    charges2assign -= 1
                    minAssignableN -= 1
                elif valence[ chains[ chid ][-1] ] == 1 :
                    self.q[ chains[ chid ][-1] ] = charge
                    self.type[ chains[ chid ][-1] ] = charge_type
                    charges2assign -= 1
                    minAssignableN -= 1
                # selection of the segment
                segment = []
                while len(segment) < minAssignableN :
                    imon = int( random.random() * len( chains[ chid ] ) )
                    segment = self.select_neutral_cluster( segment , imon , minAssignableN , IDXbondslist )
                for imon in segment :
                    if self.q[ chains[ chid ][imon] ] == 0 :
                        self.q[ chains[ chid ][imon] ] = charge
                        self.type[ chains[ chid ][imon] ] = charge_type
                        charges2assign -= 1
                    else :
                        print( "*** ERROR in method configuration.generate_charge_distribution(), charge distribution: surf_chains_blocks" )
                        exit()
                ich += 1

        elif dist == "radial" :
            # creation of the distribution
            ## changed
            distance = np.sqrt( np.add.reduce( (self.pos - self.COM)**2 , 1 ) )
            max_distance = np.amax( distance )
            Nbins = int( sqrt( len(SELECTION["IDX_LIST"]) ) )
            delta_r = max_distance / ( Nbins-1. )
            parts_in_shell = []
            radial_distribution = []
            for i in range( Nbins ) :
                parts_in_shell.append( 0 )
                lista = []
                radial_distribution.append( lista )
            for i in SELECTION["IDX_LIST"] :
                ## changed but change
                distance = np.sqrt( np.add.reduce( (self.pos[i] - self.COM)**2 , 1 ) )
                index = 0
                while (index+0.5)*delta_r < distance :
                    index += 1
                parts_in_shell[ index ] += 1
                radial_distribution[ index ].append( i )
            # assignment of charges
            max_randn = pow( max_distance , radial_power + 3. )
            i = 0
            while i < charges_number :
                rand_radius = random.random() * max_randn
                rand_radius = pow( rand_radius , 1./(radial_power+3.) )
                index = 0
                while (index+0.5)*delta_r < rand_radius :
                    index += 1
                rand_part = int( random.random()*parts_in_shell[index] )
                part_idx = radial_distribution[index][rand_part]
                if ( self.q[ part_idx ] == 0 ) and ( valence[ part_idx ] < 3 ) :
                  self.q[ part_idx ] = charge
                  self.type[ part_idx ] = charge_type
                  i += 1


    def create_cell_list_2D( self , max_cutoff ) :
        """
        Compute cell list.
        """
        box_side = self.box_sup - self.box_inf
        box_side = box_side[0:2]
        ncells = ( box_side / max_cutoff ).astype(int)
        lcells = box_side / ncells
        cells = [[[] for j in range(ncells[1])] for i in range(ncells[0])]

        wpos = self.pos[:,0:2] - ( self.periodic_image( self.pos[:,0:2] , ((self.box_sup + self.box_inf) * 0.5)[0:2] , box_side ) * box_side )
        cx = (wpos[:,0] / lcells[0] ).astype(int)
        cy = (wpos[:,1] / lcells[1] ).astype(int)
        cx[ cx == ncells[0] ] = 0
        cy[ cy == ncells[1] ] = 0
        for idxs in np.column_stack(( np.arange(self.N) , cx , cy )) :
            cells[idxs[1]][idxs[2]].append( idxs[0] )
        self.cells = cells
        return self.cells

