from .action import *
import numpy as np
import time
import itertools  ## definisce product per combinazioni liste

class MESH :
    def __init__( self ) :
        self.vertices = np.array([])
        self.faces = np.array([])
        self.normals = np.array([])

    def get_face_normals( self ) :
        self.normals = np.cross( self.vertices[ self.faces[:,1] ] - self.vertices[ self.faces[:,0] ] , self.vertices[ self.faces[:,2] ] - self.vertices[ self.faces[:,0] ] )
        norms = np.sqrt( np.add.reduce( self.normals**2 , 1 ) )
        self.normals = self.normals / np.column_stack(( norms , norms , norms ))
        return self.normals

def NN_periodic_image( points , ref , box_sides ) :
    periodic_images = np.sign(points-ref) * (  ( np.fabs(points-ref) // ( box_sides*0.5 ) + 1 ) // 2  )
    return periodic_images

def nearest_contour_point( point , triangle ) :
    EDGE = 1
    VERTEX = 2
    # here triangle is ONE face, no more
    verts_pos = triangle-point
    verts_dist = np.sqrt(np.add.reduce( (triangle-point)**2 , 1 ))
    v2 = -1
    nearest_point = 0
    nearest_pos = np.array( [ 0.0 , 0.0 , 0.0 ] )
    nearest_dist = 0.0
    v1 = np.argmin( verts_dist , axis=0 )
    for i in range(3) :
        if i != v1 :
            t = -1.0 * np.add.reduce( verts_pos[v1] * (verts_pos[i]-verts_pos[v1]) , 0 ) / np.add.reduce( (verts_pos[i]-verts_pos[v1])**2 , 0 )
            if nearest_point == 0 :
                if t <= 0.0 :
                    nearest_pos = verts_pos[v1]
                    nearest_dist = verts_dist[v1]
                    nearest_point = VERTEX
                else :
                    nearest_pos = verts_pos[v1] + t * (verts_pos[i]-verts_pos[v1])
                    nearest_dist = sqrt( np.add.reduce( nearest_pos**2 , 0 ) )
                    nearest_point = EDGE
                    v2 = i
            else :
                if t <= 0.0 :
                    pos = verts_pos[v1]
                    dist = verts_dist[v1]
                    if dist < nearest_dist :
                        nearest_pos = pos
                        nearest_dist = dist
                        nearest_point = VERTEX
                        v2 = i
                else :
                    pos = verts_pos[v1] + t * (verts_pos[i]-verts_pos[v1])
                    dist = sqrt( np.add.reduce( nearest_pos**2 , 0 ) )
                    if dist < nearest_dist :
                        nearest_pos = pos
                        nearest_dist = dist
                        nearest_point = EDGE
                        v2 = i
    return np.array( [ nearest_point , v1 , v2 ] ) , nearest_pos , nearest_dist

def point_in_triangle( point , v0 , v1 , v2 ) :
    # only if point is on the same plane of vi, i=1,2,3
    # controls are disabled to speed up the computation
    # a = [ v1[0]-v0[0] , v1[1]-v0[1] , v1[2]-v0[2] ]
    # b = [ v2[0]-v1[0] , v2[1]-v1[1] , v2[2]-v1[2] ]
    # c = [ v0[0]-v2[0] , v0[1]-v2[1] , v0[2]-v2[2] ]
    a = v1 - v0
    b = v2 - v1
    c = v0 - v2
    # p_v0 = [ point[0]-v0[0] , point[1]-v0[1] , point[2]-v0[2] ]
    # p_v1 = [ point[0]-v1[0] , point[1]-v1[1] , point[2]-v1[2] ]
    # p_v2 = [ point[0]-v2[0] , point[1]-v2[1] , point[2]-v2[2] ]
    p_v0 = point - v0
    p_v1 = point - v1
    p_v2 = point - v2
    # co-planarity condition
    # if a[0]*( b[1]*p_v0[2]-b[2]*p_v0[1] ) + a[1]*( b[2]*p_v0[0]-b[0]*p_v0[2] ) + a[2]*( b[0]*p_v0[1]-b[1]*p_v0[0] ) != 0.0 :
    triangle_ortogonal = np.cross( a , b )  # only the sign matters, for b,c and c,a it is the same
    D0 = np.add.reduce( np.cross( a , p_v0 ) * triangle_ortogonal , 1 )
    D1 = np.add.reduce( np.cross( b , p_v1 ) * triangle_ortogonal , 1 )
    D2 = np.add.reduce( np.cross( c , p_v2 ) * triangle_ortogonal , 1 )
    return np.column_stack(( D0 , D1 , D2 ))

def dump_mesh( mesh , name ) :
    outfile = open( name+".v" , "w" )
    for v in mesh.vertices :
        outfile.write( str(v[0])+" "+str(v[1])+" "+str(v[2])+"\n" )
    outfile.close()
    outfile = open( name+".f" , "w" )
    for f in mesh.faces :
        outfile.write( str(f[0])+" "+str(f[1])+" "+str(f[2])+"\n" )
    outfile.close()

def dump_points( atoms , idxs , name ) :
    outfile = open( name , "w" )
    for i in idxs :
        outfile.write( str(atoms[i,0])+" "+str(atoms[i,1])+" "+str(atoms[i,2])+"\n" )
    outfile.close()

def dump_vector( vector , name ) :
    outfile = open( name , "w" )
    for i in range(len(vector)) :
        outfile.write( str(vector[i,0])+" "+str(vector[i,1])+" "+str(vector[i,2])+"\n" )
    outfile.close()

class SHAPE_OVERLAP( ACTION ) :
    def __init__( self , argv , i ) :
        ACTION.__init__( self )
        if len(argv) == 0 :
            self.print_help()
        else :
            self.operation = "shape_overlap"
            self.mode = argv[i+1].strip()
            self.value = []
            self.method = "chull"
            self.interpenetration_style = "overlap_vol"
            self.core_radius = 0.0
            self.probing_radius = 5.0
            self.local_cutoff = 4.5
            self.density_cutoff = 0.0
            self.DUMP_MESHES = False
            if self.mode != "mol" :
                print( "*** Shape - overlap parameters can be computed only with mol mode !" )
                exit()
            for interv in argv[i+2].strip().split(",") :
                molids = interv.split(":")
                if len( molids ) == 1 :
                    self.value.append( int(molids[0]) )
                elif len( molids ) == 2 :
                    for l in range( int(molids[0]) , int(molids[1])+1 ) :
                        self.value.append( l )
                else :
                    print( "*** molecules are not correctly indicated !" )
                    exit()
            filename = "shape_overlap.dat"
            j=3
            for dummy in range(5) :
                if len( argv ) > i+j :
                    if argv[i+j].strip() == "file" :
                        j+=1
                        filename = argv[i+j].strip()
                    elif argv[i+j].strip() == "method" :
                        j+=1
                        self.method = argv[i+j].strip()
                        if self.method == "smesh" :
                            j+=1
                            self.probing_radius = float( argv[i+j].strip() )
                    elif argv[i+j].strip() == "interpenetration" :
                        j+=1
                        self.interpenetration_style = argv[i+j].strip()
                        if self.interpenetration_style not in [ "overlap_vol" , "atoms_fraction" , "core_radius" , "density_cutoff" ] :
                            print( "*** ERROR : please, specify a valid way to compute interpenetrated atoms ! " )
                        if self.interpenetration_style == "core_radius" :
                            j+=1
                            self.core_radius = float( argv[i+j].strip() )
                        elif self.interpenetration_style == "density_cutoff" :
                            j+=1
                            self.density_cutoff = float( argv[i+j].strip() )
                    elif argv[i+j].strip() == "local_cut" :
                        j+=1
                        self.local_cutoff = float( argv[i+j].strip() )
                    elif argv[i+j].strip() == "dump_meshes" :
                        self.DUMP_MESHES = True
                    j+=1
            self.outfile = open( filename , "w" )
            self.outfile.write( "# timestep <single_mgel_surface> <S_p=6*pi^0.5*V/S^1.5> <single_mgel_volume> system_overlap_volume(tot_Vov/tot_Vmgels) <interpenetrated_atoms/mgel_atoms> total_mgel-meshes_volume\n" )
            self.outfile.flush()
            self.molfiles = {}
            for l in self.value :
                self.molfiles[l] = open( filename.rsplit(".",1)[0]+"_mol"+str(l)+"."+filename.rsplit(".",1)[1] , "w" )
                self.molfiles[l].write( "# timestep surface S_p=6*pi^0.5*V/S^1.5 interpenetrated_atoms/mgel_atoms volume\n" )
                self.molfiles[l].flush()
            os.system( "mkdir log_files" )

    def point_is_inside_convex_mesh( self , point , box_sides , face_normals , triangles , NNcalc = False ) :
        is_inside = True
        FACE = 0
        EDGE = 1
        VERTEX = 2
        if NNcalc :
            periodic_images = NN_periodic_image( face_centers , point , box_sides )
            triangles = triangles - (periodic_images*box_sides)
            # the column indexes of triangles are 1: face idx , 2: vertex idx (0,1,2) , 3: xyz components (0,1,2) of the vertex
            # t0 = np.column_stack(( vertices[faces[:,0],0]-periodic_images[:,0]*box_sides[0] , vertices[faces[:,1],0]-periodic_images[:,0]*box_sides[0] , vertices[faces[:,2],0]-periodic_images[:,0]*box_sides[0] ))
            # t1 = np.column_stack(( vertices[faces[:,0],1]-periodic_images[:,1]*box_sides[1] , vertices[faces[:,1],1]-periodic_images[:,1]*box_sides[1] , vertices[faces[:,2],1]-periodic_images[:,1]*box_sides[1] ))
            # t2 = np.column_stack(( vertices[faces[:,0],2]-periodic_images[:,2]*box_sides[2] , vertices[faces[:,1],2]-periodic_images[:,2]*box_sides[2] , vertices[faces[:,2],2]-periodic_images[:,2]*box_sides[2] ))
            # triangles = np.dstack(( t0 , t1 , t2 ))
        # the displacement from point to the plane of the face is the projection of the position of one of the 3 vertexes on the face normal
        projections = np.add.reduce( ( triangles[:,0] - point ) * face_normals , 1 )
        point2faces = face_normals * np.column_stack(( projections , projections , projections ))
        point2faces_dist = np.sqrt(np.add.reduce( point2faces**2 , 1 ))
        verify_vector = point_in_triangle( point + point2faces , triangles[:,0] , triangles[:,1] , triangles[:,2] )
        # NpointLocation is a vector containing information about the nearest point of every triangles (if it is inside, on the border, a vertex, and eventually its distance vector)
        NpointLocation = np.empty( [ len(face_normals) , 3 ] )
        for i in range( len(face_normals) ) :
            if not ( verify_vector[i,0] > 0 and verify_vector[i,1] > 0 and verify_vector[i,2] > 0 ) :
                #### triangles deve essere uguale a prima ! controlla nella funzione nearest contour
                NpointLocation[i] , point2faces[i] , point2faces_dist[i] = nearest_contour_point( point , triangles[i] )
            else :
                NpointLocation[i,0] = FACE
        NNface = np.argmin( point2faces_dist , axis=0 )
        if NpointLocation[NNface,0] != FACE :
            is_inside = False
        elif np.add.reduce( point2faces[NNface] * face_normals[NNface] , 0 ) < 0 :
            is_inside = False
        else :
            is_inside = True
        return is_inside
        # mesh.get_face_adjacent_faces(fi)

    def mesh_from_AlphaSurface( self , surface_mesh ) :
        mesh = MESH()
        mesh.vertices = surface_mesh.get_vertices()
        mesh.faces = surface_mesh.get_faces()
        return mesh

    def mesh_from_hull( self , hull ) :
        meshpoints = np.empty( [ len(hull.vertices) , 3] )
        new_idxs = {}
        for v in range( len(hull.vertices) ) :
            new_idxs[ hull.vertices[v] ] = v
            meshpoints[ v,0 ] = hull.points[ hull.vertices[v] , 0 ]
            meshpoints[ v,1 ] = hull.points[ hull.vertices[v] , 1 ]
            meshpoints[ v,2 ] = hull.points[ hull.vertices[v] , 2 ]
        meshfaces = np.empty( [ len(hull.simplices) , 3] )
        for v in range( len(hull.simplices) ) :
            meshfaces[ v,0 ] = new_idxs[ hull.simplices[ v , 0 ] ]
            meshfaces[ v,1 ] = new_idxs[ hull.simplices[ v , 1 ] ]
            meshfaces[ v,2 ] = new_idxs[ hull.simplices[ v , 2 ] ]
        mesh = MESH()
        mesh.vertices = meshpoints
        mesh.faces = meshfaces.astype('int')
        return mesh

    def compute_molecule_meshes( self , equil_config , box_sides , logfile ) :
        min_box_side = min( box_sides )
        for l in self.value :
            if self.method == "chull" :
                if self.interpenetration_style == "overlap_vol" :
                    equil_config.molecules[l].properties["hull"] , equil_config.molecules[l].properties["facets_com"] , equil_config.molecules[l].properties["inertia_tensor"] = equil_config.convexHull( "indexes" , equil_config.molecules[l].atom_idxs )
                elif self.interpenetration_style in [ "atoms_fraction" , "core_radius" , "density_cutoff" ] :
                    equil_config.molecules[l].properties["hull"] , equil_config.molecules[l].properties["facets_com"] , equil_config.molecules[l].properties["inertia_tensor"] = equil_config.convexHull( "indexes" , equil_config.molecules[l].properties["bulk_atom_idxs"] )
                equil_config.molecules[l].properties["mesh"] = self.mesh_from_hull( equil_config.molecules[l].properties["hull"] )
                # equil_config.molecules[l].properties["mesh"].enable_connectivity()
                # equil_config.molecules[l].properties["mesh"].add_attribute("face_normal")
                face_normals = equil_config.molecules[l].properties["mesh"].get_face_normals()
                vertices = equil_config.molecules[l].properties["mesh"].vertices
                faces = equil_config.molecules[l].properties["mesh"].faces
                com = np.add.reduce( vertices , 0 ) / len( vertices )
                face_centers = ( vertices[ faces[:,0] ] + vertices[ faces[:,1] ] + vertices[ faces[:,2] ] ) / 3.
                normals_orientation = np.sign( np.add.reduce( ( face_centers - com ) * face_normals , 1 ) )
                face_normals = face_normals * np.column_stack(( normals_orientation , normals_orientation , normals_orientation ))
                equil_config.molecules[l].properties["mesh_face_normals"] = face_normals
            elif self.method == "smesh" :
                if self.interpenetration_style == "overlap_vol" :
                    equil_config.molecules[l].properties["AlphaSurface"] , equil_config.molecules[l].properties["facets_com"] , equil_config.molecules[l].properties["inertia_tensor"] = equil_config.surface_mesh( "indexes" , equil_config.molecules[l].atom_idxs , self.probing_radius )
                elif self.interpenetration_style in [ "atoms_fraction" , "core_radius" , "density_cutoff" ] :
                    equil_config.molecules[l].properties["AlphaSurface"] , equil_config.molecules[l].properties["facets_com"] , equil_config.molecules[l].properties["inertia_tensor"] = equil_config.surface_mesh( "indexes" , equil_config.molecules[l].properties["bulk_atom_idxs"] , self.probing_radius )
                surface_mesh = equil_config.molecules[l].properties["AlphaSurface"].surfaces['surface']
                if surface_mesh.regions.count > 2 :
                    self.lock.acquire()
                    logfile.write( "*** WARNING : the surface mesh could enclose empty regions !\n" )
                    logfile.write( "    molecule : " + str(l) + ", regions : " + str(surface_mesh.regions.count) + "\n" )
                    logfile.write( "  Filled  : " + str(surface_mesh.regions['Filled']) + "\n" )
                    logfile.write( "  Volumes : " + str(surface_mesh.regions['Volume']) + "\n" )
                    logfile.write( "  Areas   : " + str(surface_mesh.regions['Surface Area']) + "\n" )
                    logfile.flush()
                    self.lock.release()
                vertices = surface_mesh.get_vertices()
                faces = surface_mesh.get_faces()
                com = np.add.reduce( vertices , 0 ) / surface_mesh.vertices.count
            else :
                self.lock.acquire()
                print( "*** ERROR : method in shape_overlap analysis is not valid !" )
                self.lock.release()
                exit()
            # computation of the max distance among a vertex and the center of mass of the mesh
            vtx_com_dist = np.sqrt( np.add.reduce( ( vertices - com )**2 , 1 ) )
            farthest_com_vtx = np.argmax( vtx_com_dist , 0 )
            max_com_vtx_dist = vtx_com_dist[ farthest_com_vtx ]
            if max_com_vtx_dist > 0.5 * min_box_side :
                self.lock.acquire()
                logfile.write( "WARNING: The radius of the sphere enclosing the mesh exceeds the half side of the box !!!  " + str(max_com_vtx_dist) + "\n" )
                logfile.write( "         It is possible that periodic images could not be computed properly" + "\n" )
                logfile.flush()
                self.lock.release()
            # triangles contains the vertices of all faces
            # the column indexes of triangles are 1: face idx , 2: vertex idx (0,1,2) , 3: xyz components (0,1,2) of the vertex
            t0 = np.column_stack(( vertices[faces[:,0],0] , vertices[faces[:,1],0] , vertices[faces[:,2],0] ))
            t1 = np.column_stack(( vertices[faces[:,0],1] , vertices[faces[:,1],1] , vertices[faces[:,2],1] ))
            t2 = np.column_stack(( vertices[faces[:,0],2] , vertices[faces[:,1],2] , vertices[faces[:,2],2] ))
            triangles = np.dstack(( t0 , t1 , t2 ))
            equil_config.molecules[l].properties["mesh_triangles"] = triangles
            equil_config.molecules[l].properties["mesh-max_com_vtx_dist"] = max_com_vtx_dist
            equil_config.molecules[l].properties["mesh_com"] = com

    def execute( self , equil_config , col ) :
        if self.mode == "mol" :
            if col.mol == -1 :
                self.lock.acquire()
                print( " *** ERROR: I do not have information about molecules ID !" )
                self.lock.release()
                exit()
        else :
            self.lock.acquire()
            print( " *** mode not recognized !" )
            self.lock.release()
            exit()
        equil_config.unwrapped_coordinates( True )
        if len( equil_config.molecules ) == 0 :
            equil_config.create_molecules()
            for mol in equil_config.molecules.keys() :
                equil_config.molecules[mol].PERMANENT = 1
        # logfile = open( "log_files/log."+str(equil_config.time) , "w", buffering=0 ) # buffering off is allowed only with binary files
        logfile = open( "log_files/log."+str(equil_config.time) , "w" )
        box_sides = equil_config.box_sup - equil_config.box_inf
        min_box_side = min( box_sides )
        # computation of the molecules' surface meshes
        # calculation of the bulk atoms
        if self.interpenetration_style == "atoms_fraction" :
            self.compute_interpenetrated_atoms_LocalCut( equil_config , box_sides , equil_config.box_inf , logfile )
            self.compute_molecule_meshes( equil_config , box_sides , logfile )
        elif self.interpenetration_style == "core_radius" :
            self.compute_interpenetrated_atoms_OutCore( equil_config , box_sides , equil_config.box_inf , logfile )
            self.compute_molecule_meshes( equil_config , box_sides , logfile )
        elif self.interpenetration_style == "density_cutoff" :
            self.compute_interpenetrated_atoms_DensityCutoff( equil_config , box_sides , equil_config.box_inf , logfile )
        elif self.interpenetration_style == "overlap_vol" :
            self.compute_molecule_meshes( equil_config , box_sides , logfile )
            self.compute_interpenetrated_atoms_MeshVolumes( equil_config , box_sides , equil_config.box_inf , logfile )
        # calculation of the interpenetration and deformation
        N_probing_points = int( box_sides[0] * box_sides[1] * box_sides[2] )
        probing_points = np.random.rand( N_probing_points , 3 )
        probing_points = probing_points * box_sides + equil_config.box_inf
        internal_probing_points = np.zeros( [ N_probing_points , 1 ] )
        probing_idxs = np.arange(N_probing_points)
        avg_vol = 0.0
        avg_surf = 0.0
        avg_Sp = 0.0
        avg_interp_monom_frac = 0

        # stochastic computation of the meshes and intersections volume
        for l in self.value :
            equil_config.molecules[l].properties["inner_points"] = []
            periodic_com_images = NN_periodic_image( equil_config.molecules[l].properties["mesh_com"] , probing_points , box_sides )
            nearest_NNpoints = probing_points + periodic_com_images*box_sides
            points_com_dist = np.sqrt( np.add.reduce( ( equil_config.molecules[l].properties["mesh_com"] - nearest_NNpoints )**2 , 1 ) )
            select_vector = np.less_equal( points_com_dist , equil_config.molecules[l].properties["mesh-max_com_vtx_dist"] )
            nearest_NNpoints = nearest_NNpoints[ select_vector , : ]
            nearest_NNpoints_idxs = probing_idxs[ select_vector ]
            Npts2analyse = len( nearest_NNpoints )
            tv0 = time.time()
            if self.method == "chull" :
                for i in range(Npts2analyse) :
                    if i % 5000 == 0 :
                        self.lock.acquire()
                        logfile.write( str(i) + " of " + str(Npts2analyse) + " points analysed, after " + str(time.time()-tv0) + "\n" )
                        logfile.flush()
                        self.lock.release()
                    if self.point_is_inside_convex_mesh( nearest_NNpoints[i] , box_sides , equil_config.molecules[l].properties["mesh_face_normals"] , equil_config.molecules[l].properties["mesh_triangles"] , NNcalc=False ) :
                        equil_config.molecules[l].properties["inner_points"].append( nearest_NNpoints_idxs[i] )
                        internal_probing_points[ nearest_NNpoints_idxs[i] , 0 ] += 1
                volmesh = equil_config.molecules[l].properties["hull"].volume
                areamesh = equil_config.molecules[l].properties["hull"].area
            elif self.method == "smesh" :
                inner_regions = []
                surface_mesh = equil_config.molecules[l].properties["AlphaSurface"].surfaces['surface']
                for r in range(len( surface_mesh.regions['Filled'] )) :
                    if surface_mesh.regions['Filled'][r] != 0 :
                        inner_regions.append( r )
                for i in range(Npts2analyse) :
                    if i % 5000 == 0 :
                        self.lock.acquire()
                        logfile.write( str(i) + " of " + str(Npts2analyse) + " points analysed, after " + str(time.time()-tv0) + "\n" )
                        logfile.flush()
                        self.lock.release()
                    if surface_mesh.locate_point( nearest_NNpoints[i] , eps=1e-6 ) in inner_regions :
                        equil_config.molecules[l].properties["inner_points"].append( nearest_NNpoints_idxs[i] )
                        internal_probing_points[ nearest_NNpoints_idxs[i] , 0 ] += 1
                volmesh = 0.0
                areamesh = 0.0
                for r in inner_regions :
                    volmesh += surface_mesh.regions['Volume'][r]
                    areamesh += surface_mesh.regions['Surface Area'][r]
            equil_config.molecules[l].properties["S_p"] = 6 * sqrt(np.pi) * volmesh / areamesh**1.5
            equil_config.molecules[l].properties["interpenetrated_frac"] = len(equil_config.molecules[l].properties["interpenetrated_atom_idxs"]) / len(equil_config.molecules[l].atom_idxs)
            avg_Sp += equil_config.molecules[l].properties["S_p"]
            avg_interp_monom_frac += equil_config.molecules[l].properties["interpenetrated_frac"]
            avg_vol += volmesh
            avg_surf += areamesh
            volstat = float(len( equil_config.molecules[l].properties["inner_points"] )) / N_probing_points * box_sides[0] * box_sides[1] * box_sides[2]
            self.lock.acquire()
            logfile.write( "mol : " + str(l) + "\tvol stat : " + str(volstat) + "\tvol mesh : " + str(volmesh) + "\n" )
            logfile.flush()
            self.molfiles[l].write( str( equil_config.time ) + " " + repr( areamesh ) + " " + repr( equil_config.molecules[l].properties["S_p"] ) + " " + repr( equil_config.molecules[l].properties["interpenetrated_frac"] ) + " " + repr( volmesh ) + "\n" )
            if self.DUMP_MESHES :
                if self.method == "smesh" :
                    dump_mesh( self.mesh_from_AlphaSurface( equil_config.molecules[l].properties["AlphaSurface"].surfaces['surface'] ) , "mesh_"+str(l) )
                elif "mesh" in equil_config.molecules[l].properties.keys() :
                    dump_mesh( equil_config.molecules[l].properties["mesh"] , "mesh_"+str(l) )
                # dump_vector( nearest_NNpoints , "nearestp_"+str(l)+".dat" )
                dump_points( equil_config.pos , equil_config.molecules[l].atom_idxs , "atoms_"+str(l)+".dat" )
                if len(equil_config.molecules[l].properties["interpenetrated_atom_idxs"]) > 0 :
                    dump_points( equil_config.pos , equil_config.molecules[l].properties["interpenetrated_atom_idxs"] , "int_atoms_"+str(l)+".dat" )
            self.lock.release()
        avg_surf /= len( self.value )
        avg_Sp /= len( self.value )
        avg_vol /= len( self.value )
        avg_interp_monom_frac /= len( self.value )

        tot_mgels_vol = 0.0
        tot_overlap_vol = 0.0
        for i in range(N_probing_points) :
            if internal_probing_points[i,0] > 0 :
                tot_mgels_vol += 1.0
                if internal_probing_points[i,0] > 1 :
                    tot_overlap_vol += 1.0
        tot_mgels_vol = tot_mgels_vol / N_probing_points * box_sides[0] * box_sides[1] * box_sides[2]
        tot_overlap_vol = tot_overlap_vol / N_probing_points * box_sides[0] * box_sides[1] * box_sides[2]
        interp_vol_frac = tot_overlap_vol / tot_mgels_vol

        self.lock.acquire()
        logfile.write( str( tot_overlap_vol ) + " " + repr( tot_mgels_vol ) + " " + repr( interp_vol_frac ) + " " + repr( avg_interp_monom_frac ) + " " + repr( np.multiply.reduce( equil_config.box_sup - equil_config.box_inf ) ) + "\n" )
        logfile.flush()
        self.outfile.write( str( equil_config.time ) + " " + repr( avg_surf ) + " " + repr( avg_Sp ) + " " + repr( avg_vol ) + " " + repr( interp_vol_frac ) + " " + repr( avg_interp_monom_frac ) + " " + repr( tot_mgels_vol ) + "\n" )
        self.lock.release()
        logfile.close()

    def explore_cell_pairs( self , cells , cx , cy , cz , n_cells ) :
        parts = cells[ (cx+1)%n_cells[0] ][ (cy-1)%n_cells[1] ][ (cz-1)%n_cells[2] ] + cells[ (cx+1)%n_cells[0] ][ cy ][ (cz-1)%n_cells[2] ] + cells[ (cx+1)%n_cells[0] ][ (cy+1)%n_cells[1] ][ (cz-1)%n_cells[2] ]   + cells[ (cx+1)%n_cells[0] ][ (cy-1)%n_cells[1] ][ cz ] + cells[ (cx+1)%n_cells[0] ][ cy ][ cz ] + cells[ (cx+1)%n_cells[0] ][ (cy+1)%n_cells[1] ][ cz ]    + cells[ (cx+1)%n_cells[0] ][ (cy-1)%n_cells[1] ][ (cz+1)%n_cells[2] ] + cells[ (cx+1)%n_cells[0] ][ cy ][ (cz+1)%n_cells[2] ] + cells[ (cx+1)%n_cells[0] ][ (cy+1)%n_cells[1] ][ (cz+1)%n_cells[2] ]     + cells[ cx ][ (cy-1)%n_cells[1] ][ cz ]      + cells[ cx ][ (cy-1)%n_cells[1] ][ (cz+1)%n_cells[2] ] + cells[ cx ][ cy ][ (cz+1)%n_cells[2] ] + cells[ cx ][ (cy+1)%n_cells[1] ][ (cz+1)%n_cells[2] ]
        pairs = []
        for i in range( len( cells[ cx ][ cy ][ cz ] ) ) :
            for j in range( i+1 , len( cells[ cx ][ cy ][ cz ] ) ) :
                pairs.append( ( cells[ cx ][ cy ][ cz ][i] , cells[ cx ][ cy ][ cz ][j] ) )
            for j in range( len(parts) ) :
                pairs.append( ( cells[ cx ][ cy ][ cz ][i] , parts[j] ) )
        return pairs

    def local_cellist( self , parts_array_wrapped , box_sides ) :
        # building the cell list
        ###### ALTERNATIVA NUMPY :
        if 0 :
            n_cells = ( box_sides // cell_sides ).astype(int)
            cell_sides = box_sides / n_cells
            # the cell for each particle is calculated
            parts_cells = ( parts_array_wrapped // cell_sides ).astype(int)
        ###########################################################
        cell_sides = np.array( [ 4.0 , 4.0 , 4.0 ] )
        n_cells = box_sides // cell_sides
        cell_sides = box_sides / n_cells
        n_cells = n_cells.astype(int)
        # the cell for each particle is calculated
        parts_cells = ( parts_array_wrapped // cell_sides ).astype(int)
        cells = []
        for x in range(n_cells[0]) :
            ly = []
            for y in range(n_cells[1]) :
                lz = []
                for z in range(n_cells[2]) :
                    lst = []
                    lz.append( lst )
                ly.append( lz )
            cells.append( ly )
        for i in range( len(parts_cells) ) :
            cells[ parts_cells[i,0] ][ parts_cells[i,1] ][ parts_cells[i,2] ].append( i )
        parts_cells = np.array([])
        return cells , cell_sides , n_cells

    def compute_interpenetrated_atoms_LocalCut( self , equil_config , box_sides , box_infs , logfile ) :
        # wrapping the particles positions in the box, centered in the origin
        parts_array = equil_config.pos
        parts_features = np.column_stack(( equil_config.id , equil_config.type , equil_config.mol ))
        parts_features = np.insert( parts_features , 3 , np.arange(len(parts_features)) , axis=1 )
        parts_features = parts_features.astype(int)
        # all the counterions are excluded from the arrays
        parts_array = parts_array[ np.not_equal( parts_features[:,2] , 0 ) , : ]
        parts_features = parts_features[ np.not_equal( parts_features[:,2] , 0 ) , : ]
        # the box is centered in the origin and all the atoms wrapped
        parts_array_wrapped = parts_array - box_infs - ( box_sides*0.5 )
        periodic_images = np.sign(parts_array) * (  ( np.fabs(parts_array) // ( box_sides*0.5 ) + 1 ) // 2  )
        parts_array_wrapped = parts_array - ( periodic_images * box_sides )
        parts_array = np.array([])
        # the box is placed in the region x,y,z > 0
        parts_array_wrapped = parts_array_wrapped + box_sides*0.5
        cells , cell_sides , n_cells = self.local_cellist( parts_array_wrapped , box_sides )
        # calculation of the neighbours and of the lists of inner and interpenetrated monomers
        self.lock.acquire()
        tl0 = time.time()
        logfile.write( "starting local neighbours calculation, " + str(tl0) + "\n" )
        logfile.flush()
        self.lock.release()
        local_neighbours = np.array( [0]*len(parts_array_wrapped) )
        local_overlap = np.array( [0]*len(parts_array_wrapped) )
        curr_cell = 0
        tot_cells = n_cells[0] * n_cells[1] * n_cells[2]
        for cx,cy,cz in itertools.product( range(n_cells[0]) , range(n_cells[1]) , range(n_cells[2]) ) :
            curr_cell += 1
            if curr_cell % int(tot_cells//10) == 0 :
                self.lock.acquire()
                logfile.write( str(curr_cell/tot_cells*100) + "% of the cells have been explored" + "\n" )
                logfile.flush()
                self.lock.release()
            for (pi,pj) in self.explore_cell_pairs( cells , cx , cy , cz , n_cells ) :
                if sqrt(np.add.reduce( ( parts_array_wrapped[pi] - parts_array_wrapped[pj] + ( NN_periodic_image( parts_array_wrapped[pj] , parts_array_wrapped[pi] , box_sides )*box_sides ) )**2 , 0 )) < self.local_cutoff :
                    local_neighbours[pi] += 1
                    local_neighbours[pj] += 1
                    if equil_config.mol[pi] == equil_config.mol[pj] :
                        local_overlap[pi] -= 1
                        local_overlap[pj] -= 1
                    else :
                        local_overlap[pi] += 1
                        local_overlap[pj] += 1
        self.lock.acquire()
        logfile.write( "finished local neighbours calculation, after " + str(time.time()-tl0) + "\n" )
        logfile.flush()
        self.lock.release()
        for l in self.value :
            equil_config.molecules[l].properties["interpenetrated_atom_idxs"] = []
            equil_config.molecules[l].properties["bulk_atom_idxs"] = []
        for i in range( len(local_neighbours) ) :
            true_idx = parts_features[i,3] # general particle index
            l = parts_features[i,2] # molecule ID
            if local_overlap[i] >= 0 :
                equil_config.molecules[l].properties["interpenetrated_atom_idxs"].append( true_idx )
            else :
                equil_config.molecules[l].properties["bulk_atom_idxs"].append( true_idx )
        parts_array_wrapped = np.array([])
        parts_features = np.array([])
        cell_sides = np.array([])
        n_cells = np.array([])
        cells.clear()
        self.lock.acquire()
        logfile.write( "## mol_idx" + "\t" + "interpenetrated_atoms" + "\t" + "bulk_atoms" + "\n" )
        for l in self.value :
            logfile.write( str(l) + " " + str(len(equil_config.molecules[l].properties["interpenetrated_atom_idxs"])) + " " + str(len(equil_config.molecules[l].properties["bulk_atom_idxs"])) + "\n" )
        logfile.flush()
        self.lock.release()

    def compute_interpenetrated_atoms_OutCore( self , equil_config , box_sides , box_infs , logfile ) :
        for l in self.value :
            equil_config.molecules[l].properties["interpenetrated_atom_idxs"] = []
            equil_config.molecules[l].properties["bulk_atom_idxs"] = []
        for l in self.value :
            com = equil_config.com( "mol" , l )
            ## changed
            distances = np.sqrt( np.add.reduce( ( equil_config.pos[ equil_config.molecules[l].atom_idxs ] - com )**2 , 1 ) )
            equil_config.molecules[l].properties["interpenetrated_atom_idxs"] = equil_config.molecules[l].atom_idxs[ np.where( distances > self.core_radius ) ]
            equil_config.molecules[l].properties["bulk_atom_idxs"] = equil_config.molecules[l].atom_idxs[ np.where( distances <= self.core_radius ) ]
        self.lock.acquire()
        logfile.write( "## mol_idx" + "\t" + "interpenetrated_atoms" + "\t" + "bulk_atoms" + "\n" )
        for l in self.value :
            logfile.write( str(l) + " " + str(len(equil_config.molecules[l].properties["interpenetrated_atom_idxs"])) + " " + str(len(equil_config.molecules[l].properties["bulk_atom_idxs"])) + "\n" )
        logfile.flush()
        self.lock.release()

    def compute_interpenetrated_atoms_DensityCutoff( self , equil_config , box_sides , box_infs , logfile ) :
        local_cutoff = 4.0
        # wrapping the particles positions in the box, centered in the origin
        parts_array = equil_config.pos
        parts_features = np.column_stack(( equil_config.id , equil_config.type , equil_config.mol ))
        parts_features = np.insert( parts_features , 3 , np.arange(len(parts_features)) , axis=1 )
        parts_features = parts_features.astype(int)
        # all the counterions are excluded from the arrays
        parts_array = parts_array[ np.not_equal( parts_features[:,2] , 0 ) , : ]
        parts_features = parts_features[ np.not_equal( parts_features[:,2] , 0 ) , : ]
        # the box is centered in the origin and all the atoms wrapped
        parts_array_wrapped = parts_array - box_infs - ( box_sides*0.5 )
        periodic_images = np.sign(parts_array) * (  ( np.fabs(parts_array) // ( box_sides*0.5 ) + 1 ) // 2  )
        parts_array_wrapped = parts_array - ( periodic_images * box_sides )
        parts_array = np.array([])
        # the box is placed in the region x,y,z > 0
        parts_array_wrapped = parts_array_wrapped + box_sides*0.5
        cells , cell_sides , n_cells = self.local_cellist( parts_array_wrapped , box_sides )
        # calculation of the neighbours and of the lists of inner and interpenetrated monomers
        self.lock.acquire()
        tl0 = time.time()
        logfile.write( "starting local neighbours calculation, " + str(tl0) + "\n" )
        logfile.flush()
        self.lock.release()
        local_selfneighbours = np.array( [0]*len(parts_array_wrapped) )
        curr_cell = 0
        tot_cells = n_cells[0] * n_cells[1] * n_cells[2]
        for cx,cy,cz in itertools.product( range(n_cells[0]) , range(n_cells[1]) , range(n_cells[2]) ) :
            curr_cell += 1
            if curr_cell % int(tot_cells//10) == 0 :
                self.lock.acquire()
                logfile.write( str(curr_cell/tot_cells*100) + "% of the cells have been explored" + "\n" )
                logfile.flush()
                self.lock.release()
            for (pi,pj) in self.explore_cell_pairs( cells , cx , cy , cz , n_cells ) :
                if np.sqrt(np.add.reduce( ( parts_array_wrapped[pi] - parts_array_wrapped[pj] + ( NN_periodic_image( parts_array_wrapped[pj] , parts_array_wrapped[pi] , box_sides )*box_sides ) )**2 , 0 )) < local_cutoff :
                    if equil_config.mol[pi] == equil_config.mol[pj] :
                        local_selfneighbours[pi] += 1
                        local_selfneighbours[pj] += 1
        self.lock.acquire()
        logfile.write( "finished local neighbours calculation, after " + str(time.time()-tl0) + "\n" )
        logfile.flush()
        self.lock.release()
        for l in self.value :
            equil_config.molecules[l].properties["interpenetrated_atom_idxs"] = np.int32( [] )
            equil_config.molecules[l].properties["bulk_atom_idxs"] = np.int32( [] )
        local_selfdensity = local_selfneighbours / ( 4./3*np.pi*local_cutoff**3 )
        true_idx_V = parts_features[:,3] # general particle indexes
        l_V = parts_features[:,2] # molecule IDs
        for l in self.value :
            mol_i = np.argwhere( l_V == l )[:,0]
            bulk_i = mol_i[ (local_selfdensity[mol_i] > self.density_cutoff) ]
            if len(bulk_i) == 0 :
                print( "*** ERROR: The chosen value for the cutoff density is too large!!!" )
                exit(EXIT_FAILURE)
            else :
                equil_config.molecules[l].properties["bulk_atom_idxs"] = np.append( equil_config.molecules[l].properties["bulk_atom_idxs"] , true_idx_V[ bulk_i ] )
        parts_array_wrapped = np.array([])
        parts_features = np.array([])
        cell_sides = np.array([])
        n_cells = np.array([])
        cells.clear()
        ###   MODIFICAREEEEEEEEEEEEEE      e         MERGEREEEEEEEEEEEEEEEEE
        self.compute_molecule_meshes( equil_config , box_sides , logfile )
        all_atoms = np.int32( [] )
        for l in self.value :
            all_atoms = np.append( all_atoms , equil_config.molecules[l].atom_idxs , axis=0 )
        self.lock.acquire()
        logfile.write( ">> Computing the interpenetrated beads for mesh intersection volumes . . .\n" )
        logfile.flush()
        self.lock.release()
        probing_points = equil_config.pos[ all_atoms ]
        N_probing_points = len( probing_points )
        internal_probing_points = np.zeros( [ N_probing_points , 1 ] )
        probing_idxs = np.arange(N_probing_points)
        for l in self.value :
            periodic_com_images = NN_periodic_image( equil_config.molecules[l].properties["mesh_com"] , probing_points , box_sides )
            nearest_NNpoints = probing_points + periodic_com_images*box_sides
            points_com_dist = np.sqrt( np.add.reduce( ( equil_config.molecules[l].properties["mesh_com"] - nearest_NNpoints )**2 , 1 ) )
            select_vector = np.less_equal( points_com_dist , equil_config.molecules[l].properties["mesh-max_com_vtx_dist"] )
            nearest_NNpoints = nearest_NNpoints[ select_vector , : ]
            nearest_NNpoints_idxs = probing_idxs[ select_vector ]
            Npts2analyse = len( nearest_NNpoints )
            tv0 = time.time()
            if self.method == "chull" :
                for i in range(Npts2analyse) :
                    if i % 5000 == 0 :
                        self.lock.acquire()
                        logfile.write( str(i) + " of " + str(Npts2analyse) + " beads analysed, after " + str(time.time()-tv0) + "\n" )
                        logfile.flush()
                        self.lock.release()
                    if equil_config.mol[ all_atoms[ nearest_NNpoints_idxs[i] ] ] != l :
                        if self.point_is_inside_convex_mesh( nearest_NNpoints[i] , box_sides , equil_config.molecules[l].properties["mesh_face_normals"] , equil_config.molecules[l].properties["mesh_triangles"] , NNcalc=False ) :
                            internal_probing_points[ nearest_NNpoints_idxs[i] , 0 ] += 1
            elif self.method == "smesh" :
                inner_regions = []
                surface_mesh = equil_config.molecules[l].properties["AlphaSurface"].surfaces['surface']
                for r in range(len( surface_mesh.regions['Filled'] )) :
                    if surface_mesh.regions['Filled'][r] != 0 :
                        inner_regions.append( r )
                for i in range(Npts2analyse) :
                    if i % 5000 == 0 :
                        self.lock.acquire()
                        logfile.write( str(i) + " of " + str(Npts2analyse) + " beads analysed, after " + str(time.time()-tv0) + "\n" )
                        logfile.flush()
                        self.lock.release()
                    if equil_config.mol[ all_atoms[ nearest_NNpoints_idxs[i] ] ] != l :
                        if surface_mesh.locate_point( nearest_NNpoints[i] , eps=1e-6 ) in inner_regions :
                            internal_probing_points[ nearest_NNpoints_idxs[i] , 0 ] += 1
        ## changed
        mol_ids = equil_config.mol[ all_atoms ]
        for l in self.value :
            AAidxs = np.argwhere( mol_ids == l )[:,0]
            equil_config.molecules[l].properties["interpenetrated_atom_idxs"] = np.append( equil_config.molecules[l].properties["interpenetrated_atom_idxs"] , all_atoms[  AAidxs[ (internal_probing_points[AAidxs] > 0)[:,0] ]  ] )
        self.lock.acquire()
        logfile.write( "## mol_idx" + "\t" + "interpenetrated_atoms" + "\t" + "bulk_atoms" + "\n" )
        for l in self.value :
            logfile.write( str(l) + " " + str(len(equil_config.molecules[l].properties["interpenetrated_atom_idxs"])) + " " + str(len(equil_config.molecules[l].properties["bulk_atom_idxs"])) + "\n" )
        logfile.flush()
        self.lock.release()

    def compute_interpenetrated_atoms_MeshVolumes( self , equil_config , box_sides , box_infs , logfile ) :
        all_atoms = np.int32( [] )
        for l in self.value :
            equil_config.molecules[l].properties["interpenetrated_atom_idxs"] = np.int32( [] )
            equil_config.molecules[l].properties["bulk_atom_idxs"] = np.int32( [] )
            all_atoms = np.append( all_atoms , equil_config.molecules[l].atom_idxs , axis=0 )
        self.lock.acquire()
        logfile.write( ">> Computing the interpenetrated beads for mesh intersection volumes . . .\n" )
        logfile.flush()
        self.lock.release()
        probing_points = equil_config.pos[ all_atoms ]
        N_probing_points = len( probing_points )
        internal_probing_points = np.zeros( [ N_probing_points , 1 ] )
        probing_idxs = np.arange(N_probing_points)
        for l in self.value :
            periodic_com_images = NN_periodic_image( equil_config.molecules[l].properties["mesh_com"] , probing_points , box_sides )
            nearest_NNpoints = probing_points + periodic_com_images*box_sides
            points_com_dist = np.sqrt( np.add.reduce( ( equil_config.molecules[l].properties["mesh_com"] - nearest_NNpoints )**2 , 1 ) )
            select_vector = np.less_equal( points_com_dist , equil_config.molecules[l].properties["mesh-max_com_vtx_dist"] )
            nearest_NNpoints = nearest_NNpoints[ select_vector , : ]
            nearest_NNpoints_idxs = probing_idxs[ select_vector ]
            Npts2analyse = len( nearest_NNpoints )
            tv0 = time.time()
            if self.method == "chull" :
                for i in range(Npts2analyse) :
                    if i % 5000 == 0 :
                        self.lock.acquire()
                        logfile.write( str(i) + " of " + str(Npts2analyse) + " beads analysed, after " + str(time.time()-tv0) + "\n" )
                        logfile.flush()
                        self.lock.release()
                    if equil_config.mol[ all_atoms[ nearest_NNpoints_idxs[i] ] ] != l :
                        if self.point_is_inside_convex_mesh( nearest_NNpoints[i] , box_sides , equil_config.molecules[l].properties["mesh_face_normals"] , equil_config.molecules[l].properties["mesh_triangles"] , NNcalc=False ) :
                            internal_probing_points[ nearest_NNpoints_idxs[i] , 0 ] += 1
            elif self.method == "smesh" :
                inner_regions = []
                surface_mesh = equil_config.molecules[l].properties["AlphaSurface"].surfaces['surface']
                for r in range(len( surface_mesh.regions['Filled'] )) :
                    if surface_mesh.regions['Filled'][r] != 0 :
                        inner_regions.append( r )
                for i in range(Npts2analyse) :
                    if i % 5000 == 0 :
                        self.lock.acquire()
                        logfile.write( str(i) + " of " + str(Npts2analyse) + " beads analysed, after " + str(time.time()-tv0) + "\n" )
                        logfile.flush()
                        self.lock.release()
                    if equil_config.mol[ all_atoms[ nearest_NNpoints_idxs[i] ] ] != l :
                        if surface_mesh.locate_point( nearest_NNpoints[i] , eps=1e-6 ) in inner_regions :
                            internal_probing_points[ nearest_NNpoints_idxs[i] , 0 ] += 1
        ## changed
        mol_ids = equil_config.mol[ all_atoms ]
        for l in self.value :
            AAidxs = np.argwhere( mol_ids == l )[:,0]
            equil_config.molecules[l].properties["interpenetrated_atom_idxs"] = np.append( equil_config.molecules[l].properties["interpenetrated_atom_idxs"] , all_atoms[  AAidxs[ (internal_probing_points[AAidxs] > 0)[:,0] ]  ] )
            equil_config.molecules[l].properties["bulk_atom_idxs"] = np.append( equil_config.molecules[l].properties["bulk_atom_idxs"] , all_atoms[  AAidxs[ (internal_probing_points[AAidxs] == 0)[:,0] ]  ] )
        self.lock.acquire()
        logfile.write( "## mol_idx" + "\t" + "interpenetrated_atoms" + "\t" + "bulk_atoms" + "\n" )
        for l in self.value :
            logfile.write( str(l) + " " + str(len(equil_config.molecules[l].properties["interpenetrated_atom_idxs"])) + " " + str(len(equil_config.molecules[l].properties["bulk_atom_idxs"])) + "\n" )
        logfile.flush()
        self.lock.release()


    def terminate( self ) :
        self.outfile.close()
        for l in self.molfiles.keys() :
            self.molfiles[l].close()

    def return_values( self ) :
        self.lock.acquire()
        self.outfile.flush()
        for l in self.molfiles.keys() :
            self.molfiles[l].flush()
        self.lock.release()
        return {}

    def print_help( self ) :
        print( "\t\t\t{ \033[1m-shape_overlap\033[0m <mol> <mol_ids(1,2:5,...)> { file <output_fname> } { method {chull | smesh <probing_radius>} }" )
        print( "\t\t\t\t{ local_cut <interpenetration_cutoff> }                   { dump_meshes }" )
        print( "\t\t\t\t{ interpenetration overlap_vol|atoms_fraction|core_radius <core_radius>|density_cutoff <dens_cut> }" )
        print( "\t\t\t\t  [output_fname(s)=shape_overlap(_mol?).dat , method=chull , interpenetration_cutoff=4.5 , interpenetration=overlap_vol] } " )
