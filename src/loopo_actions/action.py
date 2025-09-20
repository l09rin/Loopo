import os
# MULTI THREADING CALCULATIONS
# import threading ## NOT WORKING! BECAUSE OF THE GLOBAL INTERPRETER LOCK OF PYTHON, THAT PREVENTS MORE THAN ONE THREAD AT A TIME MODIFIES AN OBJECT
import multiprocessing as mp

# MATH TOOLS
from math import *
import cmath
import random
import numpy as np    # numpy arrays are useful to arrange data to plot
# import pandas as pd
EXIT_FAILURE = -1

class COLUMNS() :
    def __init__( self ) :
        self.id = -1
        self.type = -1
        self.mol = -1
        self.x = 2
        self.y = 3
        self.z = 4
        self.q = -1
        self.vx = -1
        self.vy = -1
        self.vz = -1

    def set_id( self , val ) :
        self.id = val
    def set_type( self , val ) :
        self.type = val
    def set_mol( self , val ) :
        self.mol = val
    def set_x( self , val ) :
        self.x = val
    def set_y( self , val ) :
        self.y = val
    def set_z( self , val ) :
        self.z = val
    def set_q( self , val ) :
        self.q = val
    def set_vx( self , val ) :
        self.vx = val
    def set_vy( self , val ) :
        self.vy = val
    def set_vz( self , val ) :
        self.vz = val

class DATA_FILES() :
    def __init__( self ) :
        self.read_confs = 0
        self.files = []
        self.current_file = open(os.devnull,"w")
        self.current_file_name = ""
        self.NOT_EMPTY = True
        self.format = ""

    def initialize( self ) :
        if self.current_file_name != "" and len(self.files) != 0 :
            # input file opening
            self.current_file_name = self.files.pop(0)
            print( " >>> Opening the file : " + self.current_file_name )
            self.current_file = open( self.current_file_name , "r" )

    def update( self ) :
        if len(self.current_file_name) == 0 :
            # input file opening
            self.current_file_name = self.files.pop(0)
            print( " >>> Opening the file : " + self.current_file_name )
            self.current_file = open( self.current_file_name , "r" )
        if self.current_file.closed :
            self.current_file.close()
            if len(self.files) > 0 :
                self.current_file_name = self.files.pop(0)
                print( " >>> Opening the file : " + self.current_file_name )
                self.current_file = open( self.current_file_name , "r" )
            else :
                self.NOT_EMPTY = False

    def isnotempty( self ) :
        return self.NOT_EMPTY

    def get_read_confs( self ) :
        return self.read_confs

    def get_format( self ) :
        return self.format

    def get_file( self ) :
        return self.current_file

    def file_isclosed( self ) :
        return self.current_file.closed

    def set_format( self , fmt ) :
        self.format = fmt

    def set_files( self , lst ) :
        self.files = lst

    def increment( self ) :
        self.read_confs += 1

    def readline( self ) :
        return self.current_file.readline()

    def readlines( self , nbytes = 1 ) :
        return self.current_file.readlines( nbytes )

    def loadtxt( self , nlines = 1 ) :
        # return np.array( pd.read_csv( self.current_file , delimiter=' ', header=None, nrows=nlines) ).astype(str)
        return np.loadtxt( self.current_file , dtype = 'str' , max_rows = nlines )

    def savetxt( self , nlines = 1 ) :
        np.savetxt( self.current_file , particles_lines , fmt='%s')

    def close( self ) :
        self.current_file.close()

class ACTION() :
    def __init__( self ) :
        self.operation = ""
        self.lock = mp.Lock()
        self.MAINlock = mp.Lock()
        self.bondslist_file = ""
        self.bonded = {}
        self.xlinkers_ids = []

    def generate_bondslist( self ) :
    # creation of the bonds list, if needed
    # the bonds list is returned as a dictionary containing the list of neighbours for each particle, indicated with atom IDs
        bonded = {}
        if self.bondslist_file != "" :
            try :
                bonds_file = open( self.bondslist_file , "r" )
            except :
                self.MAINlock.acquire()
                print( " *** ERROR: The file containing the list of bonds cannot be opened ! " )
                self.MAINlock.release()
            # insertion of bonds
            self.MAINlock.acquire()
            print( "I am reading the bonds list . . ." )
            self.MAINlock.release()
            line = bonds_file.readline()
            while ( not line.isspace() ) and line != "" :
                line.strip()
                words = line.split()
                if len(words) > 1 :
                    ip = int( words[2] )
                    jp = int( words[3] )
                    if ip not in bonded :
                        bonded[ip] = []
                    if jp not in bonded :
                        bonded[jp] = []
                    bonded[ip].append( jp )
                    bonded[jp].append( ip )
                line = bonds_file.readline()
            bonds_file.close()
            return bonded

    def generate_XL_list( self ) :
        if len(self.xlinkers_ids) == 0 :
            if not self.bonded :
                for atomid in self.bonded.keys() :
                    if len(self.bonded[atomid]) > 2 :
                        self.xlinkers_ids.append(atomid)

    def terminate( self ) :
        return 0

    def return_values( self ) :
        return {}

    def merge_return_values( self , values_list ) :
        pass

    def print_help( self ) :
        pass
