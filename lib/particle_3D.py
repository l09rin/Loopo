import math

# definition of the class particle in D=3
class particle :

    def __init__( self , id0=0 , type0=0 , x0=-1.0 , y0=-1.0 , z0=-1.0 , q0=0.0 , mol0=0 ) :
        self.id = id0
        self.type = type0
        self.x = x0
        self.y = y0
        self.z = z0
        self.q = q0
        self.mol = mol0

    def __add__( self , other ) :
        sum = particle( self.id , self.type , self.x+other.x , self.y+other.y , self.z+other.z , self.q , self.mol )
        return sum

    def __iadd__( self , other ) :
        self.x = self.x + other.x
        self.y = self.y + other.y
        self.z = self.z + other.z
        return self

    def __sub__( self , other ) :
        diff = particle( self.id , self.type , self.x-other.x , self.y-other.y , self.z-other.z , self.q , self.mol )
        return diff

    def __isub__( self , other ) :
        self.x = self.x - other.x
        self.y = self.y - other.y
        self.z = self.z - other.z
        return self

    def __mul__( self , scal ) :
        scaled = particle( self.id , self.type , self.x*scal , self.y*scal , self.z*scal , self.q , self.mol )
        return scaled

    def __imul__( self , scal ) :
        self.x = self.x * scal
        self.y = self.y * scal
        self.z = self.z * scal
        return self

    def __truediv__( self , scal ) :
        scaled = particle( self.id , self.type , self.x/scal , self.y/scal , self.z/scal , self.q , self.mol )
        return scaled

    def __floordiv__( self , other ) :
        floored = particle( self.id , self.type , self.x//other.x , self.y//other.y , self.z//other.z , self.q , self.mol )
        return floored

    def __itruediv__( self , scal ) :
        self.x = self.x / scal
        self.y = self.y / scal
        self.z = self.z / scal
        return self

    def __neg__( self ) :
        opp = particle( self.id , self.type , -self.x , -self.y , -self.z , self.q , self.mol )
        return opp

    def __eq__( self , other ) :
        bit = 1
        if self.x == other.x and self.y == other.y and self.z == other.z and self.id == other.id and self.q == other.q and self.type == other.type and self.mol == other.mol :
            bit = 1
        else :
            bit = 0
        return bit

    def __ne__( self ) :
        bit = 0
        if self.x == other.x and self.y == other.y and self.z == other.z and self.id == other.id and self.q == other.q and self.type == other.type and self.mol == other.mol :
            bit = 0
        else :
            bit = 1
        return bit

    def copy( self , other ) :
        self.id = other.id
        self.type = other.type
        self.x = other.x
        self.y = other.y
        self.z = other.z
        self.q = other.q
        self.mol = other.mol

    def nearest_copy( self , other , box_edges ) :
        apparent = particle()
        apparent.copy( other )
        if abs( other.x - self.x ) < box_edges.x/2.0 :
            apparent.x = other.x
        elif other.x - self.x > 0 :
            apparent.x = other.x - box_edges.x
            while abs( apparent.x - self.x ) > box_edges.x/2.0 :
                apparent.x -= box_edges.x
        else :
            apparent.x = other.x + box_edges.x
            while abs( apparent.x - self.x ) > box_edges.x/2.0 :
                apparent.x += box_edges.x
        
        if abs( other.y - self.y ) < box_edges.y/2.0 :
            apparent.y = other.y
        elif other.y - self.y > 0 :
            apparent.y = other.y - box_edges.y
            while abs( apparent.y - self.y ) > box_edges.y/2.0 :
                apparent.y -= box_edges.y
        else :
            apparent.y = other.y + box_edges.y
            while abs( apparent.y - self.y ) > box_edges.y/2.0 :
                apparent.y += box_edges.y
        
        if abs( other.z - self.z ) < box_edges.z/2.0 :
            apparent.z = other.z
        elif other.z - self.z > 0 :
            apparent.z = other.z - box_edges.z
            while abs( apparent.z - self.z ) > box_edges.z/2.0 :
                apparent.z -= box_edges.z
        else :
            apparent.z = other.z + box_edges.z
            while abs( apparent.z - self.z ) > box_edges.z/2.0 :
                apparent.z += box_edges.z

        return apparent


    def distance( self , other , box_edges ) :
        part2 = self.nearest_copy( other , box_edges )
        return math.sqrt( (part2.x-self.x)**2 + (part2.y-self.y)**2 + (part2.z-self.z)**2 )

    def norm( self ) :
        return math.sqrt( self.x**2 + self.y**2 + self.z**2 )

    def sq_norm( self ) :
        return self.x**2 + self.y**2 + self.z**2

    def scalar_prod( self , other ) :
        return self.x*other.x + self.y*other.y + self.z*other.z

    def comp2comp_prod( self , other ) :
        prod = particle( self.id , self.type , self.x*other.x , self.y*other.y , self.z*other.z , self.q , self.mol )
        return prod

    def vol( self ) :
        return self.x * self.y * self.z

    def min_comp( self ) :
        return min( self.x , self.y , self.z )
