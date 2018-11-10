#-------------------------------------------------------------------------------
import numpy as numpy
__all__=['rotating_to_heliocentric_array','heliocentric_to_rotating_array',
'rotating_to_heliocentric_array_precalc','heliocentric_to_rotating_array_precalc']
#-------------------------------------------------------------------------------
cos = numpy.cos
arccos = numpy.arccos
sin = numpy.sin
arcsin = numpy.arcsin
tan = numpy.tan
arctan = numpy.arctan
pi = numpy.pi
#-------------------------------------------------------------------------------
# replace these constants with astropy units
G=6.67428e-11 # m3 kg-1 s-2
M_sun=1.98855e30 # kg
AU=1.496e11 # m
#-------------------------------------------------------------------------------
def rotating_to_heliocentric_array(r,v,a,t):
    '''Function to transform from rotating frame (circular orbit of radius a, at R(t=0)=(a,0,0)) to the heliocentric frame'''
    Om_k=numpy.sqrt(G*(M_sun)/(a**3.0)) #angular velocity at semimajor axis a_k
    a_vec=a*numpy.array([cos(Om_k*t),sin(Om_k*t),0.0]) #time gives location of frame
    Om_vec=numpy.array([0.0,0.0,Om_k]) # define angular velocity vector
    theta=Om_k*t # angle of rotationa s a function of time
    rot_mat=numpy.array([cos(theta),-sin(theta),0,sin(theta),cos(theta),0,0,0,1]).reshape((3,3)) # z axis rotation vector
    Om_skew=numpy.array([0,-Om_k,0,Om_k,0,0,0,0,0]).reshape((3,3)) #define skew symmetric matrix for angular velocity (0,0,Om_k)
    rot_mat_dot=numpy.dot(Om_skew,rot_mat) #time derivative of the rotation matrix

    if len(r.shape)>1:
        N=len(r)
        # R=numpy.zeros((N,3))
        # V=numpy.zeros((N,3))
        for i in range(N):
            # R[i,:]=a_vec+numpy.dot(rot_mat,r[i,:]) # transform position
            # V[i,:]=numpy.cross(Om_vec,a_vec)+numpy.dot(rot_mat_dot,r[i,:])+numpy.dot(rot_mat,v[i,:]) # transform velocity
            R=a_vec+numpy.dot(rot_mat,r.T).T # transform position
            V=numpy.cross(Om_vec,a_vec)+numpy.dot(rot_mat_dot,r.T).T+numpy.dot(rot_mat,v.T).T # transform velocity
    else:
        R=a_vec+numpy.dot(rot_mat,r) # transform position
        V=numpy.cross(Om_vec,a_vec)+numpy.dot(rot_mat_dot,r)+numpy.dot(rot_mat,v) # transform velocity
    return R,V
#-------------------------------------------------------------------------------
def heliocentric_to_rotating_array(R,V,a,t):
    '''Function to transform from heliocentric frame to rotating frame (circular orbit of radius a, at R(t=0)=(a,0,0))'''
    Om_k=numpy.sqrt(G*(M_sun)/(a**3.0)) #angular velocity at semimajor axis a_k
    a_vec=a*numpy.array([cos(Om_k*t),sin(Om_k*t),0.0])
    Om_vec=numpy.array([0.0,0.0,Om_k])
    theta=-Om_k*t #reverse angle
    rot_mat=numpy.array([cos(theta),-sin(theta),0,sin(theta),cos(theta),0,0,0,1]).reshape((3,3))
    Om_skew=numpy.array([0,-Om_k,0,Om_k,0,0,0,0,0]).reshape((3,3)) #define skew symmetric matrix for angular velocity (0,0,Om_k)
    rot_mat_dot=-numpy.dot(Om_skew,rot_mat) #time derivative of the rotation matrix

    if len(R.shape)>1:
        N=len(R)
        # r=numpy.zeros((N,3))
        # v=numpy.zeros((N,3))
        for i in range(N):
            # r[i,:]=numpy.dot(rot_mat,R[i,:])-numpy.dot(rot_mat,a_vec)
            # v[i,:]=numpy.dot(rot_mat_dot,R[i,:])+numpy.dot(rot_mat,V[i,:])-numpy.dot(rot_mat_dot,a_vec)-numpy.dot(rot_mat,numpy.cross(Om_vec,a_vec))
            r=numpy.dot(rot_mat,R.T).T-numpy.dot(rot_mat,a_vec)
            v=numpy.dot(rot_mat_dot,R.T).T+numpy.dot(rot_mat,V.T).T-numpy.dot(rot_mat_dot,a_vec)-numpy.dot(rot_mat,numpy.cross(Om_vec,a_vec))
    else:
        r=numpy.dot(rot_mat,R)-numpy.dot(rot_mat,a_vec)
        v=numpy.dot(rot_mat_dot,R)+numpy.dot(rot_mat,V)-numpy.dot(rot_mat_dot,a_vec)-numpy.dot(rot_mat,numpy.cross(Om_vec,a_vec))
    return r,v
#-------------------------------------------------------------------------------

def rotating_to_heliocentric_array_precalc(r,v,a,t):
    '''Function to transform from rotating frame (circular orbit of radius a, at R(t=0)=(a,0,0)) to the heliocentric frame
    Here we precalculate all the vectors and matrices, then loop over all particles in the array'''
    Om_k=numpy.sqrt(G*(M_sun)/(a**3.0)) #angular velocity at semimajor axis a_k
    a_vec=a*numpy.array([cos(Om_k*t),sin(Om_k*t),0.0]) #time gives location of frame
    Om_vec=numpy.array([0.0,0.0,Om_k]) # define angular velocity vector
    theta=Om_k*t # angle of rotationa s a function of time
    rot_mat=numpy.array([cos(theta),-sin(theta),0,sin(theta),cos(theta),0,0,0,1]).reshape((3,3)) # z axis rotation vector
    Om_skew=numpy.array([0,-Om_k,0,Om_k,0,0,0,0,0]).reshape((3,3)) #define skew symmetric matrix for angular velocity (0,0,Om_k)
    rot_mat_dot=numpy.dot(Om_skew,rot_mat) #time derivative of the rotation matrix
    # Several particle case
    if len(r.shape)>1:
        N=len(r[:,0])
        R=numpy.zeros((N,3))
        V=numpy.zeros((N,3))
        for i in range(N):
            #print i
            R[i,:]=a_vec+numpy.dot(rot_mat,r[i,:]) # transform position
            V[i,:]=numpy.cross(Om_vec,a_vec)+numpy.dot(rot_mat_dot,r[i,:])+numpy.dot(rot_mat,v[i,:]) # transform velocity
    # Single particle case
    else:
        R=a_vec+numpy.dot(rot_mat,r) # transform position
        V=numpy.cross(Om_vec,a_vec)+numpy.dot(rot_mat_dot,r)+numpy.dot(rot_mat,v) # transform velocity

    return R,V

def heliocentric_to_rotating_array_precalc(R,V,a,t):
    '''Function to transform from heliocentric frame to rotating frame (circular orbit of radius a, at R(t=0)=(a,0,0))
    Here we precalculate all the vectors and matrices, then loop over all particles in the array'''
    Om_k=numpy.sqrt(G*(M_sun)/(a**3.0)) #angular velocity at semimajor axis a_k
    a_vec=a*numpy.array([cos(Om_k*t),sin(Om_k*t),0.0])
    Om_vec=numpy.array([0.0,0.0,Om_k])
    theta=-Om_k*t #reverse angle
    rot_mat=numpy.array([cos(theta),-sin(theta),0,sin(theta),cos(theta),0,0,0,1]).reshape((3,3))
    Om_skew=numpy.array([0,-Om_k,0,Om_k,0,0,0,0,0]).reshape((3,3)) #define skew symmetric matrix for angular velocity (0,0,Om_k)
    rot_mat_dot=-numpy.dot(Om_skew,rot_mat) #time derivative of the rotation matrix
    # Several particle case
    if len(R.shape)>1:
        N=len(R[:,0])
        r=numpy.zeros((N,3))
        v=numpy.zeros((N,3))
        for i in range(N):
            #print i
            r[i,:]=numpy.dot(rot_mat,R[i,:])-numpy.dot(rot_mat,a_vec)
            v[i,:]=numpy.dot(rot_mat_dot,R[i,:])+numpy.dot(rot_mat,V[i,:])-numpy.dot(rot_mat_dot,a_vec)-numpy.dot(rot_mat,numpy.cross(Om_vec,a_vec))
    # Single particle case
    else:
        r=numpy.dot(rot_mat,R)-numpy.dot(rot_mat,a_vec)
        v=numpy.dot(rot_mat_dot,R)+numpy.dot(rot_mat,V)-numpy.dot(rot_mat_dot,a_vec)-numpy.dot(rot_mat,numpy.cross(Om_vec,a_vec))

    return r,v

#-------------------------------------------------------------------------------
