'''
Useful functions!

Notes:
- Units, generally units are SI (meters, seconds...) and angles are in radians, unless otherwise stated
'''

import numpy as numpy #Imports Python mathematical functions library
import math
#import scipy
#from scipy.optimize import minimize
#try:
#    import matplotlib.pyplot as pyplot
#except:
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec
import rebound
import subprocess
import datetime
import os
try:
    import pandas as pd
except:
    print "proceed without pandas"
#-------------------------------------------------------------------------------
cos = numpy.cos
arccos = numpy.arccos
sin = numpy.sin
arcsin = numpy.arcsin
tan = numpy.tan
arctan = numpy.arctan
pi = numpy.pi
#-------------------------------------------------------------------------------
'''
p=numpy.load('o_params.npy')
#assign the parameters
G=p[0]
AU=p[1]
M_sun=p[2]
mu=p[3]
rho_E=p[4]
d_E=p[5]
rho_S=p[6]
q=p[7]
m=p[8]
V_inf=p[9]
d_min=p[10]
d_max=p[11]
imp_param=p[12]
limit=p[13]
'''
G=6.67428e-11 # m3 kg-1 s-2
M_sun=1.98855e30 # kg
AU=1.496e11 # m
limit=1e-12 # fractional cut off for find_E function
#a_0=30*AU
day_s=24.0*60.0*60.0 # 1 day in seconds
year_s=365.25*day_s # 1 year in seconds

pyplot_colours=pyplot.rcParams['axes.prop_cycle'].by_key()['color'] # list of matplotlib default colours

#run_params
#cwd,rp->Ntot,rp->a,rp->R_eq,rp->rho,rp->M_tot,rp->R_c,rp->OMEGA,rp->X,rp->OM_circ,rp->f,rp->dt,run_time,np,ins

#-------------------------------------------------------------------------------
def r_elliptical(a,e,theta):
    '''
    Function to find distance, r, at true anomaly, theta, around a keplerian orbit

    Parameters
    ----------
    a
        semi major axis (m)
    e
        orbital eccentricity
    theta
        true anomaly (rad)
    '''
    r = (a*(1-e**2))/(1+e*cos(theta))
    return r
#-------------------------------------------------------------------------------
def find_e_p(OMEGA,omega,I):
    '''Function to find the normalised Lenz vector along the apsidal line, uses PSS eq11.39a'''
    e_p=numpy.zeros(3)
    e_p[0]=(cos(OMEGA)*cos(omega))-(cos(I)*sin(OMEGA)*sin(omega)) #x
    e_p[1]=(sin(OMEGA)*cos(omega))+(cos(I)*cos(OMEGA)*sin(omega)) #y
    e_p[2]=sin(I)*sin(omega) #z
    return e_p
#-------------------------------------------------------------------------------
def find_e_Q(OMEGA,omega,I):
    '''Function to find the normalised e_Q vector which is perp to h and e_a and
    lies in the orbital plane, uses PSS eq11.39b'''
    e_Q=numpy.zeros(3)
    e_Q[0]=-(cos(OMEGA)*sin(omega))-(cos(I)*sin(OMEGA)*cos(omega)) #x
    e_Q[1]=-(sin(OMEGA)*sin(omega))+(cos(I)*cos(OMEGA)*cos(omega)) #y
    e_Q[2]=sin(I)*cos(omega) #z
    return e_Q
#-------------------------------------------------------------------------------
def planet_orbit(orb,n):
    '''function to find the xyz orbit data relative to the reference point'''
    #orb is the 6 orbital elements: a,e,I,OM,om,f and n is number of points to plot with
    a=orb[0]
    e=orb[1]
    I=orb[2]
    OMEGA=orb[3]
    omega=orb[4]
    #specify theta from 0 to 2pi radians: i.e. number of points on orbit, THE TRUE ANOMALY
    theta = numpy.linspace(0,2*pi, n)
    #-------------------------------------------------------------------------------
    r = numpy.zeros(len(theta))
    for i in range(len(theta)):
        r[i]=r_elliptical(a,e,theta[i])
    #-------------------------------------------------------------------------------
    #determine orbital basis vectors
    e_p=find_e_p(OMEGA,omega,I)
    e_Q=find_e_Q(OMEGA,omega,I)
    #-------------------------------------------------------------------------------
    #define the r(x,y,z) position array
    pos=numpy.zeros((len(theta),3))
    for n in range(len(pos[:,0])):
        pos[n,:]=r[n]*((cos(theta[n])*e_p)+(sin(theta[n])*e_Q))#PSS eq 11.36a
        pos[n,:]=pos[n,:]
    return pos
#-------------------------------------------------------------------------------
def pos_vel_from_orbit(orb,mu):
    '''function to find xyz and vxvyvz from the a,e,I,OMEGA,omega,f_true orbit data, where mu=G(M+m)'''
    a=orb[0]
    e=orb[1]
    I=orb[2]
    OMEGA=orb[3]
    omega=orb[4]
    f_true = orb[5]

    n=numpy.sqrt(mu/(a**3))
    eta=numpy.sqrt(1.0-(e**2))

    r=r_elliptical(a,e,f_true)
    #-------------------------------------------------------------------------------
    #determine orbital basis vectors
    e_p=find_e_p(OMEGA,omega,I)
    e_Q=find_e_Q(OMEGA,omega,I)
    #-------------------------------------------------------------------------------
    #define the r(x,y,z) position array
    pos=numpy.zeros(3)
    pos=r*((cos(f_true)*e_p)+(sin(f_true)*e_Q))#PSS eq 11.36a
    vel=numpy.zeros(3)
    vel=(n*a/eta)*((-sin(f_true)*e_p)+((e+cos(f_true))*e_Q))
    # print (n*a/eta)
    # print numpy.linalg.norm(((-sin(f_true)*e_p)+((e+cos(f_true))*e_Q)))
    # print 2.0*pi/n
    return pos,vel

#-------------------------------------------------------------------------------
def find_E(M,e):
    '''Solve Kepler's Equation.
    Function to find E from M, to within a certain limit.
    This is an iterative method and DIVERGES for e>0.662743
    See Murray and Dermott Solar System Dynamics p.35'''
    E_0=M
    deltaE=2.0*limit #set deltaE to some value, greater than the limit
    while deltaE>limit:
        E_1=M+e*sin(E_0)
        deltaE=numpy.absolute(E_1-E_0)
        E_0=E_1
    return E_0
#-------------------------------------------------------------------------------
def find_f_true(E,e):
    '''Function to find f from E'''
    f_true=2*arctan(numpy.sqrt((1+e)/(1-e))*tan(E/2))
    return f_true
#-------------------------------------------------------------------------------
def orbit_by_time(orb,mu,t0,t):
    '''find the position and velocity of a body on a heliocentric keplerian orbit,
    as a function of time

    Parameters
    ----------
    orb
        list of 5 orbital parameters [a,e,I,OMEGA,omega],f_true is found from time, t

    mu
        gravitational parameter = G * total mass
    t0
        reference time (i.e. time when f_true = 0 ?)
    t
        the time, from t0, at which we want to fidn the position and velocity
    '''
    a=orb[0]
    e=orb[1]
    I=orb[2]
    OMEGA=orb[3]
    omega=orb[4]
    eta=numpy.sqrt(1.0-(e**2))
    n=numpy.sqrt(mu/(a**3.0)) # mean motion
    M=n*(t-t0) # mean anomaly
    E=find_E(M,e) # eccentric anomaly
    f_true=find_f_true(E,e) # true anomaly
    r=r_elliptical(a,e,f_true)
    #determine orbital basis vectors
    e_p=find_e_p(OMEGA,omega,I)
    e_Q=find_e_Q(OMEGA,omega,I)
    pos=r*((cos(f_true)*e_p)+(sin(f_true)*e_Q))#PSS eq 11.36a
    vel=(n*a/eta)*((-sin(f_true)*e_p)+((e+cos(f_true))*e_Q))
    return pos,vel
#-------------------------------------------------------------------------------
def centre_of_mass(pos,m):
    '''Function to find CoM of N particles from the position and mass arrays.
    Note: you can also find the centre of mass velocity by passing a velocity array instead.

    Parameters
    ----------
    pos
        2d numpy array, containing x,y,z position for each particle, shape (N,3)
    m
        1d numpy array of particle masses, length (N)
    '''
    xm=pos[:,0]*m
    ym=pos[:,1]*m
    zm=pos[:,2]*m
    M=numpy.sum(m)
    CoM=numpy.zeros(3)
    CoM[0]=numpy.sum(xm)/M
    CoM[1]=numpy.sum(ym)/M
    CoM[2]=numpy.sum(zm)/M
    return CoM
#-------------------------------------------------------------------------------
def N_largest(N,i,fname):
    '''Function to find the N largest particles from the data array'''
    t=numpy.loadtxt(fname,usecols=(0,)) # read only first column to get time
    dat=numpy.loadtxt(fname,skiprows=1) # read in the data skipping the first row
    t=t[0] # do some wizardry to extract the time value from the whole column
    #create structured array
    dat2 = numpy.core.records.fromarrays(dat.transpose(),names='x,y,z,vx,vy,vz,m,r',formats = 'f16,f16,f16,f16,f16,f16,f16,f16')
    dat2=numpy.sort(dat2, order='m')
    dat2=numpy.flipud(dat2)
    new_dat=numpy.zeros((N,8))
    for i in range(N):
        new_dat[i,0]=dat2[i][6] #mass
        new_dat[i,1]=dat2[i][7] #radius
        for j in range(3):
            new_dat[i,j+2]=dat2[i][j] #positions
            new_dat[i,j+5]=dat2[i][j+3] #velocities
    return new_dat,t
#-------------------------------------------------------------------------------
def N_largest2(m_lim,i,fname):
    '''Function to find the N largest particles from the data array, down to a particle mass limit'''
    t=numpy.loadtxt(fname,usecols=(0,)) # read only first column to get time
    if len(t)<1:
        print( "FILE {} IS EMPTY!".format(i))
        new_dat=numpy.zeros(0)
        t="skip"
        return new_dat,t
    dat=numpy.loadtxt(fname,skiprows=1) # read in the data skipping the first row
    t=t[0] # do some wizardry to extract the time value from the whole column
    #check that all particles are at least of mass m_lim?
    if numpy.amax(dat[:,6])<m_lim:
        print( "frame: {}, particles are too small".format(i))
        new_dat=numpy.zeros(0)
        t="skip"
        return new_dat,t
    #create structured array
    dat2 = numpy.core.records.fromarrays(dat.transpose(),names='x,y,z,vx,vy,vz,m,r',formats = 'f16,f16,f16,f16,f16,f16,f16,f16')
    dat2=numpy.sort(dat2, order='m')
    dat2=numpy.flipud(dat2)
    new_dat=numpy.zeros(8)
    new_dat[0]=dat2[0][6] #mass
    new_dat[1]=dat2[0][7] #radius
    for j in range(3):
        new_dat[j+2]=dat2[0][j] #positions
        new_dat[j+5]=dat2[0][j+3] #velocities
    temp=numpy.zeros(8)
    for i in range(1,len(dat2)):
        if dat2[i][6]<m_lim:
            break
        else:
            temp[0]=dat2[i][6] #mass
            temp[1]=dat2[i][7] #radius
            for j in range(3):
                temp[j+2]=dat2[i][j] #positions
                temp[j+5]=dat2[i][j+3] #velocities
            new_dat=numpy.vstack((new_dat,temp))
    if i==1:
        print( "only one large particle!")
        new_dat=numpy.zeros(0)
        t="skip"
        return new_dat,t
    return new_dat,t
#-------------------------------------------------------------------------------
def N_largest_mlim(m_ratio_lim,fname):
    '''Function to find the N largest particles from the data array, down to a particle mass limit.
    The mass limit is set by the minimum mass ratio of interest
    returns a dataframe'''

    #load from data file
    arrays = [numpy.array(map(float,line.split())) for line in open(fname)] #reads in line by line
    t=float(arrays[0]) #first line is time
    dat=numpy.array(arrays[1:]) #rest is the 8xN particle data array

    #create dataframe and sort
    df_dat=pd.DataFrame(dat,columns=['x(m)','y(m)','z(m)','vx(ms^-1)','vy(ms^-1)','vz(ms^-1)','m(kg)','r(m)'])
    df_dat=df_dat.sort_values('m(kg)',axis=0,ascending=False) #sort in descending order of mass
    m_lim=m_ratio_lim*numpy.amax(df_dat['m(kg)'])
    df_dat=df_dat[df_dat['m(kg)']>m_lim] #keep only masses greater than the mass limit

    return df_dat,t
#-------------------------------------------------------------------------------

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#FIX THESE TRANSFORMS!

def rotating_to_heliocentric(r,v,a,t):
    '''Function to transform from rotating frame (circular orbit of radius a, at R(t=0)=(a,0,0)) to the heliocentric frame'''
    Om_k=numpy.sqrt(G*(M_sun)/(a**3.0)) #angular velocity at semimajor axis a_k
    a_vec=a*numpy.array([cos(Om_k*t),sin(Om_k*t),0.0]) #time gives location of frame
    Om_vec=numpy.array([0.0,0.0,Om_k]) # define angular velocity vector
    theta=Om_k*t # angle of rotationa s a function of time
    rot_mat=numpy.array([cos(theta),-sin(theta),0,sin(theta),cos(theta),0,0,0,1]).reshape((3,3)) # z axis rotation vector
    R=a_vec+numpy.dot(rot_mat,r) # transform position
    #R=r+a_vec
    Om_skew=numpy.array([0,-Om_k,0,Om_k,0,0,0,0,0]).reshape((3,3)) #define skew symmetric matrix for angular velocity (0,0,Om_k)
    rot_mat_dot=numpy.dot(Om_skew,rot_mat) #time derivative of the rotation matrix
    V=numpy.cross(Om_vec,a_vec)+numpy.dot(rot_mat_dot,r)+numpy.dot(rot_mat,v) # transform velocity
    # V=v+numpy.cross(Om_vec,a_vec)#+numpy.cross(Om_vec,r) #velocity in rot frame + circular orbital velocity of rot frame + rotational velocity of rot frame
    #dot_a_vec=a*Om_k*numpy.array([-sin(Om_k*t),cos(Om_k*t),0.0])
    #dot_rot_mat=Om_k*numpy.array([-sin(theta),cos(theta),0,-cos(theta),-sin(theta),0,0,0,0]).reshape((3,3))
    #V_check=dot_a_vec+(numpy.dot(dot_rot_mat,r))+(numpy.dot(rot_mat,v))
    #V_check=dot_a_vec+numpy.cross(Om_vec,numpy.dot(rot_mat,r))
    # V_check=numpy.cross(Om_vec,a_vec)+numpy.cross(Om_vec,numpy.dot(rot_mat,r))+numpy.dot(rot_mat,v)
    #print "check V : {}".format((V-V_check)/numpy.linalg.norm(V))
    return R,V#_check

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def heliocentric_to_rotating(R,V,a,t):
    '''Function to transform from heliocentric frame to rotating frame (circular orbit of radius a, at R(t=0)=(a,0,0))'''
    Om_k=numpy.sqrt(G*(M_sun)/(a**3.0)) #angular velocity at semimajor axis a_k
    a_vec=a*numpy.array([cos(Om_k*t),sin(Om_k*t),0.0])
    Om_vec=numpy.array([0.0,0.0,Om_k])
    theta=-Om_k*t #reverse angle
    rot_mat=numpy.array([cos(theta),-sin(theta),0,sin(theta),cos(theta),0,0,0,1]).reshape((3,3))
    r=numpy.dot(rot_mat,R)-numpy.dot(rot_mat,a_vec)
    # r=numpy.dot(rot_mat,R)-a_vec
    # r=R-a_vec
    Om_skew=numpy.array([0,-Om_k,0,Om_k,0,0,0,0,0]).reshape((3,3)) #define skew symmetric matrix for angular velocity (0,0,Om_k)
    rot_mat_dot=-numpy.dot(Om_skew,rot_mat) #time derivative of the rotation matrix
    #v=V-numpy.cross(Om_vec,a_vec)#-numpy.cross(Om_vec,r)
    v=numpy.dot(rot_mat_dot,R)+numpy.dot(rot_mat,V)-numpy.dot(rot_mat_dot,a_vec)-numpy.dot(rot_mat,numpy.cross(Om_vec,a_vec))
    #dot_a_vec=a*Om_k*numpy.array([-sin(Om_k*t),cos(Om_k*t),0.0])
    #dot_rot_mat=Om_k*numpy.array([-sin(theta),cos(theta),0,-cos(theta),-sin(theta),0,0,0,0]).reshape((3,3))
    #v_check=numpy.dot(dot_rot_mat,(R-a_vec))+numpy.dot(rot_mat,(V-dot_a_vec))
    #v_check=(numpy.dot(rot_mat,V))-(numpy.dot(rot_mat,a*Om_k*numpy.array([-sin(-theta),cos(-theta),0])))-(numpy.dot(rot_mat,numpy.dot(dot_rot_mat,r)))
    # v_check=-numpy.cross(Om_vec,numpy.dot(rot_mat,(R-a_vec)))+numpy.dot(rot_mat,(V-numpy.cross(Om_vec,a_vec)))
    #v_check=numpy.cross(Om_vec,numpy.dot(rot_mat,(-R+a_vec)))-numpy.dot(rot_mat,(-V+numpy.cross(Om_vec,a_vec)))
    #print "check v : {}".format(v-v_check/numpy.linalg.norm(v))
    return r,v#_check

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    if len(r.shape)>1:
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
def find_orbits3(df_dat,i_o,inds,a_0,t,orb_fname,i,N):
    '''Function to find the orbits relative to a particle i_o. Orbit is an array: a, e, I, OMEGA, omega, f, E_J'''

    orbf = open(orb_fname, 'a') #open file to store orbits
    print "open file {}".format(orb_fname)
    pos=numpy.array(df_dat.loc[:,['x(m)','y(m)','z(m)']])
    vel=numpy.array(df_dat.loc[:,['vx(ms^-1)','vy(ms^-1)','vz(ms^-1)']])
    m=numpy.array(df_dat.loc[:,'m(kg)'])
    r=numpy.array(df_dat.loc[:,'r(m)'])
    #N=len(m)
    #sim = add_to_sim(dat)
    #sim.G=G
    oc=0
    for j in range(N):
        duplicate=0 #DOUBLE CHECK WHERE THIS SHOULD BE TO CATCH ALL ORBITS?
        if j!=i_o:
            #RELATIVE COORDINATES
            R=pos[j,:]-pos[i_o,:]
            V=vel[j,:]-vel[i_o,:]
            #V=(vel[j,:]+(1.5*OMEGA*pos[j,0]))-(vel[i_o,:]+(1.5*OMEGA*pos[i_o,0])) #NEED TO REMOVE SHEAR??!?!?!
            mu=G*(m[j]+m[i_o])
            orbit,f_true=find_elements(R,V,mu)
            orbit=numpy.append(orbit,f_true)

            #FIND BINARY COM ORBIT FOR E_J
            #rotating frame CoM position and velocity
            r_com=((m[j]*pos[j,:])+(m[i_o]*pos[i_o,:]))/(m[j]+m[i_o])
            v_com=((m[j]*vel[j,:])+(m[i_o]*vel[i_o,:]))/(m[j]+m[i_o])
            #transform to heliocentric frame CoM position and velocity
            R_com,V_com=rotating_to_heliocentric(r_com,v_com,a_0,t)
            #get a from orbit, give to function below so that we use the correct heliocentric angular velocity
            mu_com=G*(M_sun+m[j]+m[i_o])
            orbit_com,f_true_com=find_elements(R_com,V_com,mu_com)
            orbit_com=numpy.append(orbit_com,f_true_com)

            #print "com orbit position, velocity and mu: ",R_com,V_com,mu_com
            a_c=orbit_com[0]
            E_J=jacobi_rotating(pos[i_o,:],pos[j,:],vel[i_o,:],vel[j,:],m[i_o],m[j],a_c) #IS THIS FUNCTION CORRECT?

            orbit=numpy.append(orbit,E_J)
            #also output m1 and m2
            orbit=numpy.append(orbit,m[i_o])
            orbit=numpy.append(orbit,m[j])
            if (orbit[1]<1.0) and (orbit[1]>=0.0): # 0<e<1, an elliptical orbit exists
                oc+=1
                inds.append([i_o,j]) #store indices or orbital particles

                for k in range(len(inds)):
                    if inds[k][0]==inds[-1][1] and inds[k][1]==inds[-1][0]:
                        duplicate+=1
                        print "duplicate orbit"
                        break

                if duplicate<1:
                    #orbits.append(orbit) #Don't store duplicate orbits
                    '''
                    print "orbit: t={}\ti_o={}\ta={}".format(t,i_o,orbit[0])
                    orbf.write("{}\t{}\t".format(i_o,i)) #primary particle index and file number
                    for k in range(len(orbit)):
                        orbf.write("{}\t".format(orbit[k])) #orbits contains the 6 orbital elements, E_J and m1, m2
                    orbf.write("{}\t".format(t)) #record time
                    #also add the binary CoM orbit
                    for k in range(len(orbit_com)):
                        orbf.write("{}\t".format(orbit_com[k])) #orbits contains the 6 orbital elements, E_J and m1, m2
                    orbf.write("{}\t".format(f_true_com)) #record binary com heliocentric orbit
                    orbf.write("{}\t".format(j)) #record secondary particle index

                    orbf.write("\n")'''

                    orbf.write("{}\t".format(t)) #record time
                    orbf.write("{}\t{}\t{}\t".format(i_o,j,i)) #primary particle index and file number
                    for k in range(len(orbit)):
                        orbf.write("{}\t".format(orbit[k])) #orbits contains the 6 orbital elements, E_J and m1, m2
                    #also add the binary CoM orbit
                    for k in range(len(orbit_com)):
                        orbf.write("{}\t".format(orbit_com[k])) #helio orbits contains the 6 orbital elements

                    orbf.write("\n")

                    print "suitable orbit found between {} and {}".format(i_o,j)
                    '''print "t = {}".format(t)
                    print "m1 = {}, m2 = {}".format(m[j],m[i_o])
                    print "pos_j = {}\npos_io = {}".format(pos[j,:],pos[i_o,:])
                    print "vel_j = {}\nvel_io = {}".format(vel[j,:],vel[i_o,:])
                    print "r_com = {}\nR_com = {}\nv_com = {}\nV_com = {}".format(r_com,R_com,v_com,V_com)
                    print "com orbit: {}, E_J={}\n".format(orbit_com,E_J)'''
            else:
                print "no suitable orbits for {} and {}".format(i_o,j)

    orbf.close()
    return
#-------------------------------------------------------------------------------

def orbit_search(m_ratio_lim,files,run_path,dirname,t_quit,N_lim):
    '''Function to search the data files for orbits, outputs to a file'''
    directory=run_path+"/"+dirname

    #Read in from run_params
    #run_params=numpy.genfromtxt(directory+"/"+"run_params_0.txt",skip_header=1,dtype=None) #genfromtxt to avoid problems with string 'new' in run_params
    run_params=numpy.genfromtxt(directory+"/"+"run_params_0.txt",dtype=None) #genfromtxt to avoid problems with string 'new' in run_params
    print "run parameters:\n {}".format(run_params)
    a_0=float(run_params[2])
    OMEGA=float(run_params[7])
    Req=float(run_params[3])
    rho=float(run_params[4])
    X=float(run_params[8])
    Mtot=(4.0/3.0)*pi*rho*(Req**3.0)
    print "a0={},OM={},Req={},X={},Mtot={}".format(a_0,OMEGA,Req,X,Mtot)

    #DEFINE MLIM RELATIVE TO MTOT?

    print "analyse {}".format(directory)
    orb_fname = directory+"/"+dirname+'_orbits.txt'

    orbf1 = open(orb_fname, 'w') #open file to store orbits

    files=list(reversed(files))
    for fi in files:
        i=int(fi[3:10])

        fi2=directory+"/"+fi
        df_dat,t=N_largest_mlim(m_ratio_lim,fi2)

        if len(df_dat)<=1:
            print "the particles are too small or there is only one large particle"
            break
        if t<t_quit:
            print "{} analysed out of total time {}".format(t,t_quit)
            break

        N=len(df_dat)
        print "N={}".format(N)
        if N>N_lim:
            print "N = {} > N_lim = {}".format(N,N_lim)
            N=N_lim
        print "N={}".format(N)
        i_M=0 #the central particle is taken to be the most massive particle, index 0 after sorting

        print("\nrun {}, frame: {}, time: {}, number of particles: {}\n".format(dirname,fi,t,N))

        inds=[]

        #go through the particle list to find all orbits relative to that particle
        for i_o in range(N):
            print("find orbit relative to particle: {} (out of {})".format(i_o,N))
            find_orbits3(df_dat,i_o,inds,a_0,t,orb_fname,i,N)

    #close files
    orbf1.close()
#-------------------------------------------------------------------------------

def angle_fixer180(x):
    '''Function to take an angle x radians and ensure it lies between 0 and pi'''
    x=angle_fixer360(x) #first ensure it lies between 0 and 2pi
    if x>pi:
        while x>pi:
            print( 'too big')
            x=(2*pi)-x
    return x
#-------------------------------------------------------------------------------
def angle_fixer360(x):
    '''Function to take an angle x radians and ensure it lies between 0 and 2pi'''
    if x<0:
        while x<0:
            x+=(2*pi)
    if x>(2*pi):
        while x>(2*pi):
            x-=(2*pi)
    return x

#-------------------------------------------------------------------------------
def find_elements(R,V,mu):
    '''Function to find orbital elements (a,e,I,OMEGA,omega,ftrue) from position and velocity, R in m, V in m/s, angles in radians
    This is better implemented within rebound.
    Solution from Astromechanics.pdf, see also Murray and Dermott SSD sec.2.8 '''
    #print mu
    dat=numpy.zeros(5)
    #basis vectors
    ex=numpy.array([1.0,0,0])
    ey=numpy.array([0,1.0,0])
    ez=numpy.array([0,0,1.0])
    #position and velocity magnitude
    r=numpy.linalg.norm(R)
    v=numpy.linalg.norm(V)
    #Energy and semimajor axis
    e=(v**2/2)-(mu/r) #redefine this for clarity?
    a=-mu/(2*e)
    dat[0]=a
    #angular momentum
    H=numpy.cross(R,V)
    h=numpy.linalg.norm(H)
    #eccentricity
    E=(1/mu)*(numpy.cross(V,H)-(mu*(R/r)))
    e=numpy.linalg.norm(E)
    dat[1]=e
    #Inclination
    A=numpy.dot(H,ez)/h
    I=arccos(A)
    I=angle_fixer180(I)
    if A>0 and I>(pi/2):
        I=(2*pi)-I
    if A<0 and I<(pi/2):
        I=(2*pi)-I
    dat[2]=I #store as radians
    #nodal vector
    N=numpy.cross(ez,H)
    n=numpy.linalg.norm(N)
    #longitude of ascending node
    B=numpy.dot(N,ex)/n
    OMEGA=arccos(B)
    OMEGA=angle_fixer360(OMEGA)
    ny=numpy.dot(N,ey)
    if ny>0 and OMEGA>pi:
        OMEGA=(2*pi)-OMEGA
    if ny<0 and OMEGA<pi:
        OMEGA=(2*pi)-OMEGA
    dat[3]=OMEGA
    #argument of pericentre
    C=numpy.dot(N,E)/(n*e)
    omega=arccos(C)
    omega=angle_fixer360(omega)
    e_z=numpy.dot(E,ez)
    if e_z>0 and omega>pi:
        omega=(2*pi)-omega
    if e_z<0 and omega<pi:
        omega=(2*pi)-omega
    dat[4]=omega
    #true anomaly
    D=numpy.dot(R,E)/(r*e)
    f_true=arccos(D)
    f_true=angle_fixer360(f_true)
    RV=numpy.dot(R,V)
    if RV>0 and f_true>pi:
        f_true=(2*pi)-f_true
    if RV<0 and f_true<pi:
        f_true=(2*pi)-f_true
    return dat,f_true

#-------------------------------------------------------------------------------
def jacobi_rotating(r1,r2,v1,v2,m1,m2,a_c):
    '''find jacobi integral of orbit of particle 2, wrt particle 1 (Hill approxiamtion, rotating frame: Nakazawa 1988 15.6.1)
    CHECK THIS'''
    p_rel=r2-r1 #relative position
    v_rel=v2-v1 #relative velocity
    r=numpy.linalg.norm(p_rel)
    OM_k=numpy.sqrt(G*(M_sun)/(a_c**3.0)) #angular velocity at CoM semimajor axis a_c
    M=m1+m2
    r_H=pow(M/(3.0*M_sun),1.0/3.0)*(a_c) #mutual Hill radius
    A=0.5*((v_rel[0]**2.0)+(v_rel[1]**2.0)+(v_rel[2]**2.0))
    B=-(1.5*((OM_k*p_rel[0])**2.0))+(0.5*((OM_k*p_rel[2])**2.0))
    C=-(G*(M)/r)+(4.5*((r_H*OM_k)**2.0))
    E_J=A+B+C #Jacobi Energy Integral
    #print "r_H = {}, sep = {}, sep/r_H = {}, E_J = {}".format(r_H,r,r/r_H,E_J)
    return E_J

#-------------------------------------------------------------------------------
def add_to_sim(dat):
    '''creates a simulation and adds the particles'''
    sim = rebound.Simulation()
    pos=dat[:,2:5]
    vel=dat[:,5:]
    m=dat[:,0]
    r=dat[:,1]
    for i in range(len(dat[:,0])):
        sim.add(m=m[i],r=r[i],x=pos[i,0],y=pos[i,1],z=pos[i,2],vx=vel[i,0],vy=vel[i,1],vz=vel[i,2])
    return sim

def add_file_to_sim(sim,fname):
    '''loads particles from a file as a pandas dataframe, ands adds them to an existing rebound sim'''
    t,df_dat=load_dat_file(fname)
    pos=numpy.array(df_dat.loc[:,['x(m)','y(m)','z(m)']])
    vel=numpy.array(df_dat.loc[:,['vx(ms^-1)','vy(ms^-1)','vz(ms^-1)']])
    m=numpy.array(df_dat.loc[:,'m(kg)'])
    r=numpy.array(df_dat.loc[:,'r(m)'])
    N=len(m)
    for i in range(N):
        sim.add(m=m[i],r=r[i],x=pos[i,0],y=pos[i,1],z=pos[i,2],vx=vel[i,0],vy=vel[i,1],vz=vel[i,2])
    #return sim

#-------------------------------------------------------------------------------
def find_orbits(N,i,fname2,a_0):
    orbits=[]
    dat,t=N_largest(N,i,fname2)
    i_M=0 #largest particle
    #find the closest particle to most massive particle 0
    dist=[]
    for l in range(N-1):
        n=l+1
        R=numpy.linalg.norm(dat[n,2:5]-dat[0,2:5])
        dist.append(R)
    i_R=numpy.argmin(dist)+1
    dat_keep=numpy.vstack((dat[0,:],dat[i_R,:]))
    #dat=dat_keep #use only these two particles

    pos=dat[:,2:5]
    vel=dat[:,5:]
    m=dat[:,0]
    r=dat[:,1]

    print("closest particle: ",i_R)
    print( "t: ",t)
    print( dat)
    print( 'plot ',i)
    print( 'time ',t)

    #find CoM
    xm=pos[:,0]*m
    ym=pos[:,1]*m
    zm=pos[:,2]*m
    M=numpy.sum(m)
    CoM=numpy.zeros(3)
    CoM[0]=numpy.sum(xm)/M
    CoM[1]=numpy.sum(ym)/M
    CoM[2]=numpy.sum(zm)/M

    print( "masses: ",m)
    print( "mass ratio: ",m[1]/m[0])
    print( "radii: ",r)
    print( "radii ratio: ",r[1]/r[0])
    print( "CoM: ",CoM)

    #find orbits using rebound
    sim = rebound.Simulation()

    '''
    numpy.savetxt("reb_orb_dat.txt",dat)
    with open('reb_orb_dat.txt', 'r') as op:
        s = op.read()
    sim.add_particles_ascii(s)
    '''
    for i in range(N):
        sim.add(m=m[i],r=r[i],x=pos[i,0],y=pos[i,1],z=pos[i,2],vx=vel[i,0],vy=vel[i,1],vz=vel[i,2])
        if i>0:
            E_J=jacobi_rotating(pos[i_M,:],pos[i,:],vel[i_M,:],vel[i,:],m[i_M],m[i],a_0)
            print( "E_J={}".format(E_J))
    N=sim.N
    sim.G=G
    #print "Gravitational constant, G = ",sim.G
    #print "\n"
    orbs = sim.calculate_orbits(heliocentric=True)
    for orbit in orbs:
        if orbit.e<1.0:
            print(orbit)
            orbits.append([orbit.a,orbit.e,orbit.inc,orbit.Omega,orbit.omega,orbit.f])
        else:
            print( "not elliptical (0<e<1)")

    #record orbital elements of closest particle
    i_R-=1 #reindex as orbital elements do not contain particle 0
    a=(orbits[i_R][0])
    e=(orbits[i_R][1])
    I=(orbits[i_R][2])
    time=(t)
    #print "elements: \n",a,time

    return orbits,a,e,I,time
#-------------------------------------------------------------------------------
def find_orbits2(dat,i_o,inds,a_0):
    '''Function to find the orbits relative to a particle i_o. Orbit is an array: a, e, I, OMEGA, omega, f, E_J'''
    orbits=[]
    orb=numpy.zeros(6)
    N=len(dat[:,0])
    pos=dat[:,2:5]
    vel=dat[:,5:]
    m=dat[:,0]
    r=dat[:,1]
    sim = add_to_sim(dat)
    sim.G=G
    oc=0
    for j in range(N):
        duplicate=0 #DOUBLE CHECK WHERE THIS SHOULD BE TO CATCH ALL ORBITS?
        if j!=i_o:
            R=pos[j,:]-pos[i_o,:]
            V=vel[j,:]-vel[i_o,:]
            #V=(vel[j,:]+(1.5*OMEGA*pos[j,0]))-(vel[i_o,:]+(1.5*OMEGA*pos[i_o,0])) #NEED TO REMOVE SHEAR??!?!?!
            mu=G*(m[j]+m[i_o])
            orbit,f_true=find_elements(R,V,mu)
            orbit=numpy.append(orbit,f_true)
            E_J=jacobi_rotating(pos[i_o,:],pos[j,:],vel[i_o,:],vel[j,:],m[i_o],m[j],a_0)
            orbit=numpy.append(orbit,E_J)
            #also output m1 and m2
            orbit=numpy.append(orbit,m[i_o])
            orbit=numpy.append(orbit,m[j])
            if orbit[1]<1.0: #an elliptical orbit exists
                #print "particle: {}, {:.11e} {:.11e} {:.11e} E_J={}".format(j,orbit[0],orbit[1],orbit[2],E_J)
                oc+=1
                inds.append([i_o,j]) #store indices or orbital particles
                for k in range(len(inds)):
                    if inds[k][0]==inds[-1][1] and inds[k][1]==inds[-1][0]:
                        duplicate+=1
                        #print "DUPLICATE ORBIT: {} {}".format(i_o,j)
                        break
                if duplicate<1:
                    #print "APPEND"
                    orbits.append(orbit) #Don't store duplicate orbits
            #else:
                #print "particle: {}, not elliptical".format(j)

    #print "number of orbits e<1: {}".format(oc)
    return orbits
#-------------------------------------------------------------------------------
def orbital_plotting(fig,runname,fname,dat_file,i,t,pos,CoM,lim,lim2,i_M,orbits):
    '''Function containing the plotting routines for orbits.py (figure already open)'''
    pyplot.rc('text', usetex=True)
    pyplot.rc('font', family='serif')
    fig.suptitle('run: {}, frame: {}, particle: {}, time: {:.2e} s'.format(runname,dat_file,i_M,t))
    plot_option=3
    if plot_option==0: # xy plane
        #plot only the N largest particles as simple markers
        print( "plot option: ",plot_option)
        fig.set_size_inches(15.5, 10.5)
        gs = gridspec.GridSpec(1, 2,width_ratios=[1,1],height_ratios=[1,1])
        ax1 = pyplot.subplot(gs[0,0])
        ax2 = pyplot.subplot(gs[0,1])
        ax1.set_aspect("equal")
        ax1.set_xlabel('x /m')
        ax1.set_ylabel('y /m')
        ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax1.set_xlim(-lim,lim)
        ax1.set_ylim(-lim,lim)
        ax2.set_aspect("equal")
        ax2.set_xlabel('x /m')
        ax2.set_ylabel('y /m')
        ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax2.set_xlim(-lim2,lim2)
        ax2.set_ylim(-lim2,lim2)

        ax1.scatter(pos[:,0],pos[:,1])
        ax1.scatter(pos[i_M,0],pos[i_M,1],marker='x',c='r')
        ax1.scatter(CoM[0],CoM[1],marker='x',c='k')
        ax2.scatter(pos[:,0]-pos[i_M,0],pos[:,1]-pos[i_M,1])

    if plot_option==2: # xy plane
        #simple N particles and orbit plot. Single Panel
        print( "plot option: ",plot_option)
        ax1 = fig.add_subplot(111)
        ax1.set_aspect("equal")
        lim2=1e9
        ax1.set_xlim(-lim2,lim2)
        ax1.set_ylim(-lim2,lim2)
        ax1.scatter(pos[:,0]-pos[i_M,0],pos[:,1]-pos[i_M,1])
        ax1.scatter(CoM[0]-pos[i_M,0],CoM[1]-pos[i_M,1],marker='x',c='k')
        for l in range(len(orbits[:,0])):
            pos_orb=planet_orbit(orbits[l][:],1000)
            ax1.plot(pos_orb[:,0],pos_orb[:,1])

    if plot_option==3: # xy plane
        #simple N particles and orbit plot. Double panel
        print( "plot option: ",plot_option)
        fig.set_size_inches(15.5, 10.5)
        gs = gridspec.GridSpec(1, 2,width_ratios=[1,1],height_ratios=[1,1])
        ax1 = pyplot.subplot(gs[0,0])
        ax2 = pyplot.subplot(gs[0,1])
        ax1.set_aspect("equal")
        ax2.set_aspect("equal")
        ax1.set_xlabel('x /m')
        ax1.set_ylabel('y /m')
        ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax2.set_xlabel('x /m')
        ax2.set_ylabel('y /m')
        ax2.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        lim1=1e9
        ax1.set_xlim(-lim1,lim1)
        ax1.set_ylim(-lim1,lim1)

        lim2=1e7
        ax2.set_xlim(-lim2,lim2)
        ax2.set_ylim(-lim2,lim2)

        ax1.scatter(pos[:,0]-pos[i_M,0],pos[:,1]-pos[i_M,1])
        ax2.scatter(pos[:,0]-pos[i_M,0],pos[:,1]-pos[i_M,1])

        for l in range(len(orbits[:])):
            pos_orb=planet_orbit(orbits[l],1000)
            ax1.plot(pos_orb[:,0],pos_orb[:,1])
            ax2.plot(pos_orb[:,0],pos_orb[:,1])

    #save the figure
    picname="./"+fname+str(plot_option)+"/"+fname+str(plot_option)+"min_%07d.png"%i
    pyplot.savefig(picname, bbox_inches='tight')
#-------------------------------------------------------------------------------
def plot_elements(a,e,I,time):
    '''Fucntion to plot orbital elements with time'''
    fig = pyplot.figure()
    gs = gridspec.GridSpec(1, 3)
    ax1 = pyplot.subplot(gs[0,0])
    ax2 = pyplot.subplot(gs[0,1])
    ax3 = pyplot.subplot(gs[0,2])
    ax1.set_xlabel('t /s')
    ax1.set_ylabel('a /m')
    ax2.set_xlabel('t /s')
    ax2.set_ylabel('e')
    ax3.set_xlabel('t /s')
    ax3.set_ylabel('I /rad')
    ax1.plot(time,a)
    ax2.plot(time,e)
    ax3.plot(time,I)
    pyplot.show()
#-------------------------------------------------------------------------------
def orbit_plot3d(fig,pos,CoM,lim,orbits):
    '''Function to plot positions and orbits in 3d'''
    ax1 = fig.add_subplot(111,projection='3d')
    ax1.scatter(0,0,0)
    ax1.axis('equal') # set aspect ratio to equal

    lim=1e7
    ax1.set_xlim(-lim,lim)
    ax1.set_ylim(-lim,lim)
    ax1.set_zlim(-lim,lim)

    ax1.scatter(pos[:,0]-pos[i_M,0],pos[:,1]-pos[i_M,1],pos[:,2]-pos[i_M,2])
    ax1.scatter(CoM[0]-pos[i_M,0],CoM[1]-pos[i_M,1],CoM[2]-pos[i_M,2],marker='x',c='k')
    ax1.plot(pos[:,0]-pos[i_M,0],pos[:,1]-pos[i_M,1],pos[:,2]-pos[i_M,2])

    for l in range(len(orbits[:,0])):
        pos_orb=planet_orbit(orbits[l][:],1000)
        ax1.plot(pos_orb[:,0],pos_orb[:,1])
#-------------------------------------------------------------------------------
def orbit_plot2d(fig,ax1,pos,orbits,i_M,i_o,lim):
        '''simple N particles and orbit plot. Single Panel. i_M is the central particle, i_o is the orbit reference particle'''
        #print "i_M={}, i_o={}".format(i_M,i_o)
        #ax1 = fig.add_subplot(111)
        ax1.set_aspect("equal")
        ax1.set_xlim(-lim,lim)
        ax1.set_ylim(-lim,lim)

        ax1.set_xlabel('x /m')
        ax1.set_ylabel('y /m')
        ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

        pos_rel=pos-pos[i_M,:]
        ax1.scatter(pos_rel[:,0],pos_rel[:,1]) #plot positions relative to central particle
        for l in range(len(orbits)):
            pos_orb=planet_orbit(orbits[l][:],1000)
            pos_orb=pos_orb+pos_rel[i_o,:] #move orbits to relative to particle i_o
            if i_o==i_M:
                col="r"
            else:
                col="b"
            if orbits[l][6]<0.0: #if E_J is negative
                ls="-"
            else:
                ls=":"
            ax1.plot(pos_orb[:,0],pos_orb[:,1],linestyle=ls,c=col)

#-------------------------------------------------------------------------------
def all_dat_sort():
    path='../completed_runs/'

    all_runs=numpy.genfromtxt(path+"all_runs_essential.txt",dtype=None,skip_header=1)
    all_orbits=numpy.genfromtxt(path+"all_orbits.txt",dtype=None,skip_header=1)
    for i in range(len(all_runs)):
        run=all_runs[i][0]
        run=int(run[:3])
        all_runs[i][0]=run

    inds1=[int(item[0]) for item in all_runs]
    inds2=[int(item[0]) for item in all_orbits]

    dat=numpy.zeros((len(all_orbits),10))
    initial=[]

    for i in range(len(all_runs)):
        if int(all_runs[i][0]) in inds2:
            j=inds2.index(int(all_runs[i][0]))
            dat[j,0]=int(all_runs[i][0]) #run number
            dat[j,1]=float(all_runs[i][1]) #N
            dat[j,2]=all_runs[i][2] #X
            dat[j,3]=all_runs[i][3] #f
            dat[j,4]=all_runs[i][4] #dt
            dat[j,5]=all_orbits[j][1] #m1
            dat[j,6]=all_orbits[j][2] #m2m1
            dat[j,7]=all_orbits[j][3] #a
            dat[j,8]=all_orbits[j][4] #e
            dat[j,9]=all_orbits[j][5] #I
            initial.append(all_runs[i][7]) #initial

    #create structured data array for sorting
    dat_sort = numpy.core.records.fromarrays(dat.transpose(),names='run,N,X,f,dt,m1,m2m1,a,e,I',formats = 'i8,f8,f8,f8,f8,f8,f8,f8,f8,f8')
    #print( "dat_sort:\n",dat_sort)

    #print 'initial conditions',initial
    dat_sort=numpy.lib.recfunctions.append_fields(dat_sort,'int',data=initial)

    return all_runs,all_orbits,dat_sort
#-------------------------------------------------------------------------------
def all_dat_sort2():
    #no path
    all_runs=numpy.genfromtxt("all_runs_essential.txt",dtype=None,skip_header=1)
    all_orbits=numpy.genfromtxt("all_orbits.txt",dtype=None,skip_header=1)
    for i in range(len(all_runs)):
        run=all_runs[i][0]
        run=int(run[:3])
        all_runs[i][0]=run

    inds1=[int(item[0]) for item in all_runs]
    inds2=[int(item[0]) for item in all_orbits]

    dat=numpy.zeros((len(all_orbits),10))
    initial=[]

    for i in range(len(all_runs)):
        if int(all_runs[i][0]) in inds2:
            j=inds2.index(int(all_runs[i][0]))
            dat[j,0]=int(all_runs[i][0]) #run number
            dat[j,1]=float(all_runs[i][1]) #N
            dat[j,2]=all_runs[i][2] #X
            dat[j,3]=all_runs[i][3] #f
            dat[j,4]=all_runs[i][4] #dt
            dat[j,5]=all_orbits[j][1] #m1
            dat[j,6]=all_orbits[j][2] #m2m1
            dat[j,7]=all_orbits[j][3] #a
            dat[j,8]=all_orbits[j][4] #e
            dat[j,9]=all_orbits[j][5] #I
            initial.append(all_runs[i][7]) #initial

    #create structured data array for sorting
    dat_sort = numpy.core.records.fromarrays(dat.transpose(),names='run,N,X,f,dt,m1,m2m1,a,e,I',formats = 'i8,f8,f8,f8,f8,f8,f8,f8,f8,f8')
    #print( "dat_sort:\n",dat_sort)

    #print 'initial conditions',initial
    dat_sort=numpy.lib.recfunctions.append_fields(dat_sort,'int',data=initial)

    return all_runs,all_orbits,dat_sort
#-------------------------------------------------------------------------------
def hist_dist(x,n_bins):
    '''Function to find the arrays containing histogram data

    CHECK LENGTHS

    Parameters
    ----------
    x
        1d array of the data set
    n_bins
        number of bins

    Returns
    -------
    center
        array of the center positions of the bins, length=n_bins
    hist
        the number of points in each bin, length=n_bins
    width
        array of the bin widths, length=n_bins-1'''
    hist,bins=numpy.histogram(x,bins=n_bins)
    center = (bins[:-1] + bins[1:]) / 2
    width = numpy.diff(bins)
    return center, hist, width
#-------------------------------------------------------------------------------
def find_time_from_file(dirpath):
    '''Function to find the most recent time stamp in a run directory'''
    file_times=subprocess.check_output(['ls -alT {}/*bin'.format(dirpath)], shell=True)
    file_lines=file_times.split("\n")
    file_lines=file_lines[:-1]
    line=file_lines[-1]
    line2=line.split(" ")
    fname=(line2[-1][:-3])+"txt"
    day=line2[-5]
    month=line2[-4]
    time=line2[-3]
    year=line2[-2]
    ftime=datetime.datetime.strptime("{}/{}/{} {}".format(month,day,year,time),"%b/%d/%Y %H:%M:%S")
    #get time from file
    t=numpy.loadtxt(fname,usecols=(0,)) # read only first column to get time
    dat=numpy.loadtxt(fname,skiprows=1) # read in the data skipping the first row
    t=t[0]
    return ftime,t_progress
#-------------------------------------------------------------------------------
def find_time_from_file_year(dirpath,year):
    '''Function to find the most recent time stamp in a run directory'''
    file_times=subprocess.check_output(['ls -alt {}/*bin'.format(dirpath)], shell=True)
    file_lines=file_times.split("\n")
    file_lines=file_lines[:-1]
    line=file_lines[0]
    print line
    line2=line.split(" ")
    print line2
    fname=(line2[-1][:-3])+"txt"
    day=line2[-4]
    month=line2[-3]
    time=line2[-2]
    print day,month,time
    try:
        ftime=datetime.datetime.strptime("{}/{}/{} {}".format(month,day,year,time),"%b/%d/%Y %H:%M")
    except:
        ftime=datetime.datetime.strptime("{}/{}/{} {}".format(month,day,year,time),"%d/%b/%Y %H:%M")

    #get time from file
    #print fname
    #print ftime.year
    t=numpy.loadtxt(fname,usecols=(0,)) # read only first column to get time
    dat=numpy.loadtxt(fname,skiprows=1) # read in the data skipping the first row
    t=t[0]
    return ftime,t
#-------------------------------------------------------------------------------
def load_collision_file(f):
    '''Function to load the collisions_0.txt file as a dataframe'''
    arrays = [numpy.array(map(float,line.split())) for line in open(f)] #reads in line by line
    coll_dat=numpy.array(arrays) #rest is the 8xN particle data array

    try: #this will load for reb_funcs.h: reb_collision_resolve_merge_out_damp_cut
        df_coll=pd.DataFrame(coll_dat,columns=['t(s)','i','j','mi(kg)','mj(kg)','ri(m)','rj(m)',
        'xi(m)','yi(m)','zi(m)','xj(m)','yj(m)','zj(m)','vxi(ms^-1)','vyi(ms^-1)',
        'vzi(ms^-1)','vxj(ms^-1)','vyj(ms^-1)','vzj(ms^-1)','m(kg)','r(m)','x(m)','y(m)',
        'z(m)','vx(ms^-1)','vy(ms^-1)','vz(ms^-1)','N_boulder'])
    except: #this will load for reb_funcs.h: reb_collision_resolve_merge_out
        df_coll=pd.DataFrame(coll_dat,columns=['t(s)','i','j','mi(kg)','mj(kg)',
        'xi(m)','yi(m)','zi(m)','xj(m)','yj(m)','zj(m)','vxi(ms^-1)','vyi(ms^-1)',
        'vzi(ms^-1)','vxj(ms^-1)','vyj(ms^-1)','vzj(ms^-1)','v_rel(ms^-1)'])

    return df_coll
#-------------------------------------------------------------------------------
def create_file_list(d):
    '''Function to create the list of dat.txt files'''
    files=next(os.walk(d))[2] #retrieve the files in the run directory
    files = [ fi for fi in files if fi.endswith(".txt") and fi.startswith("dat") ] #keep only dat*.txt files
    files.sort() #ensure that the files are always sorted the same?
    return files
#-------------------------------------------------------------------------------

def create_dir_list(d):
    '''Function to create the list of directories at location d'''
    dirs=next(os.walk(d))[1]
    dirs.sort() #ensure that the files are always sorted the same?
    return dirs

#-------------------------------------------------------------------------------
def load_dat_file(f):
    '''Load dat.txt file and return a time and a dataframe containing all particle data

    Parameters
    ----------
    f
        path to data file, formatted with a time stamp on the first line, and then x,y,z,vx,vy,vz,m,r; one line per particle
    '''
    arrays = [numpy.array(map(float,line.split())) for line in open(f)] #reads in line by line
    t=float(arrays[0]) #first line is time
    dat=numpy.array(arrays[1:])

    '''N=len(dat[:,0])
    pos=dat[:,:3]
    vel=dat[:,3:6]
    m=dat[:,6]
    r=dat[:,7]
    return pos'''

    df_dat=pd.DataFrame(dat,columns=['x(m)','y(m)','z(m)','vx(ms^-1)','vy(ms^-1)','vz(ms^-1)','m(kg)','r(m)'])
    return t,df_dat

#-------------------------------------------------------------------------------
def load_orb_file(f):
    '''Function to load orbits file as a dataframe'''
    arrays = [numpy.array(map(float,line.split())) for line in open(f)] #reads in line by line
    orb_dat=numpy.array(arrays)

    try:
        df_orb=pd.DataFrame(orb_dat,columns=['t(s)','i_o','j','file','a(m)','e','I(rad)',
        'omega(rad)','OMEGA(rad)','f_true(rad)','E_J(J)','m1(kg)','m2(kg)','a_hel(m)',
        'e_hel','I_hel(rad)','omega_hel(rad)','OMEGA_hel(rad)','f_true_hel(rad)'])
    except:
        df_orb=pd.DataFrame(orb_dat,columns=['t(s)','file','i','j','a(m)','e','I(rad)',
        'omega(rad)','OMEGA(rad)','f_true(rad)','m1(kg)','m2(kg)'])

    return df_orb
#-------------------------------------------------------------------------------
def load_run_params(file_path):
    '''
    Load the run_params file as a pandas dataframe
    ----------
    file_path
        path to run_params txt file
    '''

    params_columns=['run','run_dir','N_tot','a_orb(m)','R_eq(m)','rho(kgm-3)',
    'M_tot(kg)','R_c(m)','OM_orb(s-1)','X','OM_circ(s-1)','f','dt(s)','t_run(s)',
    'n_cores','initial']
    df_params= pd.DataFrame([],columns=params_columns)
    # params=numpy.genfromtxt(file_path,dtype=None) # this lines causes warnings
    with open(file_path) as f:
        content = f.readlines()
    params = [x.strip() for x in content]
    try:
        d=int(params[0].split('/')[-1].split('_')[0]) #load a run with run number at front of dirname
    except:
        d=int(params[0].split('/')[-1].split('_')[-1]) #load a run with run number at end of dirname
    params=numpy.insert(params,0,0) #insert placeholder for run number
    df_params_add=pd.DataFrame([params],columns=params_columns)
    df_params_add['run']=d
    df_params=df_params.append(df_params_add,ignore_index=True)

    return df_params

#-------------------------------------------------------------------------------
def file_stats(d,files,n_samp):
    '''
    Function to find the file stats of a particular run. Returns a dataframe
    containing file number, number of particles, file size (MB) and simulation time (s)

    Parameters
    ----------
    d
        path to run directory where files are kept, str
    files
        list of files to analyse in directory d, str
    n_samp
        Sampling interval to use when loading files, int
    '''
    f_num=[]
    f_N=[]
    f_s_MB=[]
    f_t=[]
    for i in range(0,len(files),n_samp):
        f=files[i]
        f_path='{}/{}'.format(d,f)
        num=int(f.split('_')[0][3:])
        s_b=float(os.path.getsize(f_path))
        s_MB=s_b/(1024.0**2.0)
        N=int(subprocess.check_output(['wc -l {}'.format(f_path)], shell=True).split()[0])-1 #subtract 1 to deal with header line
        f_num.append(num)
        f_N.append(N)
        f_s_MB.append(s_MB)
        with open(f_path, 'r') as _f:
            first_line = _f.readline().strip()
        t=(float(first_line))
        f_t.append(t)
        print '{} {} {} {}'.format(num,N,s_MB,t)
    fstats=numpy.stack((numpy.array(f_num),numpy.array(f_N),numpy.array(f_s_MB),numpy.array(f_t)),axis=1)
    print fstats.shape
    cols=['num','N','s_MB','t(s)']
    df_stats=pd.DataFrame(fstats,columns=cols)

    return df_stats
#-------------------------------------------------------------------------------

def read_dir_sizes(f_size):
    '''
    Function to read run directory sizes, created in the main directory with:
    $ du -sk ./*/ > run_dir_sizes_K.txt
    Returns a dataframe of size in bytes (binary 1024 scale) and the run number

    Parameters
    ----------
    f_size
        The path to the file run_dir_sizes.txt
    '''
    df_sizes=pd.read_csv(f_size,sep='\t',header=None,names=['size(B)','run']) #load txt file as a dataframe
    for i in range(len(df_sizes)): #for each entry in dataframe
        size=df_sizes.iloc[i,0]
        run=df_sizes.iloc[i,1].split('/')[1].split('_')[0]
        df_sizes.iloc[i,0]=size*1024 #replace the dataframe entries with the cleaned values
        df_sizes.iloc[i,1]=int(run)
    return df_sizes
#-------------------------------------------------------------------------------
def read_problem_file(run_path):
    '''
    create a dataframe of the run details by reading problem.c files. Extracts:
    run, seed, R_eq(m), X, f, R_b(m)

    If the problem.c file does not contain the 'reb_collision_resolve_merge_out_damp_cut' routine,
    the value of R_b(m) will be NaN

    Parameters
    ----------
    run_path
        The path to the directory containing each run as a sub directory
    '''

    jobs=next(os.walk(run_path))[1] #retrieve the files in the run directory
    jobs.sort() #ensure that the files are always sorted the same?
    dat,_num,_seed,_R_eq,_X,_f,_R_b,_N_tot,_rho,_t_max,_a_orb=[],[],[],[],[],[],[],[],[],[],[]
    for j in jobs:
        num=j[:3]
        damp_check=0
        for line in open('{}/{}/problem.c'.format(run_path,j)):
            #test for certain strings in line, then use this to read value
            if 'double X=' in line:
                X=line.split('=')[-1][:-2]
            if 'double _f=' in line:
                f=line.split('=')[-1][:-2]
            if 'rp->R_eq=' in line:
                R_eq=line.split('=')[-1][:-2]
            if 'char ins[64]=' in line:
                seed=line.split('=')[-1][:-2][1:-1] #cut the quotation marks
            if 'double r_boulder=' in line:
                R_b=line.split('=')[-1][:-2]
            if 'reb_collision_resolve_merge_out_damp_cut' in line: #check which collision routine is used
                damp_check=1
            if 'double Ntot=' in line:
                N_tot=line.split('=')[-1].split(';')[0]
            if 'rp->rho=' in line:
                rho=line.split('=')[-1].split(';')[0]
            if 'static double t_max=' in line:
                t_max=line.split('=')[-1].split(';')[0]
            if 'rp->a=' in line:
                a_orb=line.split('=')[-1].split(';')[0]

        if damp_check==0: #if the damped routine is not used we return a value of NaN for R_b
            R_b='NaN'
        _num.append(num)
        _seed.append(seed)
        _R_eq.append(float(R_eq))
        _X.append(float(X))
        _f.append(float(f))
        _R_b.append(float(R_b))
        _N_tot.append(float(N_tot))
        _rho.append(float(rho))
        _t_max.append(t_max)
        _a_orb.append(a_orb)
        # dat.append([num,seed,R_eq,X,f])

    n_jobs=len(_t_max)
    _M_tot=numpy.zeros(n_jobs)
    for k in range(n_jobs):
        # find t_max
        t=_t_max[k].split('*')
        tm=float(t[0])
        for t in t[1:]:
            tm*=float(t)
        _t_max[k]=tm
        # find a_orb
        _a_orb[k]=float(_a_orb[k].split('*')[0])*AU
        # find M_tot
        _M_tot[k]=(4.0/3.0)*pi*_rho[k]*(_R_eq[k]**3.0)

    #create dataframe
    df_deets=pd.DataFrame({'run':pd.Series(_num, dtype='int'),
                       'initial':pd.Series(_seed, dtype='str'),
                       'R_eq(m)':pd.Series(_R_eq, dtype='float'),
                       'X':pd.Series(_X, dtype='float'),
                       'f':pd.Series(_f, dtype='float'),
                       'N_tot':pd.Series(_N_tot, dtype='float'),
                       'rho(kgm-3)':pd.Series(_rho, dtype='float'),
                       't_run(s)':pd.Series(_t_max, dtype='float'),
                       'a_orb(m)':pd.Series(_a_orb, dtype='float'),
                       'M_tot(kg)':pd.Series(_M_tot, dtype='float'),
                       'R_b(m)':pd.Series(_R_b, dtype='float')})
    # df_deets=df_deets[['run','initial','R_eq(m)','X','f','N_tot','R_b(m)']] #set order of dataframe
    return df_deets
#-------------------------------------------------------------------------------

def read_problem_file_cores(prob_path):
    '''
    Returns the number of cores as stated in the problem.c file
    '''
    for line in open(prob_path):
        #test for certain strings in line, then use this to read value
        if 'omp_set_num_threads' in line:
            return line

#-------------------------------------------------------------------------------
def ang_mom_com(pos,vel,m):
    '''
    Calculates the angular momentum of a group of N particles, relative to centre of mass

    Parameters
    ----------
    pos
        x, y, z positions
    vel
        vx, vy, vz velocities
    m
        N particle masses
    '''

    pos_com=centre_of_mass(pos,m)
    vel_com=centre_of_mass(vel,m)
    mass=numpy.column_stack((m,m,m))
    L=mass*numpy.cross(pos-pos_com,vel-vel_com)
    return L
#-------------------------------------------------------------------------------
def ang_mom(pos,vel,m):
    '''
    Calculates the angular momentum of a group of N particles, relative to origin

    Parameters
    ----------
    pos
        x, y, z positions
    vel
        vx, vy, vz velocities
    m
        N particle masses
    '''
    mass=numpy.column_stack((m,m,m))
    L=mass*numpy.cross(pos,vel)
    return L

#-------------------------------------------------------------------------------
def kinetic_energy(vel,m):
    '''
    Calculates the kinetic energy of a group of N particles

    Parameters
    ----------
    vel
        vx, vy, vz velocities
    m
        N particle masses
    '''
    KE=0.5*m*(numpy.linalg.norm(vel,axis=1)**2.0)
    return KE
#-------------------------------------------------------------------------------
def grav_pot_energy_direct1(pos,m):
    '''
    Calculates the gravitational potential energy of each particle, in a group of N particles, relative to the origin
    Uses for loops, so is sloooooowwwww

    Parameters
    ----------
    pos
        x, y, z positions
    m
        N particle masses
    '''
    N=len(m)
    PE=[]
    for i in range(N):
        #print "particle {}\r".format(i),
        _PE=0
        for j in range(N):
            if j!=i:
                sep=numpy.linalg.norm(pos[j,:]-pos[i,:])
                _PE+=G*m[i]*m[j]/sep
        PE.append(_PE)
    return PE
#-------------------------------------------------------------------------------
def grav_pot_energy_direct2(pos,m):
    '''
    Calculates the gravitational potential energy of each particle, in a group of N particles, relative to the origin
    Uses a matrix of separations for speed, works for up to 1e4 particles

    Parameters
    ----------
    pos
        x, y, z positions
    m
        N particle masses
    '''
    N=len(m)
    print "define sep matrix"
    # if N<=1e4:
    #     print "big ol' matrix"
    #     rj_mat=numpy.tile(pos,(N,1,1))
    #     sep=numpy.linalg.norm(rj_mat-numpy.rot90(rj_mat,3,axes=(0,1)),axis=2) #separations of particle j to particle i (matrix: symmetric, diagonal = 0)
    # elif N>1e4:
    #     print "big N, split things up"
    #     print "x component"
    #     rjx_mat=numpy.tile(pos[:,0],(N,1))
    #     sepx=rjx_mat-rjx_mat.T
    #     print "y component"
    #     rjy_mat=numpy.tile(pos[:,1],(N,1))
    #     sepy=rjy_mat-rjy_mat.T
    #     print "z component"
    #     rjz_mat=numpy.tile(pos[:,2],(N,1))
    #     sepz=rjz_mat-rjz_mat.T
    #     print "magnitude"
    #     sep=numpy.sqrt((sepx*sepx)+(sepy*sepy)+(sepz*sepz))
    # else:
    #     print "problem with N"

    rj_mat=numpy.tile(pos,(N,1,1))
    sep=numpy.linalg.norm(rj_mat-numpy.rot90(rj_mat,3,axes=(0,1)),axis=2) #separations of particle j to particle i (matrix: symmetric, diagonal = 0)

    print "define mass matrix"
    mj=numpy.tile(m,(N,1))
    mimj=mj*mj.T #matrix of all values m_j*m_i (symmetric, with diagonal values = m_j*m_j)
    #create matrix of potential energy due to each particle, sum to get total PE on each particle
    print "calculate PE"
    PE=numpy.sum(numpy.divide(G*mimj,sep,where=sep!=0),axis=1) #note that we ignore/set equal to zero cases of division by zero, this drops diagonal terms
    return PE

#-------------------------------------------------------------------------------

def N_body_sphere(N,M,R,rho,X,fname):
    '''
    Creates a file containing time on the first line, then particle properties (x,y,z,vx,vy,vz,m,r) for N particles in a uniform sphercial distribution
    Parameters
    ----------
    N
        number of particles
    M
        total cloud mass
    R
        cloud radius
    rho
        particle material density
    X
        factor of circular rotation of cloud about the z axis
    fname
        name of file to save properties in
    '''
    t=0.0
    m = M/N
    r = (3.0*m/(4.0*pi*rho))**(1.0/3.0)
    OM_circ=X*numpy.sqrt(G*M/(R**3))
    f=open(fname,'w')
    f.write('{:.18e}\n'.format(t))
    for i in range(N):
        #evenly distribute the particles in the rotating frame of radius r2
        x1 = numpy.random.rand()
        x2 = numpy.random.rand()
        x3 = numpy.random.rand()
        z = R*(1.0-(2.0*x2))
        x = numpy.sqrt((R*R)-(z*z))*cos(2.0*pi*x3)
        y = numpy.sqrt((R*R)-(z*z))*sin(2*pi*x3)
        pos=numpy.array([x,y,z])
        Om_circ=numpy.array([0,0,OM_circ])		#solid body rotation
        vel=numpy.cross(Om_circ,pos)
        vx = vel[0]
        vy = vel[1]
        vz = vel[2]
        f.write("{:.18e}\t{:.18e}\t{:.18e}\t{:.18e}\t{:.18e}\t{:.18e}\t{:.18e}\t{:.18e}\n".format(x,y,z,vx,vy,vz,m,r))
    f.close()
    return
#-------------------------------------------------------------------------------

def calc_particle_energies(fname):
    '''
    Returns an array of total energy of each particle. Total energy is the sum of kinetic energy and potential energy of a particle (due to all other particles)

    Parameters
    ----------
    fname
        txt file containing particle data to load
    '''
    sim = rebound.Simulation() #create rebound sim
    sim.G=G #set grav constant
    add_file_to_sim(sim,fname) #load particles from file to simualtion
    sim.calculate_energy_out() #generate a txt file 'energy.txt' containing KE and PE of each particle
    E=numpy.genfromtxt('energy.txt') #load energy.txt
    KE=(E[:,0])
    PE=(E[:,1])
    E_tot = KE+PE #find total energy of each particle. E<0 is bound
    return E_tot

def bound_ang_mom(fname):
    '''
    Returns the angular momentum vector of the cloud, relative to the origin, which is the sum of angular momentum of only the bound particles in a system
    Also returns the number of bound particles, and the simulation timestamp of the data file
    Parameters
    ----------
    fname
        txt file containing particle data to load
    '''
    # find bound particles by calculating the sum of KE and PE for each particle
    t,df_dat=load_dat_file(fname)
    # N=len(df_dat)
    # print N
    df_dat['E(J)']=numpy.nan
    df_dat['E(J)']=calc_particle_energies(fname)
    #calculate angular momentum, only for bound particles
    df_dat=df_dat[df_dat['E(J)']<0]
    pos=numpy.array(df_dat[['x(m)','y(m)','z(m)']])
    vel=numpy.array(df_dat[['vx(ms^-1)','vy(ms^-1)','vz(ms^-1)']])
    m=numpy.array(df_dat['m(kg)'])
    N=len(df_dat)
    L=numpy.sum(ang_mom(pos,vel,m),axis=0)
    # print "number of bound particles = {}".format(N)
    # print "angular momentum of cloud = {}".format(L)
    # print sum(L)
    return L,N,t

def bound_ang_mom2(fname,a):

    # find bound particles by calculating the sum of KE and PE for each particle
    t,df_dat=load_dat_file(fname)
    df_dat['E(J)']=numpy.nan
    df_dat['E(J)']=calc_particle_energies(fname)
    df_dat=df_dat[df_dat['E(J)']<0]     #only bound particles
    # Transform to fixed reference frame
    pos=numpy.array(df_dat[['x(m)','y(m)','z(m)']])
    vel=numpy.array(df_dat[['vx(ms^-1)','vy(ms^-1)','vz(ms^-1)']])
    m=numpy.array(df_dat['m(kg)'])
    N=len(df_dat)
    pos,vel=rotating_to_heliocentric_array_precalc(pos,vel,a,t)
    # Transform to CoM
    pos=centre_of_mass(pos,m)
    vel=centre_of_mass(vel,m)
    # find angular momentum
    L=numpy.sum(ang_mom(pos,vel,m),axis=0)
    return L,N,t

def bound_ang_mom3(df_dat,a,t):

    # find bound particles by calculating the sum of KE and PE for each particle
    # t,df_dat=load_dat_file(fname)
    # df_dat['E(J)']=numpy.nan
    # df_dat['E(J)']=calc_particle_energies(fname)
    # df_dat=df_dat[df_dat['E(J)']<0]     #only bound particles
    # Transform to fixed reference frame
    pos=numpy.array(df_dat[['x(m)','y(m)','z(m)']])
    vel=numpy.array(df_dat[['vx(ms^-1)','vy(ms^-1)','vz(ms^-1)']])
    m=numpy.array(df_dat['m(kg)'])
    N=len(df_dat)
    pos,vel=rotating_to_heliocentric_array_precalc(pos,vel,a,t)
    # Transform to CoM
    pos=centre_of_mass(pos,m)
    vel=centre_of_mass(vel,m)
    # find angular momentum
    L=numpy.sum(ang_mom(pos,vel,m),axis=0)
    return L,N,t


def line_fit(x,m,c):
    y=m*x+c
    return y

def find_multiple_system_mass(df_orb):
    '''
    Returns the total mass of all components that make up a binary/multiple system

    Parameters
    ----------
    df_orb
        dataframe containing the orbital data (at a single timestep?)
    '''

    # if len(numpy.unique(df_orb['t(s)']))!=1:
    #     print 'error, only pass orbits for a single time?'
    #     M_tot=-1
    #     return M_tot

    for k in range(len(df_orb)):
        orb=df_orb.iloc[k]
        primary_int=int(orb['i'])
        secondary_int=int(orb['j'])
        m2=float(orb['m2(kg)'])
        m1=float(orb['m1(kg)'])

        # Always add the first orbit, and initialise the lists
        if k==0:
            primary_list=[primary_int]
            secondary_list=[secondary_int]
            M_tot=m1+m2
            continue

        # Check if subsequent orbits are around a previous primary
        if (primary_int in primary_list): #the main binary orbits
            primary_list.append(primary_int)
            secondary_list.append(secondary_int)
            M_tot+=(m2)
        elif (primary_int in secondary_list): # this primary was a secondary in a previous orbit
            primary_list.append(primary_int)
            secondary_list.append(secondary_int)
            M_tot+=(m2)
        else: # this binary is unconnected to the system
            continue

    return M_tot
#-------------------------------------------------------------------------------
def find_bin_num(x,bins):
    '''
    This function finds what bin number (described by array bins) a value, from the array x, resides in

    Parameters
    ----------
    x
        an array of values to be tested
    bins
        the bins to sort values of x into (must be an ordered array, whose bounds contain all x values)

    Returns
    -------
    x_bin_num
        an array of ints for the bin number for each value of x: first bin = 1, second bin = 2...
    '''
    x_bin_num=[]
    for j in range(len(x)):
        _x=x[j]
        bin_num=-1 # variable to detect problem values
        for k in range(1,len(bins)):
            if k==len(bins)-1: # special case for the point that defines the last bin
                if _x>bins[k-1] and _x<=bins[k]:
                    x_bin_num.append(k)
                    bin_num=k
                    break
            else:
                if _x>=bins[k-1] and _x<bins[k]:
                    x_bin_num.append(k)
                    bin_num=k
                    break
        if bin_num==-1:
            x_bin_num.append('NaN')
            print 'ERROR'
    return x_bin_num
