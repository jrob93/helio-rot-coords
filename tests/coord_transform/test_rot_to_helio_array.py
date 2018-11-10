'''This script tests the coordinate transforms between heliocentric and rotating reference frames
Calculate position and velocity for concentric circular orbits, in the rotating frame,
as a function of time according to the eqns of motion for non interacting particles.
These coordinates are then transformed into the heliocentric frame'''

import numpy
import matplotlib.pyplot as pyplot
import sys
sys.path.insert(0, '..') #use this line in other completed runs places
import py_func as pf
import matplotlib.gridspec as gridspec

#set up rotating frame parameters
a0=30*pf.AU
mu=pf.G*pf.M_sun
disp=1e6
n_samp=10
T0=((a0/pf.AU)**(3.0/2.0))*(365.0*24*60*60)
print T0
Omega=numpy.sqrt(mu/(a0**3.0))
t=numpy.linspace(0,T0,n_samp)
print t

# test with system of concentric circular orbits, position at t=0
r0=numpy.array([0.0,0.0,0.0])
r1=numpy.array([disp,0.0,0.0])
r2=numpy.array([-disp,0.0,0.0])

# define figure
fig = pyplot.figure() #open figure once
gs = gridspec.GridSpec(1,2)
ax1 = pyplot.subplot(gs[0,0])
ax2 = pyplot.subplot(gs[0,1])

ax2.set_xlabel('x /m')
ax2.set_ylabel('y /m')
# ax1.set_xlim(-2*disp,2*disp)
# ax1.set_ylim(-2*disp,2*disp)
# ax1.set_aspect("equal")
ax1.axvline(0)
ax1.axvline(disp)
ax1.axvline(-disp)
ax2.set_aspect("equal")
ax2.set_xlabel('x /m')
ax2.set_ylabel('y /m')
ax2.set_xlim(-2*a0,2*a0)
ax2.set_ylim(-2*a0,2*a0)
ax2.set_title("rotating frame")
ax2.set_title("heliocentric frame")

N=3
r_pos=[r0,r1,r2]
r_pos=numpy.array(r_pos).reshape((N,3))
print r_pos,len(r_pos),len(r_pos.shape)
# c_range=['b','r','g']
for i in range(len(t)):
    print t[i]
    # for j in range(len(r_pos)):
    r=numpy.zeros((N,3))
    v=numpy.zeros((N,3))
    r[:,0]=r_pos[:,0]
    r[:,1]=-1.5*r_pos[:,0]*Omega*t[i]
    v[:,1]=-1.5*Omega*r_pos[:,0]
    ax1.scatter(r[:,0],r[:,1])#,c=c_range[j])

    R,V=pf.rotating_to_heliocentric_array(r,v,a0,t[i])    # transfrom back from rotating to helio frame

    ax2.scatter(R[:,0],R[:,1])#,c=c_range[j])
    print "rotating",r,v
    print "helio",R,V
    # V_circ=numpy.sqrt(mu/(a0+r[0]))
    # V_mag=numpy.linalg.norm(V)
    # print "v_circ={},mag(V)={},V_circ/V_mag={}".format(V_circ,V_mag,V_circ/V_mag)

    r,v=pf.heliocentric_to_rotating_array(R,V,a0,t[i])    # transform from helio to rotating frame
    ax1.scatter(r[:,0],r[:,1],marker ='x')#,c=c_range[j])

    # print "rotating",r,v
    # orb=[a0+r[0],0,0,0,0,2.0*numpy.pi*t[i]/T0]        # define heliocentric orbit
    # R,V=pf.orbit_by_time(orb,mu,0.0,t[i])   # find heliocntric position and velocity as a function of time
    # ax2.scatter(R[0],R[1],marker='x',c=c_range[j])
    print "helio",R,V

pyplot.show()
