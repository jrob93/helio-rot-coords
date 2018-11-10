'''This script tests the coordinate transforms between heliocentric and rotating reference frames
We calculate position and velocity for particles on concentric heliocentric circular orbits,
from their orbits, as a function of time.
we then transfrom in the rotating frame, and transform back to check'''

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

# test with system of concentric circular orbits
a1=a0+disp
a2=a0+-disp

# define figure
fig = pyplot.figure() #open figure once
gs = gridspec.GridSpec(1,2)
ax1 = pyplot.subplot(gs[0,0])
ax2 = pyplot.subplot(gs[0,1])

ax2.set_xlabel('x /m')
ax2.set_ylabel('y /m')
# ax2.set_xlim(-2*disp,2*disp)
# ax2.set_ylim(-2*disp,2*disp)
# ax2.set_aspect("equal")
ax2.axvline(0)
ax2.axvline(disp)
ax2.axvline(-disp)
ax1.set_aspect("equal")
ax1.set_xlabel('x /m')
ax1.set_ylabel('y /m')
ax1.set_xlim(-2*a0,2*a0)
ax1.set_ylim(-2*a0,2*a0)
ax2.set_title("rotating frame")
ax1.set_title("heliocentric frame")

a_range=[a0,a1,a2]
c_range=['b','r','g']
for i in range(len(t)):
    print t[i]
    for j in range(len(a_range)):
        orb=[a_range[j],0,0,0,0,2.0*numpy.pi*t[i]/T0]        # define heliocentric orbit
        R,V=pf.orbit_by_time(orb,mu,0.0,t[i])   # find heliocntric position and velocity as a function of time
        r,v=pf.heliocentric_to_rotating(R,V,a0,t[i])    # transform from helio to rotating frame
        print "helio",R,V
        print "rotating",r,v
        ax1.scatter(R[0],R[1],c=c_range[j])
        ax2.scatter(r[0],r[1],c=c_range[j])
        R,V=pf.rotating_to_heliocentric(r,v,a0,t[i])    # transfrom back from rotating to helio frame
        print "helio",R,V
        ax1.scatter(R[0],R[1],marker='x',c=c_range[j])
        v_shear=-1.5*Omega*r[0]     # compare rotating velocity to expected shear velocity
        v_mag=numpy.linalg.norm(v)
        print "v_shear={},mag(v)={},v_shear/v_mag={}\n".format(v_shear,v_mag,v_shear/v_mag)

pyplot.show()
