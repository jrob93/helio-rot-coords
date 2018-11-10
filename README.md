# helio-rot-coords
Here we provide functions to transform position and velocity coordinates for particles from the heliocentric reference frame to the rotating reference frame.

NOTES
--------

* Define a orbital radius for the rotating frame
* Define the position of the frame using a time variable. At t=0 the frame origin is at X=a, Y=0 in heliocentric coordinates
* The assumed units are length in metres and time in seconds, therefore velocities are ms^{-1} 
* gravitational constant G=6.67428e-11 m3 kg-1 s-2 and we assume a solar mass central body at the heliocentric origin, M=1.98855e30 # kg

Installation
-----------------------

Install like so ::
```
pip install transforms
```
Here we have a simple example::

```python
import numpy
import coord_transforms

pos_rot=numpy.array([1,0,0])
vel_rot=numpy.array([0,0,0])
a=30*1.496e+11   
t=0

print(coord_transforms.rotating_to_heliocentric_array(pos,vel,a,t))
```

Documentation
-------------
I'm working on it!

Please do not hesitate to get in touch with any questions
