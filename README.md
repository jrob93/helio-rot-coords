# helio-rot-coords
Here we provide functions to transform position and velocity coordinates for particles from the heliocentric reference frame to the rotating reference frame.

![alt text](https://github.com/jrob93/helio-rot-coords/blob/master/2d_Hill_fig_2.png)

NOTES
--------

* Define a orbital radius for the rotating frame
* Define the position of the frame using a time variable. At t=0 the frame origin is at X=a, Y=0 in heliocentric coordinates
* The assumed units are length in metres and time in seconds, therefore velocities are ms^{-1} 
* gravitational constant G=6.67428e-11 m^{3} kg^{-1} s^{-2} and we assume a solar mass central body at the heliocentric origin, M=1.98855e30 kg
* The functions accept a numpy array of (x,y,z) and (vx,vy,vz) for position and velocity respectively.
* Multiple particles can be transformed if you pass arrays of shape (N,3)
* The transformed position and velocity arrays will also be shape (N,3)

Installation
-----------------------
Try install with pip:
```
pip install coord_transforms
```
\* pip install incoming, need to add setup.py scripts to appropriate places\*

Or you should be able to install manually like so:
```
git clone git@github.com:jrob93/helio-rot-coords.git
cd helio-rot-coords
python setup.py build
python setup.py install
```
(This means that you will need to update the package manually by running the setup.py commands again)

Note that you will then need to uninstall manually, by removing the install files. The install files can be found be reinstalling and recording where everything goes:
```
python setup.py install --record files.txt
```

Example
-----------------------
Here we have a simple example:

```python
import numpy
import coord_transforms

pos_rot=numpy.array([1,0,0])
vel_rot=numpy.array([0,0,0])
a=30*1.496e+11   
t=0

print(coord_transforms.rotating_to_heliocentric_array(pos_rot,vel_rot,a,t))
```

Documentation
-------------
I'm working on it! I will put here the full list of fucntions available and how they work.

Please do not hesitate to get in touch with any questions
