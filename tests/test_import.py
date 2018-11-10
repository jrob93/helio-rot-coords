import numpy
import coord_transforms

pos=numpy.array([1,0,0])
vel=numpy.array([0,0,0])
a=1
t=0

print(coord_transforms.rotating_to_heliocentric_array(pos,vel,a,t))
print(coord_transforms.heliocentric_to_rotating_array(pos,vel,a,t))
print(coord_transforms.rotating_to_heliocentric_array_precalc(pos,vel,a,t))
print(coord_transforms.heliocentric_to_rotating_array_precalc(pos,vel,a,t))
