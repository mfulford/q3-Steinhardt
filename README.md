# q3-Steinhardt
Fortran code to calculate the 3rd order Steinhardt parameter, q3. 

Calculate q3 for atom i by projecting vectors connecting atom i to its neighbours j onto the l=3 spherical harmonics.
Include real and imaginary component of spherical harmonics for completeness. 

q3 parameter quantifies order/disorder and allows us to assign molecules as liquid-like and ice-like. 


### Compile code:

```
gfortran -fbounds-check q_3_4_6.f90 -o q346.x90
```

### Usage:

```
.`q346.x90 [some .dcd trajectory file] [cell dimensions file] [cutoff]
```

e.g

```
.`q346.x90 Trajectory_ice.dcd cell_ice.dat 3.5
```

where cell_ice.dat contains cell dimensions in Angstroms

```shell
cat cell_ice.dat
44.33853912,76.00162125,43.67976761
```


