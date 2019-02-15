# Steinhardt Parameters 
Fortran code to calculate the 3rd, 4th & 6th order Steinhardt parameters. 

Calculate q3, q4, & q6 for atom i by projecting vectors connecting atom i to its neighbours j onto the l=3, l=4 and l=6 spherical harmonics.
Include real and imaginary component of spherical harmonics for completeness. 

Steinhardt parameter quantifies the local atomic structure enabling the molecular phase to be determined. 


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
./q346.x90 Trajectory_ice.dcd cell_ice.dat 3.5
```

where cell_ice.dat contains xyz cell dimensions in Angstroms

```shell
$ cat cell_ice.dat
44.33853912,76.00162125,43.67976761
```

### Output:

Steinhardt parameters for each molecule at each timestep of trajectory. 
