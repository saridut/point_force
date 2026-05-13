# point_force
Calculates the velocity field due to a point force in Stokes flow.

Compile the code as follows
```
gfortran parameters.F90 main.F90
```
and run it
```
a.out xmin xmax zmin zmax dx dz
```
where `[xmin, xmax]` and `[zmin, zmax]` gives the domain over which the calculate the velocity field. `dx` and `dz` are the grid spacings along `x` and `z`. The velocity field is calculated at `y = 0`.
To change the point force and its location directly change the variables `pf` and `r0` in the file *main.F90*. The output files will contain the `x` and `z` components of the velocity field at the grid points.
