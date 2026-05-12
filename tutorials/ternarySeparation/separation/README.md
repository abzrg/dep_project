# Separation

After the fluid velocity and the electric field strength is solved in the
previous step, namely the initial state step, now we can run the Lagrangian
solver to obtain the path of particle over a steady fluid flow at different
times.

To run the case simply run the `./Allrun` case.

```console
#! Serial/single-core run
$ ./Allrun

#! Parallel run (4 cores by default)
$ ./Allrun -parallel

#! Parallel run (8 cores)
$ NUMBER_OF_SUBDOMAINS=8 ./Allrun -parallel

#! Only reconstruct the last time-step
$ ./Allrun -prallel -latestTime
```

## TODO: Visualization HOWTO

```console
$ foamToVTK
$ foamToVTK -overwrite
$ foamToVTK -parallel -overwrite
```

Paraview steps...

1. slice the domain
1. ...


