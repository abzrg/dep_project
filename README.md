# Di-Electro-Phoresis (DEP) Project

The code in this project is, mainly, an extension of the `lagrangian` library
that includes two
[Dielectrophoretic](https://en.wikipedia.org/wiki/Dielectrophoresis) particle
forces.

- [`DEPForce`](./src/lagrangian/intermediate/submodels/Kinematic/ParticleForces/DEPForce/)
- [`DEPRepulsiveForce`](./src/lagrangian/intermediate/submodels/Kinematic/ParticleForces/DEPRepulsiveForce/)

The first calculates the DEPForce particle force in terms of particle's radius,
the strength of the electric field $E \cdot \nabla E$, relative permittivity of
the medium $\epsilon_m$, and the so-called Clausius-Mossotti Factor,
$\mathrm{CM}$. The second one calculates a specialized Body Force that prevents
particles to hit the channel wall&mdash;as they do not in reality&mdash;when
they experience a Dielectrophoretic forces in the channel.

Along with the library comes two solvers that are used to simulate particle
separation, in microfluidics devices. They are:

- [`electricPotentialFoam`](./applications/solvers/basic/electricPotentialFoam/)
- [`icoUncoupledKinematicParcelDEPFoam`](./applications/solvers/lagrangian/icoUncoupledKinematicParcelDEPFoam/)

The former is just a simple Laplacian solver that solves for the electric
potential. It also computes the electric field and its strength. The latter
however is a Lagrangian particle tracking solver that is linked against the
Lagrangian library in this project.


## Usage

First load the OpenFOAM environment, and then build the project:

```console
$ ./Allwmake -l -j8
```

Finally, go to the tutorial directory and have a look at the two cases and run
them. Read the steps in the corresponding tutorials' `README.md` file.


## Reference

The corresponding paper can be found here:

> Reza Derakhshan, Ali Bozorgzadeh, and Abas Ramiar. **Numerical Investigation of
> Ternary Particle Separation in a Microchannel with a Wall-Mounted Obstacle
> Using Dielectrophoresis**. *Journal of Chromatography A* (2023) doi:
> [10.1016/j.chroma.2023.464079](https://doi.org/10.1016/j.chroma.2023.464079)
