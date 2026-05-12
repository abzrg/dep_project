# Ternary Separation

The simulation comprises of two steps, which manifests itself in two separate
cases:

- [`initialState`](./initialState/): uses `simpleFoam` to solve a steady-state
  fluid flow problem and obtains the continuum phase velocity field. Also, uses
  `electricPotentialFoam` to solve for electric field strength.
- [`separation`](./separation/): uses the computed fields from the
  `initialState` case as initial condition to track the particles throughout
  the domain at different times.
