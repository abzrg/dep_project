# The Lagrangian Solver: `icoUncoupledKinematicParcelDEPFoam`

The source of this solver is the same as the `icoUncoupledKinematicParcelFoam`
solver. It only differs in the libraries it links against. While
`icoUncoupledKinematicParcelFoam` links against OpenFOAM's
`lagrangian/intermediate` solver, `icoUncoupledKinematicParcelDEPFoam` does so
against the `lagrangian/intermediate` library present in this project.
