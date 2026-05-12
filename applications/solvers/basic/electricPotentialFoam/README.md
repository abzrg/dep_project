# Electric Potential Solver

It solves the Laplace equation to obtain the electric potential distribution
throughout the domain.

```math
\begin{equation}
\nabla^2 V = 0
\end{equation}
```

After solving for $V$, we, specific to this project, compute electric field,
$E$, and electric field strength $E \cdot \nabla E$. These two fields are
required by the Lagrangian solver to calculate the path of particles based on
various forces that act on them, in particular the `DEPForce` that is
implemented in the `lagrangian/intermediate` library of this project.
