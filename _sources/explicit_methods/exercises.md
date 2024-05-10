(exercises_expl)=
# Exercises for explicit methods

## Exercise 1
Estimate the CFL condition of the Leap-Frog timestepping for the DG system. To this end estimate the largest eigenvalue of $\mathbf M^{-1}\mathbf K$ by $\mu_N$ for large enough $N$ given by a power iteration
```{math}
\mathbf x_{n+1}=\mathbf M^{-1}\mathbf K\mathbf x_n,\quad \mu_n:=\frac{\mathbf x_n^\top\mathbf M^{-1}\mathbf K\mathbf x_n}{\|\mathbf x_n\|}
```
with a random starting vector $\mathbf x_n$. 


## Exercise 2
Implement the DG-Method for the acoustic wave equation on the `unit_square` using suitable initial and boundary conditions and the Leap-Frog time-stepping. Choose the time-step as large as possible to satisfy the CFL-condition. 


## Exercise 3
Implement the first order mass-lumping method by supplying the according `IntegrationRule` to the differential symbol `dx`. Explore the sparsity pattern of the resulting mass matrix. Experiment with the argument `diagonal = True` for the mass `BilinearForm`.
Replace the `H1` space by `H1LumpingFESpace` (only for orders $1,2$). The `IntegrationRule`s can by obtained by `H1LumpingFESpace.GetIntegrationRules()`.

## Exercise 4
Use a suitable example to compare the efficiency of the explicit methods (Mass-Lumping, DG) to the efficiency of the conforming method (second order system, implicit time-stepping). Plot the respective errors against the computation times.