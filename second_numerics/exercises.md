(exercises_th)=
# Exercises for time-harmonic waves

## Exercise 1

Solve the scattering problem with a first order absorbing boundary condition on a domain $[0,1]^2\setminus \bar\Omega_0$ for a compact domain $\Omega_0\subset[0,1]^2$ with a plane-wave incident field. Visualize the incident, scattered and total field.

(exercise_2_2)=
## Exercise 2

Solve the Helmholtz problem on the geometry from [square_inf.ipynb](https://markuswess.github.io/waves/_sources/second_numerics/square_inf.ipynb) for right hand side $\hat f(x)=\exp(-30(x^2+y^2))$ with homogeneous Dirichlet boundary conditions on the boundary of the square-shaped scatterer and first-order absorbing boundary conditions on the outer boundary. Plot the $L^2$-norm of the solution depending on $k^2/\pi^2\in[0,11]$. How do you interpret the results? What happens if the right hand side is a Gaussian peak which is slightly off-centered?

## Exercise 3

Solve the time-domain problem corresponding to the Helmholtz equation from {ref}`exercise_2_2` for time harmonic right hand sides $f(t,x) = \cos(\omega t)\exp(-30((x-0.2)^2+y^2))$ for $\omega=\sqrt 2\pi,\sqrt 3\pi,\sqrt {5}\pi$. Plot the $L^2$-norm of the solution with respect to time.


## Exercise 4

Argue that the discrete resonance problem corresponding to {ref}`exercise_2_2` is of the form: find $\mathbf x\in\mathbb C^N,\omega\in\mathbb C$
```{math}
:label: quad_evp
M_0 \mathbf x+\omega M_1\mathbf x+\omega^2 M_2\mathbf x = 0,
```
for certain matrices $\mathbf M_j$.
Rewrite {eq}`quad_evp` as a (twice as large) generalized linear eigenvalue problem in $\omega$ by introducing $\mathbf y = \omega \mathbf x$. Solve the resulting linear generalized eigenvalue problem and plot the resulting eigenvalues and eigenfunctions.
