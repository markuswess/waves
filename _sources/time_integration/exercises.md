(exercises_th)=
# Exercises for time-stepping schemes

## Exercise 1
Implement the Newmark time-stepping for general $\beta,\gamma$. Choose the second order acoustic wave problem on the unit-square with initial data 
```{math}
u_0 = \exp(-30(x-1/2)^2),\quad v_0=0
```
and homogeneous Neumann boundary conditions. Argue that the exact solution for $t=1$ is identical to the initial values. Compute the error at $t=1$ in the $L^2$ norm for $\beta = \gamma/2$ and varying $\gamma$. What order of convergence with respect to $\tau$ do you observe? (Pay particular interest to the case $\gamma=1/2$.)

## Exercise 2
Use the example of the previous exercise to check conservation of energies $\mathbf E$ and $E_j^{\tau,\beta,\gamma}$ for the parameters $\gamma_0=1/2,\beta_0=1/4$, $\gamma_1=1/2,\beta_1 =0$ and $\gamma_2 = 1/3, \beta_2 = 1/6$.

## Exercise 3
Use the example of the previous exercises to verify the CFL condition: Compute the largest magnitude eigenvalue $\lambda_{\mathrm {max}}$ of the resulting problem $\mathbf K \mathbf u = \lambda \mathbf M \mathbf u$. Then plot the energies $\mathbf E, E^{\tau,\beta,\gamma}$ over time for the explicit method $\beta=0,\gamma=1/2$ and timesteps $\tau_{\pm}^2=\frac{4}{\lambda_{\mathrm{max}}}(1\pm\varepsilon)$ and small $\varepsilon>0$.


## Exercise 4
Implement the Crank-Nicholson time-stepping for the $H^1-L^2$ space pairing {eq}`wave_eq_h1l2` for a similar problem as in the previous exercises. Fix the orders of spaces and time-step and check the convergence of the error with respect to the mesh-size $h$. How do the orders of the two spaces have to be related to obtain an efficient method?

## Exercise 5
Implement the Leap-Frog time-stepping for the first order equation using an $L^2-H(\mathrm{div}))$ pairing. Preserving the skew-symmetric structure of the system, what options do you find for imposing homogeneous boundary conditions?

## Exercise 6
For Galerkin approximations using a discrete space $V_h$, the error of the spacial discretization is usually governed by the error of the spacial approximation $I-\Pi_h$ where the projection $\Pi_h:L^2\to V_h$ is defined by
```{math}
\int_\Omega \Pi_h f v = \int_\Omega f v
```
for all $v\in V_h$.
Compare the projection errors for $f = \exp\left(-30((x-1/2)^2+(y-1/2)^2)\right)$ on the `unit_square` for the spaces `H1` and `L2` for different mesh-sizes $h$ with respect to the degrees of freedom (`H1.ndof`, `L2.ndof`).

 Assume you were able to use the space `L2` for a discretization of the wave equation. Compare the computational effort for factorizing the mass matrices for the `H1` and `L2` spaces and $200$ applications of the inverse mass matrix with respect to accuracy (measured by the projection error) on different meshes. Use also `L2.Mass().Inverse()` instead of assembling the matrix. 

Consider also a three-dimensional example.
