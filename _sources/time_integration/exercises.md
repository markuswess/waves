(exercises_th)=
# Exercises for time-stepping schemes

## Exercise 1
Implement the Newmark time-stepping for general $\alpha,\beta$. Choose the second order acoustic wave problem on the unit-square with initial data 
```{math}
u_0 = \exp(-30(x-1/2)^2),\quad v_0=0
```
and homogeneous Neumann boundary conditions. Argue that the exact solution for $t=1$ is identical to the initial values. Compute the error at $t=1$ in the $L^2$ norm for $\beta = \gamma/2$ and varying $\gamma$. What order of convergence with respect to $\tau$ do you observe? (Pay particular interest to the case $\gamma=1/2$.)

## Exercise 2
Use the example of the previous exercise to check conservation of energies $\mathbf E$ and $E_j^{\tau,\alpha,\beta}$ for the parameters $\gamma_0=1/2,\beta_0=1/4$, $\gamma_1=1/2,\beta_1 =0$ and $\gamma_2 = 1/3, \beta_2 = 1/6$
