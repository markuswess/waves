(elastic_td)=
# Elastic waves

Similar to {numref}`acoustic_td` the 
**geometry** is a  given open domain $\Omega\subset\mathbb R^3$, and time interval $[0,T], T>0$

The **constituents** are the displacement $u:[0,T]\times\Omega\to\mathbb R^3$ (unit $m$), the particle velocity $v=\partial_t u$ (unit $m/s$), the strain tensor 
```{math}
\varepsilon(u) = \frac{1}{2}\left(\nabla u+ (\nabla u)^\top\right)
```
 the stress tensor $\sigma:[0,T]\times\Omega\to \mathbb R^{3\times 3}_{sym},$ (unit $kg/(ms^2)$) where $\mathbb R^{3\times 3}_{sym}$ denotes the space of symmetric, real three-by-three matrices.

We either choose $u$ or $v,\sigma$ as our primal unknowns.

**Material parameters** here are the mass density $\rho:\Omega\to (0,\infty)$ and the material stiffness given by Hooke's tensor (also known as elastic tensor) $C:\Omega\to \mathcal L(\mathbb R^{3\times 3}_{sym},\mathbb R^{3\times 3}_{sym})$ (unit $kg/(ms^2)). Note that here, different from linear acoustics we already assumed that the density $\rho$ is uniform in time.

As **balancing law** we impose linear conservation of momentum
```{math}
:label: momentum_el
\rho \partial_t v = \nabla\cdot \sigma
```
where the divergence $\nabla\cdot\sigma$ has to be understood row wise, i.e. $(\nabla\cdot\sigma)_i=\sum_j \partial_{x_j}\sigma_{i,j}$.
The **material law** to link stress and strain is Hooke's law given by
```{math}
:label: hooke
\sigma = C(\varepsilon(u))
```

Provided boundary and initial data this leads to the following first and second order formulation.
````{card}
```{math}
\begin{aligned}
\rho\partial_t^2 u-\nabla\cdot C(\varepsilon(u))&=f,&\text{on }& (0,T)\times \Omega,
\\
u(0,\cdot)&=u_0,&\text{on }& (0,T)\times\Omega,
\\
\partial_t u(0,\cdot)&=v_0,&\text{on }& \Omega,
\\
u &= u_V,&\text{on }&\partial \Omega,
\end{aligned}
```
````
````{card}
```{math}
\begin{aligned}
\rho\partial_tv-\nabla\cdot \sigma &= f,&\text{on }& (0,T)\times\Omega,\\
\partial_t \sigma-C(\varepsilon(v))&=0,&\text{on }& (0,T)\times\Omega,\\
v(0)&=v_0,&\text{on }& \Omega,\\
\sigma(0)&=C\varepsilon(u_0),&\text{on }& \Omega,\\
v(t) &= \partial_t u_V(t),&\text{on }& \partial\Omega,
\end{aligned}
```
````
where we assumed dynamic boundary conditions. Alternatively one may assume static boundary conditions, namely conditions on the normal component of $\sigma$.

Hooke's tensor $C$, also called the 4-index elastic tensor has 81 components. In *isotropic elasitcity* it is given by only two parameters, namely the so-called *Lamé parameters* $\mu,\lambda$
```{math}
C(\varepsilon) = 2\mu\varepsilon+\lambda\mathrm{tr}(\varepsilon)I=2\mu(\varepsilon-\frac{1}{3}\mathrm{tr}(\varepsilon I))+\kappa\mathrm{tr}(\varepsilon)I.
```
Here the two parameters $\mu,\lambda$ correspond to the decomposition of waves into shear waves, depending on the shear modulus $\mu$ and compressional waves (pressure waves) depending on the compression modulus $\kappa=\frac{2}{3}+\lambda$.
Then the linear second order elastic wave equation in isotropic and homogeneous media takes the form
````{card}
```{math}
\rho\partial_t^2 u+\mu \nabla\times\nabla\times u-3\kappa\nabla(\nabla\cdot u)=f
```
````
Assuming a vanishing shear modulus $\mu\to 0$ leads to compressional waves only, and thus to the linear acoustic wave equation for the hydrostatic pressure $p=\frac{1}{3}\mathrm{tr}\sigma$ already introduced in  {numref}`acoustic_td`. Note however, that historically the sign conventions for pressure and stress are flipped in fluid and solid mechanics.
