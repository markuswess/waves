# Time-domain

We derive examples of linear wave models.

## Acoustic waves

We consider the propagation of sound waves of small amplitude in a homogeneous
isotropic medium viewed as an inviscid fluid.

- **Geometry:** A given open domain $\Omega\subset\mathbb R^3$, and time interval $[0,T], T>0$

- **Constituents:** The (vectorial) velocity field (unit $m/s$), the pressure ($kg/(ms^2)$) and the density ($kg/m^3$) of the fluid given as functions $v:\Omega\times [0,T]\to\mathbb R^d$, $p,\rho:\Omega\times [0,T]\to\mathbb R$ which are assumed to be small perturbations of the static state $v_0=0,p_0,\rho_0=const.$

- **Balance relations:** 
We impose the linearized Euler equation of momentum
\begin{align*}
\rho_0\partial_t v+\nabla p = 0,
\end{align*}
and the equation of continuity
\begin{align*}
\partial_t\rho+\rho_0\nabla\cdot v=0.
\end{align*}
For our model we choose $p,v$ to be the primal quantities.

To close our system of equations we need another equation which is provided by the 
- **Material law,** a linearized equation of state, relating $p$ and $\rho$:
\begin{align*}
\frac{1}{c^2}\partial_t p = \partial_t \rho,
\end{align*}
where the **parameter** $c:=\sqrt{\frac{\kappa}{\rho_0}}$ is the speed of sound in the respective medium (unit $m/s$).

Combining the balance relations and material laws above yields the first order (mixed) form of the acoustic wave equation
````{card}
\begin{align*}
\rho_0\partial_tv+\nabla p &= 0,&\text{in }&\Omega,\\
\frac{1}{c^2}\partial_t p+\rho_0 \nabla\cdot v &=f,&\text{in }&\Omega,\\
v(0)&=v_0\\
p(0)&=p_0\\
p(t)&=0,&\text{on }&\partial\Omega,\\
\end{align*}
````

where we collect all **external forces** as the right-hand-side term $f$ and added homogeneous Dirichlet boundary conditions, and initial data to close the system.


Defining a velocity potential $\phi$ we obtain
\begin{align*}
\partial_t \phi &= -p\\
\rho_0 v &= \nabla \phi
\end{align*}
and thus the scalar wave equation for $\phi$
````{card}
\begin{align}
  \frac{1}{c^2}
\partial_t^2\phi&=\Delta\phi-f\\
&+\text{i.c., b.c.}
\end{align}
````
For $f=0$ the pressure $p$ fulfills the scalar wave equation as well i.e., 
````{card}
\begin{align*}
  \frac{1}{c^2}
\partial_t^2 p&=\Delta p\\
&+\text{i.c., b.c.}
\end{align*}
````

### Boundary conditions
The homogeneous Dirichlet boundary conditions used above corresponds to a sound-soft boundary, the pressure of the total wave vanishes at the boundary. Alternatively one could use sound-hard boundary conditions which are in the acoustic case modeled by the Neumann boundary condition
\begin{align*}
\partial_n p &= 0,&\text{on }\partial\Omega.
\end{align*}

Another important type of boundary conditions for acoustic problems are so-called *impedance boundary conditions*
given by 
\begin{align*}
\partial_n p &= \lambda \partial_t p,&\text{on }\partial\Omega.
\end{align*}
for a suitable function $\lambda$.

```{card} first order absorbing boundary condition
 Example
^^^^
An important example of impedance boundary conditions are first order absorbing boundary conditions which can be motivated as follows:
D'Alembert's solution of the one-dimensional wave-equation with initial data $p(0)=p_0, v(0)=v_0$ is given by
\begin{align*}
p(t,x) = \frac{1}{2}\left(p_0(x-ct)+p_0(x+ct)+\frac{1}{c}\right)
\end{align*} 

````
