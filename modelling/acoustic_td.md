(acoustic_td)=
# Acoustic waves


We consider the propagation of sound waves in a homogeneous
isotropic medium viewed as an inviscid fluid.

- **Geometry:** A given open domain $\Omega\subset\mathbb R^3$, and time interval $[0,T], T>0$

- **Constituents:** The (vectorial) velocity field (unit $m/s$), the pressure ($kg/(ms^2)$) and the density ($kg/m^3$) of the fluid given as functions $v:\Omega\times [0,T]\to\mathbb R^3$, $p,\rho:\Omega\times [0,T]\to\mathbb R$. 

- **Balance relations:** 
We impose the Euler equation of momentum
```{math}
:label: momentum
\rho\left(\partial_t v+v\cdot\nabla v\right)+\nabla p = 0,
```
and the equation of continuity
```{math}
:label: continuity
\partial_t\rho+\nabla\cdot (\rho v)=0,
```
where $\partial_t$ denotes the time-derivative and $\nabla\cdot$ the divergence.
Note that these balance relations ad-hoc hold only in integral form (integration over each volume and time) and can be stated pointwise only for sufficiently regular functions.


For our model we choose $p,v$ to be the primal quantities.

To close our system of equations we need another equation which is provided by the 
- **Material law,** an  equation of state, relating $p$ and $\rho$:
```{math}
:label: constitutive
p = h(\rho)
```
for some scalar function $h$, which is usually assumed to be sufficiently regular and strictly increasing.

Next we assume that the quantities $p,\rho,v$ are small perturbations of an equilibrium state, i.e., 
```{math}
\begin{aligned}
\rho&=\rho_0+\rho_1, 
&p&=p_0+p_1,
&v&=v_0+v_1,
\end{aligned}
```
thus to derive a linear model we may neglect the quadratic quantities of $\rho_1,p_1,v_1$.
We impose that the medium is in a static state, i.e. $p_0,\rho_0$ are independent of $t$ and $v_0=0$ and obtain the linearized equation of momentum

```{math}
:label: lin_momentum
\rho_0\partial_t v_1+\nabla p_0+\nabla p_1 = 0.
```
Since the static state $p_1,\rho_1,v_1=0$ also satisfies the (linearized) balance of momentum we  immediately obtain that $p_0(x)=const.$ From the constitutive relation {eq}`constitutive` we also obtain that $\rho_0(x)=const.$
Thus the linearized conservation of mass and momentum become
```{math}
:label: lin_conservation
\begin{aligned}
\partial_t\rho_1+\rho_0\nabla\cdot v_1&=0,\\
\rho_0\partial_t v_1+\nabla p_1 &= 0.
\end{aligned}
```

It remains to use material law {eq}`constitutive`. Taking the derivative with respect to time and linearizing we obtain
```{math}
\partial_t \rho_1= h'(\rho_0)\partial_t\rho_1=c^2\partial_t\rho_1,
```
where the **parameter** $c^2:=h'(\rho_0)$ is the speed of sound in the respective medium (unit $m/s$). As another **parameter** we introduce the bulk modulus $\kappa_0$ of the material as $\kappa_0=\rho_0 c^2$, (unit $kg/(ms^2)$).

Combining the balance relations and material laws above, and ommiting the index $\cdot_1$ yields the first order (mixed) form of the acoustic wave equation
````{card}
```{math}
:label: wave_eq_1o
\begin{aligned}
\partial_tv+\frac{1}{\rho_0}\nabla p &= f,&\text{in }&[0,T]\times \Omega,\\
\partial_t p+\kappa_0\nabla\cdot v &=0,&\text{in }&[0,T]\times \Omega,\\
v(0,\cdot)&=v_0,&\text{in }&\Omega,\\
p(0,\cdot)&=p_0,&\text{in }&\Omega,\\
p&=0,&\text{in }&[0,T]\times \partial\Omega,\\
\end{aligned}
```
````
where we also add an **external force** as the right-hand-side term (acceleration, unit $m/s^2$) $f$, in the vectorial equation and added homogeneous Dirichlet boundary conditions, and initial data to close the system.

````{prf:Remark} Variable materials
For non-constant material-parameters $\rho_0,\kappa_0$ the constitutive relation {eq}`constitutive` needs to be replaced by
```{math}
p = f(\rho,s)
``` 
where $s$ is the entropy (unit $kg\,m^2/(s^2\,K)$) which satisfies the conservation law (adiabatic hypothesis)
```{math}
\partial_t s +v\cdot\nabla s = 0.
```
A similar reasoning as above then leads to {eq}`wave_eq_1o` with variable in space $\rho_0,\kappa_0$ (but still uniform initial pressure $p_0$).
````

The well-known second order form of the acoustic wave equation is then obtained by taking the time-derivative of the scalar equation and inserting the vectorial equation to obtain
````{card}
```{math}
:label: wave_eq_1o
\begin{aligned}
\partial_t^2 p-\kappa_0\nabla\cdot\frac{1}{\rho_0}\nabla p &= -\kappa_0 \nabla\cdot f,&\text{in }&[0,t]\times \Omega,\\
p(0,\cdot)&=p_0,&\text{in }&\Omega,\\
\partial_t p(0,\cdot)&=-\kappa_0\nabla\cdot v_0,&\text{in }&\Omega,\\
p&=0,&\text{on }&[0,t]\times \partial\Omega,\\
\end{aligned}
```
````
Defining a velocity potential $\phi$ we obtain
\begin{align*}
\partial_t \phi &= -\frac{1}{\rho_0}p\\
v &= \nabla \phi
\end{align*}
and thus the scalar wave equation for $\phi$
````{card}
\begin{align}
  \frac{1}{c^2}
\partial_t^2\phi&=\Delta\phi-f\\
&+\text{i.c., b.c.}
\end{align}
````
In the mathematical treatment of the acoustic wave equation the physical meaning of the modeled quantities is often neglected and the letter $u$ is used for the scalar unknown.
## Boundary conditions
The homogeneous Dirichlet boundary conditions used above corresponds to a sound-soft boundary, the pressure of the total wave vanishes at the boundary. Alternatively one could use sound-hard boundary conditions which are in the acoustic case modeled by the Neumann boundary condition
\begin{align*}
\partial_n p &= 0,&\text{on }\partial\Omega.
\end{align*}
Also a mix of boundary conditions is possible.

Another important type of boundary conditions for acoustic problems are so-called *impedance boundary conditions*
given by 
\begin{align*}
\partial_n p &= \lambda \partial_t p,&\text{on }\partial\Omega.
\end{align*}
for a suitable function $\lambda$.

````{prf:Example} First order absorbing boundary condition
:label: first_order
An important example of impedance boundary conditions are first order absorbing boundary conditions which can be motivated as follows:
D'Alembert's solution of the one-dimensional (homogeneous) wave-equation with initial data $p(0)=p_0, v(0)=v_0$ is given by
```{math}
:label: dalembert
\begin{aligned}
p(t,x) = \frac{1}{2}\left(p_0(x-ct)+p_0(x+ct)+\frac{1}{c}\int_{x-ct}^{x+ct} v_0(\xi)d\xi\right).
\end{aligned}
```
The first two terms constitute waves travelling to the right and left respectively.
Thus, if we assume compactly supported initial data on a finite interval $(a,b)$ and want to impose boundary conditions at $a,b$ that do not change the solution we may use
```{math}
:label: 1st_order_1d
\begin{aligned}
  c\partial_x p(t,a)&=\,\partial_t p(t,a),\\
  c\partial_x p(t,b)&=-\,\partial_t p(t,b),
\end{aligned}
```
for all times $t>0$.

In higher dimensions it is not so straightforward any more to differentiate incoming from outgoing waves. Still the higher dimensional counterpart of {eq}`1st_order_1d` can be stated as
```{math}
:label: 1st_order_xd
c\nabla p(t,x)\cdot n =\partial_t p(t,x)
```
for all $x\in\partial\Omega$ and times $t>0$, where $n$ is the outward normal vector of the boundary $\partial \Omega$.
The condition {eq}`1st_order_xd` can be interpreted as a prescribed impedance at the interface $\partial\Omega$ or, alternativly as a rough approximation to the exact absorbing boundary conditions.

````

