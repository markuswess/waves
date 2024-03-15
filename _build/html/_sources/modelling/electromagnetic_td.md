# Electromagnetic waves


Similar to {numref}`acoustic_td` and {numref}`elastic_td` the 
**geometry** is a  given open domain $\Omega\subset\mathbb R^3$, and time interval $[0,T], T>0$.

The **constituents** are the electric field $E:[0,T]\times \Omega\to\mathbb R^3$ (unit $m\,kg/(s^3A)$) and the magnetic field intensity $H:[0,T]\times \Omega\to\mathbb R^3$ (unit $A/m$), as well as the electric flux density $D:[0,T]\times \Omega\to\mathbb R^3$ (also called electric displacement field, unit $As/m^2$) and the magnetic induction $B:[0,T]\times \Omega\to\mathbb R^3$ (also called magnetic flux density unit $kg/(s^2 A)$). Further  quantities, which act as **forcings** are the electric current density $J:[0,T]\times \Omega\to\mathbb R^3$ (unit $A/m^2$) and the electric charge density $\rho:[0,T]\times\Omega\to\mathbb R$ (unit $As/m^3$).

As **balance relations** we use Faraday's law
```{math}
-\partial_t B = \nabla\times E,
```
Ampère’s law with Maxwell’s correction

```{math}
\partial_t D = \nabla\times H -J
```
and Gauß' laws
```{math}
\begin{aligned}
\nabla\cdot D& = \rho,&\nabla\cdot B &= 0.
\end{aligned}
```
Note that conservation of charge
```{math}
\partial_t\rho+\nabla\cdot J = 0
```
follows from Ampere's law and Gauß' law for the electric field.
One possible viewpoint is to take $E,H$ as the primary quantities leading to the so called $EH$-formulation of Maxwell's equations. To this end we also need the linear
**material laws**
```{math}
\begin{aligned}
D&=\varepsilon E,&B&=\mu H,
\end{aligned}
```
where the **material parameters** $\varepsilon,\mu$ are the permittivity (unit $A^2s^4/(kg\,m^3)$) and permeability (unit $kg\,m/(s^2A^2)$) of the materials in question.
The above leads to 
````{card}
```{math}
\partial_t \mu H + \nabla\times E &=0,\\
\partial_t \varepsilon E-\nabla\times H&=-J,\\
E(0)&=E_0,\\
H(0)&= H_0,\\
E\times n = 0
```
````
Similar to acoustic and elastic waves for constant $\mu,\varepsilon$ we also may obtain the second order equation for $E$

````{card}
```{math}
:label: maxwell_2o
\partial_t^2 E+c^2\nabla\times\nabla\times E = 0,
```
````
where we assumed the absence of electric currents and charges and defined the speed of light by $c^2=\frac{1}{\varepsilon\mu}$.
