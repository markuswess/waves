(mol)=
# Space-time methods vs. the method of lines




The first decision that has to be made when deciding on a numerical method for approximately solving the PDEs from {numref}`modelling` is whether to discretize the whole space-time cylinder $[0,T]\times \Omega\subset \mathbb R^{d+1}$ or apply the space and time discretization one after another.

In the following we choose the latter approach which is known in the literature as method of lines and can be traced back to the 1960's where it was used with finite-difference schemes as spacial discretization.

We present this approach in the context of Galerkin approximations.

Suppose we want to solve the homogeneous acoustic wave problem initial value problem
```{math}
:label: wave_2o
\begin{aligned}
\partial_t^2 u-\Delta u&=0,& \text{in }&[0,T]\times\Omega,\\
u(0,\cdot) &= u_0,& \text{in }&\Omega,\\
\partial t(0,\cdot) &= 0,& \text{in }&\Omega,\\
\nabla u\cdot n &= 0,& \text{on }&[0,T]\times \partial\Omega.
\end{aligned}
```

Up to this point we kept the space where we actually look for the solution $u$ vague. Multiplying the first line of {eq}`wave_2o` by a test-function $u'\in H^1(\Omega)$ integrating over the domain $\Omega$ and using integration by parts we obtain the *spacial weak formulation*
```{math}
:label: weak_space
\int_\Omega \partial_t^2 u(t,x)u'(x)dx+\int_\Omega\nabla u(t,x)\cdot\nabla u'(x)dx-\int_{\partial\Omega}\nabla u(t,x)\cdot n(x) u'(x) dS(x)=0.
```
By the boundary conditions we assume the boundary term above to vanish.

Next we pick a finite dimensional subspace $V\subset H^1(\Omega)$ with a suitable basis $b_0,\ldots, b_{N}$, i.e.,
```{math}
V = \mathrm{span}\{b_0,\ldots,b_{N}\}
``` 
and approximate the time-dependent solution $u$ by 
```{math}
:label: disc_sol
\tilde u(t,x) = \sum_{j=0}^N \tilde u_j(t)b_j(x),
```
where $\tilde u_j\in C^2([0,T])$.
Inserting {eq}`disc_sol` and replacing $u'$ by each basis function $b_j$ respectively (the equality has to hold for every choice of test function) yields the system of $N$ ordinary differential equations (ODEs)

```{math}
:label: odes
\sum_{j=0}^N\frac{d^2}{dt^2} \tilde u_j(t)
\int_\Omega  b_j(x)b_i(x) dx+\sum_{j=0}^N\tilde u_j(t)\int_\Omega\nabla b_j(x) \cdot\nabla b_i(x)=0,\quad i=0,\ldots,N.
```
If $u_0$ is not by chance an element of $V$ the initial condition $u(0,\cdot)=u_0$ cannot be satisfied by our approximate solution $\tilde u$ thus we have to use an approximate initial condition $\tilde u(0,\cdot) = \tilde u_0\in V$ where $\tilde u_0=\sum_{j=0}^{N}\tilde u^0_j b_j\approx u_0$ is usually chosen as an orthogonal projection or interpolation.

Using the vectors $\mathbf u=(\tilde u_0,\ldots, \tilde u_{N})^\top$ and $\mathbf u_0=(\tilde u^0_0,\ldots, \tilde u^0_{N})^\top$ and the matrices $\mathbf M, \mathbf S$
defined by
```{math}
\begin{aligned}
(\mathbf M)_{i,j}&:=\int_\Omega b_j(x)b_i(x)dx,&
(\mathbf S)_{i,j}&:=\int_\Omega \nabla b_j(x)\cdot \nabla b_i(x)dx
\end{aligned}
```
we may rewrite the system {eq}`odes` including initial conditions in the compact form

```{math}
:label: odes_mat
\begin{aligned}
\frac{d^2}{dt^2}\mathbf M\mathbf u+\mathbf S\mathbf u &= 0,&\text{on } [0,T],\\
\mathbf u(0) &= \mathbf u_0.
\end{aligned}
```
The solution of the  semi-discrete system {eq}`odes_mat` may now be approximated using appropriate time-integration techniques for ordinary differential equations.

Thus, to obtain a full method it remains to pick the discrete basis $b_j$ and the time-stepping scheme.
