---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

(scattering_numerics)=
# FEM for scattering problems

[download as jupyter notebook](https://markuswess.github.io/waves/_sources/second_numerics/helmholtz.ipynb)

For a bounded domain $\Omega\subset\mathbb R^2$ we want to approximate the solution of the Helmholtz problem

````{card}
```{math}
:label: helmholtz_bounded
\begin{aligned}
-k^2 u-\Delta u &= 0,&\text{in }\Omega,\\
u &= u_D,&\text{on } \Gamma_D,\\
\nabla u\cdot n &= u_N,&\text{on } \Gamma_N,
\end{aligned}
```
````
where this time we prescribe Dirichlet boundary conditions on a part of the boundary $\Gamma_D\subset\partial\Omega$ and Neumann boundary conditions on another part $\Gamma_N\subset\partial\Omega,\Gamma_N\cap\Gamma_D=\emptyset$


We choose $\Omega = [0,1]\times[0,1]$ and $\Gamma_N=\{0\}\times[0,1]\cup [0,1]\times\{0\}$. and set $u_0=\exp(-20(y-1/3)^2)$ and neglect the Dirichlet boundary conditions for a start, i.e., $\Gamma_D = \emptyset$.

Similar to the time-domain problems we generate the geometry, the mesh and the discrete space $V$. 

```{code-cell} ipython

from ngsolve import *
from ngsolve.webgui import Draw
from netgen.occ import *

geo = OCCGeometry(unit_square.shape, dim = 2) 
mesh = Mesh(geo.GenerateMesh(maxh = 0.05))
Draw(mesh);
V = H1(mesh,order = 4, complex = True)
```

Note that, since we approximate a complex solution, we have to add the flag `complex = True`. Thus the vector of a `GridFunction` has complex entries (stored as tuples):

```{code-cell} ipython
:tags: [scroll-output]
gfu = GridFunction(V)

print(gfu.vec)
```
As in {numref}`mol` and {numref}`basic_fe_wave` we derive the (discrete weak formulation by multiplication by test functions and partial integration to obtain

```{math}
:label: weak_semidisc_hh
\begin{aligned}
-k^2\int_\Omega u_hu_h'+\int_\Omega\nabla u_h\cdot\nabla u'_h&=\int_{\Gamma_N} u_{h,N} u_h',
\end{aligned}
```
for all $u_h'\in V$
where $u_{h,N}\in V$ is an approximation of $u_N$.
We assemble the matrices choosing $u_N=\exp(-20(x^2+y^2))\exp(20\pi(x+y))$ and solve the problem
```{code-cell} ipython
u,v = V.TnT()
k = 15

a = BilinearForm(V)
a += (grad(u)*grad(v)-k**2*u*v)*dx
a.Assemble()

u_N = exp(-20*(x**2+y**2))
f = LinearForm(u_N*v*ds).Assemble()

gfu = GridFunction(V)
gfu.vec.data = a.mat.Inverse() * f.vec

# draw the solution
Draw (gfu, mesh, deformation=True);
# draw the gradient
Draw (grad(gfu).real, mesh, vectors = True);
```
Note that by default the real part of the scalar function is drawn, for the vectorial function we have to add `.real` to draw the vectors.

## Essential boundary conditions
Implementing dirichlet Boundary conditions works a little differently.

They are called *essential boundary conditions*: We have to set the solution field to the given Dirichlet values, and restrict the test-functions to 0 on the Dirichlet boundary:

$$
\text{find } u \in H^1, u = u_D \text{ on } \Gamma_D \text{ s.t. } A(u,v) = f(v) \quad \forall \, v \in H^1, v = 0 \text{ on } \Gamma_D
$$
where $A(u,v)$ is the bilinear form from the weak formulation.

We split the solution vector into two parts: The given coefficients on the Dirichlet boundary, and all other including internal coefficients and coefficients on the natural boundaries:

$$
u = \left( \begin{array}{c} u_D \\ u_f \end{array} \right)
$$

(f like free). Accordingly, the matrix and the right hand side are split as

$$
A = \left( \begin{array}{cc} A_{DD} & A_{Df} \\ A_{fD} & A_{ff} \end{array} \right)
\qquad
f = \left( \begin{array}{c} f_D \\ f_f \end{array} \right)
$$

The test functions are reduced to the free nodes. Thus, the equations to solve are

$$
A_{fD} u_D + A_{ff} u_f = f_f
$$

The given $u_D$ is moved to the right hand side, and thus

$$
u_f = {A_{ff}}^{-1} (f_f - A_{fD} u_D).
$$

It is convenient to assemble the complete matrix, invert a part of the matrix, end extend by zero:

$$
A^{-1,ff} := 
\left( \begin{array}{cc} 0 & 0 \\ 0 & A_{ff}^{-1} \end{array} \right)
$$


The case of homogeneous Dirichlet boundary condtions $u_D = 0$ is simple: We just say

$$
u = A^{-1,ff} f
$$
In NGSolve, the finite element space maintains a BitArray marking the free degrees of freedom. It is used to invert the sub-matrix:


```{code-cell} ipython
V = H1(mesh, dirichlet="top|right", order = 5)
u,v = V.TnT()

a = BilinearForm((grad(u)*grad(v)-k**2*u*v)*dx).Assemble()
f = LinearForm(u_N*v*ds).Assemble()

gfu = GridFunction(V)
gfu.vec.data = a.mat.Inverse(V.FreeDofs()) * f.vec
Draw (gfu, mesh, deformation=True);
```
Non-homogeneous Dirichlet b.c. are reduced to homogeneous once as follows: Choose some function $\widetilde u$ such that 

$$
\widetilde u = u_D \quad \text{on } \Gamma_D
$$

and set $u = \widetilde u + w$. So we have to find a correction $w$ with $w = 0$ on $\Gamma_D$ and

$$
A(w,v) = f(v) - A(\widetilde u , v) \qquad \forall \, v \text{ with }  v = 0 \text{ on } \Gamma_D
$$

In matrix notation, the correction $w$ is

$$
w = A^{-1,ff} (f - A \tilde u)
$$

Now we set the Dirichlet boundary values and don't worry about the rest. The Set - function of a GridFunction does some kind of interpolation.

```{code-cell} ipython
uD = 0.1*exp(-20*((x-1)**2+(y-3/4)**2))
gfu.Set (uD, definedon=mesh.Boundaries("top|right"))
Draw(gfu);
```

Now we compute the correction
```{code-cell} ipython
r = f.vec - a.mat * gfu.vec
gfu.vec.data += a.mat.Inverse(V.FreeDofs()) * r

Draw(gfu, mesh, deformation = True);
Draw (grad(gfu).real, mesh, vectors = True);
```
