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

(resonance_numerics)=
# FEM for resonance problems

[download as jupyter notebook](https://markuswess.github.io/waves/_sources/second_numerics/resonance.ipynb)

Our goal is to solve the Neumann Helmholtz resonance problem on a bounded domain $\Omega$ find $u\in H^1,\omega^2\in \mathbb R$ such that

```{math}
\begin{aligned}
\omega^2 u &=-\Delta u,&\text{in }\Omega.
\nabla u\cdot n &= 0&\text{on }\partial\Omega
\end{aligned}
```

For a discrete (finite element) space $V$ the weak discrete formulation is to find $u_h\in V$, $\omega^2>0$
```{math}
\begin{aligned}
\omega^2\int_{\Omega} u_hu'_h &=\int_{\Omega}\nabla u_h\cdot\nabla u'_h,&\text{in }\Omega.
\end{aligned}
```
for all $u_h'\in V$.
By expanding the solution $u_h$ into a suitable basis we obtain the linear (in $\omega^2$) generalized matrix eigenvalue problem

```{math}
\omega^2\mathbf M \mathbf u = \mathbf S\mathbf u
```
In `NGSOlve` we assemble the matrices as usual:

```{code-cell} ipython
from ngsolve import *
from ngsolve.webgui import Draw
from netgen.occ import *

geo = OCCGeometry(unit_square.shape, dim = 2) 
mesh = Mesh(geo.GenerateMesh(maxh = 0.1))

V = H1(mesh,order = 3)

u,v = V.TnT()
S = BilinearForm(grad(u)*grad(v)*dx).Assemble()
M = BilinearForm(u*v*dx).Assemble()
```

We may solve the matrix eigenvalue problem using standard libraries from `numpy` or `scipy`:

```{code-cell} ipython
import scipy as sp
lam,vecs = sp.linalg.eigh(S.mat.ToDense(),M.mat.ToDense())
print(lam[:10]/sp.pi**2)
```

Lastly we draw the resulting eigenfunctions

```{code-cell} ipython
gfu = GridFunction(V)
for i in range(1,5):
  gfu.vec[:] = vecs[:,i]
  Draw(gfu);
```
