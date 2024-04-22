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

(mixed_methods)=
# Mixed formulations for wave equations

To obtain a semi- and fully discrete system from the first order systems from {numref}`modelling` we apply a similar reasoning as in {numref}`mol` to the first order acoustic wave system
````{card}
```{math}
:label: wave_eq_1o
\begin{aligned}
\partial_tv+\nabla p &= f,&\text{in }&[0,T]\times \Omega,\\
\partial_t p+\nabla\cdot v &=0,&\text{in }&[0,T]\times \Omega,\\
p(0,\cdot)&=p_0,&\text{in }&\Omega,\\
v(0,\cdot)&=v_0,&\text{in }&\Omega,\\
v\cdot n&=0,&\text{in }&[0,T]\times \partial\Omega,\\
\end{aligned}
```
````
```{prf:remark}
In the strong sense above we have to look for pairs of solutions $p,v$ where the gradient and divergence as well as the time derivative are well-defined. By using the following weak formulation we may choose larger spaces to look for $p,v$.
```

Multiplying {eq}`wave_eq_1o` by smooth test functions $q,w$, integrating over the domain $\Omega$ and applying integration by parts yields
```{math}
:label: wave_eq_h1l2
\begin{aligned}
\partial_t\int_{\Omega}vw+\int_\Omega\nabla p\cdot w &= \int_\Omega fw,\\
\partial_t \int_\Omega pq-\int_\Omega v\cdot \nabla q &=0,\\
&\text{+ i.c}.
\end{aligned}
```
In the formulation above we look for $p\in C^1([0,T];H^1(\Omega))$, $v\in C^2([0,T];L^2(\Omega)^2)$, and the equations have to hold for all $q\in H^1(\Omega)$, $v\in L^2(\Omega)^2$.
```{prf:remark}
Alternatively one may apply integration by parts in the first equation. Then the derivative is shifted to the unknown $v$ which has to be in $C^1([0,T];H(\mathrm{div}(\Omega))$. The unknown $p$ is sufficient to be in $C^1([0,1];L^2(\Omega))$. The test functions have to be chosen accordingly.
```

## Mixed formulations in NGSolve
[download as jupyter notebook](https://markuswess.github.io/waves/_sources/time-integration/mixed_methods_py.ipynb)

As usual we start by importing and generating a mesh

```{code-cell} ipython
# import ngsolve and webgui
from ngsolve import *
from ngsolve.webgui import Draw

mesh = Mesh(unit_square.GenerateMesh(maxh=0.3));
```
### Product spaces
In `NGSolve` spaces maybe tensorized by using the `*` symbol

```{code-cell} ipython
order = 2
fesh1 = H1(mesh, order = order)
print("h1 dofs: ", fesh1.ndof)
fesl2 = VectorL2(mesh, order = order-1)
print("l2 dofs: ", fesl2.ndof)

fes = fesh1*fesl2
print("total dofs: ", fes.ndof," = ", fesh1.ndof+fesl2.ndof)
```
To assemble bilinear forms we need test and trial functions. From product spaces test and trial functions are returned as tuples:

```{code-cell} ipython
(p,v), (q,w) = fes.TnT()
bfgrad = BilinearForm(grad(p)*w*dx).Assemble()
```
The matrix is still an operator on the large space!
```{code-cell} ipython
import matplotlib.pyplot as pl
pl.spy(bfgrad.mat.ToDense())
```
Gridfunctions on product spaces consist of components
```{code-cell} ipython
gfu = GridFunction(fes)
print("length of whole vector: ",len(gfu.vec))
gfup = gfu.components[0]
print("length of first component: ",len(gfup.vec))
gfuv = gfu.components[1]
print("length of second component: ",len(gfuv.vec))
```
Restriction and embedding operators are available:
```{code-cell} ipython
emb_p = fes.Embedding(0)
emb_v = fes.Embedding(1)

print(emb_p.shape)
print(emb_v.shape)

res_p = fes.Restriction(0)
res_v = fes.Restriction(1)

print(res_p.shape)
print(res_v.shape)
```
### Mixed operators

Alternatively one may define mixed operators. Note that in this case the test and trial functions must be obtained from the base spaces.
```{code-cell} ipython
fesh1 = H1(mesh,order=order)
feshdiv = HDiv(mesh,order=order)

p_,q_ = fesh1.TnT()
v_,w_ = feshdiv.TnT()

bfgrad = BilinearForm(grad(p_)*w_*dx).Assemble()
pl.spy(bfgrad.mat.ToDense())
```
Matrices may also be transposed:
```{code-cell} ipython
n = specialcf.normal(2)
bfdiv = BilinearForm(-div(v_)*q_*dx+v_.Trace()*n*q_*ds).Assemble()
print((bfdiv.mat.T.ToDense()-bfgrad.mat.ToDense()).Norm())
```
