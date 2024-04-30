(dg)=
# Discontinuous Galerkin methods for the wave equation
We want to solve the first order system 
````{card}
```{math}
:label: wave_eq_1o
\begin{aligned}
\partial_tv+\nabla p &= f,&\text{in }&[0,T_0]\times \Omega,\\
\partial_t p+\nabla\cdot v &=0,&\text{in }&[0,T_0]\times \Omega,\\
p(0,\cdot)&=p_0,&\text{in }&\Omega,\\
v(0,\cdot)&=v_0,&\text{in }&\Omega,\\
v\cdot n&=0,&\text{in }&[0,T_0]\times \partial\Omega,\\
\end{aligned}
```
````
numerically,  using {ref}`leap_frog`.
Since the computational effort depends mainly on the efficiency of the factorization of the mass matrices and the application of their inverses it would be benificial to use piecewise polynomial spaces, with basis functions supported only on single elements.

Let $\mathcal T_h$ be a decomposition of the domain $\Omega$ into triangles/tetrahedra such that
```{math}
\bigcup_{T\in\mathcal T_h}\bar T = \bar \Omega.
```
Then we want to use the discrete spaces
```{math}
:label: dg_spaces
\begin{aligned}
\mathcal W_h&:=\left\{p\in L^2(\Omega): p|_T\in\mathcal P^k(T)\, \forall T\in\mathcal T_h\right\},\\
\mathcal V_h&:=\left\{v\in L^2(\Omega)^2: v|_T\in\mathcal P^{k-1}(T)^2\, \forall T\in\mathcal T_h\right\},
\end{aligned}
```
where $\mathcal P^k(T)=\mathrm{span}\{p:T\to\mathbb R: (x,y)\mapsto x^ly^m, l+m\leq k\}$ for a fixed polynomial order $k\in\mathbb R$.

However, since $\mathcal W_h\subsetneq H^1(\Omega)$ and $\mathcal V_h\subsetneq H(\mathrm{div})(\Omega)$ the usual derivation of weak formulations fails. A naive attempt to just insert a basis of the discrete spaces {eq}`dg_spaces` into a weak formulation like {eq}`semi_disc_1o` would yield a familiy of decoupled local problems on the single elements $T$ with homogenuous boundary conditions but not an approxmation of the global problem.

This problem is fixed by modifying the weak forms $b,\tilde b$ of the gradient and divergence, which for $p\in H^1(\Omega),v\in H(\mathrm{div})(\Omega)$ and homogeneous Neumann boundary conditions should satisfy 

```{math}
:label: weak_grad
\begin{aligned}
b(p,v) &= \int_\Omega \nabla p\cdot v,\\
\tilde b(v,p) &= \int_\Omega \nabla\cdot v p=\int_\Omega v\cdot\nabla p = b(p,v).
\end{aligned}
```

To motivate the following we consider the one-dimensional case first.
## Numerical fluxes in one dimension

Consider a one-dimensional mesh $\mathcal T = \{T_j = (a_{j-1},a_{j}), j = 1,\ldots N\}$ for $a_0<a_1<\ldots<a_N$.
Since $p\in \mathcal W_h$  is still piecewise smooth we may apply the gradient on each element. To make up for the lost information across the element boundaries we add a linear combination of the inner and outer boundary values at $a_j,a_{j-1}$ with respect to each element, i.e.
```{math}
:label: dg_grad
b(p,v) := \sum_{j=1}^N\int_{T_j} \nabla p\cdot v- \left(\alpha p_+(a_{j-1})+\beta p_-(a_{j-1})\right) v(a_{j-1})+\left(\alpha p_-(a_{j})+\beta p_+(a_{j})\right) v(a_j)
```
with $p_\pm(x)=\lim_{\varepsilon\to 0} p(x\pm\varepsilon)$, where we test with a smooth function $v$.
We want a consistent bilinear form, i.e., if we plug in $p\in H^1(\Omega)$ we want to obtain the weak gradient from {eq}`weak_grad`. It immediately follows that for this hold we need to set $\beta = -\alpha$.
Using the notation $\{p\}:=\frac{1}{2}(p_++p_-)$ we may rewrite {eq}`dg_grad` with $\beta=-\alpha$ by

```{math}
:label: dg_grad_mean
b(p,v) = \sum_{j=1}^N\int_{T_j} \nabla p\cdot v- 2\alpha \left(\{p\}(a_{j-1})-p_+(a_{j-1})\right) v(a_{j-1})+ 2\alpha \left(\{p\}(a_{j})-p_-(a_{j})\right) v(a_{j})
```
Now note that the definition of $b$ in {eq}`dg_grad` is also valid for $v\in L^2(\Omega)^2$ (i.e., also for $v\in\mathcal V_h$) if we define it locally on each element. To obtain an energy preserving method we need $b(p,v)=\tilde b(v,p)$. Moreover we also want consistency of the $\tilde b$ bilinear form, i.e., if we insert $v\in H^1(\Omega)$ to obtain the weak divergence from {eq}`weak_grad`. Applying integration by parts in {eq}`dg_grad_mean`

```{math}
:label: dg_grad_mean
\begin{aligned}
b(p,v) &= \sum_{j=1}^N\int_{T_j} \nabla p\cdot v- 2\alpha \left(\{p\}(a_{j-1})-p_+(a_{j-1})\right) v(a_{j-1})+ 2\alpha \left(\{p\}(a_{j})-p_-(a_{j})\right) v(a_{j})\\
&= \sum_{j=1}^N-\int_{T_j}  p\nabla\cdot v- \left(2\alpha \{p\}(a_{j-1})+(1-2\alpha)p_+(a_{j-1})\right)v(a_{j-1})+ \left(2\alpha \{p\}(a_{j})-(2\alpha-1)p_-(a_{j})\right) v(a_{j})
\end{aligned}
```
which will be consistent for $\alpha=1/2$.

## Jumps and averages

For the discontinuous Galerkin (DG) spaces {eq}`dg_spaces` in higher dimensions we may define the average and jump operators as follows:

```{math}
\begin{aligned}
\{\cdot\}&:\begin{cases} W_h&\to L^2(\mathcal S)\\
p&\mapsto \left(x\mapsto \frac{1}{2}\lim_{\varepsilon\to 0}p(x+\varepsilon n)+p(x-\varepsilon n)\right)
\end{cases},\\
[\cdot]&:\begin{cases} V_h&\to L^2(\mathcal S)\\
v&\mapsto \left(x\mapsto \frac{1}{2}\lim_{\varepsilon\to 0}v(x+\varepsilon n)\cdot n-v(x-\varepsilon n)\cdot n\right)
\end{cases},
\end{aligned}
```
where $n$ is the outward normal and the skeleton $\mathcal S$ is defined by 
```{math}
\mathcal S:=\bigcup_{T\in\mathcal T}\partial T.
```
````{prf:remark}
The average and jump are independent of the choice of the normal vector $n$.

On the domain boundary the averaging operator can be adapted to fit boundary conditions
````

## DG in higher dimensions

Following the reasoning for the 1d case we obtain the general DG formulation


````{card}
```{math}
:label: wave_eq_1o
\begin{aligned}
\int_\Omega\partial_tv w+b(p,w) &= \int_\Omega f w,\\
\partial_t \int_\Omega pq-b(q,v) &=0,\\
p(0,\cdot)&=p_0,\\
v(0,\cdot)&=v_0,
\end{aligned}
```
for all $w\in\mathcal V_h,q\in\mathcal W_h$ with
```{math}
b(p,w) = \sum_{T\in\mathcal T} \int_T \nabla p\cdot w + \int_{\partial T}(\{p\}-p)w
```
````
