(time_domain_pmls)=
# Perfectly Matched Layers in the time domain

Again we start by motivating the construction in one dimension

## Time domain PMLs in one dimension

We have seen in {numref}`fq_domain_pmls` that the one-dimensional, untruncated complex scaled problem {eq}`cs_problem` is solved (outside of the support of inhomogeneities or sources) by
```{math}
\tilde u(x)=
\tilde c\exp(-\alpha k(x-a))\exp(ikx).
```
Now since we impose the condition that the solution is square integrable this only holds for $k>0$, since for negative $k$ this function is exponentially increasing.

In the frequency domain this does not pose a problem, since we may just use positive wavenumbers $k$. To go back to time-domain however, one needs to use the inverse Fourier transform, which involves an integral with respect to $k$ over **the whole real axis**. Thus, simply replacing the factors $-ik$ by time-derivatives and solving the resulting time-domain problem will result in exploding solutions since the PML actually imposes the *wrong* radiation condition for $k<0$.

To avoid this problem usually a frequency (or wavenumber) dependent damping parameter $\alpha(k):=\alpha_0/k$ for some $\alpha_0>0$ is chosen.

Since it leads to a more convenient formulation we state the equivalent first order system to find $p\in H^1_0(0,\infty), v\in L^2(0,\infty)$


```{math}
:label: weak_scaled
\begin{aligned}
-ik\int_0^a vw -ik(1+i\alpha/k)\int_a^\infty vw +\int_0^\infty p'w &= \int_0^\infty fq,\\
-ik\int_0^a pq -ik(1+i\alpha/k)\int_a^\infty pq -\int_0^\infty vq' &= 0,
\end{aligned}
```
for all $w\in H^1_0(0,\infty), q\in L^2(0,\infty)$. 
This can be conveniently rewritten by

```{math}
:label: weak_scaled
\begin{aligned}
-ik\int_0^\infty vw  &=-\alpha\int_a^\infty vw -\int_0^\infty p'w +\int_0^\infty fq,\\
-ik\int_0^\infty pq  &=-\alpha\int_a^\infty pq +\int_0^\infty vq'.
\end{aligned}
```
As the scaling by $k$ now flips the sign of the damping this equation yields correct outgoing solutions for positive and negative values of $k$.
Thus we might transform the equation above back to time-domain to obtain the damped wave equation in first order form (where we use the same symbols for the time domain functions as for the frequency domain ones) 

```{math}
:label: weak_scaled
\begin{aligned}
\partial_t \int_0^\infty vw  &=-\alpha\int_a^\infty vw -\int_0^\infty \partial_x pw +\int_0^\infty fq,\\
\partial_t \int_0^\infty pq  &=-\alpha\int_a^\infty pq +\int_0^\infty v\partial_x q.
\end{aligned}
```
Now again the domain may be truncated and a finite element discretization may be applied.


## Time domain PMLs in higher dimensions

In higher dimensions the same principles as in 1d apply: We have to use a similar frequency/wavenumber dependent scaling parameter to make the complex scaling work for **all** frequencies. Subsequently the resulting equation is transformed back to time domain. However in higher dimensions there arise a few more complications to consider. 

Again, we start with the general first-order complex-scaled Helmholtz equation
```{math}
\begin{aligned}
-ik\int_{\Omega}v\cdot w+\int_\Omega\nabla p\cdot w &= \int_\Omega fw,\\
-ik \int_\Omega pq-\int_\Omega v\cdot \nabla q &=0.
\end{aligned}
```

As in {numref}`fq_domain_pmls` we decompose $\Omega = \Omega_{\mathrm{int}}\cup\Omega_{\mathrm{ext}}$ and apply a complex scaling $\mathbf x\mapsto\tilde{\mathbf x}$ in $\Omega_{\mathrm{ext}}$ with Jacobian matrix $\mathbf J$. This leads to 
```{math}
\begin{aligned}
-ik\int_{\Omega}v\cdot w\det\mathbf J+\int_\Omega\mathbf J^{-\top}\nabla p\cdot w\det\mathbf J &= \int_{\Omega_{\mathrm{int}}} fw,\\
-ik \int_\Omega pq\det\mathbf J-\int_\Omega v\cdot \mathbf J^{-\top}\nabla q\det\mathbf J &=0.
\end{aligned}
```
Due to the fact that $\mathbf J$ depends on $k$ this equation might have rational dependencies on $-ik$ and thus lead to convolution operators in time, when transforming back to time domain. Usually these rational dependencies are resolved by introducing additional unknowns.


To make the formulation feasible for explicit methods we need to seperate the differentials in time (factors of $k$) and space. Replacing the vectorial functions $v,w$ by their Piola transforms $v =\frac{1}{\det\mathbf J}\mathbf J \tilde v$ leads to scaling independent mixed bilinear forms:
```{math}
\begin{aligned}
-ik\int_{\Omega}\tilde v\cdot \mathbf J^{\top}\mathbf J \tilde w\frac{1}{\det\mathbf J}+\int_\Omega\nabla p\cdot \tilde w &= \int_{\Omega_{\mathrm{int}}} fw,\\
-ik \int_\Omega pq\det\mathbf J-\int_\Omega \tilde v\cdot \nabla q &=0.
\end{aligned}
```

### Cartesian time domain PMLs

Using the cartesian scaling as in the frequency domain (but with scaling parameter $i\alpha/k$ we have
```{math}
\begin{aligned}
\mathbf J(\mathbf x) &=
\begin{cases}
\begin{pmatrix}1-\frac{\alpha}{ik}&0\\0&1\end{pmatrix},&\mathbf x\in\Omega_{lr}\\
\begin{pmatrix}1&0\\0&1-\frac{\alpha}{ik}\end{pmatrix},&\mathbf x\in\Omega_{tb}\\
\begin{pmatrix}1-\frac{\alpha}{ik}&0\\0&1-\frac{\alpha}{ik}\end{pmatrix},&\mathbf x\in\Omega_{c}\\
I,&\mathbf x\in\Omega_{\mathrm{int}}
\end{cases}\\
\det\mathbf J(\mathbf x) &=
\begin{cases}
(1-\frac{\alpha}{ik})^2,&\mathbf x\in\Omega_{c}\\
(1-\frac{\alpha}{ik}),&\mathbf x\in\Omega_{lr}\cup\Omega_{tb}\\
1,&\mathbf x\in\Omega_{\mathrm{int}}
\end{cases}
\end{aligned}
```
with 
```{math}
\begin{aligned}
\Omega_{\mathrm{int}}&=\{(x,y):x\in (a_0,a_1),y\in (b_0,b_1)\},\\
\Omega_{lr}&=\{(x,y):x\in (a_0,a_1)^c,y\in (b_0,b_1)\},\\
\Omega_{tb}&=\{(x,y):x\in (a_0,a_1),y\in (b_0,b_1)^c\},\\
\Omega_{c}&=\mathbb R^2\setminus\left(\Omega_{\mathrm{int}}\cup\Omega_{tb}\cup\Omega_{lr}\right).
\end{aligned}
```
To obtain a suitable time-domain equation we have to look at the expressions $\frac{-ik}{\det\mathbf J}\mathbf J^\top\mathbf J$ and $-ik\det\mathbf J$. 
