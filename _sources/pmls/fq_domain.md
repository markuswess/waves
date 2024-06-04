(fq_domain_pmls)= 
# Perfectly Matched Layers in the Frequency Domain

We start by introducing PMLs in the simple one-dimensional setting.

## PMLs in one dimension

Consider the one dimensional Helmholtz problem for given $k>0$ and compactly in $(0,a)$, $a>0$ supported $f$ to find $u: (0,\infty)\to\mathbb R$ such that 
```{math}
:label: helmholtz_1d
\begin{aligned}
-u''-k^2 u &= f,\\
u(0) &= 0,\\
\lim_{x\to\infty} u'(x)-iku(x) &= 0.
\end{aligned}
```
The last condition is the Sommerfeld radiation condition which accounts for the missing boundary condition at the right.

The general solution to {eq}`helmholtz_1d` (without boundary conditions) is given by
```{math}
:label: gen_sol
u(x) = \frac{1}{2ik}\left(c_+-\int_0^x\exp(-ik\xi)f(\xi) d\xi\right) \exp(ikx)+\left(c_-+\frac{1}{2ik}\int_0^x\exp(ik\xi)f(\xi)d\xi\right) \exp(-ik x).
```


Inserting the boundary conditions we obtain that the solution of {eq}`helmholtz_1d` is given in closed form by

```{math}
u(x) = \frac{1}{2ik}\left(\int_0^a\exp(ik\xi)f(\xi) d\xi-\int_0^x\exp(-ik\xi)f(\xi) d\xi\right) \exp(ikx)-\frac{1}{2ik}\int_x^a\exp(ik\xi)f(\xi)d\xi \exp(-ik x).
```
Note that for $x\geq a$ this boils down to
```{math}
u(x)=\tilde c \exp(ikx),\quad \tilde c = \frac{1}{k}\int_0^a\sin(k\xi)f(\xi) d\xi.
```
Now the crucial observation here is that $u$ for $x>a$ is an analytic function which allows a holomorphic extension to $\mathbb C$, which we also denote by $u$. Moreover if we bend the path in the complex plane from the real axis upwards the value of $u$ decreases exponentially with increasing imaginary part of the argument, e.g., if we choose for any $\alpha>0$
```{math}
:label: scaling_1d
\tilde x(x):=\begin{cases}x,&x<a\\x+i\alpha(x-a),&x\geq a\end{cases}
```
we have that
```{math}
:label: u_scaled
\tilde u(x):=u(\tilde x(x))=
\tilde c\exp(-\alpha k(x-a))\exp(ikx),\quad x\geq a.
```
Note that the *unwanted* term $\exp(-ikx)$ in the general solution {eq}`gen_sol` on the other hand increases exponentially along $\tilde x$.
We have found a way to alter the solution of the original problem in a way that it remains unchanged in the interval of interest $(0,a)$ and the Sommerfeld radiation condition reduces to the condition that the altered solution is square integrable on the real axis. This condition fits nicely into the desired variational framework.

To obtain a numerical method one now alters the Helmholtz equation such that {eq}`u_scaled` is the (unique) solution. Applying the chain rule one obtains easily:
````{prf:Theorem}
For $k>0$ the function $\tilde u$ {eq}`u_scaled` is the (unique) solution to the problem to find $\tilde u\in C^2([0,a)\cup (a,\infty))\cap L^2((0,\infty))$ such that
```{math}
:label: cs_problem
\begin{aligned}
-\tilde u''-k^2 \tilde u &= f,&\text{in }(0,a),\\
-\frac{1}{1+i\alpha}\tilde u''-(1+i\alpha)k^2 \tilde u &= 0,&\text{in }(a,\infty),\\
\tilde u(0) &= 0,\\
\lim _{\varepsilon\to 0}\tilde u(x+\varepsilon)-\tilde u(x-\varepsilon) &= 0,\\
\lim _{\varepsilon\to 0}\tilde u'(x+\varepsilon)-(1+i\alpha)\tilde u'(x-\varepsilon) &= 0.
\end{aligned}
```
````
As usual we also may derive the weak formulation to find $\tilde u\in H^1_0((0,\infty))$ s.t.
```{math}
:label: weak_scaled
\int_0^a \tilde u' v'+\frac{1}{1+i\alpha}\int_a^\infty \tilde u'v'-k^2\int_0^a \tilde uv - k^2(1+i\alpha)\int_a^\infty \tilde uv=\int_0^a fv
```
for all $v\in H^1_0((0,\infty))$.
Note that the weak formulation (also known as **complex scaled** equation) {eq}`weak_scaled` is still posed on an unbounded domain and thus not feasible for a FEM discretization. Thus we truncate the domain to $(0,T)$ for some $T>a$ and impose homogeneous Dirichlet or Neumann boundary conditions at $x=T$: Namely we look for $\tilde u\in H^1_0((0,T))$ such that

```{math}
:label: weak_pml
\int_0^a \tilde u' v'+\frac{1}{1+i\alpha}\int_a^T \tilde u'v'-k^2\int_0^a \tilde uv - k^2(1+i\alpha)\int_a^T \tilde uv=\int_0^a fv
```
for all $v\in H^1_0((0,T))$.
 Since the solution decays exponentially for $x\to\infty$ we expect an exponentially small (with respect to the truncation $T$) error.


 ## PMLs for waveguides


### Decomposition into modes and radiation condition 
 The ideas from 1d can be easily carried over to waveguide geometries: Consider the Helmholtz equation on the semi infinite waveguide $\Omega=(0,\infty)\times (0,\pi)$ with homogeneous Dirichlet boundary conditions.
 A decomposition of the solution $u(x,y)$ into transversal modes yields
 ```{math}
:label: waveguide_modes
 u(x,y)=\sum_{n=1}^\infty \sin(ny) u_n(x),\quad f(x,y)=\sum_{n=1}^\infty \sin(ny)f_n(x),
 ```
 and the functions $u_n$ have to solve
 ```{math}
 -u_n''-(k^2-n^2)u_n = f_n.
 ```
 Again we assume that $f_n$ vanishes for $x>a$ and $k\notin\mathbb N$. Then the family of general solutions $u_n$ for $x>a$ is given by
 ```{math}
 u_n(x) = \begin{cases}
 c_+\exp(i\sqrt{k^2-n^2}x)+c_-\exp(-i\sqrt{k^2-n^2}x),&n<k\\
 c_+\exp(\sqrt{n^2-k^2}x)+c_-\exp(-\sqrt{n^2-k^2}x),&n>k.
 \end{cases}
 ```
 The radiation condition is now formulated based on the decomposition above: A solution is called outgoing if it is of the form {eq}`waveguide_modes` with 

```{math}
:label: rad_cond_waveguide
 u_n(x) = \begin{cases}
 c_n\exp(i\sqrt{k^2-n^2}x),&n<k\\
 c_n\exp(-\sqrt{n^2-k^2}x),&n>k.
 \end{cases}
```
Modes with $n<k$ are called *propagating* waves, while the ones with $n>k$ are called *evanescent*.


### Complex scaling

Propagating solutions as defined above allow again a holomorphic continuation in $x$-direction. Applying the change of coordinates {eq}`scaling_1d` to the outgoing modes {eq}`rad_cond_waveguide` leads to 
```{math}
 \tilde u_n(x) = \begin{cases}
 c_n\exp(i\sqrt{k^2-n^2}x)\exp(-\alpha\sqrt{k^2-n^2}x),&n<k\\
 c_n\exp(-\sqrt{n^2-k^2}x)\exp(-i\alpha\sqrt{n^2-k^2}x),&n>k.
 \end{cases}
```
which also decay exponentially for $\alpha>0$ and $x\to\infty$.
Thus again we may truncate the domain to $\tilde\Omega = \Omega_{\mathrm {int}}\cup\Omega_{\mathrm{PML}}=(0,a)\times (0,\pi)\cup (a,T)\times(0,\pi)$ to obtain the PML formulation to find $\tilde u\in H^1_0(\tilde\Omega)$ such that

```{math}
:label: weak_pml_waveguide
\int_{\Omega_{\mathrm{int}}} \tilde u' v'+\frac{1}{1+i\alpha}\int_{\Omega_{\mathrm{PML}}} \tilde u'v'-k^2\int_{\Omega_{\mathrm{int}}} \tilde uv - k^2(1+i\alpha)\int_{\Omega_{\mathrm{PML}}} \tilde uv=\int_{\Omega_{\mathrm{int}}} fv
```
for all $v\in H^1_0((0,T))$.

````{prf:remark}
Similar arguments can be applied for two-sided or three-dimensional waveguides.
````



## Free space

In free space $\Omega = \mathbb R^{2,3}$ the situation is a bit more evolved, since there is no unique direction of propagation as in one dimension or waveguides. However we still may decompose the solution of the Helmholtz homogeneous equation into modes. This time we use polar coordinates to obtain for  $r$ large enough such that no sources or obstacles are present

```{math}
u(r\hat x)=
\sum_{n=-\infty}^\infty\left( c_{1,n} H_{|n|}^{(1)}(kr)+c_{2,n} H_{|n|}^{(2)}(kr)\right)\Phi_n(\hat x),
```
where $H_m^{(1,2)}$ are the (cylindrical) [Hankel functions](https://en.wikipedia.org/wiki/Bessel_function#Hankel_functions) of the first and second kind and $\Phi_n$ are the cylindrical harmonics.

The Hankel functions in the decomposition above now play the role of the exponentials in one dimension and the waveguide case: The Hankel functions of the first kind are the sought outgoing waves and the ones of the second kind are the unwanted incoming waves.
Thus, the radiation condition can again be stated as
```{math}
u(r\hat x)=
\sum_{n=-\infty}^\infty c_{n} H_{|n|}^{(1)}(kr)\Phi_n(\hat x),
```
for $r>r_0$ such that all sources and obstacles are contained in the ball $B_{r_0}(0)$.

The Hankel functions also have a suitable asymptotic behavior for the idea of complex scaling to work, namely

```{math}
H_n^{(1,2)}(z)\approx \sqrt{\frac{2}{\pi u}}\exp\left(\pm i(z-(2n+1)\pi/4)\right),\quad -\pi<\arg z <\pi
```

Thus, as long as the complex scaling $x\mapsto \tilde {\mathbf x}$ has the property that $\Im\left(\tilde{\mathbf x}(\mathbf x)\right)\to \infty$ for $|\mathbf x|\to\infty$ we obtain that the (wanted) Hankel function of the first kind is exponentially damped, while the (unwanted) Hankel function of the second kind increases.

Subsequently the complex scaled formulation is derived via the transformation rule: 

````{prf:Theorem}
Let $J(\mathbf x)$ denote the Jacobian matrix of the (piecewise smooth) scaling $\tilde {\mathbf x}(\mathbf x)\in\mathbb C^d$ such that $\Im\left(\tilde{\mathbf x}(\mathbf x)\right)\to \infty$ for $|\mathbf x|\to\infty$ and let $\Omega_{\mathrm{int}}$ be the interior domain where $\tilde{\mathbf x}=\mathbf x$. Then the complex scaled, radiating solution $\tilde u\in H^1(\mathbb R^d)$ of the Helmholtz equation fulfills 

```{math}
\int_{\Omega_{\mathrm{int}}}\nabla \tilde u\cdot \nabla v+
\int_{\Omega_{\mathrm{ext}}}\left(J^{-\top}\nabla \tilde u\right)\cdot\left(J^{-\top}\nabla v\right)\det J
-k^2\int_{\Omega_{\mathrm{int}}} \tilde u v 
-k^2\int_{\Omega_{\mathrm{ext}}} \tilde u v \det J=\int_{\Omega_{\mathrm{int}}}f v
```
for all $v\in H^1(\mathbb R^d)$ for right hand sides $f$ supported in $\Omega_{\mathrm{int}}$.
````
Subsequently the complex scaled equation is again truncated to a finite domain and homogeneous (Dirichlet) boundary conditions are imposed.


There are several variants to choose the change of variables $\tilde {\mathbf x}$. Two of them are presented in the following.


### Cartesian scalings

The simple idea of cartesian scalings is to superimpose scalings in $x$ and $y$ directions to obtain the $\tilde x(x)$. For given values $a_0,a_1,b_0,b_1$ such that all the sources and inhomogeneities are contained in the box $\Omega_{\mathrm{int}}=(a_0,a_1)\times(b_0,b_1)$ we define the scaling

```{math}
\tilde{\mathbf x}(\mathbf x)=
\begin{pmatrix}
\tilde x\\\tilde y
\end{pmatrix}
((x,y)^\top):=\begin{pmatrix}
x\\y\end{pmatrix}
+i\alpha\begin{pmatrix}
 (x-a_1) \chi_{x>a_1}(x)+(x-a_0) \chi_{x<a_0}(x)\\
 (y-b_1) \chi_{y>b_1}(y)+(y-b_0) \chi_{y<a_0}(y)
\end{pmatrix}
```
for some $\alpha>0$ and where
```{math}
\chi_{x>a}(x)=\begin{cases}1,&x>a\\0,&x\leq a\end{cases}
```
and similar.

To implement the PML formulation one needs to derive the Jacobian matrix of the mapping $\tilde{\mathbf x}$ which is given by

```{math}
J((x,y)^\top) =\mathbf I+i\alpha\begin{pmatrix}\chi_{x<a_0}(x)+\chi_{x>a_1}(x)&0\\0&\chi_{y<b_0}(y)+\chi_{y>b_1}(y)\end{pmatrix},
```
where $\mathbf I$ denotes the identity matrix.

Cartesian scalings are very popular due to their easy implementation.

### Radial scalings

Another idea is to scale the variable in radial direction, i.e., to use 

```{math}
\tilde{\mathbf x}(\mathbf x)=\left(r_0+i\alpha (\|\mathbf x\|-r_0)\right) \frac{\mathbf x}{\|\mathbf x\|}=r_0\left(1-i\alpha\right)\frac{\mathbf x}{\|\mathbf x\|}+i\alpha\mathbf x
```
for $\|x\|>r_0$.

Again we may directly compute the Jacobian matrix by

```{math}
J(\mathbf x)=i\alpha \mathbf I-r_0(1-i\alpha)\frac{1}{\|\mathbf x\|^3}\mathbf x\mathbf x^\top
```

