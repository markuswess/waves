# Time-harmonic waves

````{prf:Example} 
:label: example_1d
Consider the following simple one-dimensional wave example to find $u:[0,T]\times[0,\infty)\to\mathbb R$ such that
```{math}
:label: example_1d_td
\begin{aligned}
\partial_t^2 u-c^2\partial_x^2 u &= 0,&& \text{in }(0,T)\times (0,\infty),\\
u(\cdot,0) &= a,&&\text{in } (0,T],\\
u(0,\cdot ) &=0,&&\text{in } (0,\infty),\\
\partial_t u(0,\cdot ) &=0,&&\text{in } (0,\infty),
\end{aligned}
```
for a smooth function $a:[0,\infty)\to\mathbb R$ with $a(0)=a'(0)=a''(0)=0$. Then the (unique, strong) solution is given by 
```{math}
:label: example_1d_sol
u(t,x) = \begin{cases}a(t-x/c),&ct> x,\\ 0,&ct\leq x\end{cases}.
```
More specifically consider $a$ to be a "smoothened" cosine, such that
$a(t) = \cos(-\omega t)$ for $t>\varepsilon$  and some $\omega>0$ and $a(0)=a'(0)=a''(0)=0$. Then for any bounded interval $[0,b]$ for some $b>0$ we have that for $t>b/c+\varepsilon$ the solution $u$ can be written as 
```{math}
u(t,x) = \cos(-\omega(t-x/c))=\Re(\exp(-i\omega t)\exp(i\omega/c x))
```
````


In applications often scenarios are of interest where the solutions are (approximately) time-harmonic waves, i.e., are(approximately) of the form 
```{math}
:label: time_harm_wave
u(t,x) = \Re\left(\exp(-i\omega t) \hat u(x)\right),
```
where $\omega >0$ is a given angular frequency and $u(x)$ the amplitude as function in space.


Thus it is not necessary to simulate the full time-domain systems from {numref}`acoustic_td` {numref}`elastic_td` {numref}`electromagnetic_td` which in their second-order form can all be stated as
```{math}
:label: gen_second_order_wave
\begin{aligned}
\frac{1}{c^2}\partial_t^2 u-Du &= f
\end{aligned}
```
where $D$ is one of the second-order differential operators in space from the previous sections. For simplicity we assume homogeneous boundary conditions.
Inserting the time harmonic wave {eq}`time_harm_wave` into {eq}`gen_second_order_wave` for a time harmonic right hand side
```{math}
f(t,x)=\exp(-i\omega t)\hat f(x)
```
leads to the equation
```{math}
:label: gen_helmholtz
-\frac{\omega^2}{c^2} \hat u-D\hat u = \hat f,
```
````{prf:remark}
again with homogeneous boundary conditions.
Note that we neglected the initial conditions which are necessary to close the time-domain system {eq}`gen_second_order_wave`. They have to be chosen fitting to a time-harmonic solution.
````

In the case of acoustic waves (i.e., $D=\Delta$) with a constant wave speed $c$ this leads to the inhomogenous Helmholtz equation
````{card}
```{math}
:label: helmholtz
-k^2\hat u-\Delta\hat u = \hat f,
```
````
where the *wave number* $k>$ is defined by $k=\omega/c$.
When exclusively the time-harmonic regime is of interest the $\hat\cdot$ is usually ommitted.

````{prf:remark}
A more general way to derive the time-harmonic counterparts to the time-domain equations is to introduce the Fourier transform of a function $h\in L^1(\mathbb R)$ by
```{math}
\hat h(\omega) = \int_{-\infty}^\infty h(t)\exp(-i\omega t) dt,
```
and apply it to the space-time solutions of the time-domain wave problems (with respect to the time variable).
````

````{prf:Remark}
Note that the sign convention to use $\exp(-i\omega t)$ in the Fourier transform and the time harmonic ansatz {eq}`time_harm_wave` is arbitrary. Flipping the sign does **not** change the resulting time-harmonic equation. However, it corresponds to a time-reversal, thus whenever time-domain solutions, first-order time derivatives or (non-squared) factors $\omega$ are used in any reasoning this has to be taken into account.
````


## (Non)-locality

````{prf:Example}  ctd.
Considering again the one-dimensional example {prf:ref}`example_1d` we immediately obtain from {eq}`example_1d_sol` that for fixed $t_0>0$ we have $\mathrm{supp}(u(t_0,\cdot)) \subset [0,ct_0)$.
```` 
Looking e.g., at D'Alembert's solution for the one-dimensional wave equation or the example above it is comprehensible that time-domain waves have a finite speed of propagation, e.g., for initial conditions with local support inside of a compact set $\Omega$ in space and a finite time interval $[0,T]$ there exists another compact set $\Omega'$ in space such that for every $t\in [0,T]$ the support of the solution is contained in $\Omega'$.

This can be made mathematically rigorous by looking at fundamental solutions of the respective time-domain wave equations. Similar properties hold for a compactly supported (in space) forcing term.

Thus, mathematically, the problem of solving the wave-equation in the free space for finite times can always be re-stated as solving the wave-equation on a large enough bounded domain with homogeneous boundary conditions. This heuristic argument already points out that for mathematical formulations (and numerical approximations) of time-domain waves one can always avoid the difficulties that unbounded domains cause (see example below).

````{prf:Example} ctd.
:label: example_1d_fqd
The time-harmonic counterpart to {eq}`example_1d_td` can be stated as
```{math}
-\omega^2/c^2 u-u'' &= 0,&&\text{in }(0,\infty),\\
u(0)&=1,
```
However, although this equation is solved by $u(x)=\exp(i\omega/c x)$ the solution is not unique since also $u_-(x)=\exp(-i\omega/c x)$ (and every linear combination of $u,u_-$) is also a solution. Thus the time-harmonic equation needs to be supplied with another (boundary) condition to obtain well-posedness.

Transferring the solution $u_-$ back to time-domain, i.e.,
```{math}
\Re(\exp(-i\omega t)u_-(x))=\Re(\exp(-i\omega(t+x/c)))=\cos(\omega(t+x/c))
```
shows that $u_-$ corresponds to a time-harmonic wave travelling from right to left.
````


## Sommerfeld radiation condition

To obtain unique solvability of the time-harmonic equations stated above we need to pose an additional condition. As motivated in {prf:ref}`example_1d` and {prf:ref}`example_1d_fqd` this condition should separate incoming from outgoing waves. 

% In the following we focus on the wave equation in three dimensions.

The sought condition is the so-called Sommerfeld radiation condition named after Arnold Sommerfeld which is given by 
```{math}
:label: sommerfeld
\lim_{|x|\to\infty} |x|^{(d-1)/2}\left(\frac{\partial u}{\partial|x|}-iku\right) = 0
```
where $d$ is the spacial dimension.


## Scattering by an obstacle

A typical application of {eq}`helmholtz` on an unbounded domain is the so-called *scattering problem*. Thereby incident field $u^i$ on $\mathbb R^d$ which can be e.g. a plane wave
```{math}
u^i(t,x) = \exp(ik e\cdot x-\omega t)
```
for a given vector $e\in\mathbb R^d$ with $|e|=1$ and an angular frequency $\omega>0$ and a compactly supported obstacle $\Omega\subset\mathbb R^d$ is given. The *total field* $u:\mathbb R^d\setminus \Omega$ can be decomposed into the incident field $u^i$ and a scattered field $u^s$ and should satisfy
```{math}
\begin{aligned}
-k^2 u-\Delta u&=0,&\text{on }&\mathbb R^d\setminus\bar\Omega,\\
u &= u^i+u^s,&\text{on }&\mathbb R^d\setminus\bar\Omega,\\
u&=0,&\text{on }&\partial\Omega,\\
\lim_{|x|\to\infty} |x|^{(d-1)/2}\left(\frac{\partial u^s}{\partial|x|}-iku^s\right) &= 0.
\end{aligned}
```
By linearity it is clear that the scattered field also satisfies the Helmholtz equation with inhomogeneous Dirichlet boundary conditions on $\partial\Omega$ and the Sommerfeld radiation condition for $|x|\to\infty$, i.e.,

```{math}
\begin{aligned}
-k^2 u^s-\Delta u^s&=0,&&\text{on }&\mathbb R^d\setminus\bar\Omega,\\
u&=u^i,&&\text{on }&\partial\Omega,\\
\lim_{|x|\to\infty} |x|^{(d-1)/2}\left(\frac{\partial u^s}{\partial|x|}-iku^s\right) &= 0.
\end{aligned}
```

% ## The resonance problem

% Resonance problems in general are problems where not only the function $u$ but also the frequency $\omega$ is unknown. The usually lead to eigenvalue problems.

% ### Resonances on bounded domains

% Let $\Omega\subset \mathbb R^d$ be a bounded convex Lipschitz domain. Then it is well known that the 


% ### Resonances on unbounded domains
