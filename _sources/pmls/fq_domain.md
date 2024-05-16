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
\tilde x(x):=\begin{cases}x,&x<a\\x+i\alpha(x-a),&x\geq a\end{cases}
```
we have that
```{math}
:label: u_scaled
\tilde u(x):=u(\tilde x(x))=\begin{cases}
\tilde c\exp(ikx),&x<a,\\
\tilde c\exp(-\alpha k(x-a))\exp(ikx),&x\geq a.
\end{cases}
```
Note that the *unwanted* term $\exp(-ikx)$ in the general solution {eq}`gen_sol` on the other hand increases exponentially along $\tilde x$.
We have found a way to alter the solution of the original problem in a way that it remains unchanged in the interval of interest $(0,a)$ and the Sommerfeld radiation condition reduces to the condition that the altered solution is square integrable on the real axis. This condition fits nicely into the desired variational framework.

To obtain a numerical method one now alters the Helmholtz equation such that {eq}`u_scaled` is the (unique) solution. Applying the chain rule one obtains easily:
````{prf:Theorem}
The function $\tilde u$ {eq}`u_scaled` is the (unique) solution to the problem to find $\tilde u\in C^2([0,a)\cup (a,\infty))\cap L^2((0,\infty))$ such that
```{math}
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
:label: weak_scaled
\int_0^a \tilde u' v'+\frac{1}{1+i\alpha}\int_a^T \tilde u'v'-k^2\int_0^a \tilde uv - k^2(1+i\alpha)\int_a^T \tilde uv=\int_0^a fv
```
for all $v\in H^1_0((0,T))$.
 Since the solution decays exponentially for $x\to\infty$ we expect an exponentially small (with respect to the truncation $T$) error.


% ## PMLs for waveguides
% 
% The ideas from 1d can be easily carried over to waveguide geometries: Consider the Helmholtz equation on the semi infinite waveguide $\Omega=(0,\infty)\times (0,\pi)$ with Dirichlet boundary conditions.
% A decomposition of the solution $u(x,y)$ into transversal modes yields
% ```{math}
% u(x,y)=\sum_{n=1}^\infty \sin(ny) u_n(x),\quad f(x,y)=\sum_{n=1}^\infty \sin(ny)f_n(x),
% ```
% and the functions $u_n$ have to solve
% ```{math}
% -u_n''-(k^2-n^2)u_n = f_n.
% ```
% Again we assume that $f_n$ vanishes for $x>a$ and $k\notin\mathbb N$. Then the family of general solutions $u_n$ for $x>a$ is given by
% ```{math}
% u_n(x) = \begin{cases}
% c_+\exp(i\sqrt{k^2-n^2}x)+c_-\exp(-i\sqrt{k^2-n^2}x),&n<k\\
% c_+\exp(\sqrt{n^2-k^2}x)+c_-\exp(-\sqrt{n^2-k^2}x),&n>k.
% \end{cases}
% ```
% Solutions where $n<k$ are called *propagating* waves, and the radiation condition is 
