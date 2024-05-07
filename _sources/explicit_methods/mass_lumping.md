(mass_lumping)=
# Mass lumping for the wave equation

A different approach to DG methods from {numref}`dg` to obtain easily invertible mass matrices is the so called mass lumping. We consider the one-dimensional case first

(mass_lumping_1o)=
## Mass lumping for linear finite elements 

We focus on the simple case of linear finite elements, first in one and subsequently in higher dimensions.

### Mass lumping for linear 1d elements

Suppose our spacial domain is given by the interval $\Omega=(a,b)$ decomposed into $N\in\mathbb N$ subintervals $T_j:=(a_{j-1},a_j):=(a_{j-1},a_{j-1}+d_j$ with $a=a_0<a_j<\ldots<a_N=b$ and $d_j=a_{j}-a_{j-1}$ for $j=1,\ldots, N$.

Then the classical hat function basis of a first order finite element space is given by
```{math}
V_h = \mathrm{span}\{b_j, j = 0,\ldots, N\},\quad b_j =\begin{cases}\frac{x-a_{j-1}}{d_j},&x\in T_j\\\frac{a_{j+1}-x}{d_{j+1}},&x\in T_{j+1}\\0,&\mathrm{else}\end{cases}
```
with the appropriate restrictions for $b_0,b_N$.

The entries of the mass matrix can be easily computed by hand via
```{math}
m_{i,j}=\int_a^b b_j b_i = \begin{cases}
\int_{a_{j-1}}^{a_{j+1}}b_j^2=(d_j+d_{j+1})\int_0^1 x^2 dx=\frac{d_j+d_{j+1}}{3},&i=j,\\
\int_{a_j}^{a_{j+1}}b_j b_{j+1}=d_{j+1}\int_0^1 x(1-x)dx = \frac{d_{j+1}}{6},&i=j+1,\\
\int_{a_{j-1}}^{a_{j}}b_j b_{j-1}=d_{j}\int_0^1 x(1-x)dx = \frac{d_{j}}{6},&i=j-1,\\
0,&\mathrm{else,}\end{cases}
```
where we set $d_0=d_{N+1}=0$.
Instead of computing the integrals for the entries of the mass matrix exactly we may approximate them using the trapezoidal rule
```{math}
\int_a^b f \approx I_{a,b}f:= \frac{b-a}{2}(f(a)+f(b)).
```
which leads to an approximated mass matrix $\tilde{\mathbf M}$ with entries 
```{math}
\begin{aligned}
\tilde m_{i,j}&=I_{a,b} (b_j b_i) \\
&=\sum_{k=1}^N I_{a_{k-1},a_k}(b_j,b_i)\\
&=\sum_{k=1}^N \frac{d_k}{2}(b_j(a_{k-1})b_i(a_{k-1})+b_j(a_{k})b_i(a_{k}))\\
&=\frac{d_1}{2}b_i(a_0)b_j(a_0)+\sum_{k=1}^{N-1}\frac{d_k+d_{k+1}}{2}b_j(a_k)b_i(a_k)+\frac{d_N}{2}b_j(a_N)b_i(a_N)\\
&=
\begin{cases}
\frac{d_i+d_{i+1}}{2},& i=j,\\
0,&\mathrm{else.}\end{cases}
\end{aligned}
```
% ### Mass lumping in higher dimensions
