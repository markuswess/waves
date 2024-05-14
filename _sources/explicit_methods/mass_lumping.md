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
### Mass lumping in higher dimensions

For higher dimensions the first order finite element spaces (for triangular or tetrahedral elements $T\in\mathcal T$) are again given by piecewise linear functions. A convenient basis is given by the according functions which vanish on all but one vertex points, i.e.,

```{math}
V_h = \{f\in H^1(\Omega): f|_T\in \mathcal P^1,\forall T\in\mathcal T\}
```
where the polynomial space $\mathcal P^k$ is defined by
```{math}
\mathcal P^k:=\{(x,y,z)\mapsto x^ny^mz^l:n+m+l\leq k\},
```
(and in two dimensions the exponent $l$ of the $z$ coordinate is always $0$).
If we denote the set of vertices by $V_0,\ldots,V_N$ the nodal basis is given by the condition
```{math}
b_j(V_i)=\delta_{j,i}.
```
Following the same reasoning as in one dimension we may construct integration rules locally on each element $T$ by
```{math}
I_T(f):=\frac{|T|}{d+1}\sum_{i=0}^{d+1} f(V_{T,i})
```
where $V_{T_i}$ are the vertices of the element $T$ and $d$ is the spacial dimension.
Again an approximation to the correct mass integrals is given by

```{math}
\tilde m_{i,j}:=\sum_{T\in\mathcal T} I_T(b_ib_j)=\delta_{i,j}\frac{1}{d+1}\sum_{k=1}^{K_i} |T_{i,k}|
```
where  $T_{i,k},k=1,\ldots,K_i$ are the elements sharing the vertex $V_i$ and $K_i$ is the number of these elements.

## Higher order mass lumping

Usually finite element basis functions are defined on a reference element first. Then they are appropriately transformed to the physical elements and basis functions on neighbouring elements are glued together to obtain conforming basis functions. I.e., if $T$ is a physical element (interval, triangle, quadrilateral, tetrahedron, hexahedron,...) we have a reference element $\hat T$ and a diffeomorphism  $\phi:\hat T\to T$. 
If a local basis function $\hat b:\hat T\to\mathbb R$ is given then the basis function on the phyiscal element is given by $b := \hat b \circ \phi ^{-1}$.
We use the following reference elements:
```{math}
\begin{aligned}
\hat T_{\mathrm{segm}} &= \mathrm{conv}\{0,1\}=[0,1],\\
\hat T_{\mathrm{trig}} &= \mathrm{conv}\{(0,0),(1,0),(0,1)\},\\
\hat T_{\mathrm{tet}} &= \mathrm{conv}\{(0,0,0),(1,0,0),(0,1,0),(0,0,1)\} ,\\
\hat T_{\mathrm{quad}} &= \mathrm{conv}\{(0,0),(1,0),(0,1),(1,1)\}=[0,1]^2,\\
\hat T_{\mathrm{hex}} &= \mathrm{conv}\{(0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,0,1),(1,0,1),(0,1,1),(1,1,1)\}=[0,1]^3,
\end{aligned}
```
where $\mathrm{conv}$ denotes the convex hull of a set.

For a given integration rule consisting of integration points $\hat x_j\in \hat T$ we use basis functions which vanish in all but one integration points, i.e., 
```{math}
\hat b_i(\hat x_j)=\delta_{i,j}.
```
The task is now to find suitable integration rules and spaces on the reference elements.

### Higher order mass lumping in one dimension

To be able to couple the mapped basis functions on two adjoining elements (intervals) we need to fix the endpoints $0,1$ of the reference interval $T_{\mathrm{segm}}$ as integration points. Then for a given total number of points $N+1$ the integration rule with highest order is the Gauss-Lobatto rule which is exact for polynomials of up to order $2N-1$.

As local spaces we may choose the spaces of polynomials of up to order $N$ with the Lagrangian basis 
```{math}
\hat b_j(x):= \prod_{k=0,k\neq j}^N\frac{x-\hat x_k}{\hat x_j-\hat x_k},
```
where $\hat x_j$ denote the integration points on the reference element. The basis functions are then mapped to the physical elements as described above and the first basis function of an element is glued to the last basis function of the previous element to obtain globally continuous basis functions.

### Higher order mass lumping for quadrilaterals and hexahedra

In the case of quadrilaterals and hexahedra integration rules may be obtained by tensorizing the one dimensional Gauss-Lobatto rules. I.e., if a Gauss-Lobatto rule is given by the points $\hat x_0,\ldots,\hat x_N$ and weights $\hat w_0,\ldots, \hat w_N$ we define
```{math}
\hat {\mathbf x}_{i+(N+1)j}:=(\hat x_i,\hat x_j)^\top,\quad \hat {\mathbf w}_{i+(N+1)j}:=\hat w_i\hat w_j,
```
on the unit square.
The according Lagrangian polynomial basis is given by
```{math}
\hat{\mathbf b}_{i+(N+1)j}(x,y) = \hat b_i(x)\hat b_j(y).
```
Note that here we use the polynomial space 
```{math}
\mathcal Q^N:=\{x^ky^l:k,l\leq N\}
```
which differs from $\mathcal P^N$.
For three dimensions the construction works similarly.

### Higher order mass lumping on simplexes

A naive ansatz to construct a second order mass-lumping method for triangles would be to follow the ideas from the first-order case: There we used the given basis and constructed a suitable integration rule. Since the dimension of $\mathcal P^2$ is $6$ and we need three symmetric degrees of freedom on each boundary to ensure continuity of the correctly coupled basis functions. The integration points on the reference triangle are fixed by
```{math}
\begin{aligned}
\hat{\mathbf x}_1&=(0,0),&
\hat{\mathbf x}_2&=(1/2,0),&
\hat{\mathbf x}_3&=(1,0),\\
\hat{\mathbf x}_4&=(1/2,1/2),&
\hat{\mathbf x}_5&=(0,1),&
\hat{\mathbf x}_6&=(0,1/2)\}.
\end{aligned}
```
The nodal basis functions are given by
 ```{math}
\begin{aligned}
\hat b_1 &= 2(1-x-y)(1-x-y-\frac{1}{2}),&
\hat b_2 &= 4(1-x-y)x,&
\hat b_3 &= 2x(x-\frac{1}{2})\\
\hat b_4 &= 4xy,&
\hat b_5 &= 2y(1-\frac{1}{2}),&
\hat b_6 &= 4y(1-x-y).
\end{aligned}
``` 
It remains to specify the integration weights to obtain an optimal accuracy. Due to symmetry we only have to specify two weights $w_2=w_4=w_6=w_E$ for the mid-points of the edges and $w_1=w_3=w_5=w_V$ for the vertices, i.e., our integration rule is given by
```{math}
I_{\hat T}(f)=w_V\left(f(\hat{\mathbf x}_1)+f(\hat{\mathbf x}_3)+f(\hat{\mathbf x}_5)\right)+w_E\left(f(\hat{\mathbf x}_2)+f(\hat{\mathbf x}_4)+f(\hat{\mathbf x}_6)\right).
```
The condition that constants are integrated exactly yields
```{math}
\frac{1}{2}=3w_V+3w_E.
```
Imposing that also the function $(x,y)\mapsto x$ is integrated exactly yields
```{math}
\int_{\hat T}xdx=\frac{1}{6}=w_V\left(0+1+0\right)+w_E\left(\frac{1}{2}+\frac{1}{2}+0\right)=w_V+w_E,
```
which is the same condition as before. Imposing that also the function $(x,y)=x^2$ is integrated exactly yields
```{math}
\int_{\hat T} x^2 dx = \int_0^1\int_0^{1-x} x^2 dydx = \int_0^1x^2-x^3 dx = \frac{1}{12}
=w_V\left(0+1+0\right)+w_E\left(\frac{1}{4}+\frac{1}{4}+0\right)=w_V+\frac{1}{2}w_E.
```
Inserting the first condition $w_E=1/6-w_V$ yields
```{math}
w_V = 0,\quad w_E = \frac{1}{6}.
```
Using this quadrature for approximating our mass matrix is not an option since we would end up with a singular approximated mass matrix. We may settle for an integration rule which only integrates linear functions exactly, however it turns out that this will also reduce the overall convergence order of the resulting method.

The remedy to this problem is to add an additional integration point to the integration rule. From symmetry considerations this can only be the center $\hat{\mathbf x}_7=(1/3,1/3)$ of the reference triangle.

We also have to enlargen the local space of polynomials by another basis function.

Since this has to be a basis function corresponding to the integration point in the center we have to choose a bubble function namely $\tilde b_7(x,y):= 27xy(1-x-y)\}.$
The resulting space 
```{math}
\hat V_h:=\mathcal P^2\oplus \mathrm{span}\{\tilde b_7\}=\mathrm{span}\left\{\hat b_1,\hat b_2,\hat b_3,\hat b_4,\hat b_5,\hat b_6,\tilde b_7\right\}
```
now lies in between the polynomial spaces $\mathcal P^2$ and $\mathcal P^3$. To obtain a Lagrangian basis with respect to the integration points one may simply use the Lagrangian basis of $\mathcal P^2$ with respect to the original integration points and subtract the bubble function, i.e.,
```{math}
\tilde b_j = \frac{\hat b_j+\tilde b_7}{1+\hat b_j(1/3,1/3)},\quad j=1,\ldots,6.
```
It remains to specify the correct weights for the integration rule. 

````{prf:Theorem}
Let $\tilde b_j$, $\mathbf x_j$ be given as above. Moreover we define 
```{math}
w_1=w_3=w_5=\frac{1}{40},\quad w_2=w_4=w_6=\frac{1}{15},\quad w_7 = \frac{9}{40}.
```
Then the quadrature rule given by $\mathbf x_j,\omega_j$ is exact for all functions from $\mathcal P^3$.
````

% ````{prf:proof}
% Due to symmetry it is sufficient to show that the integrals 
% ````


````{prf:remark}
It turns out that a sufficient condition for the integration rule to not spoil the convergence of the method is that it is exact for $\mathcal P^{k+\tilde k-2}$ if the space on the reference element $\hat V_h$ is given by $\mathcal P^k\oplus \tilde V_h$ with $\tilde V_h\subset \mathcal P^{\tilde k}$. This means that the integration rule defined above fo the second order mass lumping space is sufficient.
````
