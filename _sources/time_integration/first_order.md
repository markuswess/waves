(time_integration_fo)=
# Time-integration for the first order system

Applying a method of lines approach described in {numref}`mol` to the acoustic, elastic and electromagnetic time-domain wave problems (with possibly added damping or some variant of not entirely reflecting boundary) in first order form leads to a semi-discrete problem of the form (denoting the time derivative by a dot)

```{math}
:label: semi_disc_1o
\begin{aligned}
\mathbf M_u\dot{\mathbf u}&=-\mathbf D_u{\mathbf u}+\mathbf B\mathbf v +\mathbf f,\\
\mathbf M_v\dot{\mathbf v}&=-\mathbf D_v{\mathbf v}-\mathbf B^\top\mathbf u,\\
\mathbf u(0)&=\mathbf u^0,\quad {\mathbf v}(0)=\mathbf v^0
\end{aligned}
```
for some symmetric positive definite mass and damping matrices $\mathbf M_u, \mathbf D_u\in\mathbb R^{N_u\times N_u}$, $\mathbf M_v,\mathbf D_v\in\mathbb R^{N_v\times N_v}$, the discrete differential operator $\mathbf B\in\mathbb R^{N_u\times N_v}$ and the unknown vector functions $\mathbf u\in C^1([0,T];\mathbb R^{N_u}),\mathbf v\in C^1([0,T];\mathbb R^{N_v})$.

````{prf:Remark}
The system {eq}`semi_disc_1o` can also be rewritten more compactly by
```{math}
:label: semi_disc_1o_compact
\mathcal M \dot{\mathbf w} =-\mathcal D\mathbf w + \mathcal B\mathbf w +\mathcal F,
```
for positive (semi-)definite matrices $\mathcal D,\mathcal M \in\mathbb R^{N\times N}$ with $N:=N_u+N_v$ and a skew-symmetric matrix $\mathcal B\in\mathbb R^{N\times N}$, where we set $\mathbf w = (\mathbf u^\top,\mathbf v^\top)^\top$, $\mathcal F= (\mathbf f^\top,0)^\top$and
```{math}
\mathcal B = \begin{pmatrix} 0&\mathbf B\\-\mathbf B^\top&0\end{pmatrix}.
```
````

Similar to the {ref}`time_integration_so` we have energy conservation in the form

````{prf:Theorem}
Let $\mathbf u,\mathbf v$ solve {eq}`semi_disc_1o`. Then
```{math}
\frac{d}{d t} E(\mathbf u,\mathbf v) = -(\mathbf u^\top\mathbf D_u\mathbf u  +\mathbf v^\top\mathbf D_v\mathbf v)+\mathrm u^\top \mathrm f,
```
where the energy $E(\mathbf u,\mathbf v)$ is defined by
```{math}
E(\mathbf u,\mathbf v)=\frac{1}{2}\left(\mathbf u^\top\mathbf M_u{\mathbf u}+\mathbf v^\top\mathbf M_v\mathbf v\right).
```
In particular, if there is no damping ($\mathbf D_v,\mathbf D_u=0$) and external sources $\mathbf f = 0$ the energy $E$ is constant.
````
````{prf:proof}
We multiply {eq}`semi_disc_1o` from the left by $(\mathrm u,\mathrm v)^\top$ to obtain
```{math}
\mathbf u^\top\mathbf M_u\dot{\mathbf u} +\mathbf v^\top\mathbf M_v\dot{\mathbf v} =-\mathbf u^\top\mathbf D_u\mathbf u-\mathbf v^\top\mathbf D_v\mathbf v- \mathbf u^\top\mathbf B\mathbf v-\mathbf v^\top\mathbf B^\top\mathbf u + \mathbf u^\top\mathbf f
```
which yields the claim since
```{math}
\frac{d}{dt}\left(\mathbf g(t)^\top\mathbf A\mathbf g(t)\right)=
\dot{\mathbf g}(t)^\top\mathbf A\mathbf g(t)+
\mathbf g(t)^\top\mathbf A\dot{\mathbf g}(t)
=
2\dot{\mathbf g}(t)^\top\mathbf A\mathbf g(t),
```
for symmetric $\mathbf A$ and arbitrary $\mathbf g$.
````
(crank_nicholson)= 
## The Crank Nicholson time stepping


The Crank Nicholson time stepping can motivated for systems of the form {eq}`semi_disc_1o_compact` as follows:
As in the previous section for a given timestep $\tau>0$  we denote the discrete approximation at $t=j\tau$ for $j\in\mathbb N$ by $\mathbf w_j$. Then
```{math}
\mathcal M\mathbf w(t_{j+1})=\mathcal M\mathbf w(t_j)+\int_{t_j}^{t_j+1}\mathcal M\dot{\mathbf w}(s)ds=
\mathcal M\mathbf w(t_j)+\int_{t_j}^{t_j+1}(-\mathcal D+\mathcal B)\mathbf w(s)+\mathcal F(s)ds.
```
Now approximating the integral on the right hand side by the trapezoidal rule (and replacing the exact values of $\mathbf w$ by the approximations) we obtain the relation
```{math}
\mathcal M\mathbf w_{j+1}=\mathcal M\mathbf w_{j}+\frac{\tau}{2}\left(-\mathcal D+\mathcal B\right)\left(\mathbf w_{j}+\mathbf w_{j+1}\right)+\frac{1}{2}\left(\mathcal F_j+\mathcal F_{j+1}\right).
```
and thus
```{math}
\begin{aligned}
\mathbf w_{j+1}&=\left(\mathcal M+\frac{\tau}{2}(\mathcal D-\mathcal B)\right)^{-1}\left(\mathcal M-\frac{\tau}{2}(\mathcal D-\mathcal B)\right)\mathbf w_j+\frac{1}{2}\left(\mathcal M+\frac{\tau}{2}(\mathcal D-\mathcal B)\right)^{-1}\left(\mathcal F_j+\mathcal F_{j+1}\right)\\
&=\mathbf w_j-\tau\mathcal S^{-1}(\mathcal D-\mathcal B)\mathbf w_j+\frac{\tau}{2}\mathcal S^{-1}\left(\mathcal F_j+\mathcal F_{j+1}\right).
\end{aligned}
```
with $\mathcal S=\mathcal M+\frac{\tau}{2}(\mathcal D-\mathcal B)$.

The Crank-Nicholson time stepping is **implicit** since the inverse of a matrix, different from the mass matrix has to be applied in each time step.

````{prf:theorem}
Let $\mathbf w_j^\top=(\mathbf u_j^\top,\mathbf v_j^\top)$ for $j\in\mathbf N$ be the discrete Crank-Nicholson approximation of {eq}`semi_disc_1o` (or {eq}`semi_disc_1o`) with $\mathbf D, \mathbf f=0$.

Then
```{math}
E(\mathbf u_j,\mathbf v_j)=E(\mathbf u_0,\mathbf v_0),
```
for all $j\in\mathbb N$.
````
````{prf:proof}
It follows from basic linear algebra that the eigenpairs of 
```{math}
\lambda\mathcal M \phi = \mathcal B \phi
```
consist of pairs of imaginary eigenvalues $\lambda_l, \lambda_{-l}=-\lambda_l\in i\mathbb R$, $l=1,2,\ldots,N$ with conlugate complex eigenfunctions $\phi_l,\phi_{-l}=\bar\phi_l$.
We fix $\mathbf w^\top=(\mathbf u^\top,\mathbf v^\top)=\mathbf w_j^\top$ and $\tilde{\mathbf w}^\top=(\tilde{\mathbf u}^\top,\tilde{\mathbf v}^\top)=\mathbf w_{j+1}^\top$. 
Decomposing $\mathbf w,\tilde{\mathbf w}$ into the according orthonormal (with respect to $\mathcal M$) eigensystem i.e.,
```{math}
\begin{aligned}
\mathbf w&=\sum_{l=- N}^ N c_l\phi_l,\\
\tilde{\mathbf w}&=\sum_{l=- N}^ N \tilde c_l\phi_l,
\end{aligned}
```
with $c_l = \bar c_{-l}\in\mathbb C$ (and for convenience we set $c_0,\phi_0=0$) 
and inserting the expansions into the Crank-Nicholson stepping leads to
```{math}
\sum_{l=- N}^ N \left(\mathcal M-\frac{\tau}{2}\lambda_l\right)\tilde c_l\phi_l=\left(\mathcal M-\frac{\tau}{2}\mathcal B\right)\sum_{l=- N}^ N \tilde c_l\phi_l=\left(\mathcal M+\frac{\tau}{2}\mathcal B\right)\sum_{l=- N}^ N c_l\phi_l=\sum_{l=- N}^ N \left(\mathcal M+\frac{\tau}{2}\lambda_l\right)c_l\phi_l.
```
Multiplying this by $\phi_k^H$ from the left yields due to the orthonormality of the eigenbase 
```{math}
\left(1-\frac{\tau}{2}\lambda_k\right)\tilde c_k=\left(1+\frac{\tau}{2}\lambda_k\right)c_k.
```
i.e.,
```{math}
\tilde{c}_k=\frac{1+\frac{\tau}{2}\lambda_k}{1-\frac{\tau}{2}\lambda_k}c_k.
```
Since $\lambda_k$ is purely imaginary we have $|\tilde c_k|=|c_k|$ since for any complex number $z$ we have $|z/\bar z| = 1$.

Thus
```{math}
E(\mathbf u,\mathbf v) =\frac{1}{2} \mathbf w^\top\mathcal M\mathbf w  =\frac{1}{2} \left(\sum_{l=1}^ N c_l\phi_l^\top+\bar c_{l}\bar\phi_l^\top\right)\mathcal M \left(\sum_{l=1}^ N c_l\phi_l+\bar c_{l}\bar\phi_l\right)=\sum_{l=1}^ N|c_l|^2=\sum_{l=1}^ N|\tilde c_l|^2=E(\tilde{\mathbf u},\tilde{\mathbf v})
```
````

 (leap_frog)= 
 ## The Leap Frog time stepping

While the Crank-Nicholson scheme can be easily derived for general first-order in time equations the Leap Frog scheme exploits the underlying skew symmetric structure. 

The basic idea is to approximate $\mathbf v$ at half steps namely to set
```{math}
:label: leap_frog
\mathbf v_{1/2}&=\mathbf v_0 -\frac{\tau}{2} \mathbf M_v^{-1}\mathbf B^\top\mathbf u_0,\\
\mathbf u_{j+1}&=\mathbf u_{j} +\tau \mathbf M_u^{-1}\mathbf B\mathbf v_{j+1/2},\\
\mathbf v_{j+1/2}&=\mathbf v_{j-1/2} -\tau \mathbf M_v^{-1}\mathbf B^\top\mathbf u_j.
```

````{prf:remark}
To be able to evaluate also the velocity at full time-steps one may rewrite the method as

```{math}
:label: leap_frog_whole
\mathbf v_{j+1/2}&=\mathbf v_j -\frac{\tau}{2} \mathbf M_v^{-1}\mathbf B^\top\mathbf u_j,\\
\mathbf u_{j+1}&=\mathbf u_{j} +\tau \mathbf M_u^{-1}\mathbf B\mathbf v_{j+1/2},\\
\mathbf v_{j+1}&=\mathbf v_{j+1/2} -\frac{\tau}{2} \mathbf M_v^{-1}\mathbf B^\top\mathbf u_{j+1}.
```
This version is also known as kick-drift-kick formulation.
````

````{prf:theorem}
Let $\mathbf u_j,\mathbf v_j$ be the Leap Frog approximations applied to {eq}`semi_disc_1o` with $\mathbf D_u,\mathbf D_v,\mathbf f = 0$ and $\hat{\mathbf u}_j,\dot{\mathbf u}_j$ be the approximations by the Verlet time stepping applied to the equation
```{math}
:label: second_order_fo
\begin{aligned}
\mathbf M_u \ddot{\mathbf u}+\mathbf B\mathbf M_v^{-1}\mathbf B^\top\mathbf u &= 0,\\
\mathbf u(0) &= \mathbf u^0,\\
 \dot{\mathbf u}(0) &= \mathbf M_u^{-1}\mathbf B \mathbf v^0.
\end{aligned}
```
Then $\mathbf u_j=\hat{\mathbf u}_j$ for all $j\in\mathbb N$.
````

````{prf:proof}
Multiplying the first and last line of {eq}`leap_frog_whole` by $\mathbf M_u^{-1}\mathbf B$ yields
```{math}
\begin{aligned}
\mathbf M_u^{-1}\mathbf B\mathbf v_{j+1/2}&=\mathbf M_u^{-1}\mathbf B\mathbf v_j -\frac{\tau}{2} \mathbf M_u^{-1}\mathbf B\mathbf M_v^{-1}\mathbf B^\top\mathbf u_j,\\
\mathbf u_{j+1}&=\mathbf u_{j} +\tau \mathbf M_u^{-1}\mathbf B\mathbf v_{j+1/2},\\
\mathbf M_u^{-1}\mathbf B\mathbf v_{j+1}&=\mathbf M_u^{-1}\mathbf B\mathbf v_{j+1/2} -\frac{\tau}{2} \mathbf M_u^{-1}\mathbf B\mathbf M_v^{-1}\mathbf B^\top\mathbf u_{j+1}.
\end{aligned}
```

On the other hand applying the equivalent Verlet algorithm {eq}`newmark_explicit_leap` to the second order problem  {eq}`second_order_fo` yields
```{math}
\begin{aligned}
\tilde{\mathbf u}&=\dot{\mathbf u}_j-\frac{\tau}{2}\mathbf M_u^{-1}\mathbf B\mathbf M_v^{-1}\mathbf B^\top{\mathbf u}_j,\\
\mathbf u_{j+1}&={\mathbf u}_j+\tau\tilde{\mathbf u},\\
\dot{\mathbf u}_{j+1}&=\tilde{\mathbf u}-\frac{\tau}{2}\mathbf M_u^{-1}\mathbf B\mathbf M_v^{-1}\mathbf B^\top{\mathbf u}_{j+1}.
\end{aligned}
```
Thus with $\tilde {\mathbf u}:=\mathbf M^{-1}_u\mathbf B\mathbf v_{j+1/2}$ and $\dot{\mathbf u}_j:=\mathbf M_u^{-1}\mathbf B \mathbf v_{j}$ the algorithms are identical (note that also the inital data coincide).
````
By the above theorem we immediately obtain that the Leapfrog time stepping inherits all of the properties of the Verlet time stepping, in particular the preservation of the modified energy and the CFL condition.
