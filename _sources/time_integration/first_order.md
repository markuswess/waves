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
for positive (semi-)definite matrices $\mathcal D,\mathcal M \in\mathbb R^{(N_u+N_v)\times (N_u+N_v)}$ and a skew-symmetric matrix $\mathcal B\in\mathbb R^{N_u\times N_v}$, where we set $\mathbf w = (\mathbf u^\top,\mathbf v^\top)^\top$ and $\mathcal F= (\mathbf f^\top,0)^\top$.
````

Similar to the {ref}`time_integration_so` we have energy conservation in the form

````{prf:Theorem}
Let $\mathbf u,\mathbf v$ solve {eq}`semi_disc_1o`. Then
```{math}
\frac{d}{d t} E = -(\mathbf u^\top\mathbf D_u\mathbf u  +\mathbf v^\top\mathbf D_v\mathbf v)+\mathrm u^\top \mathrm f,
```
where the energy $E$ is defined by
```{math}
E=\frac{1}{2}\left(\mathbf u^\top\mathbf M_u{\mathbf u}+\mathbf v^\top\mathbf M_v\mathbf v\right).
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

% (leap_frog)= 
% ## The Leap Frog time stepping
%
% While the Crank-Nicholson scheme can be easily derived for general first-order in time equations the Leap Frog scheme exploits the underlying skew symmetric structure. Thus, we assume that no damping is present
