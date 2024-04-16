(time_integration_so)=
# Time-integration for the second order systems

Applying a method of lines approach described in {numref}`mol` to the acoustic, elastic and electromagnetic time-domain wave problems (with possibly added damping or some variant of not entirely reflecting boundary) leads to a semi-discrete problem of the form (denoting the time derivative by a dot)
```{math}
:label: semi_disc_2o
\mathbf M\ddot{\mathbf u} + \mathbf C\dot{\mathbf u}+\mathbf K\mathbf u=\mathbf f,\quad
\mathbf u(0)=\mathbf u^0,\quad \dot{\mathbf u}(0)=\mathbf v^0
```

where $\mathbf u\in C^2([0,T];\mathbb R^N)$ and $N\in\mathbb N$ is the dimension of the discrete space. The matrices $\mathbf M,\mathbf C,\mathbf K\in\mathbb R^{N\times N}$ are the discrete mass-, damping- and stiffness matrix and $\mathbf f\in C([0,T];\mathbb R^N)$


## Energy conservation
The semi-discrete system {eq}`semi_disc_2o` conserves the discrete energy:

````{prf:Theorem} Discrete energy conservation
Let $\mathbf u$ solve {eq}`semi_disc_2o` with symmetric matrices $\mathbf M, \mathbf C, \mathbf K$. Then 
```{math}
\dot {\mathbf E}(\mathbf u):=\frac{d}{dt}\frac{1}{2}\left(\dot{\mathbf u}^\top\mathbf M\dot{\mathbf u}+\mathbf u^\top\mathbf K\mathbf u\right) = -\dot{\mathbf u}^\top\mathbf C\dot{\mathbf u}+\dot{\mathbf u}^\top\mathbf f.
```
In particular, if there are no external sources ($\mathbf f =0$) and damping ($\mathbf C = 0$), the discrete energy $\mathbf E$ is conserved over time.
````
````{prf:proof}
Multiplying {eq}`semi_disc_2o` by $\dot{\mathbf  u}^\top$ from the left yields
```{math}
\dot {\mathbf u}^\top\mathbf M\ddot{\mathbf u} +\dot{\mathbf  u}^\top \mathbf C\dot{\mathbf u}+\dot{\mathbf  u}^\top\mathbf K\mathbf u=\dot{\mathbf  u}^\top\mathbf f.
```
Since for any symmetric matrix $\mathbf A$ and vector function $\mathbf g\in C^1([0,T];\mathbb R^N)$
```{math}
\frac{d}{dt}\left(\mathbf g(t)^\top\mathbf A\mathbf g(t)\right)=
\dot{\mathbf g}(t)^\top\mathbf A\mathbf g(t)+
\mathbf g(t)^\top\mathbf A\dot{\mathbf g}(t)
=
2\dot{\mathbf g}(t)^\top\mathbf A\mathbf g(t)
```
the claim follows.
````

````{prf:Remark}
Usually the matrices $\mathbf M,\mathbf C,\mathbf K$ are assume to be positive (semi-)definite. Then it follows that the energy $\mathbf E$ is always non-negative.
````

For the time integration of wave problems one usually tries to preserve energy.

## The Newmark $\beta-$method

For the given time interval $[0,T]$ for $T>0$ we want to approximate the solution $\mathbf u$ to {eq}`semi_disc_2o` at $M+1$ discrete timesteps $t_j=j\tau$ for $j=0,\ldots,M$, $\tau = T/M$. We write $\mathbf u_j$ for the approximation of $\mathbf u(t_j)$ (and use similar notation for the derivatives).

Assuming sufficient regularity we have by the mean value theorem that for a fixed $j$ there exist $\xi_1,\xi_2\in[t_j,t_{j+1}]$ such that 
```{math}
\dot{\mathbf u}(t_{j+1})=\dot{\mathbf u}(t_j)+\tau\ddot{\mathbf u}(\xi_1),\qquad
\mathbf u(t_{j+1})={\mathbf u}(t_j)+\tau\dot{\mathbf u}(t_j)+\frac{\tau^2}{2}\ddot{\mathbf u}(\xi_2).
```
Approximating $\ddot{\mathbf u}(\xi_1),\ddot{\mathbf u}(\xi_2)$ by convex combinations of the approximated values $\dot{\mathbf u}_j,\dot{\mathbf u}_{j+1}$, i.e., for fixed $\gamma,2\beta\in [0,1]$
```{math}
\ddot{\mathbf u}(\xi_1)\approx (1-\gamma)\ddot{\mathbf u}_j+\gamma\ddot{\mathbf u}_{j+1},\qquad
\ddot{\mathbf u}(\xi_2)\approx (1-2\beta)\ddot{\mathbf u}_j+2\beta\ddot{\mathbf u}_{j+1}
```
leads to the discrete approximations
````{card}
```{math}
:label: newmark_approx
\begin{aligned}
\dot{\mathbf u}_{j+1}&=\dot{\mathbf u}_j+\tau\left((1-\gamma)\ddot{\mathbf u}_j+\gamma\ddot{\mathbf u}_{j+1}\right),\\
\mathbf u_{j+1}&={\mathbf u}_j+\tau\dot{\mathbf u}_j+\tau^2\left(\left(\frac{1}{2}-\beta\right)\ddot{\mathbf u}_j+\beta\ddot{\mathbf u}_{j+1}\right).
\end{aligned}
```
````
Additionally approximating ${\mathbf u}(t_{j+1})\approx{\mathbf u}_{j+1}$, $\dot{\mathbf u}(t_{j+1})\approx\dot{\mathbf u}_{j+1}$, $\ddot{\mathbf u}(t_{j+1})\approx\ddot{\mathbf u}_{j+1}$ in {eq}`semi_disc_2o` yields the third equation

```{math}
:label: newmark_state
\mathbf M\ddot{\mathbf u}_{j+1} + \mathbf C\dot{\mathbf u}_{j+1}+\mathbf K\mathbf u_{j+1}=\mathbf f_{j+1} 
```
 to close the linear system for the unknowns $\mathbf u_{j+1},\dot{\mathbf u}_{j+1},\ddot{\mathbf u}_{j+1}$ for given values for
$\mathbf u_{j},\dot{\mathbf u}_{j},\ddot{\mathbf u}_{j}$. Note that the initial data is given by
```{math}
:label: newmark_init
\mathbf u_0 = \mathbf u^0,\quad 
\dot{\mathbf u}_0 = \mathbf v^0,\quad 
\mathbf M\ddot{\mathbf u}_0 = \mathbf f_0-\mathbf C\mathbf v^0-\mathbf K\mathbf u^0.
```

The Newmark $\beta$-method is well-defined by the discrete approximations {eq}`newmark_approx`, the underlying system of ODE's for each time step {eq}`newmark_state` and the initial conditions {eq}`newmark_init`. In this general form a $3N\times 3N$ system would have to be solved in each time-step. In practice this is never done explicitely since there exist various ways to equivalently rewrite it. In fact, the following Theorem gives one general way to compute the single steps.
````{prf:Theorem}
Let $\mathbf M$ be positive definite and $\mathbf C,\mathbf K$ be positive semi-definite. Then the Newmark method is well defined for all $\tau>0$, $\alpha,2\beta\in[0,1]$ and initial data $\mathbf u^0,\mathbf v^0\in\mathbb R^N$. Each step requires one application of the inverse of the matrix $\mathbf S:=\mathbf M+\tau\gamma\mathbf C+\tau^2\beta\mathbf K$ and vector operations (additions, multiplications by scalars).
````
````{prf:Proof}
Inserting {eq}`newmark_approx` into {eq}`newmark_state` yields and collecting all terms containing $\ddot{\mathbf u}_{n+1}$ on the left hand side yields
```{math}
(\mathbf M+\gamma\tau\mathbf C+\beta\tau^2\mathbf K)\ddot{\mathbf u}_{j+1}=\mathbf f_{j+1}-\mathbf C(\dot{\mathbf u}_j+(1-\gamma)\tau\ddot{\mathbf u}_j)-\mathbf K(\mathbf u_j+\tau\dot{\mathbf u}_j+(\frac{1}{2}-\beta)\tau^2\ddot{\mathbf u}_j). 
```

Now the method can be written in matrix form
```{math}
\mathcal R
\begin{pmatrix}
\mathbf u_{j+1}\\
\dot{\mathbf u}_{j+1}\\
\ddot{\mathbf u}_{j+1}
\end{pmatrix}
=
\mathcal L
\begin{pmatrix}
\mathbf u_{j}\\
\dot{\mathbf u}_{j}\\
\ddot{\mathbf u}_{j}
\end{pmatrix}
+\begin{pmatrix}
0\\0\\\mathbf f_{j+1}
\end{pmatrix}.
```
with 
```{math}
\mathcal R:=
\begin{pmatrix}
\mathbf I&0&-\tau^2\beta \mathbf I\\
0&\mathbf I&-\tau\gamma \mathbf I\\
0&0&\mathbf M+\tau\gamma\mathbf C+\tau^2\beta\mathbf K
\end{pmatrix},\quad
\mathcal L:=
\begin{pmatrix}
\mathbf I&\tau&\tau^2(\frac{1}{2}-\beta) \mathbf I\\
0&\mathbf I&\tau(1-\gamma) \mathbf I\\
-\mathbf K&-\mathbf C-\tau\mathbf K&-\tau(1-\gamma)\mathbf C-\tau^2(\frac{1}{2}-\beta)\mathbf K
\end{pmatrix}
```
The system above can be solved by backwards elimination where merely in the process of computing $\ddot{\mathbf u}_{j+1}$ the inverse of the matrix $\mathbf S:=\mathbf M+\tau\gamma\mathbf C+\tau^2\beta\mathbf K$ has to be applied, all other operations are vector operations.
Clearly the matrix $\mathbf S$ is positive definite for all parameters $\tau,\gamma,\beta$ and thus regular.
````

### Stability of the Newmark method

````{prf:Theorem}
:label: newmark_energy
Let $\mathbf u_j, \dot{\mathbf u}_j, \ddot{\mathbf u}_j$ for $j=0,\ldots, M$ be the discrete Newmark approximations for given $\tau,\alpha,\beta$ and $\mathbf C = 0,\mathbf f=0$. Then 
```{math}
E^{\tau,\alpha,\beta}_{j+1}-E^{\tau,\alpha,\beta}_j=-(\gamma-\frac{1}{2})(\mathbf u_{j+1}-\mathbf u_j)^\top\mathbf K_{\tau,\alpha,\beta}(\mathbf u_{j+1}-\mathbf u_j)
```
with
```{math}
E^{\tau,\alpha,\beta}_{j}=\frac{1}{2}\left(\dot{\mathbf u}_j^\top\mathbf M\dot{\mathbf u}_j+{\mathbf u}_j^\top\mathbf K_{\tau,\alpha,\beta}{\mathbf u}_j\right),\quad
\mathbf K_{\tau,\alpha,\beta}=\mathbf K+\left(\beta-\frac{1}{2}\gamma\right)\tau^2\mathbf K\mathbf M^{-1}\mathbf K.
```
````

% ````{prf:Proof}
% We start with the change of the energy $\mathbf E(\mathbf u_j)$.
% Due to symmetry of $\mathbf M$ and $\mathbf K$ we have
% ```{math}
% \mathbf E(\mathbf u_{j+1})
% -\mathbf E(\mathbf u_{j})=\frac{1}{2}(\dot{\mathbf u}_{j+1}+\dot{\mathbf u}_{j})^\top\mathbf M(\dot{\mathbf u}_{j+1}-\dot{\mathbf u}_{j})
% +\frac{1}{2}({\mathbf u}_{j+1}+{\mathbf u}_{j})^\top\mathbf K(\dot{\mathbf u}_{j+1}-{\mathbf u}_{j})
% ```
% Using the fact that 
% ````

The quantity $E^{\tau,\alpha,\beta}_j$ is always positive for non-zero $\mathbf u_j,\dot{\mathbf u}_j$ if the matrix $\mathbf K_{\tau,\alpha,\beta}$ is positive semi-definite.
In this case we call $E^{\tau,\alpha,\beta}_j$ the discrete energy at timestep $j$. The matrix $\mathbf K_{\tau,\alpha,\beta}$ is certainly positive semi-definite if $\beta\geq \gamma/2$, for any $\tau>0$.

{prf:ref}`newmark_energy` shows that the discrete energy is preserved if $\gamma=1/2$. Moreover, the discrete energy $E^{\tau,\alpha,\beta}_j$ coincides with the Energy $\mathbf E(\mathbf u_j)$ iff $\beta=\gamma/2$.


We discuss some specific choices of $\gamma,\beta$ in more detail:

### The explicit method $\beta = 0$, $\gamma\geq 1/2$

for $\mathbf C=0$ the matrix $\mathbf S$ which has to be inverted reduces to $\mathbf M$. Thus the method is called **explicit**.
Moreover, let $\phi_n,\lambda_n$ be the normalized eigenpairs of 
```{math}
:label: disc_evp
\mathbf K u=\lambda \mathbf M
```
Then we have 
```{math}
\phi_n^\top\mathbf M\phi_k=\delta_{n,k},\quad
\phi_n^\top\mathbf K\phi_k=\lambda_n\delta_{n,k}.
```
Expanding $\mathbf u_j$ with respect to the basis $\phi_n$ leads to
```{math}
\mathbf u_j=\sum_{n=1}^N c_n\phi_n,
```
and due to the orthogonality of the eigenvectors we have
```{math}
\begin{aligned}
\mathbf u_j^\top\mathbf K_{\tau,\alpha,0}\mathbf u_j&= \mathbf u_j^\top\mathbf K\mathbf u_j-\frac{\gamma\tau^2}{2}\mathbf u_j^\top\mathbf K\mathbf M^{-1}\mathbf K\mathbf u_j\\
&=\sum_{n=1}^N\lambda_n\phi_n-\frac{\gamma\tau^2}{2}\lambda_n^2\phi_n.
\end{aligned}
```
Thus
```{math}
\mathbf u_j^\top\mathbf K_{\tau,\alpha,0}\mathbf u_j\geq 0
```
for all $\mathbf u_j$ if 
```{math}
:label: newmark_cfl
\tau^2\leq \frac{2}{\gamma\lambda_{\mathrm {max}}},
```
where $\lambda_{\mathrm {max}}$ is the maximal eigenvalue of {eq}`disc_evp`.
The condition {eq}`newmark_cfl` is called CFL condition and yields conditional stability of the explicit scheme, depending on the chosen time-step $\tau$.
