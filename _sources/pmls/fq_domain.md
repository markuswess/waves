(fq_domain_pmls)= 
# Perfectly Matched Layers in the Frequency Domain

We start by introducing PMLs in the simple one-dimensional setting.

## PMLs in one dimension

Consider the one dimensional Helmholtz problem for given $k>0$ and compactly in $(0,a)$, $a>0$ supported $f$ to find $u: (0,\infty)\to\mathbb R$ such that 
```{math}
\begin{aligned}
-u''-k^2 u &= f,\\
u(0) &= 0,\\
\lim_{x\to\infty} u'(x)-iku(x) &= 0.
\end{aligned}
```
The last condition is the Sommerfeld radiation condition which accounts for the missing boundary condition at the right.



