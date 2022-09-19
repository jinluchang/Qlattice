# Note

> Luchang Jin

$$
\def \ra {\rangle}
\def \la {\langle}
\def \nn {\nonumber}
\def \ba {\begin{eqnarray}}
\def \ea {\end{eqnarray}}
\def \a {a}
\nn
$$

## Propagator

$$
G_\pi(x-y) = \int\frac{d^4 q}{(2\pi)^4} \frac{1}{q^2 + M_\pi^2} e^{i q\cdot (x-y)}
$$

Let $p= (iM_\pi , \vec 0)$, then:
$$
\ba
\int d^4 x~e^{-i q \cdot x}
{1 \over 2 M_\pi}
\la \pi(\vec 0) | J_\mu(x) J_\nu(0) | \pi(\vec 0) \ra
&=&
e^2
{
F_\pi^2(q^2)
\over
2 M_\pi
}
\Big[
2 \delta_{\mu,\nu}
-
{
(2p + q)_\mu(2p+q)_\nu
\over
(p + q)^2 + M_\pi^2
}
\nn\\
&&
\hspace{3cm}
-
{
(2p - q)_\mu(2p - q)_\nu
\over
(p - q)^2 + M_\pi^2
}

\Big]
\ea
$$
for zero spatial momentum, we set
$$
\ba
\int d^3 x
{1 \over 2 M_\pi}
\la \pi(\vec 0) | J_t(x) J_t(0) | \pi(\vec 0) \ra
&=&
e^2
\ea
$$
where
$$
F_\pi(q^2) = 1 - {r_\pi^2 \over 6}q^2 + \cdots
$$


With lattice regularization
$$
q_\mu \to 2 \sin \Big({q_\mu \over 2}\Big)
$$

$$
M_\pi^2 \to 2\big(\cosh(M_\pi) - 1\big) = 4 \sinh^2\Big({M_\pi \over 2}\Big)
$$

$$
p_t \to 2 i \sinh \Big({M_\pi \over 2}\Big).
$$

In for $q$ and $\mu$ that $q_\mu = \pi$, we will use the following special replacement $q_\mu^2 \to 4$ and $q_\mu \to 0$.

In the code, we set $e=1$.

