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
F_\pi^2(q^2)
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
where
$$
F_\pi(q^2) = 1 - {r_\pi^2 \over 6}q^2
$$


With lattice regularization
$$
q_\mu \to 2 \sin({q_\mu \over 2})
$$

$$
M_\pi^2 \to 2\big(\cosh(M_\pi) - 1\big)
$$

$$
p_t \to i M_\pi
$$

