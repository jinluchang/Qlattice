$$
\def\ba#1\ea{\begin{align}#1\end{align}}
\newcommand{\nn}{\nonumber}
\newcommand{\ra}{\rangle}
\newcommand{\la}{\langle}
\newcommand{\bra}{\big\rangle}
\newcommand{\bla}{\big\langle}
\newcommand{\Bra}{\Big\rangle}
\newcommand{\Bla}{\Big\langle}
\newcommand{\ud}{\mathrm{d}}
\newcommand{\a}{a}
\nn
$$



# HMC formulation

> Luchang Jin
>
> 2024/06/12

## Formula

### Gauge link convention

$$
\psi(x) \sim U(x,\mu) \psi(x + \hat\mu)
$$

### Gauge momentum convention

$$
\pi(x,\mu) = \pi(x,\mu)^a T^a
$$

$$
\pi^a \sim N\Big(0,\frac{1}{\sqrt{2}}\Big) \sim e^{-(\pi^a)^2}
$$

where
$$
\ba
{T^a}^\dagger &= - T^a
\\
\mathrm{tr}(T^a T^b) &= -2 \delta^{a,b}
\\
\mathrm{tr}(T^a) &= 0
\\
T^a T^a &= -\frac{16}{3}
\ea
$$

### Hamiltonian

$$
\ba
S(U)
=
- \beta
\Big(
(1 - 8 c_1) \sum_P \Big(\frac{1}{3} \mathrm{tr}(U_P) - 1\Big)
+ c_1 \sum_R \Big(\frac{1}{3} \mathrm{tr}(U_R) - 1\Big)
\Big)
\ea
$$

$$
\ba
H(\pi^a, U) =
\sum_{x,\mu} \big(\pi(x,\mu)^a\big)^2
+ S(U)
\ea
$$

Total number of plaq $P$ is $6 L^3 T$, total number of rectangular plaq $R$ is $12 L^3 T$​.

Wilson: $c_1 = 0$

Iwasaki: $c_1 = -0.331$

DBW2: $c_1 = -1.4008$

### Gauge evolve

$$
U(x,\mu) \gets e^{\pi(x,\mu) dt} U(x,\mu) = e^{\pi(x,\mu)^a T^a dt} U(x,\mu)
$$

Convention for gauge field change:
$$
U(x,\mu) \gets e^{d s (x,\mu)^a T^a} U(x,\mu)
$$

$$
d U(x,\mu) \approx d s(x,\mu)^a T^a U(x,\mu)
$$

### Gauge momentum evolve

$$
\ba
\frac{d}{dt}s(x,\mu)^a &= \pi(x,\mu)^a
\\
\frac{d}{dt}\pi(x,\mu)^a &= - \frac{1}{2} \frac{\delta H}{\delta s(x,\mu)^a}
\ea
$$

$$
\ba
\frac{d}{dt}\pi(x,\mu)
&=
-\frac{1}{2} T^a \frac{\delta S}{\delta s(x,\mu)^a}
\\
&=
\frac{\beta}{6} T^a \mathrm{Re}\, \mathrm{tr}\big( T^a U(x,\mu) C^\dagger(x,\mu) \big)
\\
&=
-\frac{\beta}{3} \mathcal P \big\{ U(x,\mu) C^\dagger(x,\mu) \big\}
\ea
$$

where
$$
\ba
C(x,\mu)
&=
(1 - 8 c_1)
\sum_{\nu\,(\nu\neq\mu)}
U(x,\nu)U(x + \hat\nu,\mu)U(x + \hat\nu + \hat\mu,-\nu)
\nn\\&\quad
+
c_1
\sum_{\nu\,(\nu\neq\mu)}
U(x,\nu)U(x+\hat\nu,\nu)U(x + 2\hat\nu,\mu)U(x + 2\hat\nu + \hat\mu,-\nu)U(x + \hat\nu + \hat\mu,-\nu)
\nn\\&\quad
+
c_1
\sum_{\nu\,(\nu\neq\mu)}
U(x,-\mu)U(x-\hat\mu,\nu)U(x-\hat\mu + \hat\nu,\mu)U(x + \hat\nu,\mu)U(x + \hat\nu + \hat\mu,-\nu)
\nn\\&\quad
+
c_1
\sum_{\nu\,(\nu\neq\mu)}
U(x,\nu)U(x + \hat\nu,\mu)U(x + \hat\mu + \hat\nu,\mu)U(x + \hat\nu + 2\hat\mu,-\nu)U(x + 2\hat\mu,-\mu)
\ea
$$

$$
\ba
\mathcal P \{m\}
=
\frac{1}{2}(m - m^\dagger) - \frac{1}{6}\mathrm{tr}(m - m^\dagger)
\ea
$$

Note $\mathcal P\{T^a\} = T^a$.

## Code

```
P_a = basis_projection_anti_hermitian_matrix(P)
make_anti_hermitian_matrix(P_a) ==> \sum_a P_a * T_a
P = make_anti_hermitian_matrix(P_a)
```

```
make_g_rand_anti_hermitian_matrix(rng_state, sigma) # default sigma=1.0
P = make_g_rand_anti_hermitian_matrix(rng_state, sigma)
```

```
matrix_evolve(U_x_mu,P_x_mu,dt) ==> make_matrix_exp(P_x_mu * dt) * U_x_mu
U_x_mu = matrix_evolve(U_x_mu,P_x_mu,dt)
```

$H$

```
double gm_hamilton_node(const GaugeMomentum& gm)
double gf_hamilton_node_no_comm(const GaugeField& gf, const GaugeAction& ga)
```

$\mathcal P\{m\}$

```
ColorMatrix make_tr_less_anti_herm_matrix(const ColorMatrix& m)
```

