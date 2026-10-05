:orphan:

# SHMC formulation

> Luchang Jin
>
> 2024/06/13

## Formula

### Gauge link convention

$$
\psi(x) \sim U(x,\mu) \psi(x + \hat\mu)
$$

Let

$$
\begin{align}
U(x,\mu) = U_1(x,\mu) U_2^\dagger(x,\mu)
\end{align}
$$

We can also treat $U_1$ and $U_2$ as the integration variable in the path integral.

### Gauge momentum convention

$$
\begin{align}
\pi_1(x,\mu) =& \pi_1(x,\mu)^a T^a
\\
\pi_2(x,\mu) =& \pi_2(x,\mu)^a T^a
\end{align}
$$

$$
\begin{align}
\pi_1^a \sim N\Big(0,\frac{1}{\sqrt{2}}\Big) \sim e^{-(\pi_1^a)^2}
\\
\pi_2^a \sim N\Big(0,\frac{1}{\sqrt{2}}\Big) \sim e^{-(\pi_2^a)^2}
\end{align}
$$

where

$$
\begin{align}
{T^a}^\dagger &= - T^a
\\
\mathrm{tr}(T^a T^b) &= -2 \delta^{a,b}
\\
\mathrm{tr}(T^a) &= 0
\\
T^a T^a &= -\frac{16}{3}
\end{align}
$$

### Hamiltonian

$$
\begin{align}
S(U)
=
- \beta
\Big(
(1 - 8 c_1) \sum_P \Big(\frac{1}{3} \mathrm{tr}(U_P) - 1\Big)
+ c_1 \sum_R \Big(\frac{1}{3} \mathrm{tr}(U_R) - 1\Big)
\Big)
\end{align}
$$

$$
\begin{align}
H(\pi^a, U) =
\sum_{x,\mu} \big(\pi_1(x,\mu)^a\big)^2
+
\sum_{x,\mu} \big(\pi_2(x,\mu)^a\big)^2
+ S(U)
\end{align}
$$

Total number of plaq $P$ is $6 L^3 T$, total number of rectangular plaq $R$ is $12 L^3 T$​.

Wilson: $c_1 = 0$

Iwasaki: $c_1 = -0.331$

DBW2: $c_1 = -1.4008$

### Gauge evolve

$$
\begin{align}
U(x,\mu)
\gets
e^{\pi_1(x,\mu) dt} U(x,\mu) e^{-\pi_2(x,\mu) dt}
=
e^{\pi_1(x,\mu)^a T^a dt} U(x,\mu) e^{-\pi_2(x,\mu)^a T^a dt}
\end{align}
$$

Conceptually:

$$
\begin{align}
U_1(x,\mu) \gets e^{\pi_1(x,\mu) dt} U(x,\mu) =& e^{\pi_1(x,\mu)^a T^a dt} U_1(x,\mu)
\\
U_2(x,\mu) \gets e^{\pi_2(x,\mu) dt} U(x,\mu) =& e^{\pi_2(x,\mu)^a T^a dt} U_2(x,\mu)
\end{align}
$$

Convention for gauge field change:

$$
\begin{align}
U(x,\mu)
\gets
e^{d s_1 (x,\mu)^a T^a}
U(x,\mu)
e^{-d s_2 (x,\mu)^a T^a}
\end{align}
$$

$$
\begin{align}
d U(x,\mu)
\approx &
d s_1(x,\mu)^a T^a
U(x,\mu)
-
d s_2(x,\mu)^a
U(x,\mu) T^a
\\
=&
\Big(
d s_1(x,\mu)^a
T^a
-
d s_2(x,\mu)^a
U(x,\mu) T^a U^{-1}(x,\mu)
\Big)
U(x,\mu)
\\
=&
\Big(
d s_1(x,\mu)^a
-
d s_2(x,\mu)^b
U^{ba}(x,\mu)
\Big)
T^a
U(x,\mu)
\end{align}
$$

where we define

$$
\begin{align}
U^{ba}(x,\mu) T^a = U(x,\mu) T^b U^{-1}(x,\mu)
\end{align}
$$

therefore

$$
\begin{align}
U^{ba}(x,\mu) T^b = U^{-1}(x,\mu) T^a U(x,\mu)
\end{align}
$$

### Gauge momentum evolve

$$
\begin{align}
\frac{d}{dt}s_1(x,\mu)^a &= \pi_1(x,\mu)^a
\\
\frac{d}{dt}s_2(x,\mu)^a &= \pi_2(x,\mu)^a
\\
\frac{d}{dt}\pi_1(x,\mu)^a &= - \frac{1}{2} \frac{\delta H}{\delta s_1(x,\mu)^a}
\\
\frac{d}{dt}\pi_2(x,\mu)^a &= - \frac{1}{2} \frac{\delta H}{\delta s_2(x,\mu)^a}
\end{align}
$$

$$
\begin{align}
\frac{d}{dt}\pi_1(x,\mu)
&=
-\frac{1}{2} T^a \frac{\delta S}{\delta s_1(x,\mu)^a}
\\
&=
\frac{\beta}{6} T^a \mathrm{Re}\, \mathrm{tr}\big( T^a U(x,\mu) C^\dagger(x,\mu) \big)
\\
&=
-\frac{\beta}{3} \mathcal P \big\{ U(x,\mu) C^\dagger(x,\mu) \big\}
\end{align}
$$

Note:

$$
\begin{align}
\frac{\delta S}{\delta s_2(x,\mu)^a}
=
-U^{ba}(x,\mu)
\frac{\delta S}{\delta s_1(x,\mu)^b}
\end{align}
$$

Therefore

$$
\begin{align}
\frac{d}{dt}\pi_2(x,\mu)
&=
-\frac{1}{2} T^b \frac{\delta S}{\delta s_2(x,\mu)^b}
\\
&=
\frac{1}{2} U^{ba}(x,\mu) T^b \frac{\delta S}{\delta s_1(x,\mu)^a}
\\
&=
\frac{1}{2} U(x,\mu) T^a U(x,\mu) \frac{\delta S}{\delta s_1(x,\mu)^a}
\\
&=
- U^{-1}(x,\mu) \frac{d}{dt}\pi_1(x,\mu) U(x,\mu)
\end{align}
$$

where

$$
\begin{align}
C(x,\mu)
&=
(1 - 8 c_1)
\sum_{\nu\,(\nu\neq\mu)}
U(x,\nu)U(x + \hat\nu,\mu)U(x + \hat\nu + \hat\mu,-\nu)
\nonumber\\&\quad
+
c_1
\sum_{\nu\,(\nu\neq\mu)}
U(x,\nu)U(x+\hat\nu,\nu)U(x + 2\hat\nu,\mu)U(x + 2\hat\nu + \hat\mu,-\nu)U(x + \hat\nu + \hat\mu,-\nu)
\nonumber\\&\quad
+
c_1
\sum_{\nu\,(\nu\neq\mu)}
U(x,-\mu)U(x-\hat\mu,\nu)U(x-\hat\mu + \hat\nu,\mu)U(x + \hat\nu,\mu)U(x + \hat\nu + \hat\mu,-\nu)
\nonumber\\&\quad
+
c_1
\sum_{\nu\,(\nu\neq\mu)}
U(x,\nu)U(x + \hat\nu,\mu)U(x + \hat\mu + \hat\nu,\mu)U(x + \hat\nu + 2\hat\mu,-\nu)U(x + 2\hat\mu,-\mu)
\end{align}
$$

$$
\begin{align}
\mathcal P \{m\}
=
\frac{1}{2}(m - m^\dagger) - \frac{1}{6}\mathrm{tr}(m - m^\dagger)
\end{align}
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
matrix_evolve(U_x_mu,P_x_mu,dt) ==> make_color_matrix_exp(P_x_mu * dt) * U_x_mu
U_x_mu = matrix_evolve(U_x_mu,P_x_mu,dt)
```

$H$

```
double gm_hamilton_node(const GaugeMomentum& gm)
double gf_hamilton_node(const GaugeField& gf, const GaugeAction& ga)
```

$\mathcal P\{m\}$

```
ColorMatrix make_tr_less_anti_herm_matrix(const ColorMatrix& m)
```

```
void gf_evolve(GaugeField& gf, const GaugeMomentum& gm, const double step_size)
void gf_evolve_dual(GaugeField& gf, const GaugeMomentum& gm_dual, const double step_size)
void set_gm_force(GaugeMomentum& gm_force, const GaugeField& gf, const GaugeAction& ga)
void set_gm_force_dual(GaugeMomentum& gm_force_dual, const GaugeField& gf, const GaugeMomentum& gm_force)
```



