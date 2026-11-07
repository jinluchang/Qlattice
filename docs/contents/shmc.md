:orphan:

# SHMC formulation

> Luchang Jin
>
> 2024/06/13

## Formula

### Gauge link convention

$$
\psi(x) \sim U(x,\mu) \psi(x + \hat\mu)
\tag{1}
$$

Let

$$
U(x,\mu) = U_1(x,\mu) U_2^\dagger(x,\mu)
\tag{2}
$$

We can also treat $U_1$ and $U_2$ as the integration variable in the path integral.

### Gauge momentum convention

$$
\pi_1(x,\mu) = \pi_1(x,\mu)^a T^a
\tag{3}
$$

$$
\pi_2(x,\mu) = \pi_2(x,\mu)^a T^a
\tag{4}
$$

$$
\pi_1^a \sim N\Big(0,\frac{1}{\sqrt{2}}\Big)
\tag{5}
$$

$$
\pi_1^a \sim e^{-(\pi_1^a)^2}
\tag{6}
$$

$$
\pi_2^a \sim N\Big(0,\frac{1}{\sqrt{2}}\Big)
\tag{7}
$$

$$
\pi_2^a \sim e^{-(\pi_2^a)^2}
\tag{8}
$$

where

$$
{T^a}^\dagger = - T^a
\tag{9}
$$

$$
\mathrm{tr}(T^a T^b) = -2 \delta^{a,b}
\tag{10}
$$

$$
\mathrm{tr}(T^a) = 0
\tag{11}
$$

$$
T^a T^a = -\frac{16}{3}
\tag{12}
$$

### Hamiltonian

$$
S(U) = - \beta \Big( (1 - 8 c_1) \sum_P \Big(\frac{1}{3} \mathrm{tr}(U_P) - 1\Big) + c_1 \sum_R \Big(\frac{1}{3} \mathrm{tr}(U_R) - 1\Big) \Big)
\tag{13}
$$

$$
H(\pi^a, U) = \sum_{x,\mu} \big(\pi_1(x,\mu)^a\big)^2 + \sum_{x,\mu} \big(\pi_2(x,\mu)^a\big)^2 + S(U)
\tag{14}
$$

Total number of plaq $P$ is

$$
6 L^3 T
\tag{15}
$$

Total number of rectangular plaq $R$ is

$$
12 L^3 T
\tag{16}
$$

Wilson:

$$
c_1 = 0
\tag{17}
$$

Iwasaki:

$$
c_1 = -0.331
\tag{18}
$$

DBW2:

$$
c_1 = -1.4008
\tag{19}
$$

### Gauge evolve

$$
U(x,\mu) \gets e^{\pi_1(x,\mu) dt} U(x,\mu) e^{-\pi_2(x,\mu) dt}
\tag{20}
$$

$$
e^{\pi_1(x,\mu) dt} U(x,\mu) e^{-\pi_2(x,\mu) dt} = e^{\pi_1(x,\mu)^a T^a dt} U(x,\mu) e^{-\pi_2(x,\mu)^a T^a dt}
\tag{21}
$$

Conceptually:

$$
U_1(x,\mu) \gets e^{\pi_1(x,\mu) dt} U(x,\mu)
\tag{22}
$$

$$
e^{\pi_1(x,\mu) dt} U(x,\mu) = e^{\pi_1(x,\mu)^a T^a dt} U_1(x,\mu)
\tag{23}
$$

$$
U_2(x,\mu) \gets e^{\pi_2(x,\mu) dt} U(x,\mu)
\tag{24}
$$

$$
e^{\pi_2(x,\mu) dt} U(x,\mu) = e^{\pi_2(x,\mu)^a T^a dt} U_2(x,\mu)
\tag{25}
$$

Convention for gauge field change:

$$
U(x,\mu) \gets e^{d s_1 (x,\mu)^a T^a} U(x,\mu) e^{-d s_2 (x,\mu)^a T^a}
\tag{26}
$$

$$
d U(x,\mu) \approx d s_1(x,\mu)^a T^a U(x,\mu) - d s_2(x,\mu)^a U(x,\mu) T^a = \Big( d s_1(x,\mu)^a T^a - d s_2(x,\mu)^a U(x,\mu) T^a U^{-1}(x,\mu) \Big) U(x,\mu) = \Big( d s_1(x,\mu)^a - d s_2(x,\mu)^b U^{ba}(x,\mu) \Big) T^a U(x,\mu)
\tag{27}
$$

where we define

$$
U^{ba}(x,\mu) T^a = U(x,\mu) T^b U^{-1}(x,\mu)
\tag{28}
$$

therefore

$$
U^{ba}(x,\mu) T^b = U^{-1}(x,\mu) T^a U(x,\mu)
\tag{29}
$$

### Gauge momentum evolve

$$
\frac{d}{dt}s_1(x,\mu)^a = \pi_1(x,\mu)^a
\tag{30}
$$

$$
\frac{d}{dt}s_2(x,\mu)^a = \pi_2(x,\mu)^a
\tag{31}
$$

$$
\frac{d}{dt}\pi_1(x,\mu)^a = - \frac{1}{2} \frac{\delta H}{\delta s_1(x,\mu)^a}
\tag{32}
$$

$$
\frac{d}{dt}\pi_2(x,\mu)^a = - \frac{1}{2} \frac{\delta H}{\delta s_2(x,\mu)^a}
\tag{33}
$$

$$
\frac{d}{dt}\pi_1(x,\mu) = -\frac{1}{2} T^a \frac{\delta S}{\delta s_1(x,\mu)^a} = \frac{\beta}{6} T^a \mathrm{Re}\, \mathrm{tr}\big( T^a U(x,\mu) C^\dagger(x,\mu) \big) = -\frac{\beta}{3} \mathcal P \big\{ U(x,\mu) C^\dagger(x,\mu) \big\}
\tag{34}
$$

Note:

$$
\frac{\delta S}{\delta s_2(x,\mu)^a} = -U^{ba}(x,\mu) \frac{\delta S}{\delta s_1(x,\mu)^b}
\tag{35}
$$

Therefore

$$
\frac{d}{dt}\pi_2(x,\mu) = -\frac{1}{2} T^b \frac{\delta S}{\delta s_2(x,\mu)^b} = \frac{1}{2} U^{ba}(x,\mu) T^b \frac{\delta S}{\delta s_1(x,\mu)^a} = \frac{1}{2} U(x,\mu) T^a U(x,\mu) \frac{\delta S}{\delta s_1(x,\mu)^a} = - U^{-1}(x,\mu) \frac{d}{dt}\pi_1(x,\mu) U(x,\mu)
\tag{36}
$$

where

$$
C(x,\mu) = (1 - 8 c_1) \sum_{\nu\,(\nu\neq\mu)} U(x,\nu)U(x + \hat\nu,\mu)U(x + \hat\nu + \hat\mu,-\nu) + c_1 \sum_{\nu\,(\nu\neq\mu)} U(x,\nu)U(x+\hat\nu,\nu)U(x + 2\hat\nu,\mu)U(x + 2\hat\nu + \hat\mu,-\nu)U(x + \hat\nu + \hat\mu,-\nu) + c_1 \sum_{\nu\,(\nu\neq\mu)} U(x,-\mu)U(x-\hat\mu,\nu)U(x-\hat\mu + \hat\nu,\mu)U(x + \hat\nu,\mu)U(x + \hat\nu + \hat\mu,-\nu) + c_1 \sum_{\nu\,(\nu\neq\mu)} U(x,\nu)U(x + \hat\nu,\mu)U(x + \hat\mu + \hat\nu,\mu)U(x + \hat\nu + 2\hat\mu,-\nu)U(x + 2\hat\mu,-\mu)
\tag{37}
$$

$$
\mathcal P \{m\} = \frac{1}{2}(m - m^\dagger) - \frac{1}{6}\mathrm{tr}(m - m^\dagger)
\tag{38}
$$

Note

$$
\mathcal P\{T^a\} = T^a
\tag{39}
$$

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


