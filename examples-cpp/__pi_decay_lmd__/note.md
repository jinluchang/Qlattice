# Note

## Code conventions

$$
\langle J_\mu(x)  J_\nu(y) \pi(w) \rangle = - \int_z F_{\mu,\nu}(x,y,z) G_{m_\pi}(w-z)
$$

$$
F_{\mu,\nu}(x,y,z) = \frac{-i}{4\pi^2 F_\pi} \epsilon_{\mu,\nu,\rho,\sigma}
(-i\partial_\rho^x)(-i\partial_\sigma^y)
F(x,y,z)
$$

For VMD model:
$$
F(x,y,z) = M_V^4 G_{M_V}(x-z) G_{M_V}(y-z)
$$

$$
\begin{eqnarray}
&&\text{set_photon_pion_photon_vertex}(\Gamma_1(x)[\mu], \pi(z), \Gamma_2(y)[\nu], M_V, F_\pi)
\\
\implies
&&\Gamma_1(x)[\mu] \leftarrow -\int_{y,z} F_{\mu,\nu}(x,y,z) \pi(z) \Gamma_2(y)[\nu]
\nonumber
\end{eqnarray}
$$
For VMD
$$
\Gamma_1(x)[\mu]
\leftarrow
(-i\partial_\rho^x) \int_z G_{M_V}(x-z)
\frac{-i M_V^4}{4\pi^2 F_\pi} \epsilon_{\mu,\nu,\rho,\sigma}
\int_{y}  \pi(z)
G_{M_V}(y-z)(-i\partial_\sigma^y)\Gamma_2(y)[\nu]
$$

## Saved fields

### pion-gg-decay

$$
\begin{eqnarray}
F(x)[4\mu+\nu]
&=&  e^{m_\pi t_\text{sep}}
\int_\vec w \langle J_\mu(x) J_\nu(0) \pi(\vec w, -t_\text{sep}) \rangle
\\
&=& - \int_z F_{\mu,\nu} (x,0,z)
\int_\vec{w} e^{m_\pi t_\text{sep}} G_{m_\pi}(z-(\vec{w},-t_\text{sep}))
\end{eqnarray}
$$

### pion-corr

$$
C(t) = \int_{\vec{w},\vec{v}}G_{m_\pi}((\vec w, t) - (\vec v, 0))
$$

### parameters (physical pion mass)

$$
\begin{eqnarray}
m_\pi &=& 134.9766~\mathrm{MeV}
\\
F_\pi &=& 92 ~\mathrm{MeV}
\\
M_V &=& 770~\mathrm{MeV}
\end{eqnarray}
$$

|      size        | $t_\text{sep}$ | $a^{-1}/\mathrm{GeV}$ |
| :-------------: | :-------------------: | :------------: |
| $24^3\times 96$ |          24          |      1.0      |
| $32^3\times 128$ | 32 | 1.0 |
| $32^3\times 128$ | 32 | 1.3333 |
| $48^3\times 192$ | 48 | 1.0 |
| $48^3\times 192$ | 48 | 2.0 |

### parameters (heavy pion mass)

$$
\begin{eqnarray}
m_\pi &=& 340 ~\mathrm{MeV}
\\
F_\pi &=& 105 ~\mathrm{MeV}
\\
M_V &=& 830 ~\mathrm{MeV}
\end{eqnarray}
$$

|       size       | $t_\text{sep}$ | $a^{-1}/\mathrm{GeV}$ |
| :--------------: | :------------: | :-------------------: |
| $24^3\times 96$  |       24       |          1.0          |
| $32^3\times 128$ |       32       |          1.0          |
| $32^3\times 128$ |       32       |        1.3333         |
| $48^3\times 192$ |       48       |          1.0          |
| $48^3\times 192$ |       48       |          2.0          |

