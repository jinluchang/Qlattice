# Heatbath for Scalar Field Theory

## Action

$$
\begin{eqnarray}
S = \sum_x \Big(
-\sum_\mu \phi(x+\mu)\phi(x)
+ \big(4+ \frac{1}{2}m^2\big) \phi^2(x)
+ \frac{1}{4!}\lambda \phi^4(x)
\Big)
\end{eqnarray}
$$

## Correlation functions

$$
\begin{eqnarray}
\phi^2 &=& \langle \phi^2(x) \rangle
\\
C_2(t) &=& \langle \phi(t) \phi(0) \rangle
\\
C_4(t) &=& \langle \phi^2(t) \phi^2(0) \rangle
\end{eqnarray}
$$

where
$$
\begin{eqnarray}
\phi(t) = \sum_{\vec x} \sum_{t'=t}^{t+\delta t-1} \phi(\vec x, t')
\end{eqnarray}
$$

## Observables

$$
R_4(t) = \frac{C_4(t) - C_2^2(0)}{C_2^2(t)}
$$

$$
\begin{eqnarray}
m_\text{eff}(t_1,t_2)
&=&
\frac{1}{t_2-t_1} \log\Bigg(\frac{C_2(t_1)}{C_2(t_2)} \Bigg)
\\
V_\text{eff}(t_1,t_2)
&=&
\frac{1}{t_2-t_1}
\log\Bigg(
\frac{R_4(t_1)}{R_4(t_2)}
\Bigg)
\end{eqnarray}
$$

## Heatbath

Action related to one site $x$ is:
$$
- \sum_\mu (\phi(x+\mu) + \phi(x-\mu))\phi(x) + \big(4+ \frac{1}{2}m^2\big) \phi^2(x)
+ \frac{1}{4!}\lambda \phi^4(x)
$$

$$
\begin{eqnarray}
k_1 &=& 4+ \frac{1}{2}m^2
\\
k_2 &=& \frac{1}{4!}\lambda
\end{eqnarray}
$$

### Sample results

#### v1

```cpp
const Coordinate total_site = Coordinate(4,4,4,256);
const double mass_sqr = 0.04;
const double lambda = 0.0;
const int t1 = 2;
const int t2 = 4;
const int dt = 1;
```

Results:

```cpp
n_traj=15175 ; m_eff=0.199134427039918 ; v_eff=-0.000256958791712.
```

#### v2

```cpp
const Coordinate total_site = Coordinate(4,4,4,256);
const double mass_sqr = 0.00;
const double lambda = 0.4;
const int t1 = 2;
const int t2 = 4;
const int dt = 1;
```

Results:

```cpp
n_traj=38629 ; m_eff=0.188785922011518 ; v_eff=0.018863497373480.
```

## HMC

$$
\begin{eqnarray}
H(\pi, \phi) = T(\pi) + S(\phi) 
\end{eqnarray}
$$

The kinetic term is
$$
\begin{eqnarray}
T(\pi) &=& \sum_p \frac{\pi(p)\pi(-p)}{2 \Big(4 \sin (\frac{p}{2})^2 + M^2\Big)}
\\
\pi(p) &=& \frac{1}{V} \sum_x \pi(x) e^{-i p \cdot x}
\end{eqnarray}
$$
Force is
$$
\begin{eqnarray}
F(x)
&=& -\frac{\delta S(\phi)}{\delta \phi(x)}
\nonumber\\
&=& \sum_\mu (\phi(x+\mu) + \phi(x-\mu)) - (8+m^2)\phi(x) - \frac{1}{6} \lambda \phi^3(x)
\nonumber
\end{eqnarray}
$$

## Running of coupling

$$
\begin{eqnarray}
\frac{d}{d \log(1/a)} \frac{1}{\lambda(a)} = -\frac{3}{16\pi^2}
\end{eqnarray}
$$

After reduce $a$ by a factor of $2$, the change on $1/\lambda$ should be
$$
\begin{eqnarray}
\Delta\frac{1}{\lambda} = - \frac{3}{16\pi^2} \log(2) = 0.0132
\end{eqnarray}
$$


 