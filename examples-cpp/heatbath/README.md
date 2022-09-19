# Heatbath for Scalar Field Theory

$$
\def \rr {\rangle}
\def \ll {\langle}
\def \nn {\nonumber}
\def \ba {\begin{eqnarray}}
\def \ea {\end{eqnarray}}
\nn
$$

## Action

$$
\ba
S = \sum_x \Big(
-\sum_\mu \phi(x+\mu)\phi(x)
+ \big(4+ \frac{1}{2}m^2\big) \phi^2(x)
+ \frac{1}{4!}\lambda \phi^4(x)
\Big)
\ea
$$

## Correlation functions

$$
\ba
\phi^2 &=& \ll \phi^2(x) \rr
\\
C_2(t) &=& \ll \phi(t) \phi(0) \rr
\\
C_4(t) &=& \ll \phi^2(t) \phi^2(0) \rr
\ea
$$

where
$$
\ba
\phi(t) = \sum_{\vec x} \sum_{t'=t}^{t+\delta t-1} \phi(\vec x, t')
\ea
$$

## Observables

$$
R_4(t) = \frac{C_4(t) - C_2^2(0)}{C_2^2(t)}
$$

$$
\ba
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
\ea
$$

## Heatbath

Action related to one site $x$ is:
$$
- \sum_\mu (\phi(x+\mu) + \phi(x-\mu))\phi(x) + \big(4+ \frac{1}{2}m^2\big) \phi^2(x)
+ \frac{1}{4!}\lambda \phi^4(x)
$$

$$
\ba
k_1 &=& 4+ \frac{1}{2}m^2
\\
k_2 &=& \frac{1}{4!}\lambda
\ea
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
\ba
H(\pi, \phi) = T(\pi) + S(\phi) 
\ea
$$

The kinetic term is
$$
\ba
T(\pi) &=& \sum_p \frac{\pi(p)\pi(-p)}{2 \Big(4 \sin (\frac{p}{2})^2 + M^2\Big)}
\\
\pi(p) &=& \frac{1}{V} \sum_x \pi(x) e^{-i p \cdot x}
\ea
$$
Force is
$$
\ba
F(x)
&=& -\frac{\delta S(\phi)}{\delta \phi(x)}
\nn\\
&=& \sum_\mu (\phi(x+\mu) + \phi(x-\mu)) - (8+m^2)\phi(x) - \frac{1}{6} \lambda \phi^3(x)
\nn
\ea
$$

## Running of coupling

$$
\ba
\frac{d}{d \log(1/a)} \frac{1}{\lambda(a)} = -\frac{3}{16\pi^2}
\ea
$$

After reduce $a$ by a factor of $2$, the change on $1/\lambda$ should be
$$
\ba
\Delta\frac{1}{\lambda} = - \frac{3}{16\pi^2} \log(2) = 0.0132
\ea
$$


 