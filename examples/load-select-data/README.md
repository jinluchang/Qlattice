# README

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

[TOC]

## Notation

### Gamma matrix

$$
\begin{eqnarray}
  \gamma^{\text{va}}_{{\mu}} & = & \left\{ \begin{array}{ll}
    \gamma_{{\mu}} & 0 \leqslant {\mu}< 4\\
    \gamma_{{\mu}} \gamma_5 & 4 \leqslant {\mu}< 8
  \end{array} \right.
\end{eqnarray}
$$

$$
\begin{eqnarray}
  \Gamma_{a + 2 b + 4 c + 8 d} & = & \gamma_x^a \gamma_y^b \gamma_z^c
  \gamma_t^d
\end{eqnarray}
$$

$$
\begin{eqnarray}
  \gamma_5 & = & \Gamma_{15} = \gamma_x \gamma_y \gamma_z \gamma_t
\end{eqnarray}
$$

### Propagator

Wall source propagator, with Coulomb gauge fixing:
$$
\begin{eqnarray}
  S (\vec{x}, t_{\text{snk}} ; t_{\text{src}}) & = & \sum_{\vec{y}} S
  (\vec{x}, t_{\text{snk}} ; \vec{y}, t_{\text{src}})
\end{eqnarray}
$$

$$
\begin{eqnarray}
  S (t_{\text{snk}} ; \vec{y}, t_{\text{src}}) & = & \sum_{\vec{x}} S
  (\vec{x}, t_{\text{snk}} ; \vec{y}, t_{\text{src}})
\end{eqnarray}
$$

$$
\begin{eqnarray}
  S (t_{\text{snk}} ; t_{\text{src}}) & = & \sum_{\vec{x}, \vec{y}} S
  (\vec{x}, t_{\text{snk}} ; \vec{y}, t_{\text{src}})
\end{eqnarray}
$$

## ```compute-two-point-func.h```

```cpp
ssprintf("analysis/lat-two-point/%s/results=%d", job_tag.c_str(), traj)
```

$$
\text{ld}
[0\le t_\text{sep} < T]
[0\le\text{op}_\text{src}<16]
[0\le\text{op}_\text{snk}<16]
$$

```cpp
ssprintf("/two-point-%d-%d.lat", type1, type2)
```

$$
\text{ld}[t_\text{sep}][\text{op}_\text{src}][\text{op}_\text{snk}]
= \text{Tr}\Big( \big(\sum_\vec x S_1(\vec x,t_\text{snk};t_\text{src}) \Gamma_{\text{op}_\text{src}} \gamma_5
S_2(\vec x,t_\text{snk};t_\text{src})^\dagger \gamma_5\big) \Gamma_{\text{op}_\text{snk}} \Big)
$$

```cpp
ssprintf("/two-point-wall-snk-sparse-corrected-%d-%d.lat", type1, type2)
```

$$
\ba
\text{ld}[t_\text{sep}][\text{op}_\text{src}][\text{op}_\text{snk}]
&=& \text{Tr}\Big( \big(\sum_\vec x S_1(\vec x,t_\text{snk};t_\text{src}) \Gamma_{\text{op}_\text{src}}
\sum_\vec y S_2(t_\text{src};\vec y,t_\text{snk})\big) \Gamma_{\text{op}_\text{snk}} \Big)
\\
&=& \text{Tr}\Big( \big(\sum_\vec x S_1(\vec x,t_\text{snk};t_\text{src}) \Gamma_{\text{op}_\text{src}} \gamma_5
\sum_\vec y S_2(\vec y,t_\text{snk};t_\text{src})^\dagger \gamma_5\big) \Gamma_{\text{op}_\text{snk}} \Big)
\ea
$$

## ```compute-three-point-func.h```

```cpp
ssprintf("analysis/lat-three-point/%s/results=%d", job_tag.c_str(), traj)
```

$$
\text{ld}
[0\le t_\text{sep} < T]
[0\le t_\text{op} < T]
[0\le \text{op} < 16]
$$

```cpp
ssprintf("/three-point-%d-%d-%d.lat", type1, type2, type3)
```

$$
\ba
\text{ld}[t_\text{sep}][t_\text{op}][\text{op}]
= \text{Tr}
\sum_\vec x
\Big(
S_1(t_\text{op},\vec x;t_\text{src})  \gamma_5 
S_3(t_\text{src};t_\text{snk})
\gamma_5
\big( \gamma_5 S_2(t_\text{snk};t_\text{op},\vec x)^\dagger \gamma_5 \big)
\Big)
\Gamma_{\text{op}} 
\ea
$$

