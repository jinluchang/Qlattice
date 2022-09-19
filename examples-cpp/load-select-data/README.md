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

## Contraction functions

See ``Qlattice/docs/contraction.md`` and ``Qlattice/docs/file-format.md``.

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

$$
\theta_\mu =  \left\{ \begin{array}{ll}
    1 & 0 \leqslant {\mu}< 4\\
   -1 & 4 \leqslant {\mu}< 8
  \end{array} \right.
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

## ``compute-two-point-func.h``

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

## ``compute-three-point-func.h``

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

## ``compute-psel-fsel-distribution.h``

```cpp
ssprintf("analysis/field-psel-fsel-distribution/%s/results=%d", job_tag.c_str(), traj)
```

Data format: ``FieldM<Complex, 1>`` with ``write_field_double``.

```cpp
ssprintf("/pos.field")
```
The expectation value is:
$$
H(x-y) = 1
$$

The data is created by summing over all selected points for $x$ and all point source locations for $y$, and then properly normalize the data.

```cpp
ssprintf("/neg.field")
```

Same as the ``pos.field`` but with $x$ and $y$ reversed.

```cpp
ssprintf("/avg.field")
```

The average of the above two data sets.

## ``compute-meson-vv.h``

<img src="figs/matrix-elements/png/fig-6.png" width=400px />

```cpp
ssprintf("analysis/field-meson-vv/%s/results=%d", job_tag.c_str(), traj)
```

Data format: ``FieldM<Complex, 8 * 8>`` with ``write_field_float_from_double``.

Use $y$ as the point source location in calculation.

```cpp
ssprintf("/decay-%d-%d-%d.field", type1, type2, type3)
```

$$
\ba
H_\text{decay-1-2-3}(x-y)[8\mu+\nu]
&=&\mathrm{Tr}
[S_3(x;y)\gamma^{\mathrm{va}}_\nu S_2(y;t_\text{src})\gamma_5 S_1(t_\text{src};x)\gamma^{\mathrm{va}}_\mu]
\ea
$$

where:
$$
\ba
t_\text{src} &=& \min(x_t,y_t) - t_\text{sep}
\\
t_\text{snk} &=& \max(x_t,y_t) + t_\text{sep}
\ea
$$
and for $t_\text{sep}$:

```cpp
inline int tsep_op_wall_src(const std::string& job_tag)
// parameter
{
  if (job_tag == "24D" or job_tag == "32D" or job_tag == "24DH") {
    return 8;
  } else if (job_tag == "32Dfine") {
    return 10;
  } else if (job_tag == "48I") {
    return 12;
  } else if (job_tag == "64I") {
    return 18;
  } else {
    qassert(false);
  }
  return 8;
}
```

## ``compute-meson-vv-meson.h``

<img src="figs/matrix-elements/png/fig-3.png" width=400px />

```cpp
ssprintf("analysis/field-meson-vv-meson/%s/results=%d", job_tag.c_str(), traj)
```

Data format: ``FieldM<Complex, 8 * 8>`` with ``write_field_float_from_double``.

Use $y$ as the point source location in calculation.

```cpp
ssprintf("/forward-%d-%d-%d-%d.field", type1, type2, type3, type4)
```

$$
\ba
H_\text{forward}(x-y)[8\mu+\nu]
&=&\mathrm{Tr}[
S_1(t_\text{snk};x)\gamma^{\mathrm{va}}_\mu
S_4(x;y)
\gamma^{\mathrm{va}}_\nu S_2(y;t_\text{src})
\gamma_5 S_3(t_\text{src};t_\text{snk})\gamma_5
]
\ea
$$

where:
$$
\ba
t_\text{src} &=& \min(x_t,y_t) - t_\text{sep}
\\
t_\text{snk} &=& \max(x_t,y_t) + t_\text{sep}
\ea
$$
and for $t_\text{sep}$:

```cpp
inline int tsep_op_wall_src(const std::string& job_tag)
// parameter
{
  if (job_tag == "24D" or job_tag == "32D" or job_tag == "24DH") {
    return 8;
  } else if (job_tag == "32Dfine") {
    return 10;
  } else if (job_tag == "48I") {
    return 12;
  } else if (job_tag == "64I") {
    return 18;
  } else {
    qassert(false);
  }
  return 8;
}
```

## ``compute-meson-snk-src.h``

```cpp
ssprintf("analysis/lat-meson-snk-src/%s/results=%d", job_tag.c_str(), traj);
```

$$
\text{ld}
[0\le t_\text{snk}< T]
[0\le t_\text{src}< T]
$$

```cpp
ssprintf("/meson-snk-src-%d-%d.lat", type1, type2);
```

$$
\text{ld-1-2}[t_\text{snk}][t_\text{src}]
=
\mathrm{Tr}[
S_1(t_\text{snk};t_\text{src})
\gamma_5 S_2(t_\text{src};t_\text{snk})\gamma_5
]
$$

## ``compute-chvp.h``

```cpp
ssprintf("analysis/field-chvp/%s/results=%d", job_tag.c_str(), traj);
```

Data format: ``FieldM<Complex, 8 * 8>`` with ``write_field_float_from_double``.

Use $y$ as the point source location in calculation.

```cpp
ssprintf("/chvp-%d-%d.field", type1, type2);
```

$$
H_\text{chvp-1-2} (x-y) [8\mu+\nu] = \mathrm{Tr}[
S_1(x;y)\gamma^{\mathrm{va}}_\nu
S_2(y;x)\gamma^{\mathrm{va}}_\mu
]
$$

## ``compute-meson-chvp.h``

<img src="figs/matrix-elements/png/fig-5.png" width=400px />

```cpp
ssprintf("analysis/lat-meson-snk-src-shift-weight/%s/results=%d", job_tag.c_str(), traj);
```

$$
\text{ld}
[0\le t_\text{snk}< T]
[0\le t_\text{src}< T]
$$

```cpp
ssprintf("/meson-snk-src-%d-%d-%d-%d.field", type1, type2, type3, type4);
```

$$
\text{ld-1-2}[t_\text{snk}][t_\text{src}]
=
\mathrm{Tr}[
S_1(t_\text{snk};t_\text{src})
\gamma_5 S_2(t_\text{src};t_\text{snk})\gamma_5
]
$$

Weighted properly for different time slice according to number of point source propagator available.

```cpp
ssprintf("analysis/field-meson-chvp/%s/results=%d", job_tag.c_str(), traj);
```

Data format: ``FieldM<Complex, 8 * 8>`` with ``write_field_float_from_double``.

Use $y$ as the point source location in calculation.

```cpp
ssprintf("/mchvp-%d-%d-%d-%d.field", type1, type2, type3, type4);
```

$$
\ba
H_\text{1-2-3-4}(x-y)[8\mu+\nu]
&\texttt{ += }&
\mathrm{Tr}[
S_1(t_\text{snk};t_\text{src})
\gamma_5 S_2(t_\text{src};t_\text{snk})\gamma_5
]
\mathrm{Tr}[
S_3(x;y)\gamma^{\mathrm{va}}_\nu
S_4(y;x)\gamma^{\mathrm{va}}_\mu
]
\ea
$$

where:
$$
\ba
t_\text{src} &=& \min(x_t,y_t) - t_\text{sep}
\\
t_\text{snk} &=& \max(x_t,y_t) + t_\text{sep}
\ea
$$
and for $t_\text{sep}$:

```cpp
inline int tsep_op_wall_src(const std::string& job_tag)
// parameter
{
  if (job_tag == "24D" or job_tag == "32D" or job_tag == "24DH") {
    return 8;
  } else if (job_tag == "32Dfine") {
    return 10;
  } else if (job_tag == "48I") {
    return 12;
  } else if (job_tag == "64I") {
    return 18;
  } else {
    qassert(false);
  }
  return 8;
}
```

