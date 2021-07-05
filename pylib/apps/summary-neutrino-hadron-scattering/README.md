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
  \gamma_5 & = & \gamma_x \gamma_y \gamma_z \gamma_t
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

The subscript of $S$ specifies the flavor of the quark. E.g. $S_1$ stand for quark propagator of $\text{type1}$, which is specified in the dataset tag, $0$ for light quark and $1$ for strange quark.

## HDF5 file content

### ``lmom-list``

List of momentum projections.
$$
\text{dataset}[l] = [ p_{l,x}, p_{l,y}, p_{l,z}, p_{l,t} ]
$$

### ``two-point-func-{type1}-{type2}``

$$
\text{dataset}[t_\text{sep}]
= \text{Tr}\Big( \big(\sum_\vec x S_1(\vec x,t_\text{snk};t_\text{src}) \gamma_5 \gamma_5
\sum_\vec y S_2(\vec y,t_\text{snk};t_\text{src})^\dagger \gamma_5\big) \gamma_5\Big)
$$

where $t_\text{sep} = t_\text{snk} - t_\text{src}$.

### ``three-point-{type1}-{type2}-{type3}``

$$
\ba
\text{dataset}[t_\text{sep}][t_\text{op}]
= \text{Tr}
\sum_\vec x
\Big(
S_1(\vec x,t_\text{op};t_\text{src})  \gamma_5 
S_3(t_\text{src};t_\text{snk})
\gamma_5
\big( \gamma_5 S_2(t_\text{snk};\vec x,t_\text{op})^\dagger \gamma_5 \big)
\gamma_t
\Big)
\ea
$$

where $t_\text{sep} = t_\text{snk} - t_\text{src}$.

### ``field-meson-vv-meson-{type1}-{type2}-{type3}-{type4}``

<img src="figs/matrix-elements/png/fig-3.png" width=400px />
$$
\ba
H_\text{forward}[l][x_t-y_t][\mu][\nu]
&=&
\int d^3\vec x\,
e^{-i p_l \cdot (x - y)}
\mathrm{Tr}[
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

```python
def get_tsep(job_tag):
    if job_tag in [ "24D", "32D", "24DH", ]:
        return 8
    elif job_tag == "32Dfine":
        return 10
    elif job_tag == "48I":
        return 12
    elif job_tag == "64I":
        return 18
    elif job_tag == "test-4nt16":
        return 2
    else:
        assert False
```

