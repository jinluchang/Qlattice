# Contraction

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

## ```contraction-pion.h```

### Only Pion or Kaon correlation function

```cpp
inline LatData mk_pion_corr_table(const Coordinate& total_site)
{
  LatData ld;
  ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  set_zero(ld);
  return ld;
}
```

$$
\text{ld}[0\le t_\text{sep} < T]
$$

```cpp
inline LatData contract_pion(const SelProp& prop, const int tslice_src,
                             const FieldSelection& fsel);

inline LatData contract_kaon(const SelProp& prop1, const SelProp& prop2,
                             const int tslice_src, const FieldSelection& fsel);
```

$$
\text{ld}[t_\text{sep}] =
\sum_{\vec x}
\text{Tr}\big(\text{prop}_1(t_\text{snk},\vec x)\text{prop}_2(t_\text{snk},\vec x)^\dagger\big)
$$

```cpp
inline LatData contract_pion_wall_snk(const SelProp& prop, const int tslice_src,
                                      const FieldSelection& fsel);

inline LatData contract_kaon_wall_snk(const SelProp& prop1,
                                      const SelProp& prop2,
                                      const int tslice_src,
                                      const FieldSelection& fsel);
```

$$
\text{ld}[t_\text{sep}] =
\text{Tr}\big(\sum_{\vec x}\text{prop}_1(t_\text{snk},\vec x)\sum_{\vec y}\text{prop}_2(t_\text{snk},\vec y)^\dagger\big)
$$

### Two point function

```cpp
inline LatData mk_two_point_table(const Coordinate& total_site)
{
  LatData ld;
  ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
  ld.info.push_back(lat_dim_number("op-src", 0, 15));
  ld.info.push_back(lat_dim_number("op-snk", 0, 15));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  set_zero(ld);
  return ld;
}
```

$$
\text{ld}
[0\le t_\text{sep} < T]
[0\le\text{op}_\text{src}<16]
[0\le\text{op}_\text{snk}<16]
$$

```cpp
inline LatData contract_two_point_function(const SelProp& prop1,
                                           const SelProp& prop2,
                                           const int tslice,
                                           const FieldSelection& fsel);

inline LatData contract_two_point_function(const WallSrcProps& wsp1,
                                           const WallSrcProps& wsp2,
                                           const FieldSelection& fsel);
```

$$
\text{ld}[t_\text{sep}][\text{op}_\text{src}][\text{op}_\text{snk}]
= \text{Tr}\Big( \big(\sum_\vec x \text{prop}_1(t_\text{snk},\vec x) \Gamma_{\text{op}_\text{src}} \gamma_5
\text{prop}_2(t_\text{snk},\vec x)^\dagger \gamma_5\big) \Gamma_{\text{op}_\text{snk}} \Big)
$$

```cpp
inline LatData contract_two_point_wall_snk_function(
    const std::vector<WilsonMatrix>& prop1,
    const std::vector<WilsonMatrix>& prop2, const int tslice,
    const Coordinate& total_site);
// need to be sparse corrected

inline LatData contract_two_point_wall_snk_function(const WallSrcProps& wsp1,
                                                    const WallSrcProps& wsp2,
                                                    const FieldSelection& fsel);
// need to be sparse corrected

inline LatData contract_two_point_wall_snk_function(
    const LatData& ld_two_point_wall_snk_func, const LatData& ld_two_point_func,
    const FieldSelection& fsel);
// perform sparse correction
```

$$
\text{ld}[t_\text{sep}][\text{op}_\text{src}][\text{op}_\text{snk}]
= \text{Tr}\Big( \big(\sum_\vec x S_1(t_\text{src};t_\text{snk},\vec x) \Gamma_{\text{op}_\text{src}} \gamma_5
\sum_\vec y S_2(t_\text{src};t_\text{snk},\vec y)^\dagger \gamma_5\big) \Gamma_{\text{op}_\text{snk}} \Big)
$$



