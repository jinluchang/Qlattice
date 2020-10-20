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

### Pion or Kaon correlation function

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
\text{Tr}\big(\text{prop}_1(\vec x,t_\text{snk})\text{prop}_2(\vec x,t_\text{snk})^\dagger\big)
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
\text{Tr}\big(\sum_{\vec x}\text{prop}_1(\vec x,t_\text{snk})\sum_{\vec y}\text{prop}_2(\vec y,t_\text{snk})^\dagger\big)
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
= \text{Tr}\Big( \big(\sum_\vec x S_1(\vec x,t_\text{snk};t_\text{src}) \Gamma_{\text{op}_\text{src}} \gamma_5
S_2(\vec x,t_\text{snk};t_\text{src})^\dagger \gamma_5\big) \Gamma_{\text{op}_\text{snk}} \Big)
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
\ba
\text{ld}[t_\text{sep}][\text{op}_\text{src}][\text{op}_\text{snk}]
&=& \text{Tr}\Big( \big(\sum_\vec x S_1(\vec x,t_\text{snk};t_\text{src}) \Gamma_{\text{op}_\text{src}}
\sum_\vec y S_2(t_\text{src};\vec y,t_\text{snk})\big) \Gamma_{\text{op}_\text{snk}} \Big)
\\
&=& \text{Tr}\Big( \big(\sum_\vec x S_1(\vec x,t_\text{snk};t_\text{src}) \Gamma_{\text{op}_\text{src}} \gamma_5
\sum_\vec y S_2(\vec y,t_\text{snk};t_\text{src})^\dagger \gamma_5\big) \Gamma_{\text{op}_\text{snk}} \Big)
\ea
$$

### Three point function

```cpp
inline LatData mk_three_point_table(const Coordinate& total_site)
{
  LatData ld;
  ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
  ld.info.push_back(lat_dim_number("top", 0, total_site[3] - 1));
  ld.info.push_back(lat_dim_number("op", 0, 15));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  set_zero(ld);
  return ld;
}
```

$$
\text{ld}
[0\le t_\text{sep} < T]
[0\le t_\text{op} < T]
[0\le \text{op} < 16]
$$

```cpp
inline LatData contract_three_point_function(const SelProp& prop_a,
                                             const SelProp& prop_b,
                                             const WilsonMatrix& wm_ab,
                                             const int ta, const int tb,
                                             const FieldSelection& fsel);
// ``wm_ab'' is prop from ``tb'' to ``ta''.
// |  ->- prop_a ->- op ->- inv prop_b ->- |
// a (gamma5)                              b (gamma5)
// |            -<- wm_ab -<-              |
//
// prop_a (type1)
// prop_b (type2)
// wm_ab (type3)

inline LatData contract_three_point_function(
    const WallSrcProps& wsp1, const WallSrcProps& wsp2,
    const WallSrcProps& wsp3, const FieldSelection& fsel,
    const int yt_measurement_sparsity = 1, const int yt_measurement_start = 0);
```

$$
\text{ld}[t_\text{sep}][t_\text{op}][\text{op}]
= \text{Tr}\Big(
\sum_\vec x
\big( \gamma_5 S_2(t_\text{snk};t_\text{op},\vec x)^\dagger \gamma_5 \big)
\Gamma_{\text{op}} S_1(t_\text{op},\vec x;t_\text{src})  \gamma_5 
S_3(t_\text{src};t_\text{snk})
\gamma_5
\Big)
$$



