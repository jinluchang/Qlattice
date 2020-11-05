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

### pion or kaon correlation function

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

### two point function

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

### three point function

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

## ```contraction-field.h```

$$
\ba
t_\text{src} &=& \min(x_t,y_t) - t_\text{sep}
\\
t_\text{snk} &=& \max(x_t,y_t) + t_\text{sep}
\ea
$$

### psel-fsel-distribution

Use ``xg_y`` as point source location, contraction for all available ``xg_x``.

Proper factor is compensated so it can treated as if ``xg_x`` is contracted for all lattice sites.

```cpp
inline void contract_psel_fsel_distribution_acc(FieldM<Complex, 1>& pos,
                                                const Coordinate& xg_y,
                                                const FieldSelection& fsel,
                                                const ShiftShufflePlan& ssp);
// ssp = make_shift_shuffle_plan(fsel, -xg_y);

template <class M>
void rescale_field_with_psel_fsel_distribution(Field<M>& f,
                                               const FieldM<Complex, 1>& pfdist);
```

$$
H(x-y) \texttt{ += } 1
$$

### meson-vv

Use ``xg_y`` as point source location, contraction for all available ``xg_x``.

Proper factor is compensated so it can treated as if ``xg_x`` is contracted for all lattice sites.

```cpp
inline void contract_meson_vv_acc(
    FieldM<Complex, 8 * 8>& meson_vv_decay,
    FieldM<Complex, 8 * 8>& meson_vv_fission, const WallSrcProps& wsp1,
    const WallSrcProps& wsp2, const SelProp& prop3_x_y, const Coordinate& xg_y,
    const long xg_y_psel_idx, const int tsep, const PointSelection& psel,
    const FieldSelection& fsel, const ShiftShufflePlan& ssp,
    const ShiftShufflePlan& ssp_reflect);
// xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
// ssp = make_shift_shuffle_plan(fsel, -xg_y);
// ssp_reflect = make_shift_shuffle_plan(fsel, -xg_y, true);
```

$$
\ba
H_\text{decay}(x-y)[8\mu+\nu]
&\texttt{ += }&
\frac{1}{2}\mathrm{Tr}
[S_3(x;y)\gamma^{\mathrm{va}}_\nu S_2(y;t_\text{src})\gamma_5 S_1(t_\text{src};x)\gamma^{\mathrm{va}}_\mu]
\\
H_\text{decay}(y-x)[8\mu+\nu]
&\texttt{ += }&
\frac{1}{2}\mathrm{Tr}
[S_3(y;x) \gamma^{\mathrm{va}}_\nu S_2(x;t_\text{src})\gamma_5 S_1(t_\text{src};y)\gamma^{\mathrm{va}}_\mu]
\ea
$$

$$
\ba
H_\text{fission}(x-y)[8\mu+\nu]
&\texttt{ += }&
\frac{1}{2}\mathrm{Tr}
[S_3(x;y)\gamma^{\mathrm{va}}_\nu S_2(y;t_\text{snk})\gamma_5 S_1(t_\text{snk};x)\gamma^{\mathrm{va}}_\mu]
\\
H_\text{fission}(y-x)[8\mu+\nu]
&\texttt{ += }&
\frac{1}{2}\mathrm{Tr}
[S_3(y;x) \gamma^{\mathrm{va}}_\nu S_2(x;t_\text{snk})\gamma_5 S_1(t_\text{snk};y)\gamma^{\mathrm{va}}_\mu]
\ea
$$

### meson-vv-meson

Use ``xg_y`` as point source location, contraction for all available ``xg_x``.

Proper factor is compensated so it can treated as if ``xg_x`` is contracted for all lattice sites.

```cpp
inline void contract_meson_vv_meson_acc(
    FieldM<Complex, 8 * 8>& meson_vv_meson_forward,
    FieldM<Complex, 8 * 8>& meson_vv_meson_backward, const WallSrcProps& wsp1,
    const WallSrcProps& wsp2, const WallSrcProps& wsp3,
    const SelProp& prop4_x_y, const Coordinate& xg_y, const long xg_y_psel_idx,
    const int tsep, const PointSelection& psel, const FieldSelection& fsel,
    const ShiftShufflePlan& ssp, const ShiftShufflePlan& ssp_reflect);
// xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
// ssp = make_shift_shuffle_plan(fsel, -xg_y);
// ssp_reflect = make_shift_shuffle_plan(fsel, -xg_y, true);
```

$$
\ba
H_\text{forward}(x-y)[8\mu+\nu]
&\texttt{ += }&
\frac{1}{2}\mathrm{Tr}[
S_1(t_\text{snk};x)\gamma^{\mathrm{va}}_\mu
S_4(x;y)
\gamma^{\mathrm{va}}_\nu S_2(y;t_\text{src})
\gamma_5 S_3(t_\text{src};t_\text{snk})\gamma_5
]
\\
H_\text{forward}(y-x)[8\mu+\nu]
&\texttt{ += }&
\frac{1}{2}\mathrm{Tr}[
S_1(t_\text{snk};y)\gamma^{\mathrm{va}}_\mu
S_4(y;x)
\gamma^{\mathrm{va}}_\nu S_2(x;t_\text{src})
\gamma_5 S_3(t_\text{src};t_\text{snk})\gamma_5
]
\ea
$$

$$
\ba
H_\text{backward}(x-y)[8\mu+\nu]
&\texttt{ += }&
\frac{1}{2}\mathrm{Tr}[
S_1(t_\text{src};x)\gamma^{\mathrm{va}}_\mu
S_4(x;y)
\gamma^{\mathrm{va}}_\nu S_2(y;t_\text{snk})
\gamma_5 S_3(t_\text{snk};t_\text{src})\gamma_5
]
\\
H_\text{backward}(y-x)[8\mu+\nu]
&\texttt{ += }&
\frac{1}{2}\mathrm{Tr}[
S_1(t_\text{src};y)\gamma^{\mathrm{va}}_\mu
S_4(y;x)
\gamma^{\mathrm{va}}_\nu S_2(x;t_\text{snk})
\gamma_5 S_3(t_\text{snk};t_\text{src})\gamma_5
]
\ea
$$

### chvp

Use ``xg_y`` as point source location, contraction for all available ``xg_x``.

Proper factor is compensated so it can treated as if ``xg_x`` is contracted for all lattice sites.

```cpp
inline void contract_chvp(
    FieldM<Complex, 8 * 8>& chvp,
    FieldM<Complex, 8 * 8>& meson_vv_meson_backward, const WallSrcProps& wsp1,
    const WallSrcProps& wsp2, const WallSrcProps& wsp3,
    const SelProp& prop4_x_y, const Coordinate& xg_y, const long xg_y_psel_idx,
    const int tsep, const PointSelection& psel, const FieldSelection& fsel,
    const ShiftShufflePlan& ssp, const ShiftShufflePlan& ssp_reflect);
// xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
// ssp = make_shift_shuffle_plan(fsel, -xg_y);
// ssp_reflect = make_shift_shuffle_plan(fsel, -xg_y, true);
```

