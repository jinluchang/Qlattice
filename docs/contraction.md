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

### utils

```cpp
template <class M>
void reflect_field(Field<M>& f);
```

From ``field-shuffle.h``
$$
H(x) \to H(-x)
$$

```cpp
inline void field_permute_mu_nu(FieldM<Complex, 8 * 8>& f);
```

$$
H(x)[8\mu+\nu] = H(x)[8\nu+\mu]
$$

```cpp
inline void field_conjugate_mu_nu(FieldM<Complex, 8 * 8>& f);
```

$$
H(x)[8\mu+\nu] = \theta_\mu \theta_\nu H(x)[8\mu+\nu]
$$

```cpp
inline void field_complex_conjugate(Field<Complex>& f);
```

$$
H(x) = H(x)^*
$$

### meson-vv

Use ``xg_y`` as point source location, contraction for all available ``xg_x``.

Proper factor is compensated so it can treated as if ``xg_x`` is contracted for all lattice sites.

```cpp
inline void contract_meson_vv_acc(
    FieldM<Complex, 8 * 8>& decay, FieldM<Complex, 8 * 8>& fission,
    const WallSrcProps& wsp1, const WallSrcProps& wsp2,
    const SelProp& prop3_x_y, const Coordinate& xg_y, const long xg_y_psel_idx,
    const int tsep, const PointSelection& psel, const FieldSelection& fsel,
    const ShiftShufflePlan& ssp);
// xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
// ssp = make_shift_shuffle_plan(fsel, -xg_y);
```

$$
\ba
H_\text{decay-1-2-3}(x-y)[8\mu+\nu]
&\texttt{ += }&
\mathrm{Tr}
[S_3(x;y)\gamma^{\mathrm{va}}_\nu S_2(y;t_\text{src})\gamma_5 S_1(t_\text{src};x)\gamma^{\mathrm{va}}_\mu]
\\
H_\text{fission-1-2-3}(x-y)[8\mu+\nu]
&\texttt{ += }&
\mathrm{Tr}
[S_3(x;y)\gamma^{\mathrm{va}}_\nu S_2(y;t_\text{snk})\gamma_5 S_1(t_\text{snk};x)\gamma^{\mathrm{va}}_\mu]
\ea
$$

Some properties:
$$
\ba
\big(H_\text{decay-1-2-3}(x-y)[8\mu+\nu]\big)^\dagger
&=&
\mathrm{Tr}
[S_3(x;y)\gamma^{\mathrm{va}}_\nu S_2(y;t_\text{src})\gamma_5 S_1(t_\text{src};x)\gamma^{\mathrm{va}}_\mu]^\dagger
\\
&=&
\mathrm{Tr}
[{\gamma^{\mathrm{va}}_\mu}^\dagger S_1(x;t_\text{src}) \gamma_5 S_2(t_\text{src};y){\gamma^{\mathrm{va}}_\nu}^\dagger S_3(y;x)]
\\
&=&
\mathrm{Tr}
[ S_3(y;x) {\gamma^{\mathrm{va}}_\mu}^\dagger S_1(x;t_\text{src}) \gamma_5 S_2(t_\text{src};y){\gamma^{\mathrm{va}}_\nu}^\dagger]
\\
&=&
\theta_\mu \theta_\nu H_\text{decay-2-1-3}(y-x)[8\nu+\mu]
\ea
$$

$$
\big(H_\text{fission-1-2-3}(x-y)[8\mu+\nu]\big)^\dagger
=
\theta_\mu \theta_\nu  H_\text{fission-2-1-3}(y-x)[8\nu+\mu]
$$

$$
\ba
H_\text{fission-1-2-3}(x-y)[8\mu+\nu]
&\iff&
\mathrm{Tr}
[S_3(-x;-y)\gamma^{\mathrm{va}}_\nu S_2(-y;-t_\text{snk})\gamma_5 S_1(-t_\text{snk};-x)\gamma^{\mathrm{va}}_\mu]
\\
&=&
\mathrm{Tr}
[S_3(y;x)\gamma^{\mathrm{va}}_\nu S_2(x;t_\text{src})\gamma_5 S_1(t_\text{src};y)\gamma^{\mathrm{va}}_\mu]
\\
&=&
H_\text{decay-1-2-3}(y-x)[8\mu+\nu]
\ea
$$

Possible post processing:

```cpp
reflect_field(fission);
decay += fission;
decay *= 0.5;
```

and

```cpp
FieldM<Complex, 8 * 8> avg;
avg = decay_2_1_3;
reflect_field(avg);
field_permute_mu_nu(avg);
field_conjugate_mu_nu(avg);
field_complex_conjugate(avg);
avg += decay_1_2_3;
avg *= 0.5;
```

### meson-vv-meson

Use ``xg_y`` as point source location, contraction for all available ``xg_x``.

Proper factor is compensated so it can treated as if ``xg_x`` is contracted for all lattice sites.

```cpp
inline void contract_meson_vv_meson_acc(
    FieldM<Complex, 8 * 8>& forward, FieldM<Complex, 8 * 8>& backward,
    const WallSrcProps& wsp1, const WallSrcProps& wsp2,
    const WallSrcProps& wsp3, const SelProp& prop4_x_y, const Coordinate& xg_y,
    const long xg_y_psel_idx, const int tsep, const PointSelection& psel,
    const FieldSelection& fsel, const ShiftShufflePlan& ssp);
// xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
// ssp = make_shift_shuffle_plan(fsel, -xg_y);
```

$$
\ba
H_\text{forward-1-2-3-4}(x-y)[8\mu+\nu]
&\texttt{ += }&
\mathrm{Tr}[
S_1(t_\text{snk};x)\gamma^{\mathrm{va}}_\mu
S_4(x;y)
\gamma^{\mathrm{va}}_\nu S_2(y;t_\text{src})
\gamma_5 S_3(t_\text{src};t_\text{snk})\gamma_5
]
\\
H_\text{backward-1-2-3-4}(x-y)[8\mu+\nu]
&\texttt{ += }&\mathrm{Tr}[
S_1(t_\text{src};x)\gamma^{\mathrm{va}}_\mu
S_4(x;y)
\gamma^{\mathrm{va}}_\nu S_2(y;t_\text{snk})
\gamma_5 S_3(t_\text{snk};t_\text{src})\gamma_5
]
\ea
$$

Some properties:
$$
\ba
&&\hspace{-2cm}
\big(H_\text{forward-1-2-3-4}(x-y)[8\mu+\nu]\big)^\dagger
\nn\\
&=&
\mathrm{Tr}[
S_1(t_\text{snk};x)\gamma^{\mathrm{va}}_\mu
S_4(x;y)
\gamma^{\mathrm{va}}_\nu S_2(y;t_\text{src})
\gamma_5 S_3(t_\text{src};t_\text{snk})\gamma_5
]^\dagger
\\
&=&
\mathrm{Tr}[
\gamma_5 S_3(t_\text{snk};t_\text{src})\gamma_5
S_2(t_\text{src};y) {\gamma^{\mathrm{va}}_\nu}^\dagger
S_4(y;x)
{\gamma^{\mathrm{va}}_\mu}^\dagger S_1(x;t_\text{snk})
]
\\
&=&
\mathrm{Tr}[
S_2(t_\text{src};y) {\gamma^{\mathrm{va}}_\nu}^\dagger
S_4(y;x)
{\gamma^{\mathrm{va}}_\mu}^\dagger S_1(x;t_\text{snk})
\gamma_5 S_3(t_\text{snk};t_\text{src})\gamma_5
]
\\
&=&
\theta_\mu \theta_\nu H_\text{backward-2-1-3-4}(y-x)[8\nu+\mu]
\ea
$$

$$
\ba
\big(H_\text{backward-1-2-3-4}(x-y)[8\mu+\nu]\big)^\dagger
&=&
\theta_\mu \theta_\nu H_\text{forward-2-1-3-4}(y-x)[8\nu+\mu]
\ea
$$

$$
\ba
&&\hspace{-2cm}H_\text{backward-1-2-3-4}(x-y)[8\mu+\nu]
\nn\\
&\iff&
\mathrm{Tr}[
S_1(-t_\text{src};-x)\gamma^{\mathrm{va}}_\mu
S_4(-x;-y)
\gamma^{\mathrm{va}}_\nu S_2(-y;-t_\text{snk})
\gamma_5 S_3(-t_\text{snk};-t_\text{src})\gamma_5
]
\\
&=&
\mathrm{Tr}[
S_1(t_\text{snk};y)\gamma^{\mathrm{va}}_\mu
S_4(y;x)
\gamma^{\mathrm{va}}_\nu S_2(x;t_\text{src})
\gamma_5 S_3(t_\text{src};t_\text{snk})\gamma_5
]
\\
&=&
H_\text{forward-1-2-3-4}(y-x)[8\mu+\nu]
\ea
$$

Possible post processing:

```cpp
reflect_field(backward);
forward += backward;
forward *= 0.5;
```

and

```cpp
FieldM<Complex, 8 * 8> tmp;
tmp = forward_2_1_3_4;
field_permute_mu_nu(tmp);
field_conjugate_mu_nu(tmp);
field_complex_conjugate(tmp);
forward_1_2_3_4 += tmp;
forward_1_2_3_4 *= 0.5;
```

### meson-snk-src

```cpp
inline LatData contract_meson(const WallSrcProps& wsp1,
                              const WallSrcProps& wsp2)
```

$$
\text{ld-1-2}[t_\text{snk}][t_\text{src}]
=
\mathrm{Tr}[
S_1(t_\text{snk};t_\text{src})
\gamma_5 S_2(t_\text{src};t_\text{snk})\gamma_5
]
$$

Some properties:
$$
\text{ld-1-2}[t_\text{snk}][t_\text{src}]^\dagger
=
\text{ld-2-1}[t_\text{snk}][t_\text{src}]
=
\text{ld-1-2}[t_\text{src}][t_\text{snk}]
$$

## chvp

Use ``xg_y`` as point source location, contraction for all available ``xg_x``.

Proper factor is compensated so it can treated as if ``xg_x`` is contracted for all lattice sites.


```cpp
inline void contract_chvp(SelectedField<Complex>& chvp,
                          const SelProp& prop1_x_y, const SelProp& prop2_x_y,
                          const FieldSelection& fsel);
```

$$
H_\text{chvp-1-2} (x-y) [8\mu+\nu] = \mathrm{Tr}[
S_1(x;y)\gamma^{\mathrm{va}}_\nu
S_2(y;x)\gamma^{\mathrm{va}}_\mu
]
$$

Some properties:
$$
H_\text{chvp-1-2} (x-y) [8\mu+\nu] = H_\text{chvp-2-1} (y-x)[8\nu+\mu]
$$

$$
\ba
\big(H_\text{chvp-1-2} (x-y) [8\mu+\nu] \big)^\dagger
&=&
\mathrm{Tr}[
{\gamma^{\mathrm{va}}_\mu}^\dagger
S_2(x;y){\gamma^{\mathrm{va}}_\nu}^\dagger
S_1(y;x)
]
\\
&=&
\theta_\mu \theta_\nu
H_\text{chvp-2-1}(x-y)[8\mu+\nu]
\\
&=&
\theta_\mu \theta_\nu
H_\text{chvp-1-2}(y-x)[8\nu+\mu]
\ea
$$

Possible post processing:

```cpp
FieldM<Complex, 8 * 8> tmp;
tmp = chvp;
reflect_field(tmp);
field_permute_mu_nu(tmp);
field_conjugate_mu_nu(tmp);
field_complex_conjugate(tmp);
chvp += tmp;
chvp *= 0.5;
```

### meson-chvp

Use ``xg_y`` as point source location, contraction for all available ``xg_x``.

Proper factor is compensated so it can treated as if ``xg_x`` is contracted for all lattice sites.


```cpp
inline void contract_meson_chvp_acc(FieldM<Complex, 8 * 8>& mchvp,
                                    const LatData& ld_meson_snk_src_1_2,
                                    const FieldM<Complex, 8 * 8>& chvp_3_4,
                                    const int tsep);
// ld_meson_snk_src_1_2 should already shifted to origin (t_y -> 0)
// chvp_3_4 should already shifted to origin (xg_y -> 0)
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
Some properties:
$$
H_\text{1-2-3-4}(x-y)[8\mu+\nu] = H_\text{1-2-4-3}(y-x)[8\nu+\mu]
$$

$$
\ba
\big(H_\text{1-2-3-4}(x-y)[8\mu+\nu]\big)^\dagger
&=&
\mathrm{Tr}[
S_1(t_\text{src};t_\text{snk})
\gamma_5 S_2(t_\text{snk};t_\text{src})\gamma_5
]
\mathrm{Tr}[
S_4(x;y){\gamma^{\mathrm{va}}_\nu}^\dagger
S_3(y;x){\gamma^{\mathrm{va}}_\mu}^\dagger
]
\\
&=&
\theta_\mu \theta_\nu
H_\text{2-1-4-3}(x-y)[8\mu+\nu]
\\
&=&
\theta_\mu \theta_\nu
H_\text{2-1-3-4}(y-x)[8\nu+\mu]
\ea
$$

$$
\ba
&&\hspace{-2cm}H_\text{1-2-3-4}(x-y)[8\mu+\nu]
\nn\\
&\iff&
\mathrm{Tr}[
S_1(t_\text{src};t_\text{snk})
\gamma_5 S_2(t_\text{snk};t_\text{src})\gamma_5
]
\mathrm{Tr}[
S_3(y;x)\gamma^{\mathrm{va}}_\nu
S_4(x;y)\gamma^{\mathrm{va}}_\mu
]
\\
&=&
H_\text{2-1-3-4}(y-x)[8\mu+\nu]
\\
&=&
\big(\theta_\mu \theta_\nu H_\text{1-2-3-4}(x-y)[8\nu+\mu]\big)^\dagger
\ea
$$

Possible post processing:

```cpp
FieldM<Complex, 8 * 8> tmp;
tmp = mchvp_2_1_3_4;
reflect_field(tmp);
tmp += mchvp_1_2_3_4;
tmp *= 0.5;
mchvp_1_2_3_4 = tmp;
field_permute_mu_nu(tmp);
field_conjugate_mu_nu(tmp);
field_complex_conjugate(tmp);
mchvp_1_2_3_4 += tmp;
mchvp_1_2_3_4 *= 0.5;
```

### meson-v-v-meson

Selectively contract all possible combinations. No point source propagators are needed.

Proper factor is compensated so it can treated as if ``xg_x`` is contracted for all lattice sites and ``xg_y`` is contracted only on a single site.

```cpp
inline void contract_meson_v_v_meson_acc(
    FieldM<Complex, 8 * 8>& mf_v_v,
    const WallSrcProps& wsp1, const WallSrcProps& wsp2,
    const WallSrcProps& wsp3, const SelProp& prop4_x_y, const Coordinate& xg_y,
    const long xg_y_psel_idx, const int tsep, const PointSelection& psel,
    const FieldSelection& fsel, const ShiftShufflePlan& ssp);
// xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
// ssp = make_shift_shuffle_plan(fsel, -xg_y);
```

$$
\ba
H_\text{1-2-3-4}(x-y)[8\mu+\nu]
&\texttt{ += }&
\mathrm{Tr}[
S_1(t_\text{snk};x)\gamma^{\mathrm{va}}_\mu
S_2(x;t_\text{src})\gamma_5
S_3(t_\text{src};y) \gamma^{\mathrm{va}}_\nu 
S_4(y;t_\text{snk})\gamma_5
]
\ea
$$

