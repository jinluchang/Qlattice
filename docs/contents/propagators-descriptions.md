# Summit and Oakforest propagators

## Locations

Data are stored at Long term storage provided by USQCD at BNL.

```
/sdcc/u/jluchang/qcdqedta/summit-oakforest-data
```

```
16IH2 -> qcddata/DWF/2+1f/16nt32/IWASAKI/b2.13/ls16/M1.8/ms0.04/ml0.01/Hlbl
24D -> qcddata/MDWF/2+1f/24nt64/IWASAKI+DSDR/b1.633/ls24b+c4/M1.8/ms0.0850/ml0.00107/Hlbl
24DH -> qcddata/MDWF/2+1f/24nt64/IWASAKI+DSDR/b1.633/ls24b+c4/M1.8/ms0.0985/ml0.0174/Hlbl
24IH1 -> qcddata/DWF/2+1f/24nt64/IWASAKI/b2.13/ls16/M1.8/ms0.04/ml0.005/Hlbl
24IH2 -> qcddata/DWF/2+1f/24nt64/IWASAKI/b2.13/ls16/M1.8/ms0.04/ml0.01/Hlbl
24IH3 -> qcddata/DWF/2+1f/24nt64/IWASAKI/b2.13/ls16/M1.8/ms0.04/ml0.02/Hlbl
24IH4 -> qcddata/DWF/2+1f/24nt64/IWASAKI/b2.13/ls16/M1.8/ms0.04/ml0.03/Hlbl
32D -> qcddata/MDWF/2+1f/32nt64/IWASAKI+DSDR/b1.633/ls24b+c4/M1.8/ms0.0850/ml0.00107/Hlbl
32IcoarseH1 -> qcddata/DWF/2+1f/32nt64/IWASAKI/b2.13/ls16/M1.8/ms0.04/mu0.005/Hlbl
32IfineH -> qcddata/DWF/2+1f/32nt64/IWASAKI/b2.37/ls12/M1.8/ms0.0186/mu0.0047/Hlbl
32IH1 -> qcddata/DWF/2+1f/32nt64/IWASAKI/b2.25/ls16/M1.8/ms0.03/mu0.004/Hlbl
32IH2 -> qcddata/DWF/2+1f/32nt64/IWASAKI/b2.25/ls16/M1.8/ms0.03/mu0.006/Hlbl
32IH3 -> qcddata/DWF/2+1f/32nt64/IWASAKI/b2.25/ls16/M1.8/ms0.03/mu0.008/Hlbl
48I -> qcddata/MDWF/2+1f/48nt96/IWASAKI/b2.13/ls24b+c2/M1.8/ms0.0362/mu0.00078/rhmc_H_R_G/HLbL
qcddata -> /sdcc/u/jluchang/archive-lqcd-long-term/RBC/qcddata
```

For each ensemble, typical contents include the following directories.
```
hlbl
prop-psrc-light
prop-psrc-strange
prop-smear-light
prop-smear-strange
field-selection.qar
gauge-transform.qar
point-selection.qar
point-selection-smear.qar
prop-rand-u1-charm.qar
prop-rand-u1-light.qar
prop-rand-u1-strange.qar
prop-wsrc-light.qar
prop-wsrc-strange.qar
psel-prop-psrc-light.qar
psel-prop-psrc-strange.qar
psel-prop-smear-light.qar
psel-prop-smear-strange.qar
psel-prop-wsrc-light.qar
psel-prop-wsrc-strange.qar
topo-measure.qar
wall-src-info-light.qar
wall-src-info-strange.qar
```

The `.qar` format is archive that allow easy access to its content without unfolding. It is described at [qar-format](qar-format.md).

## Propagators

Propagators are stored in sparse format. It include Coulomb gauge fixed wall source propagators:
```
prop-wsrc-light
prop-wsrc-strange
```
These are generated on Coulomb gauge fixed configurations (both source and sink are gauge fixed). The other propagators are not generated with gauge fixing.

```
prop-psrc-light
prop-psrc-strange
```

```
prop-smear-light
prop-smear-strange
```

```
prop-rand-u1-charm
prop-rand-u1-light
prop-rand-u1-strange
```

The information about the source locations, propagators types (0=light, 1=strange, 2=charm quark), inversion accuracy (0=sloppy, 1=moderate, 2=exact) can be inferred from the names or tags of the propagator. Some information are also available below.
```
point-selection.qar
point-selection-smear.qar
wall-src-info-light.qar
wall-src-info-strange.qar
```

### Sparsening

We only store randomly selected points of the propagators to save storage space. The ratio is usually `1/16`. Many complexities are associated with the sparsening.

The points which are stored include:

```
field-selection
point-selection
```
where `point-selection` stores the point source propagators point source locations, which are selected as well.

#### Different sinks:

There are also propagators with some special sinks are stored below. Including the wall sink propagators, smeared sink propagators, etc.
```
psel-prop-wsrc-light
psel-prop-wsrc-strange
psel-prop-smear-light
psel-prop-smear-strange
psel-prop-psrc-light
psel-prop-psrc-strange
```

### Relevant C++ functions:

```
read_field_selection
load_point_selection_info
```

```
const FieldSelection& fsel = get_field_selection(job_tag, traj);
const PointSelection& psel = get_point_selection(job_tag, traj);
FieldSelection fselc;
fselc.f_rank = fsel.f_rank;
add_field_selection(fselc.f_rank, psel);
update_field_selection(fselc);
update_field_selection(fselc, fsel.n_per_tslice);
ShuffledFieldsReader& sfr = get_shuffled_fields_reader(path);
sbs = mk_shuffled_bitset(fselc, sfr.new_size_node);
```

```
load_selected_points_complex
list_fields
read_field_double_from_float
set_selected_points(ps_prop, sprop, psel, sbs.fsel);
set_selected_field(s_prop, sprop, fsel, sbs.fsel);
```

### Example code that may serve as example:

[C++ code](https://github.com/waterret/Qlattice/blob/master/examples-cpp/load-select-data/data-load-base.h)

[Python code](https://github.com/waterret/Qlattice/blob/master/applications/auto-contract-bk-test/load_data.py)

