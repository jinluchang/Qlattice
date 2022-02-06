# Qlat Core



[TOC]

## Structure

```c++
const int DIMN = 4;
struct Coordinate : public array<int, DIMN> {};
```

### GeometryNode

```c++
struct GeometryNode {
  bool initialized;
  int num_node;
  int id_node;
  Coordinate size_node;
  Coordinate coor_node;
};
```

### Geometry

```c++
struct Geometry {
  bool initialized;
  GeometryNode geon;
  int eo;  // 0:full; 1:odd; 2:even
  int multiplicity;
  Coordinate node_site;
  Coordinate expansion_left;
  Coordinate expansion_right;
  Coordinate node_site_expanded;
};
```

### Field

```c++
template <class M>
struct Field {
  // Avoid copy constructor when possible
  // (it is likely not be what you think it is)
  bool initialized;
  box_acc<Geometry> geo;
  vector_acc<M> field;
  //
  void init();
  void init(const Geometry& geo_);
  // only initialize if uninitilized
  // if initialized already, then check for matching geo (including
  // multiplicity)
  // can have different geo expansion
  void init(const Geometry& geo_, const int multiplicity_);
  // only initialize if uninitilized
  // if initialized already, then check for matching geo (including
  // multiplicity)
  // can have different geo expansion
  void init(const Field<M>& f);
  // initialize to be identical to f if uninitilized
  // otherwise use assignment operator
  //
  Field() { init(); }
  Field(const Field<M>&) = default;
  Field(Field<M>&&) noexcept = default;
  //
  const Field<M>& operator=(const Field<M>& f);
  // skip if same object
  // otherwise:
  // 1. assert f is initialized
  // 2. init with geo_resize(f.geo())
  // 3. copy content
  //
  qacc M& get_elem(const long offset);
  qacc const M& get_elem(const long offset) const;
  //
  qacc M& get_elem(const Coordinate& x, const int m);
  qacc const M& get_elem(const Coordinate& x, const int m) const;
  //
  qacc M& get_elem(const Coordinate& x);
  qacc const M& get_elem(const Coordinate& x) const;
  //
  qacc Vector<M> get_elems(const Coordinate& x);
  qacc Vector<M> get_elems_const(const Coordinate& x) const;
  // Be cautious about the const property
  // 改不改靠自觉
  //
  qacc Vector<M> get_elems(const long index);
  // qassert(geo().is_only_local())
  qacc Vector<M> get_elems_const(const long index) const;
  // Be cautious about the const property
  // 改不改靠自觉
  // qassert(geo().is_only_local())
};
```

### FieldM

```c++
template <class M, int multiplicity>
struct FieldM : Field<M> {
  using Field<M>::init;
  void init(const Geometry& geo_)
  {
    Field<M>::init(geo_, multiplicity);
  }
  void init(const Geometry& geo_, const int multiplicity_)
  {
    qassert(multiplicity == multiplicity_);
    Field<M>::init(geo_, multiplicity);
  }
  void init(const Field<M>& f)
  {
    qassert(multiplicity == f.geo().multiplicity);
    Field<M>::init(f);
  }
  //
  FieldM<M, multiplicity>() { init(); }
};
```

### FieldSelection

```c++
struct FieldSelection {
  FieldM<int64_t, 1>
      f_rank;  // rank when the points being selected (-1 if not selected)
  //
  long n_per_tslice;  // num points per time slice (not enforced and should work
                      // properly if not true)
  double prob;        // (double)n_per_tslice / (double)spatial_vol
  //
  FieldM<long, 1>
      f_local_idx;  // idx of points on this node (-1 if not selected)
  long n_elems;     // num points of this node
  //
  vector_acc<int64_t> ranks;  // rank of the selected points
  vector_acc<long> indices;   // local indices of selected points
  //
  void init();
  //
  FieldSelection() { init(); }
};
```

### SelectedField

```c++
template <class M>
struct SelectedField {
  // Avoid copy constructor when possible
  // (it is likely not be what you think it is)
  //
  bool initialized;
  long n_elems;
  box_acc<Geometry> geo;
  vector_acc<M> field;
  //
  void init();
  void init(const Geometry& geo_, const long n_elems_, const int multiplicity);
  void init(const FieldSelection& fsel, const int multiplicity);
  //
  SelectedField() { init(); }
  //
  qacc M& get_elem(const long& idx);
  qacc const M& get_elem(const long& idx) const;
  //
  Vector<M> get_elems(const long idx);
  Vector<M> get_elems_const(const long idx) const;
  // Be cautious about the const property
  // 改不改靠自觉
};
```

## Functions

### Get data

```c++
template <class M>
qacc Vector<M> get_data(const Field<M>& f);

template <class M>
Vector<M> get_data(const SelectedField<M>& sf);
```

### Set zero

```c++
template <class M>
void set_zero(SelectedField<M>& sf);
```

### Qnorm

```c++
template <class M>
double qnorm(const Field<M>& f);

template <class M>
double qnorm(const SelectedField<M>& sf);
```

```c++
template <class M>
double qnorm_double(const Field<M>& f1, const Field<M>& f2);
```

### Qswap

```c++
template <class M>
void qswap(Field<M>& f1, Field<M>& f2);

template <class M>
void qswap(SelectedField<M>& f1, SelectedField<M>& f2);
```

### Field Selection

```c++
inline void add_field_selection(FieldM<int64_t, 1>& f_rank,
                                const PointSelection& psel,
                                const long rank_psel = 1024L * 1024L * 1024L *
                                                       1024L * 1024L);
// add psel points to f_rank. (only lower rank if already selected)

inline void mk_field_selection(FieldM<int64_t, 1>& f_rank,
                               const Coordinate& total_site,
                               const int64_t val = 0);
// select everything with val
// default val = 0 ; means selection everything
// val = -1 deselection everything

inline void mk_field_selection(FieldM<int64_t, 1>& f_rank,
                               const Coordinate& total_site,
                               const std::vector<Coordinate>& xgs,
                               const long rank_xgs = 1024L * 1024L * 1024L *
                                                     1024L * 1024L);

inline void select_rank_range(FieldM<int64_t, 1>& f_rank,
                              const long rank_start = 0,
                              const long rank_stop = -1);
// keep rank info if rank_start <= rank and (rank < rank_stop or rank_stop == -1)
// otherwise rank = -1
// default parameter does not change selection
// but will erase the rank information for points not selected (rank = -1)

inline void select_t_range(FieldM<int64_t, 1>& f_rank, const long t_start = 0,
                           const long t_stop = -1);
// keep rank info if t_start <= t and (t < t_stop or t_stop == -1)
// otherwise rank = -1
// default parameter does not change selection
// but will erase the rank information for points not selected (rank = -1)

inline void set_n_per_tslice(FieldM<int64_t, 1>& f_rank,
                             const long n_per_tslice);
// will erase the rank information for points not selected (rank = -1)
//
// if n_per_tslice == -1 then all points are selected regardless of rank
// (if a point was not selected before (rank < 0), then rank will be set to be
// rank = spatial_vol).
//
// otherwise: n_per_tslice is not enforced but only serve as an limit for f_rank
//
// 0 <= rank < n_per_tslice
//
// if point is not selected, rank = -1

inline void update_field_selection(FieldSelection& fsel);
// update fsel based only on f_rank
// do not touch n_per_tslice and prob at all

inline void update_field_selection(FieldSelection& fsel,
                                   const long n_per_tslice_);
// only adjust parameter, do not change contents

inline void set_field_selection(FieldSelection& fsel,
                                const FieldM<int64_t, 1>& f_rank,
                                const long n_per_tslice_ = 0,
                                const bool is_limit_on_rank = false);
// call set_n_per_tslice if is_limit_on_rank = true
// otherwise will strictly follow f_rank without constraint of n_per_tslice

inline void set_field_selection(FieldSelection& fsel,
                                const Coordinate& total_site);
// select everything with rank = 0

inline PointSelection psel_from_fsel(const FieldSelection& fsel);
```

### Selected Field

```c++
template <class M>
bool is_consistent(const SelectedField<M>& sf, const FieldSelection& fsel);

inline void set_selected_gindex(SelectedField<long>& sfgi,
                                const FieldSelection& fsel);

template <class M>
void only_keep_selected_points(Field<M>& f, const FieldSelection& fsel);

template <class M>
void set_selected_field(SelectedField<M>& sf, const Field<M>& f,
                        const FieldSelection& fsel);

template <class M>
void set_selected_field(SelectedField<M>& sf, const SelectedField<M>& sf0,
                        const FieldSelection& fsel, const FieldSelection& fsel0);
// Does not clear sf's original value if not assigned

template <class M>
void set_selected_points(SelectedPoints<M>& sp, const SelectedField<M>& sf,
                         const PointSelection& psel, const FieldSelection& fsel);

template <class M>
void set_field_selected(Field<M>& f, const SelectedField<M>& sf,
                        const FieldSelection& fsel,
                        const bool is_keeping_data = false);

template <class M>
void acc_field(Field<M>& f, const Complex coef, const SelectedField<M>& sf,
               const FieldSelection& fsel);

template <class M>
std::vector<M> field_sum_tslice(const SelectedField<M>& sf,
                                const FieldSelection& fsel);
// length = t_size * multiplicity

template <class M>
void field_glb_sum_tslice_double(SelectedPoints<M>& sp,
                                 const SelectedField<M>& sf,
                                 const FieldSelection& fsel);

template <class M>
void field_glb_sum_tslice_long(SelectedPoints<M>& sp,
                               const SelectedField<M>& sf,
                               const FieldSelection& fsel);
```

### Fields IO

```c++
inline ShuffledBitSet mk_shuffled_bitset(const FieldSelection& fsel,
                                         const Coordinate& new_size_node);
```

```c++
struct ShuffledFieldsWriter {
  ShuffledFieldsWriter() { init(); }
  ShuffledFieldsWriter(const std::string& path_,
                       const Coordinate& new_size_node_,
                       const bool is_append = false);
  ~ShuffledFieldsWriter() { close(); }
  void init();
  void init(const std::string& path_, const Coordinate& new_size_node_,
            const bool is_append = false);
  void close();
};

struct ShuffledFieldsReader {
  ShuffledFieldsReader();
  ShuffledFieldsReader(const std::string& path_,
                       const Coordinate& new_size_node_ = Coordinate());
  void init();
  void init(const std::string& path_,
            const Coordinate& new_size_node_ = Coordinate());
};
```

```c++
template <class M>
long write(ShuffledFieldsWriter& sfw, const std::string& fn,
           const Field<M>& field);

template <class M>
long write(ShuffledFieldsWriter& sfw, const std::string& fn,
           const SelectedField<M>& sf, const ShuffledBitSet& sbs);

template <class M>
long read(ShuffledFieldsReader& sfr, const std::string& fn, Field<M>& field);

template <class M>
long read(ShuffledFieldsReader& sfr, const std::string& fn, const ShuffledBitSet& sbs, SelectedField<M>& sf);
// sbs must match the actual data
// (code will verify & will fail if not match)

template <class M>
long read_field(Field<M>& field, const std::string& path, const std::string& fn);

template <class M>
long read_field(SelectedField<M>& sf, const std::string& path, const std::string& fn, const ShuffledBitSet& sbs);

template <class M>
long write_float_from_double(ShuffledFieldsWriter& sfw, const std::string& fn,
                             const Field<M>& field);

template <class M>
long write_float_from_double(ShuffledFieldsWriter& sfw, const std::string& fn,
                             const SelectedField<M>& sf, const ShuffledBitSet& sbs);

template <class M>
long read_double_from_float(ShuffledFieldsReader& sfr, const std::string& fn,
                            Field<M>& field);

template <class M>
long read_double_from_float(ShuffledFieldsReader& sfr, const std::string& fn,
                            const ShuffledBitSet& sbs, SelectedField<M>& sf);

template <class M>
long read_field_double_from_float(Field<M>& field, const std::string& path,
                                  const std::string& fn);

template <class M>
long read_field_double_from_float(SelectedField<M>& sf,
                                  const std::string& path,
                                  const std::string& fn,
                                  const ShuffledBitSet& sbs);
```

```c++
inline long flush(ShuffledFieldsWriter& sfw);

inline bool does_file_exist_sync_node(ShuffledFieldsReader& sfr,
                                      const std::string& fn);

inline bool check_file_sync_node(ShuffledFieldsReader& sfr,
                                 const std::string& fn,
                                 std::vector<long>& final_offsets);
// set final_offsets to be the files position after loading the data ``fn'' (zero if failed for that file)
// return if data is loaded successfully

inline std::vector<std::string> list_fields(ShuffledFieldsReader& sfr);

inline std::vector<std::string> properly_truncate_fields_sync_node(
    const std::string& path, const bool is_check_all = false,
    const bool is_only_check = false,
    const Coordinate& new_size_node = Coordinate());
// return available fns
```

### Field shuffle

```c++
template <class M>
void shuffle_field(std::vector<Field<M> >& fs, const Field<M>& f,
                   const Coordinate& new_size_node);

template <class M>
void shuffle_field_back(Field<M>& f, const std::vector<Field<M> >& fs,
                        const Coordinate& new_size_node);
```

```c++
inline ShufflePlan make_shuffle_plan(std::vector<FieldSelection>& fsels,
                                     const FieldSelection& fsel,
                                     const Coordinate& new_size_node);
```

```c++
template <class M>
void shuffle_field(std::vector<Field<M> >& fs, const Field<M>& f,
                   const ShufflePlan& sp);

template <class M>
void shuffle_field_back(Field<M>& f, const std::vector<Field<M> >& fs,
                        const ShufflePlan& sp);

template <class M>
void shuffle_field(std::vector<SelectedField<M> >& fs,
                   const SelectedField<M>& f, const ShufflePlan& sp);

template <class M>
void shuffle_field_back(SelectedField<M>& f,
                        const std::vector<SelectedField<M> >& fs,
                        const ShufflePlan& sp);
```

