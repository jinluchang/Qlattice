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
```

### Qnorm

```c++
template <class M>
double qnorm(const Field<M>& f);
```

```c++
template <class M>
double qnorm_double(const Field<M>& f1, const Field<M>& f2);
```

### Qswap

```c++
template <class M>
void qswap(Field<M>& f1, Field<M>& f2);
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

