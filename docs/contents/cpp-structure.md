# Structure of C++ code

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
  Coordinate node_site;
  Coordinate expansion_left;
  Coordinate expansion_right;
  Coordinate node_site_expanded;
};
```

### vector and box

```c++
template <class M>
struct vector {
  bool is_copy;  // do not free memory if is_copy=true
  bool is_acc; // if place data on cudaMallocManaged memory (default false)
  Vector<M> v;
};
template <class M>
struct vector_acc : vector<M> {
  // default is_acc = true
};
```

```c++
template <class M>
struct box {
  bool is_copy;  // do not free memory if is_copy=true
  bool is_acc; // if place data on cudaMallocManaged memory (default false)
  Handle<M> v;
};
template <class M>
struct box_acc : box<M> {
  // default is_acc = true
};
```

```c++
template <class M>
struct API Vector {
  M* p;
  Long n;
};
```

```c++
template <class M>
struct API Handle {
	M* p;
};
```

### Field

```c++
template <class M>
struct Field {
  bool initialized;
  int multiplicity;
  box_acc<Geometry> geo;
  vector_acc<M> field;
};
```

### FieldM

```c++
template <class M, int multiplicity>
struct FieldM : Field<M> {};
```

### SelectedPoints

```c++
enum struct PointsDistType {
  Global,
  Local,
  Random, // shuffle based on coordinate
};
```

```c++
struct API PointsSelection {
  bool initialized;
  PointsDistType points_dist_type;  // default PointsDistType::Global (all node has the same data)
  Coordinate total_site;
  vector_acc<Coordinate> xgs;
};
```

```c++
template <class M>
struct SelectedPoints {
  bool initialized;
  PointsDistType points_dist_type;  // default PointsDistType::Global (all node has the same data)
  int multiplicity;
  long n_points;
  vector_acc<M> points;
  // points.size() == n_points * multiplicity if initialized = true
};
```

### SelectedField

```c++
template <class M>
struct SelectedField {
  bool initialized;
  long n_elems;
  int multiplicity;
  box_acc<Geometry> geo;
  vector_acc<M> field;
};
```
