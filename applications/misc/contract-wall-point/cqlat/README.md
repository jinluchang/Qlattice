# README

New functions can be added this directory as ``*.C`` source files. 

Use the ``EXPORT`` macro to define functions that can be used in python.

Example:

```cpp
EXPORT(gf_avg_plaq, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_gf)) {
    return NULL;
  }
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const double ret = gf_avg_plaq(gf);
  return py_convert(ret);
});
```

