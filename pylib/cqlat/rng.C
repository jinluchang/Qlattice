#include "convert.h"
#include "dispatch.h"
#include "exceptions.h"

EXPORT(mk_rng, {
  using namespace qlat;
  RngState* prng = NULL;
  PyObject* p_seed = NULL;
  if (!PyArg_ParseTuple(args, "|lO", &prng, &p_seed)) {
    return NULL;
  }
  std::string seed;
  if (p_seed != NULL) {
    py_convert(seed, p_seed);
  }
  RngState* prng_new = NULL;
  if (NULL == prng) {
    prng_new = new RngState();
  } else {
    prng_new = new RngState(*prng, seed);
  }
  return PyLong_FromVoidPtr(prng_new);
});

EXPORT(free_rng, {
  using namespace qlat;
  RngState* prng = NULL;
  if (!PyArg_ParseTuple(args, "l", &prng)) {
    return NULL;
  }
  pqassert(prng != NULL);
  delete prng;
  Py_RETURN_NONE;
});
