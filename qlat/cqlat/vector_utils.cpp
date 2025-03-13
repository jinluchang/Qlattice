#include "lib.h"
#include <qlat/vector_utils/utils_io_vec.h>
#include <qlat/vector_utils/utils_check_fun.h>
#include <qlat/vector_utils/utils_construction.h>

EXPORT(diff_gauge, {
  using namespace qlat;
  PyObject* p0 = NULL;
  PyObject* p1 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p0, &p1)) {
    return NULL;
  }
  GaugeField &g0 = py_convert_type<GaugeField>(p0);
  GaugeField &g1 = py_convert_type<GaugeField>(p1);
  diff_gauge(g0,g1);
  Py_RETURN_NONE;
  ///std::string path;
  ///py_convert(path, p_path);
  ///ld.load(path);
  ///Py_RETURN_NONE;
})

EXPORT(load_gwu_link, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_ld, &p_path)) {
    return NULL;
  }
  GaugeField &g0 = py_convert_type<GaugeField>(p_ld);
  std::string path;
  py_convert(path, p_path);
  load_gwu_link(path,g0);
  Py_RETURN_NONE;
})

EXPORT(save_gwu_prop, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_ld, &p_path)) {
    return NULL;
  }
  Propagator4d& prop = py_convert_type<Propagator4d>(p_ld);
  std::string path;
  py_convert(path, p_path);
  save_gwu_prop(path, prop);
  Py_RETURN_NONE;
})

EXPORT(load_gwu_prop, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_ld, &p_path)) {
    return NULL;
  }
  Propagator4d& prop = py_convert_type<Propagator4d>(p_ld);
  std::string path;
  py_convert(path, p_path);
  load_gwu_prop(path, prop);
  Py_RETURN_NONE;
})

EXPORT(save_gwu_noiP, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_ld, &p_path)) {
    return NULL;
  }
  Propagator4d& prop = py_convert_type<Propagator4d>(p_ld);
  std::string path;
  py_convert(path, p_path);
  save_gwu_noiP(path, prop);
  Py_RETURN_NONE;
})

EXPORT(diff_prop, {
  using namespace qlat;
  PyObject* p_0 = NULL;
  PyObject* p_1 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_0, &p_1)) {
    return NULL;
  }
  Propagator4d& p0 = py_convert_type<Propagator4d>(p_0);
  Propagator4d& p1 = py_convert_type<Propagator4d>(p_1);
  diff_prop(p0,p1);
  Py_RETURN_NONE;

})

EXPORT(random_point_src, {
  using namespace qlat;
  PyObject* p_prop = NULL;
  PyObject* p_seed = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_prop, &p_seed)) {
    return NULL;
  }
  int seed=0;
  py_convert(seed, p_seed);
  Propagator4d& prop = py_convert_type<Propagator4d>(p_prop);

  random_point_src(prop, seed);

  Py_RETURN_NONE;
})

EXPORT(load_gwu_noiP, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_ld, &p_path)) {
    return NULL;
  }
  Propagator4d& prop = py_convert_type<Propagator4d>(p_ld);
  std::string path;
  py_convert(path, p_path);
  load_gwu_noiP(path, prop);
  Py_RETURN_NONE;
})

EXPORT(mk_output, {
  using namespace qlat;
  PyObject* p_0 = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_0)) {
    return NULL;
  }

  std::vector<int > key_T;
  py_convert(key_T, p_0);
  long size = 1;
  for(unsigned int i=0;i<key_T.size();i++)
  {
    int li = key_T[i];
    qassert(li > 0);
    size = size * li;
  }

  //PyObject* ret = PyList_New(size);
  //return ret;

  //std::vector<double >& write = set_obj<std::vector<double > >();
  //write.resize(size);

  ///std::vector<double > write;
  ///write.resize(size);
  std::vector<double >* write = new std::vector<double >(size);
  return py_convert((void*)write);
  //print0("size of write %d \n", int(size));

  //double* write = new double[size];
  //return py_convert((void*)write);

})

EXPORT(free_output, {
  using namespace qlat;
  PyObject* p_0 = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_0)) {
    return NULL;
  }
  std::vector<double >* write = (std::vector<double >*) PyLong_AsVoidPtr(p_0);
  delete write;
  Py_RETURN_NONE;
})

EXPORT(clear_output, {
  using namespace qlat;
  PyObject* p_0 = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_0)) {
    return NULL;
  }

  std::vector<double >& write = *((std::vector<double >*) PyLong_AsVoidPtr(p_0));
  print0("size of write %d \n", int(write.size()));
  for(unsigned int i=0;i<write.size();i++){write[i] = 0;}

  Py_RETURN_NONE;
})

EXPORT(write_output, {
  using namespace qlat;
  PyObject* p_0 = NULL;
  PyObject* p_1 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_0, &p_1)) {
    return NULL;
  }

  std::vector<double >& write = *((std::vector<double >*) PyLong_AsVoidPtr(p_0));
  /////std::vector<double >& write = py_convert_type<std::vector<double > >(p_0);

  std::string output;
  py_convert(output , p_1);

  write_data(write,output.c_str());

  Py_RETURN_NONE;

})

EXPORT(make_point_prop, {
  using namespace qlat;
  PyObject* p_v0 = NULL;
  PyObject* p_v1 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_v0, &p_v1)) {
    return NULL;
  }
  Coordinate sp;
  py_convert(sp, p_v1);

  Propagator4d& prop = py_convert_type<Propagator4d>(p_v0);
  make_point_prop(prop, sp);
  Py_RETURN_NONE;
})

EXPORT(make_volume_src, {
  using namespace qlat;
  PyObject* p_v0 = NULL;
  ///PyObject* p_v1 = NULL;
  //PyObject* p_v2 = NULL;
  //PyObject* p_v3 = NULL;

  int seed      =  0;
  int mix_color =  0;
  int mix_spin  =  0;
  int tini      = -1;
  if (!PyArg_ParseTuple(args, "Oi|iii", &p_v0, &seed, &mix_color, &mix_spin, &tini)) {
    return NULL;
  }

  //py_convert(seed, p_v1);
  //if(p_v2 != NULL){py_convert(mix_color, p_v2);}
  //if(p_v3 != NULL){py_convert(mix_spin , p_v3);}

  Propagator4d& prop = py_convert_type<Propagator4d>(p_v0);
  make_volume_src(prop, seed, mix_color, mix_spin, tini);
  Py_RETURN_NONE;
})

EXPORT(local_sequential_source, {
  using namespace qlat;
  PyObject* p_v0 = NULL;
  PyObject* p_v1 = NULL;
  PyObject* p_v2 = NULL;
  //PyObject* p_v3 = NULL;
  int gammai = -1;
  if (!PyArg_ParseTuple(args, "OOO|i", &p_v0, &p_v1, &p_v2, &gammai)) {
    return NULL;
  }
  //if(p_v3 != NULL){py_convert(gammai, p_v3);}
  Propagator4d& res = py_convert_type<Propagator4d>(p_v0);
  Propagator4d& src = py_convert_type<Propagator4d>(p_v1);

  const std::vector<int> tseqL = py_convert_data<std::vector<int> >(p_v2);
  qlat::vector_acc<int > tseq;tseq.resize(tseqL.size());
  for(long i=0;i<tseq.size();i++){tseq[i] = tseqL[i];}

  local_sequential_source(res, src, tseq, gammai);
  Py_RETURN_NONE;
})

EXPORT(meson_corr, {
  using namespace qlat;
  PyObject* p_v0 = NULL;
  PyObject* p_v1 = NULL;
  PyObject* p_v2 = NULL;
  PyObject* p_v4 = NULL;
  PyObject* p_v3 = NULL;

  int g0      =  0;
  int g1      =  0;
  int tini    =  0;
  int invmode =  1;
  int shift_end = 1;
  if (!PyArg_ParseTuple(args, "OOOii|iiOiO", &p_v0, &p_v1, &p_v2, &g0, &g1, &tini, &invmode, &p_v4, &shift_end, &p_v3)) {
    return NULL;
  }

  Propagator4d& p0 = py_convert_type<Propagator4d>(p_v0);
  Propagator4d& p1 = py_convert_type<Propagator4d>(p_v1);
  std::string filename;
  std::string info = std::string("NONE");
  py_convert(filename, p_v2);

  Coordinate mom = Coordinate(0, 0, 0, 0);
  if(p_v3 != NULL){py_convert(mom , p_v3);}
  if(p_v4 != NULL){py_convert(info, p_v4);}

  meson_corrE(p0, p1, g0, g1, filename, mom, invmode, tini, info, shift_end);
  Py_RETURN_NONE;
})

EXPORT(prop_corr, {
  using namespace qlat;
  PyObject* p_v0 = NULL;
  PyObject* p_v2 = NULL;
  PyObject* p_v4 = NULL;
  PyObject* p_v3 = NULL;

  int tini    =  0;
  int shift_end = 1;
  if (!PyArg_ParseTuple(args, "OO|iOiO", &p_v0, &p_v2, &tini, &p_v4, &shift_end, &p_v3)) {
    return NULL;
  }

  Propagator4d& p0 = py_convert_type<Propagator4d>(p_v0);
  std::string filename;
  std::string info = std::string("NONE");
  py_convert(filename, p_v2);

  Coordinate mom = Coordinate(0, 0, 0, 0);
  if(p_v3 != NULL){py_convert(mom , p_v3);}
  if(p_v4 != NULL){py_convert(info, p_v4);}

  prop_corrE(p0, filename, mom, tini, info, shift_end);
  Py_RETURN_NONE;
})

EXPORT(corr_dat_create, {
  using namespace qlat;
  PyObject* p_v0 = NULL;
  PyObject* p_v1 = NULL;
  PyObject* p_v2 = NULL;
  PyObject* p_v3 = NULL;

  if (!PyArg_ParseTuple(args, "OOO|O", &p_v0, &p_v1, &p_v2, &p_v3)) {
    return NULL;
  }

  std::string filename;
  std::string key_T;
  std::string dimN ;
  std::string info = std::string("NONE");

  py_convert(filename, p_v0);
  py_convert(key_T   , p_v1);
  py_convert(dimN    , p_v2);
  if(p_v3 != NULL){py_convert(info, p_v3);}

  corr_dat_create(filename, key_T, dimN, info);
  Py_RETURN_NONE;
})

EXPORT(corr_dat_info, {
  using namespace qlat;
  PyObject* p_v0 = NULL;
  PyObject* p_v1 = NULL;

  if (!PyArg_ParseTuple(args, "OO", &p_v0, &p_v1)) {
    return NULL;
  }

  std::string filename;
  std::string info = std::string("NONE");

  py_convert(filename, p_v0);
  py_convert(info    , p_v1);

  corr_dat_info(filename, info);
  Py_RETURN_NONE;
})

EXPORT(prop4d_conj, {
  using namespace qlat;
  PyObject* p_v0 = NULL;
  int rotate = 1;
  if (!PyArg_ParseTuple(args, "O|i", &p_v0, &rotate)) {
    return NULL;
  }
  Propagator4d& src = py_convert_type<Propagator4d>(p_v0);
  prop4d_conj(src, rotate);
  Py_RETURN_NONE;
})

EXPORT(prop4d_src_gamma, {
  using namespace qlat;
  PyObject* p_v0 = NULL;
  int g0      =  0;
  int Conj    =  0;
  if (!PyArg_ParseTuple(args, "Oi|i", &p_v0, &g0, &Conj)) {
    return NULL;
  }
  Propagator4d& src = py_convert_type<Propagator4d>(p_v0);

  ga_matrices_cps ga_cps;
  std::vector<ga_M > gL;gL.resize(16);
  {int o=0;
  for(int i=0;i<6;i++){gL[o] = ga_cps.ga[0][i];o+=1;}
  for(int i=2;i<6;i++){gL[o] = ga_cps.ga[1][i];o+=1;}
  for(int i=3;i<6;i++){gL[o] = ga_cps.ga[2][i];o+=1;}
  for(int i=4;i<6;i++){gL[o] = ga_cps.ga[3][i];o+=1;}
  for(int i=5;i<6;i++){gL[o] = ga_cps.ga[4][i];o+=1;}}

  prop4d_src_gamma(src, gL[g0], Conj);
  Py_RETURN_NONE;
})

EXPORT(prop4d_sink_gamma, {
  using namespace qlat;
  PyObject* p_v0 = NULL;
  int g0      =  0;
  int Conj    =  0;
  if (!PyArg_ParseTuple(args, "Oi|i", &p_v0, &g0, &Conj)) {
    return NULL;
  }
  Propagator4d& src = py_convert_type<Propagator4d>(p_v0);

  ga_matrices_cps ga_cps;
  std::vector<ga_M > gL;gL.resize(16);
  {int o=0;
  for(int i=0;i<6;i++){gL[o] = ga_cps.ga[0][i];o+=1;}
  for(int i=2;i<6;i++){gL[o] = ga_cps.ga[1][i];o+=1;}
  for(int i=3;i<6;i++){gL[o] = ga_cps.ga[2][i];o+=1;}
  for(int i=4;i<6;i++){gL[o] = ga_cps.ga[3][i];o+=1;}
  for(int i=5;i<6;i++){gL[o] = ga_cps.ga[4][i];o+=1;}}

  prop4d_sink_gamma(src, gL[g0], Conj);
  Py_RETURN_NONE;
})

