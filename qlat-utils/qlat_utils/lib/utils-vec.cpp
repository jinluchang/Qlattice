#include <qlat-utils/utils-vec.h>

namespace qlat
{  //

void convert_double_from_float(Vector<RealD> vd, const Vector<RealF> vf)
{
  qassert(vd.size() == vf.size());
  qacc_for(i, vd.size(), { vd[i] = vf[i]; });
}

void convert_float_from_double(Vector<RealF> vf, const Vector<RealD> vd)
{
  qassert(vf.size() == vd.size());
  qacc_for(i, vf.size(), { vf[i] = vd[i]; });
}

}  // namespace qlat
