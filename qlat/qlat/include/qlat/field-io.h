#pragma once

#include <qlat/field-base-io.h>
#include <qlat/field-dist-io.h>
#include <qlat/field-serial-io.h>

namespace qlat
{  //

// --------------------

#ifdef QLAT_INSTANTIATE_FIELD_IO
#define QLAT_EXTERN
#else
#define QLAT_EXTERN extern
#endif

#define QLAT_EXTERN_TEMPLATE(TYPENAME)                   \
                                                         \
  QLAT_EXTERN template Long write_field<TYPENAME>(       \
      const Field<TYPENAME>& f, const std::string& path, \
      const Coordinate& new_size_node);                  \
                                                         \
  QLAT_EXTERN template Long read_field<TYPENAME>(        \
      Field<TYPENAME> & f, const std::string& path,      \
      const Coordinate& new_size_node_)

QLAT_CALL_WITH_TYPES(QLAT_EXTERN_TEMPLATE);
#undef QLAT_EXTERN_TEMPLATE

#undef QLAT_EXTERN

}  // namespace qlat
