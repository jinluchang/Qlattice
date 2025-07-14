#include <qlat-utils/utils.h>

namespace qlat
{  //

Vector<char> get_data(const std::string& str)
{
  return Vector<char>(&str[0], str.length());
}

}  // namespace qlat
