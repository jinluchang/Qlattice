#include <qlat/fields-io.h>

namespace qlat
{  //

void close_all_shuffled_fields_writer()
// Force close all the ShuffledFieldsWriter.
// Only call this when quitting the program (e.g. in qquit(msg)).
{
  TIMER_VERBOSE("close_all_shuffled_fields_writer");
  ShuffledFieldsWriterMap& sfwm = get_all_shuffled_fields_writer();
  std::vector<Handle<ShuffledFieldsWriter> > sfwv;
  for (auto it = sfwm.begin(); it != sfwm.end(); ++it) {
    sfwv.push_back(it->second);
  }
  for (long i = 0; i < (long)sfwv.size(); ++i) {
    sfwv[i]().close();
  }
  qassert(sfwm.size() == 0);
}

}  // namespace qlat
