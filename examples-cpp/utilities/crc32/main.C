#include <qlat/qlat.h>

#include <complex>
#include <iostream>
#include <vector>

using namespace qlat;
using namespace std;

int main(int argc, char* argv[])
{
  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 2));
  size_node_list.push_back(Coordinate(1, 1, 1, 4));
  size_node_list.push_back(Coordinate(1, 1, 1, 8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  begin(&argc, &argv, size_node_list);
  if (argc < 2) {
    displayln_info("usage: ./crc32 location1 location2 ...");
    exit(-1);
  }
  std::vector<std::string> fns(argc - 1);
  for (long i = 0; i < (long)fns.size(); ++i) {
    fns[i] = remove_trailing_slashes(argv[1 + i]);
    displayln_info(
        ssprintf("fns[%5d/%d] = '%s'", i, fns.size(), fns[i].c_str()));
  }
  displayln_info("Start to calculating crc32...");
  for (long i = 0; i < (long)fns.size(); ++i) {
    displayln_info(
        ssprintf("fns[%5d/%d] = '%s'", i, fns.size(), fns[i].c_str()));
    if (is_d_field(fns[i])) {
      Field<char> f;
      const long total_bytes = read_field(f, fns[i]);
      if (total_bytes == 0 or is_checksum_missmatch()) {
        displayln_info(ssprintf("FAILED: fns[%5d/%d]\nFAILED FILENAME: %s", i,
                                fns.size(), fns[i].c_str()));
      }
    } else if (does_file_exist_sync_node(fns[i] + "/metadata.txt")) {
      if (check_compressed_eigen_vectors(fns[i])) {
        displayln_info(ssprintf("FAILED: fns[%5d/%d]\nFAILED FILENAME: %s", i,
                                fns.size(), fns[i].c_str()));
      }
    } else if (does_file_exist_sync_node(fns[i] + "/checkpoint")) {
      decompressed_eigen_vectors_check_crc32(fns[i]);
    } else if (does_file_exist_sync_node(fns[i])) {
      if (get_id_node() == 0) {
        const crc32_t crc = compute_crc32(fns[i]);
        displayln_info(ssprintf("crc32 = %08X ; fns[%5d/%d] = '%s'", crc, i,
                                fns.size(), fns[i].c_str()));
      }
    } else {
      displayln_info("Cannot calculate crc32 for this data: '" + fns[i] + "'.");
    }
  }
  Timer::display();
  end();
  return 0;
}
