#include <qlat-utils/qassert.h>
#include <qlat-utils/qacc.h>
#include <qlat-utils/timer.h>

namespace qlat
{  //

std::string show(const MemType mem_type)
{
  if (mem_type == MemType::Cpu) {
    return "cpu";
  } else if (mem_type == MemType::Acc) {
    return "acc";
  } else if (mem_type == MemType::Uvm) {
    return "uvm";
  } else if (mem_type == MemType::Comm) {
    return "comm";
  } else if (mem_type == MemType::CommAcc) {
    return "comm_acc";
  } else {
    qassert(false);
    return "";
  }
}

MemType read_mem_type(const std::string& mem_type_str)
{
  if (mem_type_str == "cpu") {
    return MemType::Cpu;
  } else if (mem_type_str == "acc") {
    return MemType::Acc;
  } else if (mem_type_str == "uvm") {
    return MemType::Uvm;
  } else if (mem_type_str == "comm") {
    return MemType::Comm;
  } else if (mem_type_str == "comm_acc") {
    return MemType::CommAcc;
  } else {
    qassert(false);
    return MemType::Cpu;
  }
}

}  // namespace qlat
