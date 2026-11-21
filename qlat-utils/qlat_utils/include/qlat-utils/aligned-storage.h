#pragma once

#include <cstddef>
#include <cstdint>
#include <new>
#include <type_traits>

namespace qlat
{

// ------------------------------------------------------------

// QlatAlignedStorage provides correctly aligned storage for one `T` inside a
// containing object.
//
// The allocator of the containing object only guarantees the fundamental
// alignment (alignof(std::max_align_t), 16 bytes on x86-64; this is also the
// guarantee CPython gives for the memory of an object).  A type that requires
// a larger alignment (e.g. qlat::MatrixT can require up to 256 bytes) must
// therefore be placed by rounding the address of an over-sized,
// max_align-aligned buffer up to the required alignment.  `data` is
// sizeof(T) + (alignment() - alignof(std::max_align_t)) bytes, so the rounded
// up `ptr()` (at most alignment() - 1 - alignof(std::max_align_t) bytes past
// the start) always leaves sizeof(T) bytes inside the buffer.
//
// The storage is raw memory.  Call construct() to placement-new the T (the
// containing Cython object memory is zeroed by the allocator and the default
// constructors of the matrix types write nothing, so this keeps the previous
// zero-initialised semantics), and destroy() before the containing object
// goes away.  T must be trivially destructible.

template <class T>
struct QlatAlignedStorage {
  static constexpr std::size_t max_align = alignof(std::max_align_t);
  static constexpr std::size_t align = alignof(T);
  static constexpr std::size_t slack =
      align > max_align ? align - max_align : 0;
  static constexpr std::size_t size = sizeof(T) + slack;
  //
  static_assert(align <= max_align || align % max_align == 0,
                "alignment must be a multiple of alignof(std::max_align_t)");
  static_assert(size >= sizeof(T), "storage must be large enough for T");
  static_assert(std::is_trivially_destructible<T>::value,
                "T must be trivially destructible to live in raw storage");
  //
  alignas(max_align) char data[size];
  //
  char* base() noexcept { return data; }
  //
  const char* base() const noexcept { return data; }
  //
  T* ptr() noexcept
  {
    const std::uintptr_t addr = reinterpret_cast<std::uintptr_t>(data);
    return reinterpret_cast<T*>((addr + align - 1) & ~(align - 1));
  }
  //
  const T* ptr() const noexcept
  {
    const std::uintptr_t addr = reinterpret_cast<std::uintptr_t>(data);
    return reinterpret_cast<const T*>((addr + align - 1) & ~(align - 1));
  }
  //
  T* construct() { return new (ptr()) T(); }
  //
  void destroy() { ptr()->~T(); }
  //
  std::size_t alignment() const noexcept { return align; }
  //
  std::size_t capacity() const noexcept { return size; }
};

// ------------------------------------------------------------

}  // namespace qlat
