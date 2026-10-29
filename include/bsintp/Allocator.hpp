#ifndef INTP_ALLOCATOR
#define INTP_ALLOCATOR

#ifndef INTP_CUSTOM_ALLOCATOR
#include "aligned-allocator.hpp"

namespace intp {

// Construction buffers can become retained control points or coefficients.
template <typename T>
using BufferAllocator = AlignedAllocator<T, Alignment::AVX>;

}  // namespace intp
#endif  // INTP_CUSTOM_ALLOCATOR

#endif  // INTP_ALLOCATOR
