#pragma once

#include <concepts>

namespace nd
{
	using Dimension = size_t;

	template <typename T>
	concept Scalar = _STD integral<T> || _STD floating_point<T>;
}
