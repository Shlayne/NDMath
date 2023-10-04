#pragma once

#include "Types.hpp"
#include <concepts>

namespace nd::impl
{
	// Simplify things with fold expressions: https://en.cppreference.com/w/cpp/language/fold
	// This is here so if I forget what it's called, I can get back to it.
	
	// Check if all dimensions are greater than some dimension K.

	template <Dimension K, Dimension... Ns>
	constexpr bool gt_v = (... && (Ns > K));


	// The common type of a list of scalars.

	template <Scalar S, Scalar... Scalars>
	constexpr bool all_convertible_to_v = (... && _STD is_convertible_v<Scalars, S>);

	template <Scalar... Scalars>
	requires(all_convertible_to_v<_STD common_type_t<Scalars...>, Scalars...>)
	using CT = _STD common_type_t<Scalars...>;


	// Generic function concept

	template <typename Func, typename ReturnType, typename... Args>
	concept Function = requires(const Func& func, Args&&... args)
	{
		{func(_STD forward<Args>(args)...)} -> _STD same_as<ReturnType>;
	};
}
