export module nd.impl;

import nd.types;
import <concepts>;

export namespace nd::impl
{
	// Check if a list of boolean expressions all evaluate to the same value.
	// https://stackoverflow.com/questions/31767645/how-to-check-that-all-types-in-variadic-template-are-convertible-to-size-t

	template<bool...>
	struct bools {};

	template<bool B, bool... Bools>
	inline constexpr bool all_same_bool_v = _STD is_same_v<bools<B, Bools...>, bools<Bools..., B>>;

	template<bool... Bools>
	inline constexpr bool all_true_v = all_same_bool_v<true, Bools...>;

	
	// Check if all dimensions are positive numbers.

	template<Dimension N = 1, Dimension... Ns>
	struct gt0
	{
		inline static constexpr bool value{N > 0 && gt0<Ns...>::value};
	};

	template<>
	struct gt0<1>
	{
		inline static constexpr bool value{true};
	};

	template<Dimension... Ns>
	inline constexpr bool gt0_v = gt0<Ns...>::value;


	// The common type of a list of scalars.

	template<Scalar... Scalars>
	requires(all_true_v<_STD is_convertible_v<Scalars, _STD common_type_t<Scalars...>>...>)
	using CT = _STD common_type_t<Scalars...>;
}
