export module nd.impl;

import nd.scalar;
import <concepts>;

export namespace nd::impl
{
	// https://stackoverflow.com/questions/31767645/how-to-check-that-all-types-in-variadic-template-are-convertible-to-size-t

	template<bool...>
	struct bools {};

	template<bool B, bool... Bools>
	constexpr bool all_same_bool_v = _STD is_same_v<bools<B, Bools...>, bools<Bools..., B>>;

	template<bool... Bools>
	constexpr bool all_true_v = all_same_bool_v<true, Bools...>;


	template<Scalar... Scalars>
	using CT = _STD common_type_t<Scalars...>;
}
