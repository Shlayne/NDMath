#pragma once

#include <array>
#include <type_traits>

namespace std::ext
{
	// all_true_v, all_false_v, are_convertible_v
	// https://stackoverflow.com/questions/31767645/how-to-check-that-all-types-in-variadic-template-are-convertible-to-size-t

	namespace impl
	{
		template<bool...>
		struct bools {};

		template<bool B, bool... Bools>
		inline constexpr bool all_same_bool_v = std::is_same_v<bools<B, Bools...>, bools<Bools..., B>>;
	}

	template<bool... Bools>
	inline constexpr bool all_true_v = impl::all_same_bool_v<true, Bools...>;

	template<bool... Bools>
	inline constexpr bool all_false_v = impl::all_same_bool_v<false, Bools...>;

	template<typename To, typename... From>
	inline constexpr bool are_convertible_v = all_true_v<std::is_convertible<From, To>{}...>;

	// Example:
	//template<typename T, typename... Args, std::enable_if_t<std::ext::are_convertible_v<T, Args...>, int> = 0>
	//inline constexpr void Function(Args&&... rrArgs)
	//{
	//
	//}



	// disable_if_t

	namespace impl
	{
		template<bool Condition, typename T = void>
		struct disable_if
		{
			using type = T;
		};

		template<typename T>
		struct disable_if<true, T> {};
	}

	template<bool Condition, typename T = void>
	using disable_if_t = typename impl::disable_if<Condition, T>::type;



	// integer_sequence_to_array, combine_integer_sequences, integer_sequences_to_array,
	// integer_sequence_pop_front, integer_sequence_pop_back

	// Converts an integer sequence to an array of its integers.
	template<typename Int, Int... Ints>
	inline constexpr std::array<Int, sizeof...(Ints)> integer_sequence_to_array(std::integer_sequence<Int, Ints...>)
	{
		return { Ints... };
	}

	// Combines a number of integer sequences.
	template<typename Int, Int... Ints1, Int... Ints2, typename... Args>
	inline constexpr auto combine_integer_sequences(std::integer_sequence<Int, Ints1...>, std::integer_sequence<Int, Ints2...>, Args&&... rrArgs)
	{
		if constexpr (sizeof...(Args) > 0)
			return combine_integer_sequences(std::integer_sequence<Int, Ints1..., Ints2...>(), std::forward<Args>(rrArgs)...);
		else
			return std::integer_sequence<Int, Ints1..., Ints2...>();
	}

	// Converts a number of integer sequences into an std::array of their integers.
	template<typename Int, Int... Ints1, Int... Ints2, typename... Args>
	inline constexpr auto integer_sequences_to_array(std::integer_sequence<Int, Ints1...> sequence1, std::integer_sequence<Int, Ints2...> sequence2, Args&&... rrArgs)
	{
		return integer_sequence_to_array(combine_integer_sequences(sequence1, sequence2, std::forward<Args>(rrArgs)...));
	}

	// Removes the first N integers in the sequence.
	template<uint32_t N, typename Int, Int... Ints, Int R = 0>
	inline constexpr auto integer_sequence_pop_front(std::integer_sequence<Int, R, Ints...>)
	{
		if constexpr (N > 0)
			return integer_sequence_pop_front<N - 1>(std::integer_sequence<Int, Ints...>());
		else
			return std::integer_sequence<Int, R, Ints...>();
	}

	// Removes the first N integers in the sequence.
	template<uint32_t N, typename Int, Int... Ints, Int R = 0>
	inline constexpr auto integer_sequence_pop_back(std::integer_sequence<Int, Ints..., R>)
	{
		if constexpr (N > 0)
			return integer_sequence_pop_back<N - 1>(std::integer_sequence<Int, Ints...>());
		else
			return std::integer_sequence<Int, Ints..., R>();
	}
}
