#pragma once

#include <array>
#include <type_traits>

namespace std::ext
{
	// integer_sequence_to_array, combine_integer_sequences, integer_sequences_to_array,
	// integer_sequence_pop_front, integer_sequence_pop_back

	// Converts an integer sequence to an array of its integers.
	template<typename Int, Int... Ints>
	inline constexpr auto integer_sequence_to_array(std::integer_sequence<Int, Ints...>) -> std::array<Int, sizeof...(Ints)>
	{
		return { Ints... };
	}

	// Combines a number of integer sequences.
	template<typename Int, Int... Ints1, Int... Ints2, typename... Sequences>
	inline constexpr auto combine_integer_sequences(std::integer_sequence<Int, Ints1...>, std::integer_sequence<Int, Ints2...>, Sequences&&... sequences)
	{
		if constexpr (sizeof...(Args) > 0)
			return combine_integer_sequences(std::integer_sequence<Int, Ints1..., Ints2...>(), std::forward<Sequences>(sequences)...);
		else
			return std::integer_sequence<Int, Ints1..., Ints2...>();
	}

	// Converts a number of integer sequences into an std::array of their integers.
	template<typename Int, Int... Ints1, Int... Ints2, typename... Sequences>
	inline constexpr auto integer_sequences_to_array(std::integer_sequence<Int, Ints1...> sequence1, std::integer_sequence<Int, Ints2...> sequence2, Sequences&&... sequences)
	{
		return integer_sequence_to_array(combine_integer_sequences(sequence1, sequence2, std::forward<Sequences>(sequences)...));
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

	// Removes the last N integers in the sequence.
	template<uint32_t N, typename Int, Int... Ints, Int R = 0>
	inline constexpr auto integer_sequence_pop_back(std::integer_sequence<Int, Ints..., R>)
	{
		if constexpr (N > 0)
			return integer_sequence_pop_back<N - 1>(std::integer_sequence<Int, Ints...>());
		else
			return std::integer_sequence<Int, Ints..., R>();
	}
}
