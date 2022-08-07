#pragma once

#include "ND/Vector.h"
#include <vector>

namespace nd
{
	using Index = uint32_t;

	template<Index... Inds>
	using Indices = std::integer_sequence<Index, Inds...>;

	template<typename T, uint32_t N>
	struct Vertex
	{
		using Position = Vector<T, N>;

		constexpr Vertex() noexcept = default;
		constexpr Vertex(const Position& crPosition) noexcept;

		Position position;
	};

	template<typename T, uint32_t N>
	struct Mesh
	{
		std::vector<Vertex<T, N>> vertices;
		std::vector<Index> indices;
	};
}

#include "ND/Mesh.inl"
