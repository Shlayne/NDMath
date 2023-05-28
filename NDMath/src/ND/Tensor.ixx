module;

#include <gcem.hpp>

export module nd.tensor;

import nd.dimension;
import nd.scalar;
import nd.impl;

#define _ND ::nd::
#define _IMPL ::nd::impl::

namespace nd::impl
{
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
	inline static constexpr bool gt0_v = gt0<Ns...>::value;
}

export namespace nd
{
	template<Scalar S, Dimension R, Dimension N = 0, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (R > 0 && N > 0 && sizeof...(Ns) == R - 1 && _IMPL gt0_v<Ns...>))
	struct Tensor
	{
	public:
		constexpr Tensor() noexcept = default;
	public:
		constexpr auto operator[](Dimension n) noexcept -> Tensor<S, R - 1, Ns...>&;
		constexpr auto operator[](Dimension n) const noexcept -> const Tensor<S, R - 1, Ns...>&;
	protected:
		Tensor<S, R - 1, Ns...> m_Scalars[N]{Tensor<S, R - 1, Ns...>{}};
	};

	template<Scalar S>
	struct Tensor<S, 0>
	{
	public:
		constexpr Tensor(S scalar = {}) noexcept;
	public:
		constexpr operator S&() noexcept;
		constexpr operator const S&() const noexcept;
	protected:
		S m_Scalar{};
	};

	// External Operators

	// Static Methods

	// Aliases

	//template<Scalar S, Dimension C, Dimension R = C>
	//using Matrix = Tensor<S, 2, C, R>;
}

// Implementation: Don't export.
namespace nd::impl
{

}

export namespace nd
{
	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (R > 0 && N > 0 && sizeof...(Ns) == R - 1 && _IMPL gt0_v<Ns...>))
	constexpr auto Tensor<S, R, N, Ns...>::operator[](Dimension n) noexcept -> Tensor<S, R - 1, Ns...>&
	{
		__assume(n < N);
		return m_Scalars[n];
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (R > 0 && N > 0 && sizeof...(Ns) == R - 1 && _IMPL gt0_v<Ns...>))
	constexpr auto Tensor<S, R, N, Ns...>::operator[](Dimension n) const noexcept -> const Tensor<S, R - 1, Ns...>&
	{
		__assume(n < N);
		return m_Scalars[n];
	}

	// Tensor<0, S>

	template<Scalar S>
	constexpr Tensor<S, 0>::Tensor(S scalar) noexcept
		: m_Scalar{scalar}
	{

	}

	template<Scalar S>
	constexpr Tensor<S, 0>::operator S&() noexcept
	{
		return m_Scalar;
	}

	template<Scalar S>
	constexpr Tensor<S, 0>::operator const S&() const noexcept
	{
		return m_Scalar;
	}
}
