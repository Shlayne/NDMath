module;

#include <gcem.hpp>

export module nd.tensor;

import nd.dimension;
import nd.scalar;

#define _ND ::nd::
#define _IMPL ::nd::impl::

namespace nd
{
	template<Scalar S, Dimension R, Dimension N = 0, Dimension... Ns>
	requires((!R && !N && !sizeof...(Ns)) || (R > 0 && N > 0 && sizeof...(Ns) == R - 1))
	struct Tensor;
}

template<_ND Scalar... Scalars>
using CT = _STD common_type_t<Scalars...>;

export namespace nd
{
	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((!R && !N && !sizeof...(Ns)) || (R > 0 && N > 0 && sizeof...(Ns) == R - 1))
	struct Tensor
	{
	public:
		constexpr Tensor() noexcept = default;
	public:
		constexpr auto operator[](Dimension n) noexcept -> Tensor<S, R - 1, Ns...>&;
		constexpr auto operator[](Dimension n) const noexcept -> const Tensor<S, R - 1, Ns...>&;
	private:
		Tensor<S, R - 1, Ns...> m_Scalars[N]{Tensor<S, R - 1, Ns...>{}};
	};

	template<Scalar S>
	struct Tensor<S, 0>
	{
	public:
		constexpr Tensor() noexcept = default;
	public:
		constexpr operator S() const noexcept;
	private:
		S m_Scalar{};
	};

	// External Operators

	// Static Methods

	// Aliases
	template<Scalar S, Dimension N>
	using Vector = Tensor<S, 1, N>;

	template<Dimension N>
	using VectorNf = Vector<float, N>;

	using Vector1f = VectorNf<1>;
	using Vector2f = VectorNf<2>;
	using Vector3f = VectorNf<3>;
	using Vector4f = VectorNf<4>;
	using Vector5f = VectorNf<5>;
	using Vector6f = VectorNf<6>;
	using Vector7f = VectorNf<7>;
	using Vector8f = VectorNf<8>;

	template<Scalar S, Dimension C, Dimension R = C>
	using Matrix = Tensor<S, 2, C, R>;

	template<Dimension C, Dimension R = C>
	using MatrixCRf = Matrix<float, C, R>;

	using Matrix1f = MatrixCRf<1>;
	using Matrix2f = MatrixCRf<2>;
	using Matrix3f = MatrixCRf<3>;
	using Matrix4f = MatrixCRf<4>;
	using Matrix5f = MatrixCRf<5>;
	using Matrix6f = MatrixCRf<6>;
	using Matrix7f = MatrixCRf<7>;
	using Matrix8f = MatrixCRf<8>;
}

// Implementation: Don't export.
namespace nd::impl
{

}

export namespace nd
{
	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((!R && !N && !sizeof...(Ns)) || (R > 0 && N > 0 && sizeof...(Ns) == R - 1))
	constexpr auto Tensor<S, R, N, Ns...>::operator[](Dimension n) noexcept -> Tensor<S, R - 1, Ns...>&
	{
		__assume(n < N);
		return m_Scalars[n];
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((!R && !N && !sizeof...(Ns)) || (R > 0 && N > 0 && sizeof...(Ns) == R - 1))
	constexpr auto Tensor<S, R, N, Ns...>::operator[](Dimension n) const noexcept -> const Tensor<S, R - 1, Ns...>&
	{
		__assume(n < N);
		return m_Scalars[n];
	}

	// Tensor<0, S>

	template<Scalar S>
	constexpr Tensor<S, 0>::operator S() const noexcept
	{
		return m_Scalar;
	}
}
