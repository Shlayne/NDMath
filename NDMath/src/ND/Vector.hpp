#pragma once

#include "Tensor.hpp"

#define _ND ::nd::
#define _IMPL ::nd::impl::
#define _GCEM ::gcem::

namespace nd
{
	template <Scalar S, Dimension N>
	struct Vector : public Tensor<S, N>
	{
	public:
		template <typename... Args>
		requires(sizeof...(Args) > 1)
		constexpr Vector(Args&&... args) noexcept
		{
			Fill<N>(_STD forward<Args>(args)...);
		}

		template <Scalar S2, Dimension N2>
		constexpr Vector(const Tensor<S2, N2>& tensor) noexcept
			: Tensor<S, N>{tensor}
		{

		}
	public:
		using Tensor<S, N>::Tensor;
		using Tensor<S, N>::operator=;
		using Tensor<S, N>::operator+; // also handles unary plus.
		using Tensor<S, N>::operator+=;
		using Tensor<S, N>::operator-; // also handles unary minus.
		using Tensor<S, N>::operator-=;
		using Tensor<S, N>::operator*;
		using Tensor<S, N>::operator*=;
		using Tensor<S, N>::operator/;
		using Tensor<S, N>::operator/=;
		using Tensor<S, N>::at;
		using Tensor<S, N>::operator[];
		using Tensor<S, N>::operator==;
		using Tensor<S, N>::operator!=;
	public:
		// NOTE: using a base class' templated conversion operator is expliticly prohibited by C++23.
		template <Scalar S2>
		requires(!_STD is_same_v<S, S2> && _STD is_convertible_v<S, S2>)
		constexpr operator Vector<S2, N>() const noexcept
		{
			Vector<S2, N> result;
			for (Dimension n{}; n < N; ++n)
				result[n] = static_cast<S2>(m_Scalars[n]);
			return result;
		}
	private:
		template <Dimension N2, Scalar S2, typename... Args>
		requires(N2 >= 1 && _STD is_convertible_v<S2, S>)
		constexpr void Fill(S2 scalar, Args&&... args) noexcept
		{
			m_Scalars[N - N2] = static_cast<S>(scalar);
			if constexpr (sizeof...(Args) > 0)
				Fill<N2 - 1>(_STD forward<Args>(args)...);
		}

		template <Dimension N2, Scalar S2, Dimension N3, typename... Args>
		requires(N2 >= N3 && _STD is_convertible_v<S2, S>)
		constexpr void Fill(const Vector<S2, N3>& vector, Args&&... args) noexcept
		{
			for (Dimension n{}; n < N3; ++n)
				m_Scalars[N - N2 + n] = static_cast<S>(vector[n]);
			if constexpr (sizeof...(Args) > 0)
				Fill<N2 - N3>(_STD forward<Args>(args)...);
		}
	protected:
		using Tensor<S, N>::m_Scalars;
	};

	// Vector of length 1.
	// When passed to functions that require a Normal,
	// this will not be normalized, however, Vectors
	// will be. If you already have a Normal, pass it
	// as usual, and if you have an unnormalized Vector,
	// also pass it as usual, and it will be normalized
	// automatically.
	// This struct shall be immutable to ensure it is
	// always normalized. If you want to change it,
	// cast it back to a Vector (unary +/- also work).
	template <Scalar S, Dimension N>
	struct Normal : public Vector<S, N>
	{
	public:
		// No default constructor since that would be
		// all zeros, and hence violate the definition.
		constexpr Normal() noexcept = delete;

		constexpr Normal(const Normal&) noexcept = default;
		constexpr Normal& operator=(const Normal&) noexcept = default;

		constexpr Normal(const Tensor<S, N>& vector);
		constexpr Normal& operator=(const Tensor<S, N>& vector);
	public:
		using Vector<S, N>::operator+; // also handles unary plus.
		using Vector<S, N>::operator-; // also handles unary minus.
		using Vector<S, N>::operator*;
		using Vector<S, N>::operator/;
		using Vector<S, N>::at;
		using Vector<S, N>::operator[];
		using Vector<S, N>::operator==;
		using Vector<S, N>::operator!=;

		// These could also violate the definition, so they too are not allowed.
	public:
		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N>& operator+=(S2) = delete;

		template <Scalar S2, Dimension N2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N>& operator+=(const Tensor<S2, N2>&) = delete;

		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N>& operator-=(S2) = delete;

		template <Scalar S2, Dimension N2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N>& operator-=(const Tensor<S2, N2>&) = delete;

		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N>& operator*=(S2) = delete;

		template <Scalar S2, Dimension N2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N>& operator*=(const Tensor<S2, N2>&) = delete;

		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N>& operator/=(S2) = delete;

		template <Scalar S2, Dimension N2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N>& operator/=(const Tensor<S2, N2>&) = delete;
		
		template <Dimension N2>
		requires(N2 < N)
		constexpr S& at() noexcept = delete;

		constexpr S& operator[](Dimension) noexcept = delete;
	protected:
		using Vector<S, N>::m_Scalars;
	};

	// Serialization.

	template <Scalar S, Dimension N>
	_STD ostream& operator<<(_STD ostream& ostream, const Tensor<S, N>& vector)
	{
		ostream << '<' << vector[0];
		for (Dimension n{1}; n < N; ++n)
			ostream << ',' << vector[n];
		return ostream << '>';
	}

	// Floating Methods

	namespace impl
	{
		template <Scalar S>
		constexpr S Lerp(S a, S b, S t)
		{
			return (S{1} - t) * a + t * b;
		};

		template <Scalar S, Dimension N, Function<S, S> Func>
		constexpr Vector<S, N> Apply(const Func& func, const Tensor<S, N>& vector)
		{
			Vector<S, N> result;
			for (Dimension n{}; n < N; ++n)
				result[n] = func(vector[n]);
			return result;
		}

		template <Scalar S, Dimension N, Function<S, S, S> Func>
		constexpr Vector<S, N> Apply(const Func& func, const Tensor<S, N>& vector1, const Tensor<S, N>& vector2)
		{
			Vector<S, N> result;
			for (Dimension n{}; n < N; ++n)
				result[n] = func(vector1[n], vector2[n]);
			return result;
		}

		template <Scalar S, Dimension N, Function<S, S, S, S> Func>
		constexpr Vector<S, N> Apply(const Func& func, const Tensor<S, N> vector1, const Tensor<S, N>& vector2, const Tensor<S, N>& vector3)
		{
			Vector<S, N> result;
			for (Dimension n{}; n < N; ++n)
				result[n] = func(vector1[n], vector2[n], vector3[n]);
			return result;
		}

		template <Scalar S, Dimension N, Function<S, S, S, S> Func>
		constexpr Vector<S, N> Apply(const Func& func, const Tensor<S, N> vector1, const Tensor<S, N>& vector2, S scalar)
		{
			Vector<S, N> result;
			for (Dimension n{}; n < N; ++n)
				result[n] = func(vector1[n], vector2[n], scalar);
			return result;
		}
	}

	template <Scalar S1, Scalar S2, Dimension N>
	constexpr _IMPL CT<S1, S2> Dot(const Tensor<S1, N>& vector1, const Tensor<S2, N>& vector2) noexcept
	{
		return _ND InnerProduct(vector1, vector2);
	}

	template <Scalar S, Dimension N>
	constexpr S Length2(const Tensor<S, N>& vector) noexcept
	{
		return _ND Dot(vector, vector);
	}

	template <Scalar S, Dimension N>
	constexpr S Length(const Tensor<S, N>& vector) noexcept
	{
		return S{_GCEM sqrt(_ND Length2(vector))};
	}

	template <Scalar S, Dimension N>
	constexpr Normal<S, N> Normalize(const Tensor<S, N>& vector) noexcept
	{
		return vector;
	}

	template <Scalar S1, Scalar S2, Dimension N>
	constexpr Vector<_IMPL CT<S1, S2>, N> Min(const Tensor<S1, N>& vector1, const Tensor<S2, N>& vector2) noexcept
	{
		return _IMPL Apply(_GCEM min, Vector<_IMPL CT<S1, S2>, N>{vector1}, Vector<_IMPL CT<S1, S2>, N>{vector2});
	}

	template <Scalar S1, Scalar S2, Dimension N>
	constexpr Vector<_IMPL CT<S1, S2>, N> Max(const Tensor<S1, N>& vector1, const Tensor<S2, N>& vector2) noexcept
	{
		return _IMPL Apply(_GCEM max, Vector<_IMPL CT<S1, S2>, N>{vector1}, Vector<_IMPL CT<S1, S2>, N>{vector2});
	}

	template <Scalar S, Dimension N>
	constexpr Vector<S, N> Abs(const Tensor<S, N>& vector) noexcept
	{
		return _IMPL Apply(_GCEM abs, vector);
	}

	template <Scalar S, Dimension N>
	constexpr Vector<S, N> Floor(const Tensor<S, N>& vector) noexcept
	{
		return _IMPL Apply(_GCEM floor, vector);
	}

	template <Scalar S, Dimension N>
	constexpr Vector<S, N> Ceil(const Tensor<S, N>& vector) noexcept
	{
		return _IMPL Apply(_GCEM ceil, vector);
	}

	template <Scalar S1, Scalar S2, Dimension N>
	constexpr Vector<_IMPL CT<S1, S2>, N> Mod(const Tensor<S1, N>& vector1, const Tensor<S2, N>& vector2) noexcept
	{
		return _IMPL Apply(_GCEM fmod, Vector<_IMPL CT<S1, S2>, N>{vector1}, Vector<_IMPL CT<S1, S2>, N>{vector2});
	}

	template <Scalar S1, Scalar S2, Scalar S3, Dimension N>
	constexpr Vector<_IMPL CT<S1, S2, S3>, N> Clamp(const Tensor<S1, N>& vector, const Tensor<S2, N>& min, const Tensor<S3, N>& max) noexcept
	{
		return _IMPL Apply(_STD clamp, Vector<_IMPL CT<S1, S2, S3>, N>{vector}, Vector<_IMPL CT<S1, S2, S3>, N>{min}, Vector<_IMPL CT<S1, S2, S3>, N>{max});
	}

	template <Scalar S, Dimension N>
	constexpr Vector<S, N> Trunc(const Tensor<S, N>& vector) noexcept
	{
		return _IMPL Apply(_GCEM trunc, vector);
	}

	template <Scalar S, Dimension N>
	constexpr Vector<S, N> Round(const Tensor<S, N>& vector) noexcept
	{
		return _IMPL Apply(_GCEM round, vector);
	}

	template <Scalar S1, Scalar S2, Scalar S3, Dimension N>
	constexpr Vector<_IMPL CT<S1, S2, S3>, N> Lerp(const Tensor<S1, N>& vector1, const Tensor<S2, N>& vector2, S3 t) noexcept
	{
		return _IMPL Apply(_IMPL Lerp, Vector<_IMPL CT<S1, S2, S3>, N>{vector1}, Vector<_IMPL CT<S1, S2, S3>, N>{vector2}, _IMPL CT<S1, S2, S3>{t});
	}

	template <Scalar S1, Scalar S2>
	requires(_STD is_convertible_v<S2, S1>)
	constexpr Vector<_IMPL CT<S1, S2>, 3> Cross(const Tensor<S1, 3>& vector1, const Tensor<S2, 3>& vector2) noexcept
	{
		Vector<_IMPL CT<S1, S2>, 3> v1{vector1};
		Vector<_IMPL CT<S1, S2>, 3> v2{vector2};

		return Vector<_IMPL CT<S1, S2>, 3>
		{
			v1[1] * v2[2] - v1[2] * v2[1],
			v1[2] * v2[0] - v1[0] * v2[2],
			v1[0] * v2[1] - v1[1] * v2[0]
		};
	}

	template <Scalar S1, Scalar S2>
	requires(_STD is_convertible_v<S2, S1>)
	constexpr Vector<_IMPL CT<S1, S2>, 7> Cross(const Tensor<S1, 7>& vector1, const Tensor<S2, 7>& vector2) noexcept
	{
		Vector<_IMPL CT<S1, S2>, 7> v1{vector1};
		Vector<_IMPL CT<S1, S2>, 7> v2{vector2};

		// https://en.wikipedia.org/wiki/Seven-dimensional_cross_product#:~:text=The%20result%20is
		return Vector<_IMPL CT<S1, S2>, 7>
		{
			v1[1] * v2[3] - v1[3] * v2[1] + v1[2] * v2[6] - v1[6] * v2[2] + v1[4] * v2[5] - v1[5] * v2[4],
			v1[2] * v2[4] - v1[4] * v2[2] + v1[3] * v2[0] - v1[0] * v2[3] + v1[5] * v2[6] - v1[6] * v2[5],
			v1[3] * v2[5] - v1[5] * v2[3] + v1[4] * v2[1] - v1[1] * v2[4] + v1[6] * v2[0] - v1[0] * v2[6],
			v1[4] * v2[6] - v1[6] * v2[4] + v1[5] * v2[2] - v1[2] * v2[5] + v1[0] * v2[1] - v1[1] * v2[0],
			v1[5] * v2[0] - v1[0] * v2[5] + v1[6] * v2[3] - v1[3] * v2[6] + v1[1] * v2[2] - v1[2] * v2[1],
			v1[6] * v2[1] - v1[1] * v2[6] + v1[0] * v2[4] - v1[4] * v2[0] + v1[2] * v2[3] - v1[3] * v2[2],
			v1[0] * v2[2] - v1[2] * v2[0] + v1[1] * v2[5] - v1[5] * v2[1] + v1[3] * v2[4] - v1[4] * v2[3]
		};
	}

	// Normal
	// These needs to be defined down here because they use nd::Length, which is defined after Normal.

	template <Scalar S, Dimension N>
	constexpr Normal<S, N>::Normal(const Tensor<S, N>& vector)
	{
		Tensor<S, N> normalized = vector / _ND Length(vector);
		for (Dimension n{}; n < N; ++n)
			m_Scalars[n] = normalized[n];
	}

	template <Scalar S, Dimension N>
	constexpr auto Normal<S, N>::operator=(const Tensor<S, N>& vector) -> Normal&
	{
		if (_STD addressof(*this) != _STD addressof(vector))
		{
			Tensor<S, N> normalized = vector / _ND Length(vector);
			for (Dimension n{}; n < N; ++n)
				m_Scalars[n] = normalized[n];
		}
		return *this;
	}

	// Aliases

	template <Dimension N>
	using VectorNf = Vector<float, N>;

	using Vector2f = VectorNf<2>;
	using Vector3f = VectorNf<3>;
	using Vector4f = VectorNf<4>;
	using Vector5f = VectorNf<5>;
	using Vector6f = VectorNf<6>;
	using Vector7f = VectorNf<7>;
	using Vector8f = VectorNf<8>;
}
