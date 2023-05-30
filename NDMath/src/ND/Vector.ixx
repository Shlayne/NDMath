module;

#include <gcem.hpp>

export module nd.vector;

import nd.types;
import nd.impl;
import nd.tensor;
import <algorithm>;

#define _ND ::nd::
#define _IMPL ::nd::impl::
#define _GCEM ::gcem:: 

export namespace nd
{
	template<Scalar S, Dimension N>
	struct Vector : public Tensor<S, 1, N>
	{
	public:
		template<typename... Args> requires(sizeof...(Args) > 1)
		constexpr Vector(Args&&... args) noexcept;
	public:
		using Tensor<S, 1, N>::Tensor;
		using Tensor<S, 1, N>::operator=;
		using Tensor<S, 1, N>::operator+; // also handles unary plus.
		using Tensor<S, 1, N>::operator+=;
		using Tensor<S, 1, N>::operator-; // also handles unary minus.
		using Tensor<S, 1, N>::operator-=;
		using Tensor<S, 1, N>::operator*;
		using Tensor<S, 1, N>::operator*=;
		using Tensor<S, 1, N>::operator/;
		using Tensor<S, 1, N>::operator/=;
		using Tensor<S, 1, N>::at;
		using Tensor<S, 1, N>::operator[];
		using Tensor<S, 1, N>::operator==;
		using Tensor<S, 1, N>::operator!=;
	public:
		// NOTE: using a base class' templated conversion operator is expliticly prohibited by C++23.
		template<Scalar S2, Dimension N2> requires(_STD is_convertible_v<S, S2>)
		constexpr operator Vector<S2, N2>() const noexcept;
	private:
		template<Dimension N2, Scalar S2, typename... Args> requires(N2 >= 1 && _STD is_convertible_v<S2, S>)
		constexpr auto Fill(const S2& scalar, Args&&... args) noexcept -> void;

		template<Dimension N2, Scalar S2, Dimension N3, typename... Args> requires(N2 >= N3 && _STD is_convertible_v<S2, S>)
		constexpr auto Fill(const Vector<S2, N3>& vector, Args&&... args) noexcept -> void;
	private:
		using Tensor<S, 1, N>::m_Scalars;
	};

	// Floating Methods

	template<Scalar S, Dimension N, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Dot(const Vector<S, N>& vector1, const Vector<S2, N>& vector2) noexcept -> _IMPL CT<S, S2>;

	template<Scalar S, Dimension N>
	constexpr auto Length2(const Vector<S, N>& vector) noexcept -> S;

	template<Scalar S, Dimension N>
	constexpr auto Length(const Vector<S, N>& vector) noexcept -> S;

	template<Scalar S, Dimension N>
	constexpr auto Normalize(const Vector<S, N>& vector) noexcept -> Vector<S, N>;

	template<Scalar S, Dimension N, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Min(const Vector<S, N>& vector1, const Vector<S2, N>& vector2) noexcept -> Vector<_IMPL CT<S, S2>, N>;

	template<Scalar S, Dimension N, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Max(const Vector<S, N>& vector1, const Vector<S2, N>& vector2) noexcept -> Vector<_IMPL CT<S, S2>, N>;

	template<Scalar S, Dimension N>
	constexpr auto Abs(const Vector<S, N>& vector) noexcept -> Vector<S, N>;

	template<Scalar S, Dimension N>
	constexpr auto Floor(const Vector<S, N>& vector) noexcept -> Vector<S, N>;

	template<Scalar S, Dimension N>
	constexpr auto Ceil(const Vector<S, N>& vector) noexcept -> Vector<S, N>;

	template<Scalar S, Dimension N, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Mod(const Vector<S, N>& vector1, const Vector<S2, N>& vector2) noexcept -> Vector<_IMPL CT<S, S2>, N>;

	template<Scalar S, Dimension N, Scalar S2, Scalar S3> requires(_STD is_convertible_v<S2, S> && _STD is_convertible_v<S3, S>)
	constexpr auto Clamp(const Vector<S, N>& vector, const Vector<S2, N>& min, const Vector<S3, N>& max) noexcept -> Vector<_IMPL CT<S, S2, S3>, N>;

	template<Scalar S, Dimension N>
	constexpr auto Trunc(const Vector<S, N>& vector) noexcept -> Vector<S, N>;

	template<Scalar S, Dimension N>
	constexpr auto Round(const Vector<S, N>& vector) noexcept -> Vector<S, N>;

	template<Scalar S, Dimension N, Scalar S2, Scalar S3> requires(_STD is_convertible_v<S2, S> && _STD is_convertible_v<S3, S>)
	constexpr auto Lerp(const Vector<S, N>& vector1, const Vector<S2, N>& vector2, const S3& t) noexcept -> Vector<_IMPL CT<S, S2, S3>, N>;

	template<Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Cross(const Vector<S, 3>& vector1, const Vector<S2, 3>& vector2) noexcept -> Vector<_IMPL CT<S, S2>, 3>;

	template<Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Cross(const Vector<S, 7>& vector1, const Vector<S2, 7>& vector2) noexcept -> Vector<_IMPL CT<S, S2>, 7>;

	// Aliases

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
}

// Implementation: Don't export.
namespace nd::impl
{
	template<Scalar S>
	constexpr auto Lerp(S a, S b, S t) -> S
	{
		return (S{1} - t) * a + t * b;
	};

	template<Scalar S, Dimension N>
	constexpr auto Apply(S(*func)(const S&), const Vector<S, N>& vector) -> Vector<S, N>
	{
		Vector<S, N> result;
		for (Dimension n{}; n < N; ++n)
			result[n] = func(vector[n]);
		return result;
	}

	template<Scalar S, Dimension N>
	constexpr auto Apply(S(*func)(const S&, const S&), const Vector<S, N>& vector1, const Vector<S, N>& vector2) -> Vector<S, N>
	{
		Vector<S, N> result;
		for (Dimension n{}; n < N; ++n)
			result[n] = func(vector1[n], vector2[n]);
		return result;
	}

	template<Scalar S, Dimension N>
	constexpr auto Apply(S(*func)(const S&, const S&, const S&), const Vector<S, N> vector1, const Vector<S, N>& vector2, const Vector<S, N>& vector3) -> Vector<S, N>
	{
		Vector<S, N> result;
		for (Dimension n{}; n < N; ++n)
			result[n] = func(vector1[n], vector2[n], vector3[n]);
		return result;
	}

	template<Scalar S, Dimension N>
	constexpr auto Apply(S(*func)(const S&, const S&, const S&), const Vector<S, N> vector1, const Vector<S, N>& vector2, const S& scalar) -> Vector<S, N>
	{
		Vector<S, N> result;
		for (Dimension n{}; n < N; ++n)
			result[n] = func(vector1[n], vector2[n], scalar);
		return result;
	}
}

export namespace nd
{
	template<Scalar S, Dimension N>
	template<typename... Args> requires(sizeof...(Args) > 1)
	constexpr Vector<S, N>::Vector(Args&&... args) noexcept
	{
		Fill<N>(_STD forward<Args>(args)...);
	}

	template<Scalar S, Dimension N>
	template<Scalar S2, Dimension N2> requires(_STD is_convertible_v<S, S2>)
	constexpr Vector<S, N>::operator Vector<S2, N2>() const noexcept
	{
		Vector<S2, N2> result;
		for (Dimension n{}; n < _GCEM min(N, N2); ++n)
			result[n] = static_cast<S2>(m_Scalars[n]);
		return result;
	}

	template<Scalar S, Dimension N>
	template<Dimension N2, Scalar S2, typename... Args> requires(N2 >= 1 && _STD is_convertible_v<S2, S>)
	constexpr auto Vector<S, N>::Fill(const S2& scalar, Args&&... args) noexcept -> void
	{
		m_Scalars[N - N2] = static_cast<S>(scalar);
		if constexpr (sizeof...(Args) > 0)
			Fill<N2 - 1>(_STD forward<Args>(args)...);
	}

	template<Scalar S, Dimension N>
	template<Dimension N2, Scalar S2, Dimension N3, typename... Args> requires(N2 >= N3 && _STD is_convertible_v<S2, S>)
	constexpr auto Vector<S, N>::Fill(const Vector<S2, N3>& vector, Args&&... args) noexcept -> void
	{
		for (Dimension n{}; n < N3; ++n)
			m_Scalars[N - N2 + n] = static_cast<S>(vector[n]);
		if constexpr (sizeof...(Args) > 0)
			Fill<N2 - N3>(_STD forward<Args>(args)...);
	}

	template<Scalar S, Dimension N, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Dot(const Vector<S, N>& vector1, const Vector<S2, N>& vector2) noexcept -> _IMPL CT<S, S2>
	{
		_IMPL CT<S, S2> dot{};
		Vector<_IMPL CT<S, S2>, N> product{Vector<_IMPL CT<S, S2>, N>{vector1} * Vector<_IMPL CT<S, S2>, N>{vector2}};
		for (Dimension n{}; n < N; ++n)
			dot += product[n];
		return dot;
	}
	
	template<Scalar S, Dimension N>
	constexpr auto Length2(const Vector<S, N>& vector) noexcept -> S
	{
		return Dot(vector, vector);
	}

	template<Scalar S, Dimension N>
	constexpr auto Length(const Vector<S, N>& vector) noexcept -> S
	{
		return S{_GCEM sqrt(Length2(vector))};
	}

	template<Scalar S, Dimension N>
	constexpr auto Normalize(const Vector<S, N>& vector) noexcept -> Vector<S, N>
	{
		return vector / Length(vector);
	}

	template<Scalar S, Dimension N, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Min(const Vector<S, N>& vector1, const Vector<S2, N>& vector2) noexcept -> Vector<_IMPL CT<S, S2>, N>
	{
		return _IMPL Apply(_GCEM min, Vector<_IMPL CT<S, S2>, N>{vector1}, Vector<_IMPL CT<S, S2>, N>{vector2});
	}

	template<Scalar S, Dimension N, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Max(const Vector<S, N>& vector1, const Vector<S2, N>& vector2) noexcept -> Vector<_IMPL CT<S, S2>, N>
	{
		return _IMPL Apply(_GCEM max, Vector<_IMPL CT<S, S2>, N>{vector1}, Vector<_IMPL CT<S, S2>, N>{vector2});
	}

	template<Scalar S, Dimension N>
	constexpr auto Abs(const Vector<S, N>& vector) noexcept -> Vector<S, N>
	{
		return _IMPL Apply(_GCEM abs, vector);
	}

	template<Scalar S, Dimension N>
	constexpr auto Floor(const Vector<S, N>& vector) noexcept -> Vector<S, N>
	{
		return _IMPL Apply(_GCEM floor, vector);
	}

	template<Scalar S, Dimension N>
	constexpr auto Ceil(const Vector<S, N>& vector) noexcept -> Vector<S, N>
	{
		return _IMPL Apply(_GCEM ceil, vector);
	}

	template<Scalar S, Dimension N, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Mod(const Vector<S, N>& vector1, const Vector<S2, N>& vector2) noexcept -> Vector<_IMPL CT<S, S2>, N>
	{
		return _IMPL Apply(_GCEM fmod, Vector<_IMPL CT<S, S2>, N>{vector1}, Vector<_IMPL CT<S, S2>, N>{vector2});
	}

	template<Scalar S, Dimension N, Scalar S2, Scalar S3> requires(_STD is_convertible_v<S2, S> && _STD is_convertible_v<S3, S>)
	constexpr auto Clamp(const Vector<S, N>& vector, const Vector<S2, N>& min, const Vector<S3, N>& max) noexcept -> Vector<_IMPL CT<S, S2, S3>, N>
	{
		return _IMPL Apply(_STD clamp, Vector<_IMPL CT<S, S2, S3>, N>{vector}, Vector<_IMPL CT<S, S2, S3>, N>{min}, Vector<_IMPL CT<S, S2, S3>, N>{max});
	}

	template<Scalar S, Dimension N>
	constexpr auto Trunc(const Vector<S, N>& vector) noexcept -> Vector<S, N>
	{
		return _IMPL Apply(_GCEM trunc, vector);
	}

	template<Scalar S, Dimension N>
	constexpr auto Round(const Vector<S, N>& vector) noexcept -> Vector<S, N>
	{
		return _IMPL Apply(_GCEM round, vector);
	}

	template<Scalar S, Dimension N, Scalar S2, Scalar S3> requires(_STD is_convertible_v<S2, S> && _STD is_convertible_v<S3, S>)
	constexpr auto Lerp(const Vector<S, N>& vector1, const Vector<S2, N>& vector2, const S3& t) noexcept -> Vector<_IMPL CT<S, S2, S3>, N>
	{
		return _IMPL Apply(_IMPL Lerp, Vector<_IMPL CT<S, S2, S3>, N>{vector1}, Vector<_IMPL CT<S, S2, S3>, N>{vector2}, _IMPL CT<S, S2, S3>{t});
	}

	template<Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Cross(const Vector<S, 3>& vector1, const Vector<S2, 3>& vector2) noexcept -> Vector<_IMPL CT<S, S2>, 3>
	{
		Vector<_IMPL CT<S, S2>, 3> v1{vector1};
		Vector<_IMPL CT<S, S2>, 3> v2{vector2};

		return Vector<_IMPL CT<S, S2>, 3>
		{
			v1[1] * v2[2] - v1[2] * v2[1],
			v1[2] * v2[0] - v1[0] * v2[2],
			v1[0] * v2[1] - v1[1] * v2[0]
		};
	}

	template<Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Cross(const Vector<S, 7>& vector1, const Vector<S2, 7>& vector2) noexcept -> Vector<_IMPL CT<S, S2>, 7>
	{
		Vector<_IMPL CT<S, S2>, 7> v1{vector1};
		Vector<_IMPL CT<S, S2>, 7> v2{vector2};

		// https://en.wikipedia.org/wiki/Seven-dimensional_cross_product#:~:text=The%20result%20is
		return Vector<_IMPL CT<S, S2>, 7>
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
}
