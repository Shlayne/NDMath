module;

#include <gcem.hpp>

export module nd.vector;

import nd.dimension;
import nd.scalar;
import nd.impl;
import <iostream>;
import <algorithm>;

#define _ND ::nd::
#define _IMPL ::nd::impl::

namespace nd
{
	template<Dimension N, Scalar S> requires(N > 0)
	struct Vector;
}

template<_ND Scalar... Scalars>
using CT = _STD common_type_t<Scalars...>;

template<_ND Dimension N, _ND Scalar... Scalars>
using CVT = _ND Vector<N, CT<Scalars...>>;

export namespace nd
{
	template<Dimension N, Scalar S> requires(N > 0)
	struct Vector
	{
	public:
		constexpr Vector() noexcept = default;

		template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
		constexpr Vector(const Vector<N2, S2>& vector) noexcept;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr Vector(const S2& scalar) noexcept;

		template<typename... Args> requires(sizeof...(Args) > 1)
		constexpr Vector(Args&&... args) noexcept;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator=(const S2& scalar) noexcept -> Vector<N, S>&;

		template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
		constexpr auto operator=(const Vector<N2, S2>& vector) noexcept -> Vector<N, S>&;
	public:
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator+(const S2& scalar) const noexcept -> CVT<N, S, S2>;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator+=(const S2& scalar) noexcept -> Vector<N, S>&;

		template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
		constexpr auto operator+(const Vector<N2, S2>& vector) const noexcept -> CVT<N, S, S2>;

		template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
		constexpr auto operator+=(const Vector<N2, S2>& vector) noexcept -> Vector<N, S>&;
	public:
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator-(const S2& scalar) const noexcept -> CVT<N, S, S2>;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator-=(const S2& scalar) noexcept -> Vector<N, S>&;

		template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
		constexpr auto operator-(const Vector<N2, S2>& vector) const noexcept -> CVT<N, S, S2>;

		template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
		constexpr auto operator-=(const Vector<N2, S2>& vector) noexcept -> Vector<N, S>&;
	public:
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator*(const S2& scalar) const noexcept -> CVT<N, S, S2>;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator*=(const S2& scalar) noexcept -> Vector<N, S>&;

		template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
		constexpr auto operator*(const Vector<N2, S2>& vector) const noexcept -> CVT<N, S, S2>;

		template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
		constexpr auto operator*=(const Vector<N2, S2>& vector) noexcept -> Vector<N, S>&;
	public:
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator/(const S2& scalar) const noexcept -> CVT<N, S, S2>;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator/=(const S2& scalar) noexcept -> Vector<N, S>&;

		template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
		constexpr auto operator/(const Vector<N2, S2>& vector) const noexcept -> CVT<N, S, S2>;

		template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
		constexpr auto operator/=(const Vector<N2, S2>& vector) noexcept -> Vector<N, S>&;
	public:
		constexpr auto operator+() const noexcept -> Vector<N, S>;
		constexpr auto operator-() const noexcept -> Vector<N, S>;
	public:
		template<Dimension N2> requires(N2 < N)
		constexpr auto at() noexcept -> S&;

		template<Dimension N2> requires(N2 < N)
		constexpr auto at() const noexcept -> const S&;

		constexpr auto operator[](Dimension index) noexcept -> S&;
		constexpr auto operator[](Dimension index) const noexcept -> const S&;
	public:
		template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S, S2>)
		constexpr operator Vector<N2, S2>() const noexcept;

		template<Dimension N2, Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator==(const Vector<N2, S2>& vector) const noexcept -> bool;

		template<Dimension N2, Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator!=(const Vector<N2, S2>& vector) const noexcept -> bool;
	public:
		template<typename THash = _STD hash<S>>
		struct Hash
		{
			constexpr auto operator()(const Vector<N, S>& vector) const noexcept -> size_t;
		};
	private:
		template<Dimension N2, Scalar S2, typename... Args> requires(N2 >= 1 && _STD is_convertible_v<S2, S>)
		constexpr auto Fill(const S2& scalar, Args&&... args) noexcept -> void;

		template<Dimension N2, Scalar S2, Dimension N3, typename... Args> requires(N2 >= N3 && _STD is_convertible_v<S2, S>)
		constexpr auto Fill(const Vector<N3, S2>& vector, Args&&... args) noexcept -> void;
	private:
		S m_Scalars[N]{S{}};
	};

	// External Operators

	template<Dimension N, Scalar S, Scalar S2> requires(N > 0 && _STD is_convertible_v<S2, S>)
	constexpr auto operator+(const S2& scalar, const Vector<N, S>& vector) noexcept -> CVT<N, S, S2>;

	template<Dimension N, Scalar S, Scalar S2> requires(N > 0 && _STD is_convertible_v<S2, S>)
	constexpr auto operator-(const S2& scalar, const Vector<N, S>& vector) noexcept -> CVT<N, S, S2>;

	template<Dimension N, Scalar S, Scalar S2> requires(N > 0 && _STD is_convertible_v<S2, S>)
	constexpr auto operator*(const S2& scalar, const Vector<N, S>& vector) noexcept -> CVT<N, S, S2>;

	template<Dimension N, Scalar S, Scalar S2> requires(N > 0 && _STD is_convertible_v<S2, S>)
	constexpr auto operator/(const S2& scalar, const Vector<N, S>& vector) noexcept -> CVT<N, S, S2>;

	// Floating Methods

	template<Dimension N, Scalar S, Scalar S2> requires(N > 0 && _STD is_convertible_v<S2, S>)
	constexpr auto Dot(const Vector<N, S>& vector1, const Vector<N, S2>& vector2) noexcept -> CT<S, S2>;

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Length2(const Vector<N, S>& vector) noexcept -> S;

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Length(const Vector<N, S>& vector) noexcept -> S;

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Normalize(const Vector<N, S>& vector) noexcept -> Vector<N, S>;

	template<Dimension N, Scalar S, Scalar S2> requires(N > 0 && _STD is_convertible_v<S2, S>)
	constexpr auto Min(const Vector<N, S>& vector1, const Vector<N, S2>& vector2) noexcept -> CVT<N, S, S2>;

	template<Dimension N, Scalar S, Scalar S2> requires(N > 0 && _STD is_convertible_v<S2, S>)
	constexpr auto Max(const Vector<N, S>& vector1, const Vector<N, S2>& vector2) noexcept -> CVT<N, S, S2>;

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Abs(const Vector<N, S>& vector) noexcept -> Vector<N, S>;

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Floor(const Vector<N, S>& vector) noexcept -> Vector<N, S>;

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Ceil(const Vector<N, S>& vector) noexcept -> Vector<N, S>;

	template<Dimension N, Scalar S, Scalar S2> requires(N > 0 && _STD is_convertible_v<S2, S>)
	constexpr auto Mod(const Vector<N, S>& vector1, const Vector<N, S2>& vector2) noexcept -> CVT<N, S, S2>;

	template<Dimension N, Scalar S, Scalar S2, Scalar S3> requires(N > 0 && _STD is_convertible_v<S2, S> && _STD is_convertible_v<S3, S>)
	constexpr auto Clamp(const Vector<N, S>& vector, const Vector<N, S2>& min, const Vector<N, S3>& max) noexcept -> CVT<N, S, S2, S3>;

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Trunc(const Vector<N, S>& vector) noexcept -> Vector<N, S>;

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Round(const Vector<N, S>& vector) noexcept -> Vector<N, S>;

	template<Dimension N, Scalar S, Scalar S2, Scalar S3> requires(N > 0 && _STD is_convertible_v<S2, S> && _STD is_convertible_v<S3, S>)
	constexpr auto Lerp(const Vector<N, S>& vector1, const Vector<N, S2>& vector2, const S3& t) noexcept -> CVT<N, S, S2, S3>;

	template<Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Cross(const Vector<3, S>& vector1, const Vector<3, S2>& vector2) noexcept -> CVT<3, S, S2>;

	template<Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Cross(const Vector<7, S>& vector1, const Vector<7, S2>& vector2) noexcept -> CVT<7, S, S2>;

	template<Dimension N, Scalar S> requires(N > 0)
	auto operator<<(_STD ostream& ostream, const Vector<N, S>& vector) -> _STD ostream&;

	// Aliases

	template<Dimension N>
	using VectorNf = Vector<N, float>;

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
	constexpr auto Lerp(const S& a, const S& b, const S& t) -> S
	{
		return (S{1} - t) * a + t * b;
	};

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Apply(S(*func)(S), const Vector<N, S>& vector) -> Vector<N, S>
	{
		Vector<N, S> result;
		for (Dimension n{}; n < N; n++)
			result[n] = func(vector[n]);
		return result;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Apply(S(*func)(S, S), const Vector<N, S>& vector1, const Vector<N, S>& vector2) -> Vector<N, S>
	{
		Vector<N, S> result;
		for (Dimension n{}; n < N; n++)
			result[n] = func(vector1[n], vector2[n]);
		return result;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Apply(S(*func)(const S&, const S&, const S&), const Vector<N, S>& vector1, const Vector<N, S>& vector2, const Vector<N, S>& vector3) -> Vector<N, S>
	{
		Vector<N, S> result;
		for (Dimension n{}; n < N; n++)
			result[n] = func(vector1[n], vector2[n], vector3[n]);
		return result;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Apply(S(*func)(const S&, const S&, const S&), const Vector<N, S>& vector1, const Vector<N, S>& vector2, const S& scalar) -> Vector<N, S>
	{
		Vector<N, S> result;
		for (Dimension n{}; n < N; n++)
			result[n] = func(vector1[n], vector2[n], scalar);
		return result;
	}
}

export namespace nd
{
	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
	constexpr Vector<N, S>::Vector(const Vector<N2, S2>& vector) noexcept
	{
		for (Dimension n{}; n < N2; n++)
			m_Scalars[n] = S{vector[n]};
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr Vector<N, S>::Vector(const S2& scalar) noexcept
	{
		for (S& e : m_Scalars)
			e = S{scalar};
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<typename... Args> requires(sizeof...(Args) > 1)
	constexpr Vector<N, S>::Vector(Args&&... args) noexcept
	{
		Fill<N>(_STD forward<Args>(args)...);
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator=(const S2& scalar) noexcept -> Vector<N, S>&
	{
		for (S& s : m_Scalars)
			s = S{scalar};
		return *this;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator=(const Vector<N2, S2>& vector) noexcept -> Vector<N, S>&
	{
		if (this != &vector)
		{
			Dimension n{};
			for (; n < N2; n++)
				m_Scalars[n] = S{vector[n]};
			for (; n < N; n++)
				m_Scalars[n] = S{};
		}
		return *this;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator+(const S2& scalar) const noexcept -> CVT<N, S, S2>
	{
		return CVT<N, S, S2>{*this} += CT<S, S2>{scalar};
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator+=(const S2& scalar) noexcept -> Vector<N, S>&
	{
		for (S& e : m_Scalars)
			e += S{scalar};
		return *this;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator+(const Vector<N2, S2>& vector) const noexcept -> CVT<N, S, S2>
	{
		return CVT<N, S, S2>{*this} += CVT<N, S, S2>{vector};
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator+=(const Vector<N2, S2>& vector) noexcept -> Vector<N, S>&
	{
		for (Dimension n{}; n < N2; n++)
			m_Scalars[n] += S{vector[n]};
		return *this;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator-(const S2& scalar) const noexcept -> CVT<N, S, S2>
	{
		return CVT<N, S, S2>{*this} -= CT<S, S2>{scalar};
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator-=(const S2& scalar) noexcept -> Vector<N, S>&
	{
		for (S& s : m_Scalars)
			s -= S{scalar};
		return *this;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator-(const Vector<N2, S2>& vector) const noexcept -> CVT<N, S, S2>
	{
		return CVT<N, S, S2>{*this} -= CVT<N, S, S2>{vector};
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator-=(const Vector<N2, S2>& vector) noexcept -> Vector<N, S>&
	{
		for (Dimension n{}; n < N2; n++)
			m_Scalars[n] -= S{vector[n]};
		return *this;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator*(const S2& scalar) const noexcept -> CVT<N, S, S2>
	{
		return CVT<N, S, S2>{*this} *= CT<S, S2>{scalar};
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator*=(const S2& scalar) noexcept -> Vector<N, S>&
	{
		for (S& s : m_Scalars)
			s *= S{scalar};
		return *this;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator*(const Vector<N2, S2>& vector) const noexcept -> CVT<N, S, S2>
	{
		return CVT<N, S, S2>{*this} *= CVT<N, S, S2>{vector};
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator*=(const Vector<N2, S2>& vector) noexcept -> Vector<N, S>&
	{
		for (Dimension n{}; n < N2; n++)
			m_Scalars[n] *= S{vector[n]};
		return *this;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator/(const S2& scalar) const noexcept -> CVT<N, S, S2>
	{
		return CVT<N, S, S2>{*this} /= CT<S, S2>{scalar};
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator/=(const S2& scalar) noexcept -> Vector<N, S>&
	{
		for (S& s : m_Scalars)
			s /= S{scalar};
		return *this;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator/(const Vector<N2, S2>& vector) const noexcept -> CVT<N, S, S2>
	{
		return CVT<N, S, S2>{*this} /= CVT<N, S, S2>{vector};
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator/=(const Vector<N2, S2>& vector) noexcept -> Vector<N, S>&
	{
		for (Dimension n{}; n < N2; n++)
			m_Scalars[n] /= S{vector[n]};
		return *this;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Vector<N, S>::operator+() const noexcept -> Vector<N, S>
	{
		return Vector<N, S>{*this};
	}

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Vector<N, S>::operator-() const noexcept -> Vector<N, S>
	{
		Vector<N, S> vector{+*this};
		for (S& s : m_Scalars)
			s = -s;
		return vector;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2> requires(N2 < N)
	constexpr auto Vector<N, S>::at() noexcept -> S&
	{
		return m_Scalars[N2];
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2> requires(N2 < N)
	constexpr auto Vector<N, S>::at() const noexcept -> const S&
	{
		return m_Scalars[N2];
	}

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Vector<N, S>::operator[](Dimension index) noexcept -> S&
	{
		__assume(index < N);
		return m_Scalars[index];
	}

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Vector<N, S>::operator[](Dimension index) const noexcept -> const S&
	{
		__assume(index < N);
		return m_Scalars[index];
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2, Scalar S2> requires(N2 <= N && _STD is_convertible_v<S, S2>)
	constexpr Vector<N, S>::operator Vector<N2, S2>() const noexcept
	{
		Vector<N2, S2> result;
		for (Dimension n{}; n < N2; n++)
			result[n] = S2{m_Scalars[n]};
		return result;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator==(const Vector<N2, S2>& vector) const noexcept -> bool
	{
		Dimension n{};
		for (; n < gcem::min(N, N2); n++)
			if (CT<S, S2>(m_Scalars[n]) != CT<S, S2>{vector[n]})
				return false;
		if constexpr (N < N2)
		{
			for (; n < N2; n++)
				if (CT<S, S2>{vector[n]} != CT<S, S2>{})
					return false;
		}
		else if constexpr (N2 < N)
		{
			for (; n < N; n++)
				if (CT<S, S2>{m_Scalars[n]} != CT<S, S2>{})
					return false;
		}
		return true;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::operator!=(const Vector<N2, S2>& vector) const noexcept -> bool
	{
		Dimension n{};
		for (; n < gcem::min(N, N2); n++)
			if (CT<S, S2>{m_Scalars[n]} == CT<S, S2>{vector[n]})
				return false;
		if constexpr (N < N2)
		{
			for (; n < N2; n++)
				if (CT<S, S2>{vector[n]} == CT<S, S2>{})
					return false;
		}
		else if constexpr (N2 < N)
		{
			for (; n < N; n++)
				if (CT<S, S2>{m_Scalars[n]} == CT<S, S2>{})
					return false;
		}
		return true;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<typename THash>
	constexpr auto Vector<N, S>::Hash<THash>::operator()(const Vector<N, S>& vector) const noexcept -> size_t
	{
		size_t hash{};
		for (Dimension n{}; n < N; n++)
		{
			size_t scalarHash{THash{}(vector[n])};
			hash ^= (scalarHash << ((7 * n) & ((1ull << 6) - 1))) | ((scalarHash >> 56) & 0x7F);
		}
		return hash;
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2, Scalar S2, typename... Args> requires(N2 >= 1 && _STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::Fill(const S2& scalar, Args&&... args) noexcept -> void
	{
		m_Scalars[N - N2] = S{scalar};
		if constexpr (sizeof...(Args) > 0)
			Fill<N2 - 1>(_STD forward<Args>(args)...);
	}

	template<Dimension N, Scalar S> requires(N > 0)
	template<Dimension N2, Scalar S2, Dimension N3, typename... Args> requires(N2 >= N3 && _STD is_convertible_v<S2, S>)
	constexpr auto Vector<N, S>::Fill(const Vector<N3, S2>& vector, Args&&... args) noexcept -> void
	{
		for (Dimension n{}; n < N3; n++)
			m_Scalars[N - N2 + n] = S{vector[n]};
		if constexpr (sizeof...(Args) > 0)
			Fill<N2 - N3>(_STD forward<Args>(args)...);
	}

	template<Dimension N, Scalar S, Scalar S2> requires(N > 0 && _STD is_convertible_v<S2, S>)
	constexpr auto operator+(const S2& scalar, const Vector<N, S>& vector) noexcept -> CVT<N, S, S2>
	{
		return CVT<N, S, S2>{scalar} += CVT<N, S, S2>{vector};
	}

	template<Dimension N, Scalar S, Scalar S2> requires(N > 0 && _STD is_convertible_v<S2, S>)
	constexpr auto operator-(const S2& scalar, const Vector<N, S>& vector) noexcept -> CVT<N, S, S2>
	{
		return CVT<N, S, S2>{scalar} -= CVT<N, S, S2>{vector};
	}

	template<Dimension N, Scalar S, Scalar S2> requires(N > 0 && _STD is_convertible_v<S2, S>)
	constexpr auto operator*(const S2& scalar, const Vector<N, S>& vector) noexcept -> CVT<N, S, S2>
	{
		return CVT<N, S, S2>{scalar} *= CVT<N, S, S2>{vector};
	}

	template<Dimension N, Scalar S, Scalar S2> requires(N > 0 && _STD is_convertible_v<S2, S>)
	constexpr auto operator/(const S2& scalar, const Vector<N, S>& vector) noexcept -> CVT<N, S, S2>
	{
		return CVT<N, S, S2>{scalar} /= CVT<N, S, S2>{vector};
	}

	template<Dimension N, Scalar S, Scalar S2> requires(N > 0 && _STD is_convertible_v<S2, S>)
	constexpr auto Dot(const Vector<N, S>& vector1, const Vector<N, S2>& vector2) noexcept -> CT<S, S2>
	{
		CT<S, S2> dot{};
		CVT<N, S, S2> product{CVT<N, S, S2>{vector1} * CVT<N, S, S2>{vector2}};
		for (Dimension n{}; n < N; n++)
			dot += product[n];
		return dot;
	}
	
	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Length2(const Vector<N, S>& vector) noexcept -> S
	{
		return Dot(vector, vector);
	}

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Length(const Vector<N, S>& vector) noexcept -> S
	{
		return S{gcem::sqrt(Length2(vector))};
	}

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Normalize(const Vector<N, S>& vector) noexcept -> Vector<N, S>
	{
		return vector / Length(vector);
	}

	template<Dimension N, Scalar S, Scalar S2> requires(N > 0 && _STD is_convertible_v<S2, S>)
	constexpr auto Min(const Vector<N, S>& vector1, const Vector<N, S2>& vector2) noexcept -> CVT<N, S, S2>
	{
		return _IMPL Apply(gcem::min, CVT<N, S, S2>{vector1}, CVT<N, S, S2>{vector2});
	}

	template<Dimension N, Scalar S, Scalar S2> requires(N > 0 && _STD is_convertible_v<S2, S>)
	constexpr auto Max(const Vector<N, S>& vector1, const Vector<N, S2>& vector2) noexcept -> CVT<N, S, S2>
	{
		return _IMPL Apply(gcem::max, CVT<N, S, S2>{vector1}, CVT<N, S, S2>{vector2});
	}

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Abs(const Vector<N, S>& vector) noexcept -> Vector<N, S>
	{
		return _IMPL Apply(gcem::abs, vector);
	}

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Floor(const Vector<N, S>& vector) noexcept -> Vector<N, S>
	{
		return _IMPL Apply(gcem::floor, vector);
	}

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Ceil(const Vector<N, S>& vector) noexcept -> Vector<N, S>
	{
		return _IMPL Apply(gcem::ceil, vector);
	}

	template<Dimension N, Scalar S, Scalar S2> requires(N > 0 && _STD is_convertible_v<S2, S>)
	constexpr auto Mod(const Vector<N, S>& vector1, const Vector<N, S2>& vector2) noexcept -> CVT<N, S, S2>
	{
		return _IMPL Apply(gcem::fmod, CVT<N, S, S2>{vector1}, CVT<N, S, S2>{vector2});
	}

	template<Dimension N, Scalar S, Scalar S2, Scalar S3> requires(N > 0 && _STD is_convertible_v<S2, S> && _STD is_convertible_v<S3, S>)
	constexpr auto Clamp(const Vector<N, S>& vector, const Vector<N, S2>& min, const Vector<N, S3>& max) noexcept -> CVT<N, S, S2, S3>
	{
		return _IMPL Apply(_STD clamp, CVT<N, S, S2, S3>{vector}, CVT<N, S, S2, S3>{min}, CVT<N, S, S2, S3>{max});
	}

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Trunc(const Vector<N, S>& vector) noexcept -> Vector<N, S>
	{
		return _IMPL Apply(gcem::trunc, vector);
	}

	template<Dimension N, Scalar S> requires(N > 0)
	constexpr auto Round(const Vector<N, S>& vector) noexcept -> Vector<N, S>
	{
		return _IMPL Apply(gcem::round, vector);
	}

	template<Dimension N, Scalar S, Scalar S2, Scalar S3> requires(N > 0 && _STD is_convertible_v<S2, S> && _STD is_convertible_v<S3, S>)
	constexpr auto Lerp(const Vector<N, S>& vector1, const Vector<N, S2>& vector2, const S3& t) noexcept -> CVT<N, S, S2, S3>
	{
		return _IMPL Apply(_IMPL Lerp, CVT<N, S, S2, S3>{vector1}, CVT<N, S, S2, S3>{vector2}, CT<S, S2, S3>{t});
	}

	template<Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Cross(const Vector<3, S>& vector1, const Vector<3, S2>& vector2) noexcept -> CVT<3, S, S2>
	{
		CVT<3, S, S2> v1{vector1};
		CVT<3, S, S2> v2{vector2};

		return CVT<3, S, S2>
		{
			v1[1] * v2[2] - v1[2] * v2[1],
			v1[2] * v2[0] - v1[0] * v2[2],
			v1[0] * v2[1] - v1[1] * v2[0]
		};
	}

	template<Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Cross(const Vector<7, S>& vector1, const Vector<7, S2>& vector2) noexcept -> CVT<7, S, S2>
	{
		CVT<7, S, S2> v1{vector1};
		CVT<7, S, S2> v2{vector2};

		// https://en.wikipedia.org/wiki/Seven-dimensional_cross_product#:~:text=The%20result%20is
		return CVT<7, S, S2>
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

	template<Dimension N, Scalar S> requires(N > 0)
	auto operator<<(_STD ostream& ostream, const Vector<N, S>& vector) -> _STD ostream&
	{
		ostream << '<' << vector[0];
		for (Dimension n{1}; n < N; n++)
			ostream << ',' << vector[n];
		return ostream << '>';
	}
}
