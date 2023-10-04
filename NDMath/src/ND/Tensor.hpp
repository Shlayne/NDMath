#pragma once

#include "Impl.hpp"
#include <gcem.hpp>
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>

#define _ND ::nd::
#define _IMPL ::nd::impl::
#define _GCEM ::gcem::

namespace nd
{
	template <Scalar S, Dimension N = 0, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	struct Tensor
	{
	public:
		static constexpr Dimension Rank = static_cast<Dimension>(N > 0 ? 1 : 0) + sizeof...(Ns);

		constexpr Tensor() noexcept;

		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor(S2 scalar) noexcept;

		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor(const Tensor<S2, N2, N2s...>& tensor) noexcept;

		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator=(S2 scalar) noexcept;

		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator=(const Tensor<S2, N2, N2s...>& tensor) noexcept;
	public:
		// Scalar addition.
		template <Scalar S2>
		constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator+(S2 scalar) const noexcept;

		// Scalar addition.
		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator+=(S2 scalar) noexcept;

		// Element-wise addition.
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator+(const Tensor<S2, N2, N2s...>& tensor) const noexcept;

		// Element-wise addition.
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator+=(const Tensor<S2, N2, N2s...>& tensor) noexcept;
	public:
		// Scalar subtraction.
		template <Scalar S2>
		constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator-(S2 scalar) const noexcept;

		// Scalar subtraction.
		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator-=(S2 scalar) noexcept;

		// Element-wise subtraction.
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator-(const Tensor<S2, N2, N2s...>& tensor) const noexcept;

		// Element-wise subtraction.
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator-=(const Tensor<S2, N2, N2s...>& tensor) noexcept;
	public:
		// Scalar multiplication.
		template <Scalar S2>
		constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator*(S2 scalar) const noexcept;

		// Scalar multiplication.
		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator*=(S2 scalar) noexcept;

		// Element-wise multiplication.
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator*(const Tensor<S2, N2, N2s...>& tensor) const noexcept;

		// Element-wise multiplication.
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator*=(const Tensor<S2, N2, N2s...>& tensor) noexcept;
	public:
		// Scalar division.
		template <Scalar S2>
		constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator/(S2 scalar) const noexcept;

		// Scalar division.
		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator/=(S2 scalar) noexcept;

		// Element-wise division.
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator/(const Tensor<S2, N2, N2s...>& tensor) const noexcept;

		// Element-wise division.
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator/=(const Tensor<S2, N2, N2s...>& tensor) noexcept;
	public:
		constexpr Tensor<decltype(+S{}), N, Ns... > operator+() const noexcept;
		constexpr Tensor<decltype(-S{}), N, Ns...> operator-() const noexcept;
	public:
		template <Dimension N2>
		requires(N2 < N)
		constexpr Tensor<S, Ns...>& at() noexcept;

		template <Dimension N2>
		requires(N2 < N)
		constexpr const Tensor<S, Ns...>& at() const noexcept;

		constexpr Tensor<S, Ns...>& operator[](Dimension n) noexcept;
		constexpr const Tensor<S, Ns...>& operator[](Dimension n) const noexcept;
	public:
		template <Scalar S2>
		requires(!_STD is_same_v<S, S2> && _STD is_convertible_v<S, S2>)
		constexpr operator Tensor<S2, N, Ns...>() const noexcept;

		template <Scalar S2>
		constexpr bool operator==(const Tensor<S2, N, Ns...>& tensor) const noexcept;

		template <Scalar S2>
		constexpr bool operator!=(const Tensor<S2, N, Ns...>& tensor) const noexcept;
	public:
		template <typename SHasher = _STD hash<S>>
		struct Hasher
		{
			constexpr size_t operator()(const Tensor<S, N, Ns...>& tensor) const noexcept;
		};

		constexpr size_t Hash() const noexcept;
	private:
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		requires(_STD is_convertible_v<S2, S>)
		constexpr auto Set(Dimension n, const Tensor<S2, N2, N2s...>& tensor) noexcept;
	protected:
		_STD array<Tensor<S, Ns...>, N> m_Scalars;
	};

	template <Scalar S>
	struct Tensor<S>
	{
	public:
		static constexpr Dimension Rank = 0;

		constexpr Tensor() noexcept = default;

		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor(S2 scalar) noexcept
			: m_Scalar{static_cast<S>(scalar)}
		{

		}

		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S>& operator=(S2 scalar) noexcept
		{
			m_Scalar = static_cast<S>(scalar);
			return *this;
		}
	public:
		constexpr operator S&() noexcept;
		constexpr operator const S&() const noexcept;
	private:
		S m_Scalar{};
	};

	// External Operators

	template <Scalar S, Scalar S2, Dimension N = 0, Dimension... Ns>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator+(S scalar, const Tensor<S2, N, Ns...>& tensor) noexcept;

	template <Scalar S, Scalar S2, Dimension N = 0, Dimension... Ns>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator-(S scalar, const Tensor<S2, N, Ns...>& tensor) noexcept;

	template <Scalar S, Scalar S2, Dimension N = 0, Dimension... Ns>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator*(S scalar, const Tensor<S2, N, Ns...>& tensor) noexcept;

	template <Scalar S, Scalar S2, Dimension N = 0, Dimension... Ns>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator/(S scalar, const Tensor<S2, N, Ns...>& tensor) noexcept;

	template <Scalar S>
	_STD ostream& operator<<(_STD ostream& ostream, const Tensor<S>& scalar);

	template <Scalar S, Dimension N>
	_STD ostream& operator<<(_STD ostream& ostream, const Tensor<S, N>& vector);

	template <Scalar S, Dimension C, Dimension R>
	_STD ostream& operator<<(_STD ostream& ostream, const Tensor<S, C, R>& matrix);

	// Static Methods

	template <Scalar S1, Scalar S2, Dimension N1 = 0, Dimension N2 = 0, Dimension... N1s, Dimension... N2s>
	constexpr Tensor<_IMPL CT<S1, S2>, N1, N1s..., N2, N2s...> OuterProduct(const Tensor<S1, N1, N1s...>& tensor1, const Tensor<S2, N2, N2s...>& tensor2) noexcept;

	template <Scalar S1, Scalar S2, Dimension N = 0, Dimension... Ns>
	constexpr _IMPL CT<S1, S2> InnerProduct(const Tensor<S1, N, Ns...>& tensor1, const Tensor<S2, N, Ns...>& tensor2) noexcept;
}

namespace nd::impl
{
	template <Scalar S, Dimension C, Dimension R>
	void PrintRow(_STD ostream& ostream, const Tensor<S, C, R>& matrix, _STD streamsize width, Dimension r) noexcept
	{
		ostream << '[' << _STD setw(width) << +matrix[0][r];
		for (Dimension c{1}; c < C; c++)
			ostream << _STD setw(width + 1) << +matrix[c][r];
		ostream << ']';
	};
}

namespace nd
{
	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	constexpr Tensor<S, N, Ns...>::Tensor() noexcept
	{
		m_Scalars.fill(Tensor<S, Ns...>{});
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Tensor<S, N, Ns...>::Tensor(S2 scalar) noexcept
	{
		for (Dimension n{}; n < N; ++n)
			m_Scalars[n] = static_cast<S>(scalar);
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2, Dimension N2, Dimension... N2s>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Tensor<S, N, Ns...>::Tensor(const Tensor<S2, N2, N2s...>& tensor) noexcept
	{
		for (Dimension n{}; n < _GCEM min(N, N2); ++n)
			Set(n, tensor);
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Tensor<S, N, Ns...>& Tensor<S, N, Ns...>::operator=(S2 scalar) noexcept
	{
		for (Dimension n{}; n < N; ++n)
			m_Scalars[n] = static_cast<S>(scalar);
		return *this;
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2, Dimension N2, Dimension... N2s>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Tensor<S, N, Ns...>& Tensor<S, N, Ns...>::operator=(const Tensor<S2, N2, N2s...>& tensor) noexcept
	{
		Dimension n{};
		for (; n < _GCEM min(N, N2); ++n)
			Set(n, tensor);
		if constexpr (N2 < N)
			for (; n < N; ++n)
				m_Scalars[n] = {};
		return *this;
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> Tensor<S, N, Ns...>::operator+(S2 scalar) const noexcept
	{
		return static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(*this) += static_cast<_IMPL CT<S, S2>>(scalar);
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Tensor<S, N, Ns...>& Tensor<S, N, Ns...>::operator+=(S2 scalar) noexcept
	{
		for (auto& s : m_Scalars)
			s += static_cast<S>(scalar);
		return *this;
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2, Dimension N2, Dimension... N2s>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> Tensor<S, N, Ns...>::operator+(const Tensor<S2, N2, N2s...>& tensor) const noexcept
	{
		return static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(*this) += static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(tensor);
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2, Dimension N2, Dimension... N2s>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Tensor<S, N, Ns...>& Tensor<S, N, Ns...>::operator+=(const Tensor<S2, N2, N2s...>& tensor) noexcept
	{
		Dimension n{};
		for (; n < _GCEM min(N, N2); ++n)
			m_Scalars[n] += tensor[n];
		return *this;
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> Tensor<S, N, Ns...>::operator-(S2 scalar) const noexcept
	{
		return static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(*this) -= static_cast<_IMPL CT<S, S2>>(scalar);
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Tensor<S, N, Ns...>& Tensor<S, N, Ns...>::operator-=(S2 scalar) noexcept
	{
		for (auto& s : m_Scalars)
			s -= static_cast<S>(scalar);
		return *this;
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2, Dimension N2, Dimension... N2s>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> Tensor<S, N, Ns...>::operator-(const Tensor<S2, N2, N2s...>& tensor) const noexcept
	{
		return static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(*this) -= static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(tensor);
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2, Dimension N2, Dimension... N2s>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Tensor<S, N, Ns...>& Tensor<S, N, Ns...>::operator-=(const Tensor<S2, N2, N2s...>& tensor) noexcept
	{
		Dimension n{};
		for (; n < _GCEM min(N, N2); ++n)
			m_Scalars[n] -= tensor[n];
		return *this;
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> Tensor<S, N, Ns...>::operator*(S2 scalar) const noexcept
	{
		return static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(*this) *= static_cast<_IMPL CT<S, S2>>(scalar);
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Tensor<S, N, Ns...>& Tensor<S, N, Ns...>::operator*=(S2 scalar) noexcept
	{
		for (auto& s : m_Scalars)
			s *= static_cast<S>(scalar);
		return *this;
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2, Dimension N2, Dimension... N2s>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> Tensor<S, N, Ns...>::operator*(const Tensor<S2, N2, N2s...>& tensor) const noexcept
	{
		return static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(*this) *= static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(tensor);
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2, Dimension N2, Dimension... N2s>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Tensor<S, N, Ns...>& Tensor<S, N, Ns...>::operator*=(const Tensor<S2, N2, N2s...>& tensor) noexcept
	{
		Dimension n{};
		for (; n < _GCEM min(N, N2); ++n)
			m_Scalars[n] *= tensor[n];
		return *this;
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> Tensor<S, N, Ns...>::operator/(S2 scalar) const noexcept
	{
		return static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(*this) /= static_cast<_IMPL CT<S, S2>>(scalar);
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Tensor<S, N, Ns...>& Tensor<S, N, Ns...>::operator/=(S2 scalar) noexcept
	{
		for (auto& s : m_Scalars)
			s /= static_cast<S>(scalar);
		return *this;
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2, Dimension N2, Dimension... N2s>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> Tensor<S, N, Ns...>::operator/(const Tensor<S2, N2, N2s...>& tensor) const noexcept
	{
		return static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(*this) /= static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(tensor);
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2, Dimension N2, Dimension... N2s>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Tensor<S, N, Ns...>& Tensor<S, N, Ns...>::operator/=(const Tensor<S2, N2, N2s...>& tensor) noexcept
	{
		Dimension n{};
		for (; n < _GCEM min(N, N2); ++n)
			m_Scalars[n] /= tensor[n];
		return *this;
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	constexpr Tensor<decltype(+S{}), N, Ns...> Tensor<S, N, Ns...>::operator+() const noexcept
	{
		Tensor<decltype(+S{}), N, Ns...> result;
		for (Dimension n{}; n < N; ++n)
			result[n] = +m_Scalars[n];
		return result;
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	constexpr Tensor<decltype(-S{}), N, Ns...> Tensor<S, N, Ns...>::operator-() const noexcept
	{
		Tensor<decltype(-S{}), N, Ns...> result;
		for (Dimension n{}; n < N; ++n)
			result[n] = -m_Scalars[n];
		return result;
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Dimension N2>
	requires(N2 < N)
	constexpr Tensor<S, Ns...>& Tensor<S, N, Ns...>::at() noexcept
	{
		return m_Scalars[N2];
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Dimension N2>
	requires(N2 < N)
	constexpr const Tensor<S, Ns...>& Tensor<S, N, Ns...>::at() const noexcept
	{
		return m_Scalars[N2];
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	constexpr Tensor<S, Ns...>& Tensor<S, N, Ns...>::operator[](Dimension n) noexcept
	{
		//__assume(n < N);
		return m_Scalars[n];
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	constexpr const Tensor<S, Ns...>& Tensor<S, N, Ns...>::operator[](Dimension n) const noexcept
	{
		//__assume(n < N);
		return m_Scalars[n];
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2>
	requires(!_STD is_same_v<S, S2> && _STD is_convertible_v<S, S2>)
	constexpr Tensor<S, N, Ns...>::operator Tensor<S2, N, Ns...>() const noexcept
	{
		Tensor<S2, N, Ns...> result;
		for (Dimension n{}; n < N; ++n)
			result.Set(n, *this);
		return result;
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2>
	constexpr bool Tensor<S, N, Ns...>::operator==(const Tensor<S2, N, Ns...>& tensor) const noexcept
	{
		if (this != (void*)&tensor)
			for (Dimension n{}; n < N; ++n)
				if (static_cast<Tensor<_IMPL CT<S, S2>, Ns...>>(m_Scalars[n]) != static_cast<Tensor<_IMPL CT<S, S2>, Ns...>>(tensor[n]))
					return false;
		return true;
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2>
	constexpr bool Tensor<S, N, Ns...>::operator!=(const Tensor<S2, N, Ns...>& tensor) const noexcept
	{
		return !(*this == tensor);
	}
	
	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <typename SHasher>
	constexpr size_t Tensor<S, N, Ns...>::Hasher<SHasher>::operator()(const Tensor<S, N, Ns...>& tensor) const noexcept
	{
		size_t hash{};
		for (Dimension n{}; n < N; ++n)
		{
			size_t tensorHash;
			if constexpr (Rank == 1)
				tensorHash = SHasher{}(tensor[n]);
			else
				tensorHash = typename Tensor<S, Ns...>::Hasher{}(tensor[n]);
			hash ^= (tensorHash << ((7 * n) & ((1ull << 6) - 1))) | ((tensorHash >> 56) & 0x7F);
		}
		return hash;
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	constexpr size_t Tensor<S, N, Ns...>::Hash() const noexcept
	{
		return Hasher{}(*this);
	}

	template <Scalar S, Dimension N, Dimension... Ns>
	requires((N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<N, Ns...>))
	template <Scalar S2, Dimension N2, Dimension... N2s>
	requires(_STD is_convertible_v<S2, S>)
	constexpr auto Tensor<S, N, Ns...>::Set(Dimension n, const Tensor<S2, N2, N2s...>& tensor) noexcept
	{
		if constexpr (Rank == 1)
			m_Scalars[n] = static_cast<S>(tensor[n]);
		else
			m_Scalars[n] = static_cast<Tensor<S, N2s...>>(tensor[n]);
	}

	// Tensor<S, 0>

	template <Scalar S>
	constexpr Tensor<S, 0>::operator S&() noexcept
	{
		return m_Scalar;
	}

	template <Scalar S>
	constexpr Tensor<S, 0>::operator const S&() const noexcept
	{
		return m_Scalar;
	}

	// External Operators

	template <Scalar S, Scalar S2, Dimension N, Dimension... Ns>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator+(S scalar, const Tensor<S2, N, Ns...>& tensor) noexcept
	{
		return Tensor<_IMPL CT<S, S2>, N, Ns...>{scalar} += tensor;
	}

	template <Scalar S, Scalar S2, Dimension N, Dimension... Ns>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator-(S scalar, const Tensor<S2, N, Ns...>& tensor) noexcept
	{
		return Tensor<_IMPL CT<S, S2>, N, Ns...>{scalar} -= tensor;
	}

	template <Scalar S, Scalar S2, Dimension N, Dimension... Ns>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator*(S scalar, const Tensor<S2, N, Ns...>& tensor) noexcept
	{
		return Tensor<_IMPL CT<S, S2>, N, Ns...>{scalar} *= tensor;
	}

	template <Scalar S, Scalar S2, Dimension N, Dimension... Ns>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator/(S scalar, const Tensor<S2, N, Ns...>& tensor) noexcept
	{
		return Tensor<_IMPL CT<S, S2>, N, Ns...>{scalar} /= tensor;
	}

	template <Scalar S>
	_STD ostream& operator<<(_STD ostream& ostream, const Tensor<S>& scalar)
	{
		return ostream << static_cast<S>(scalar);
	}

	template <Scalar S, Dimension N>
	_STD ostream& operator<<(_STD ostream& ostream, const Tensor<S, N>& vector)
	{
		ostream << '<' << vector[0];
		for (Dimension n{1}; n < N; ++n)
			ostream << ',' << vector[n];
		return ostream << '>';
	}

	template <_STD floating_point S, Dimension C, Dimension R>
	_STD ostream& operator<<(_STD ostream& ostream, const Tensor<S, C, R>& matrix)
	{
		ostream << _STD showpos << _STD scientific;
		_IMPL PrintRow(ostream, matrix, 13, 0);
		for (Dimension r{1}; r < R; ++r)
			_IMPL PrintRow(ostream << '\n', matrix, 13, r);
		return ostream << _STD defaultfloat << _STD noshowpos;
	}

	template <_STD integral S, Dimension C, Dimension R>
	_STD ostream& operator<<(_STD ostream& ostream, const Tensor<S, C, R>& matrix)
	{
		constexpr _STD streamsize width{_STD is_signed_v<S> + static_cast<_STD streamsize>(_GCEM max(
			S{_GCEM ceil(_GCEM log10(_STD make_unsigned_t<S>(_STD numeric_limits<S>::min())))},
			S{_GCEM ceil(_GCEM log10(_STD make_unsigned_t<S>(_STD numeric_limits<S>::max())))}
		))};

		_IMPL PrintRow(ostream, matrix, width, 0);
		for (Dimension r{1}; r < R; ++r)
			_IMPL PrintRow(ostream << '\n', matrix, width, r);
		return ostream;
	}

	template <Scalar S1, Scalar S2, Dimension N1, Dimension N2, Dimension... N1s, Dimension... N2s>
	constexpr Tensor<_IMPL CT<S1, S2>, N1, N1s..., N2, N2s...> OuterProduct(const Tensor<S1, N1, N1s...>& tensor1, const Tensor<S2, N2, N2s...>& tensor2) noexcept
	{
		Tensor<_IMPL CT<S1, S2>, N1, N1s..., N2, N2s...> result;
		[]<Dimension NR, Dimension NL1, Dimension NL2, Dimension... NRs, Dimension... NL1s, Dimension... NL2s>(this auto&& self, Tensor< _IMPL CT<S1, S2>, NR, NRs...>& resultL, const Tensor<_IMPL CT<S1, S2>, NL1, NL1s...>& tensorL1, const Tensor<_IMPL CT<S1, S2>, NL2, NL2s...>& tensorL2) -> void
		{
			if constexpr (_STD remove_cvref_t<decltype(tensorL1)>::Rank == 0)
				resultL = static_cast<_IMPL CT<S1, S2>>(tensorL1) * tensorL2;
			else for (Dimension n{}; n < NL1; ++n)
				self(resultL[n], tensorL1[n], tensorL2);
		}(result, static_cast<Tensor<_IMPL CT<S1, S2>, N1, N1s...>>(tensor1), static_cast<Tensor<_IMPL CT<S1, S2>, N2, N2s...>>(tensor2));
		return result;
	}

	template <Scalar S1, Scalar S2, Dimension N, Dimension... Ns>
	constexpr _IMPL CT<S1, S2> InnerProduct(const Tensor<S1, N, Ns...>& tensor1, const Tensor<S2, N, Ns...>& tensor2) noexcept
	{
		_IMPL CT<S1, S2> sum{};
		[&sum]<Dimension NL, Dimension... NLs>(this auto&& self, const Tensor<_IMPL CT<S1, S2>, NL, NLs...>& tensor) -> void
		{
			if constexpr (_STD remove_cvref_t<decltype(tensor)>::Rank == 0)
				sum += static_cast<_IMPL CT<S1, S2>>(tensor);
			else for (Dimension n{}; n < NL; ++n)
				self(tensor[n]);
		}(Tensor<_IMPL CT<S1, S2>, N, Ns...>{tensor1 * tensor2});
		return sum;
	}
}
