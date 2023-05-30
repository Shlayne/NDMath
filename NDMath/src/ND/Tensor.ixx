module;

#include <gcem.hpp>

export module nd.tensor;

import nd.types;
import nd.impl;
import <algorithm>;
import <array>;
import <iomanip>;
import <iostream>;

#define _ND ::nd::
#define _IMPL ::nd::impl::
#define _GCEM ::gcem::

export namespace nd
{
	template<Scalar S, Dimension R, Dimension N = 0, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	struct Tensor
	{
	public:
		constexpr Tensor() noexcept;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor(S2 scalar) noexcept;

		template<Scalar S2, Dimension N2 = 0, Dimension... N2s> requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor(const Tensor<S2, R, N2, N2s...>& tensor) noexcept;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator=(S2 scalar) noexcept -> Tensor<S, R, N, Ns...>&;

		template<Scalar S2, Dimension N2 = 0, Dimension... N2s> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator=(const Tensor<S2, R, N2, N2s...>& tensor) noexcept -> Tensor<S, R, N, Ns...>&;
	public:
		// Scalar addition.
		template<Scalar S2>
		constexpr auto operator+(S2 scalar) const noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>;

		// Scalar addition.
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator+=(S2 scalar) noexcept -> Tensor<S, R, N, Ns...>&;

		// Element-wise addition.
		template<Scalar S2, Dimension N2 = 0, Dimension... N2s>
		constexpr auto operator+(const Tensor<S2, R, N2, N2s...>& tensor) const noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>;

		// Element-wise addition.
		template<Scalar S2, Dimension N2 = 0, Dimension... N2s> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator+=(const Tensor<S2, R, N2, N2s...>& tensor) noexcept -> Tensor<S, R, N, Ns...>&;
	public:
		// Scalar subtraction.
		template<Scalar S2>
		constexpr auto operator-(S2 scalar) const noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>;

		// Scalar subtraction.
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator-=(S2 scalar) noexcept -> Tensor<S, R, N, Ns...>&;

		// Element-wise subtraction.
		template<Scalar S2, Dimension N2 = 0, Dimension... N2s>
		constexpr auto operator-(const Tensor<S2, R, N2, N2s...>& tensor) const noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>;

		// Element-wise subtraction.
		template<Scalar S2, Dimension N2 = 0, Dimension... N2s> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator-=(const Tensor<S2, R, N2, N2s...>& tensor) noexcept -> Tensor<S, R, N, Ns...>&;
	public:
		// Scalar multiplication.
		template<Scalar S2>
		constexpr auto operator*(S2 scalar) const noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>;

		// Scalar multiplication.
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator*=(S2 scalar) noexcept -> Tensor<S, R, N, Ns...>&;

		// Element-wise multiplication.
		template<Scalar S2, Dimension N2 = 0, Dimension... N2s>
		constexpr auto operator*(const Tensor<S2, R, N2, N2s...>& tensor) const noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>;

		// Element-wise multiplication.
		template<Scalar S2, Dimension N2 = 0, Dimension... N2s> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator*=(const Tensor<S2, R, N2, N2s...>& tensor) noexcept -> Tensor<S, R, N, Ns...>&;
	public:
		// Scalar division.
		template<Scalar S2>
		constexpr auto operator/(S2 scalar) const noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>;

		// Scalar division.
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator/=(S2 scalar) noexcept -> Tensor<S, R, N, Ns...>&;

		// Element-wise division.
		template<Scalar S2, Dimension N2 = 0, Dimension... N2s>
		constexpr auto operator/(const Tensor<S2, R, N2, N2s...>& tensor) const noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>;

		// Element-wise division.
		template<Scalar S2, Dimension N2 = 0, Dimension... N2s> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator/=(const Tensor<S2, R, N2, N2s...>& tensor) noexcept -> Tensor<S, R, N, Ns...>&;
	public:
		constexpr auto operator+() const noexcept -> Tensor<decltype(+S{}), R, N, Ns... > ;
		constexpr auto operator-() const noexcept -> Tensor<decltype(-S{}), R, N, Ns...>;
	public:
		template<Dimension N2> requires(N2 < N)
		constexpr auto at() noexcept -> Tensor<S, R - 1, Ns...>&;

		template<Dimension N2> requires(N2 < N)
		constexpr auto at() const noexcept -> const Tensor<S, R - 1, Ns...>&;

		constexpr auto operator[](Dimension n) noexcept -> Tensor<S, R - 1, Ns...>&;
		constexpr auto operator[](Dimension n) const noexcept -> const Tensor<S, R - 1, Ns...>&;
	public:
		template<Scalar S2> requires(_STD is_convertible_v<S, S2>)
		constexpr operator Tensor<S2, R, N, Ns...>() const noexcept;

		template<Scalar S2>
		constexpr auto operator==(const Tensor<S2, R, N, Ns...>& tensor) const noexcept -> bool;

		template<Scalar S2>
		constexpr auto operator!=(const Tensor<S2, R, N, Ns...>& tensor) const noexcept -> bool;
	public:
		template<typename SHasher = _STD hash<S>>
		struct Hasher
		{
			constexpr auto operator()(const Tensor<S, R, N, Ns...>& tensor) const noexcept -> size_t;
		};

		constexpr auto Hash() const noexcept -> size_t;
	private:
		template<Scalar S2, Dimension N2 = 0, Dimension... N2s> requires(_STD is_convertible_v<S2, S>)
		constexpr auto Set(Dimension n, const Tensor<S2, R, N2, N2s...>& tensor) noexcept;
	protected:
		_STD array<Tensor<S, R - 1, Ns...>, N> m_Scalars;
	};

	template<Scalar S>
	struct Tensor<S, 0>
	{
	public:
		constexpr Tensor() noexcept = default;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor(S2 scalar) noexcept;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator=(S2 scalar) noexcept -> Tensor<S, 0>&;
	public:
		constexpr operator S&() noexcept;
		constexpr operator const S&() const noexcept;
	private:
		S m_Scalar{};
	};

	// External Operators

	template<Scalar S, Scalar S2, Dimension R, Dimension N = 0, Dimension... Ns>
	constexpr auto operator+(S scalar, const Tensor<S2, R, N, Ns...>& tensor) noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>;

	template<Scalar S, Scalar S2, Dimension R, Dimension N = 0, Dimension... Ns>
	constexpr auto operator-(S scalar, const Tensor<S2, R, N, Ns...>& tensor) noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>;

	template<Scalar S, Scalar S2, Dimension R, Dimension N = 0, Dimension... Ns>
	constexpr auto operator*(S scalar, const Tensor<S2, R, N, Ns...>& tensor) noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>;

	template<Scalar S, Scalar S2, Dimension R, Dimension N = 0, Dimension... Ns>
	constexpr auto operator/(S scalar, const Tensor<S2, R, N, Ns...>& tensor) noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>;

	template<Scalar S, Dimension N>
	auto operator<<(_STD ostream& ostream, const Tensor<S, 1, N>& vector) -> _STD ostream&;

	template<Scalar S, Dimension C, Dimension R>
	auto operator<<(_STD ostream& ostream, const Tensor<S, 2, C, R>& matrix) -> _STD ostream&;

	// Static Methods

	template<Scalar S1, Scalar S2, Dimension R1, Dimension R2, Dimension N1 = 0, Dimension N2 = 0, Dimension... N1s, Dimension... N2s>
	constexpr auto OuterProduct(const Tensor<S1, R1, N1, N1s...>& tensor1, const Tensor<S2, R2, N2, N2s...>& tensor2) noexcept; // -> Tensor<_IMPL CT<S1, S2>, R1 + R2, N1, N1s..., N2, N2s...>;

	template<Scalar S1, Scalar S2, Dimension R, Dimension N = 0, Dimension... Ns>
	constexpr auto InnerProduct(const Tensor<S1, R, N, Ns...>& tensor1, const Tensor<S2, R, N, Ns...>& tensor2) noexcept -> _IMPL CT<S1, S2>;

	// Aliases
}

// Implementation: Don't export.
namespace nd::impl
{
	template<Scalar S, Dimension C, Dimension R>
	auto PrintRow(_STD ostream& ostream, const Tensor<S, 2, C, R>& matrix, _STD streamsize width, Dimension r) noexcept -> void
	{
		ostream << '[' << _STD setw(width) << +matrix[0][r];
		for (Dimension c{1}; c < C; c++)
			ostream << _STD setw(width + 1) << +matrix[c][r];
		ostream << ']';
	};
}

export namespace nd
{
	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	constexpr Tensor<S, R, N, Ns...>::Tensor() noexcept
	{
		m_Scalars.fill(Tensor<S, R - 1, Ns...>{});
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr Tensor<S, R, N, Ns...>::Tensor(S2 scalar) noexcept
	{
		for (Dimension n{}; n < N; ++n)
			m_Scalars[n] = static_cast<S>(scalar);
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2, Dimension N2, Dimension... N2s> requires(_STD is_convertible_v<S2, S>)
	constexpr Tensor<S, R, N, Ns...>::Tensor(const Tensor<S2, R, N2, N2s...>& tensor) noexcept
	{
		for (Dimension n{}; n < _GCEM min(N, N2); ++n)
			Set(n, tensor);
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Tensor<S, R, N, Ns...>::operator=(S2 scalar) noexcept -> Tensor<S, R, N, Ns...>&
	{
		for (Dimension n{}; n < N; ++n)
			m_Scalars[n] = static_cast<S>(scalar);
		return *this;
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2, Dimension N2, Dimension... N2s> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Tensor<S, R, N, Ns...>::operator=(const Tensor<S2, R, N2, N2s...>& tensor) noexcept -> Tensor<S, R, N, Ns...>&
	{
		Dimension n{};
		for (; n < _GCEM min(N, N2); ++n)
			Set(n, tensor);
		if constexpr (N2 < N)
			for (; n < N; ++n)
				m_Scalars[n] = {};
		return *this;
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2>
	constexpr auto Tensor<S, R, N, Ns...>::operator+(S2 scalar) const noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>
	{
		return static_cast<Tensor<_IMPL CT<S, S2>, R, N, Ns...>>(*this) += static_cast<_IMPL CT<S, S2>>(scalar);
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Tensor<S, R, N, Ns...>::operator+=(S2 scalar) noexcept -> Tensor<S, R, N, Ns...>&
	{
		for (auto& s : m_Scalars)
			s += static_cast<S>(scalar);
		return *this;
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2, Dimension N2, Dimension... N2s>
	constexpr auto Tensor<S, R, N, Ns...>::operator+(const Tensor<S2, R, N2, N2s...>& tensor) const noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>
	{
		return static_cast<Tensor<_IMPL CT<S, S2>, R, N, Ns...>>(*this) += static_cast<Tensor<_IMPL CT<S, S2>, R, N, Ns...>>(tensor);
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2, Dimension N2, Dimension... N2s> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Tensor<S, R, N, Ns...>::operator+=(const Tensor<S2, R, N2, N2s...>& tensor) noexcept -> Tensor<S, R, N, Ns...>&
	{
		Dimension n{};
		for (; n < _GCEM min(N, N2); ++n)
			m_Scalars[n] += tensor[n];
		return *this;
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2>
	constexpr auto Tensor<S, R, N, Ns...>::operator-(S2 scalar) const noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>
	{
		return static_cast<Tensor<_IMPL CT<S, S2>, R, N, Ns...>>(*this) -= static_cast<_IMPL CT<S, S2>>(scalar);
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Tensor<S, R, N, Ns...>::operator-=(S2 scalar) noexcept -> Tensor<S, R, N, Ns...>&
	{
		for (auto& s : m_Scalars)
			s -= static_cast<S>(scalar);
		return *this;
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2, Dimension N2, Dimension... N2s>
	constexpr auto Tensor<S, R, N, Ns...>::operator-(const Tensor<S2, R, N2, N2s...>& tensor) const noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>
	{
		return static_cast<Tensor<_IMPL CT<S, S2>, R, N, Ns...>>(*this) -= static_cast<Tensor<_IMPL CT<S, S2>, R, N, Ns...>>(tensor);
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2, Dimension N2, Dimension... N2s> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Tensor<S, R, N, Ns...>::operator-=(const Tensor<S2, R, N2, N2s...>& tensor) noexcept -> Tensor<S, R, N, Ns...>&
	{
		Dimension n{};
		for (; n < _GCEM min(N, N2); ++n)
			m_Scalars[n] -= tensor[n];
		return *this;
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2>
	constexpr auto Tensor<S, R, N, Ns...>::operator*(S2 scalar) const noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>
	{
		return static_cast<Tensor<_IMPL CT<S, S2>, R, N, Ns...>>(*this) *= static_cast<_IMPL CT<S, S2>>(scalar);
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Tensor<S, R, N, Ns...>::operator*=(S2 scalar) noexcept -> Tensor<S, R, N, Ns...>&
	{
		for (auto& s : m_Scalars)
			s *= static_cast<S>(scalar);
		return *this;
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2, Dimension N2, Dimension... N2s>
	constexpr auto Tensor<S, R, N, Ns...>::operator*(const Tensor<S2, R, N2, N2s...>& tensor) const noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>
	{
		return static_cast<Tensor<_IMPL CT<S, S2>, R, N, Ns...>>(*this) *= static_cast<Tensor<_IMPL CT<S, S2>, R, N, Ns...>>(tensor);
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2, Dimension N2, Dimension... N2s> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Tensor<S, R, N, Ns...>::operator*=(const Tensor<S2, R, N2, N2s...>& tensor) noexcept -> Tensor<S, R, N, Ns...>&
	{
		Dimension n{};
		for (; n < _GCEM min(N, N2); ++n)
			m_Scalars[n] *= tensor[n];
		return *this;
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2>
	constexpr auto Tensor<S, R, N, Ns...>::operator/(S2 scalar) const noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>
	{
		return static_cast<Tensor<_IMPL CT<S, S2>, R, N, Ns...>>(*this) /= static_cast<_IMPL CT<S, S2>>(scalar);
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Tensor<S, R, N, Ns...>::operator/=(S2 scalar) noexcept -> Tensor<S, R, N, Ns...>&
	{
		for (auto& s : m_Scalars)
			s /= static_cast<S>(scalar);
		return *this;
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2, Dimension N2, Dimension... N2s>
	constexpr auto Tensor<S, R, N, Ns...>::operator/(const Tensor<S2, R, N2, N2s...>& tensor) const noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>
	{
		return static_cast<Tensor<_IMPL CT<S, S2>, R, N, Ns...>>(*this) /= static_cast<Tensor<_IMPL CT<S, S2>, R, N, Ns...>>(tensor);
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2, Dimension N2, Dimension... N2s> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Tensor<S, R, N, Ns...>::operator/=(const Tensor<S2, R, N2, N2s...>& tensor) noexcept -> Tensor<S, R, N, Ns...>&
	{
		Dimension n{};
		for (; n < _GCEM min(N, N2); ++n)
			m_Scalars[n] /= tensor[n];
		return *this;
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	constexpr auto Tensor<S, R, N, Ns...>::operator+() const noexcept -> Tensor<decltype(+S{}), R, N, Ns...>
	{
		Tensor<decltype(+S{}), R, N, Ns...> result;
		for (Dimension n{}; n < N; ++n)
			result[n] = +m_Scalars[n];
		return result;
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	constexpr auto Tensor<S, R, N, Ns...>::operator-() const noexcept -> Tensor<decltype(-S{}), R, N, Ns...>
	{
		Tensor<decltype(-S{}), R, N, Ns...> result;
		for (Dimension n{}; n < N; ++n)
			result[n] = -m_Scalars[n];
		return result;
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Dimension N2> requires(N2 < N)
	constexpr auto Tensor<S, R, N, Ns...>::at() noexcept -> Tensor<S, R - 1, Ns...>&
	{
		return m_Scalars[N2];
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Dimension N2> requires(N2 < N)
	constexpr auto Tensor<S, R, N, Ns...>::at() const noexcept -> const Tensor<S, R - 1, Ns...>&
	{
		return m_Scalars[N2];
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	constexpr auto Tensor<S, R, N, Ns...>::operator[](Dimension n) noexcept -> Tensor<S, R - 1, Ns...>&
	{
		//__assume(n < N);
		return m_Scalars[n];
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	constexpr auto Tensor<S, R, N, Ns...>::operator[](Dimension n) const noexcept -> const Tensor<S, R - 1, Ns...>&
	{
		//__assume(n < N);
		return m_Scalars[n];
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2> requires(_STD is_convertible_v<S, S2>)
	constexpr Tensor<S, R, N, Ns...>::operator Tensor<S2, R, N, Ns...>() const noexcept
	{
		Tensor<S2, R, N, Ns...> result;
		for (Dimension n{}; n < N; ++n)
			result.Set(n, *this);
		return result;
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2>
	constexpr auto Tensor<S, R, N, Ns...>::operator==(const Tensor<S2, R, N, Ns...>& tensor) const noexcept -> bool
	{
		if (this != (void*)&tensor)
			for (Dimension n{}; n < N; ++n)
				if (static_cast<Tensor<_IMPL CT<S, S2>, R - 1, Ns...>>(m_Scalars[n]) != static_cast<Tensor<_IMPL CT<S, S2>, R - 1, Ns...>>(tensor[n]))
					return false;
		return true;
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2>
	constexpr auto Tensor<S, R, N, Ns...>::operator!=(const Tensor<S2, R, N, Ns...>& tensor) const noexcept -> bool
	{
		return !(*this == tensor);
	}
	
	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<typename SHasher>
	constexpr auto Tensor<S, R, N, Ns...>::Hasher<SHasher>::operator()(const Tensor<S, R, N, Ns...>& tensor) const noexcept -> size_t
	{
		size_t hash{};
		for (Dimension n{}; n < N; ++n)
		{
			size_t tensorHash;
			if constexpr (R == 1)
				tensorHash = SHasher{}(tensor[n]);
			else
				tensorHash = typename Tensor<S, R - 1, Ns...>::Hasher{}(tensor[n]);
			hash ^= (tensorHash << ((7 * n) & ((1ull << 6) - 1))) | ((tensorHash >> 56) & 0x7F);
		}
		return hash;
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	constexpr auto Tensor<S, R, N, Ns...>::Hash() const noexcept -> size_t
	{
		return Hasher{}(*this);
	}

	template<Scalar S, Dimension R, Dimension N, Dimension... Ns>
	requires((R == 0 && N == 0 && sizeof...(Ns) == 0) || (_IMPL gt0_v<R, N, Ns...> && sizeof...(Ns) == R - 1))
	template<Scalar S2, Dimension N2, Dimension... N2s> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Tensor<S, R, N, Ns...>::Set(Dimension n, const Tensor<S2, R, N2, N2s...>& tensor) noexcept
	{
		if constexpr (R == 1)
			m_Scalars[n] = static_cast<S>(tensor[n]);
		else
			m_Scalars[n] = static_cast<Tensor<S, R - 1, N2s...>>(tensor[n]);
	}

	// Tensor<S, 0>

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

	// External Operators

	template<Scalar S, Scalar S2, Dimension R, Dimension N, Dimension... Ns>
	constexpr auto operator+(S scalar, const Tensor<S2, R, N, Ns...>& tensor) noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>
	{
		return Tensor<_IMPL CT<S, S2>, R, N, Ns...>{scalar} += tensor;
	}

	template<Scalar S, Scalar S2, Dimension R, Dimension N, Dimension... Ns>
	constexpr auto operator-(S scalar, const Tensor<S2, R, N, Ns...>& tensor) noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>
	{
		return Tensor<_IMPL CT<S, S2>, R, N, Ns...>{scalar} -= tensor;
	}

	template<Scalar S, Scalar S2, Dimension R, Dimension N, Dimension... Ns>
	constexpr auto operator*(S scalar, const Tensor<S2, R, N, Ns...>& tensor) noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>
	{
		return Tensor<_IMPL CT<S, S2>, R, N, Ns...>{scalar} *= tensor;
	}

	template<Scalar S, Scalar S2, Dimension R, Dimension N, Dimension... Ns>
	constexpr auto operator/(S scalar, const Tensor<S2, R, N, Ns...>& tensor) noexcept -> Tensor<_IMPL CT<S, S2>, R, N, Ns...>
	{
		return Tensor<_IMPL CT<S, S2>, R, N, Ns...>{scalar} /= tensor;
	}

	template<Scalar S, Dimension N>
	auto operator<<(_STD ostream& ostream, const Tensor<S, 1, N>& vector) -> _STD ostream&
	{
		ostream << '<' << vector[0];
		for (Dimension n{1}; n < N; ++n)
			ostream << ',' << vector[n];
		return ostream << '>';
	}

	template<_STD floating_point S, Dimension C, Dimension R>
	auto operator<<(_STD ostream& ostream, const Tensor<S, 2, C, R>& matrix) -> _STD ostream&
	{
		ostream << _STD showpos << _STD scientific;
		_IMPL PrintRow(ostream, matrix, 13, 0);
		for (Dimension r{1}; r < R; ++r)
			_IMPL PrintRow(ostream << '\n', matrix, 13, r);
		return ostream << _STD defaultfloat << _STD noshowpos;
	}

	template<_STD integral S, Dimension C, Dimension R>
	auto operator<<(_STD ostream& ostream, const Tensor<S, 2, C, R>& matrix) -> _STD ostream&
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

	template<Scalar S1, Scalar S2, Dimension R1, Dimension R2, Dimension N1, Dimension N2, Dimension... N1s, Dimension... N2s>
	constexpr auto OuterProduct(const Tensor<S1, R1, N1, N1s...>& tensor1, const Tensor<S2, R2, N2, N2s...>& tensor2) noexcept
	{
		if constexpr (R1 == 0)
			return static_cast<_IMPL CT<S1, S2>>(tensor1) * tensor2;
		else
		{
			Tensor<_IMPL CT<S1, S2>, R1 + R2, N1, N1s..., N2, N2s...> result;
			for (Dimension n{}; n < N1; ++n)
				result[n] = OuterProduct(tensor1[n], tensor2);
			return result;
		}
	}

	template<Scalar S1, Scalar S2, Dimension R, Dimension N, Dimension... Ns>
	constexpr auto InnerProduct(const Tensor<S1, R, N, Ns...>& tensor1, const Tensor<S2, R, N, Ns...>& tensor2) noexcept -> _IMPL CT<S1, S2>
	{
		_IMPL CT<S1, S2> sum{};
		[&sum]<Dimension RL, Dimension NL, Dimension... NLs>(this auto&& self, const Tensor<_IMPL CT<S1, S2>, RL, NL, NLs...>& tensor) -> void
		{
			if constexpr (RL == 0)
				sum += static_cast<_IMPL CT<S1, S2>>(tensor);
			else for (Dimension n{}; n < NL; ++n)
				self(tensor[n]);
		}(Tensor<_IMPL CT<S1, S2>, R, N, Ns...>{tensor1 * tensor2});
		return sum;
	}

	// Tensor<S, 0>
	// NOTE: intellisense doesn't think these two are valid even though they are.
	// Keep them at the bottom of this file, even though they're out of order.

	template<Scalar S>
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr Tensor<S, 0>::Tensor(S2 scalar) noexcept
		: m_Scalar{static_cast<S>(scalar)}
	{

	}

	template<Scalar S>
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Tensor<S, 0>::operator=(S2 scalar) noexcept -> Tensor<S, 0>&
	{
		m_Scalar = static_cast<S>(scalar);
		return *this;
	}
}
