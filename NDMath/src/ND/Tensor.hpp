#pragma once

#include "Impl.hpp"
#include <gcem.hpp>
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <utility>

#define _ND ::nd::
#define _IMPL ::nd::impl::
#define _GCEM ::gcem::

namespace nd
{
	template <Scalar S, Dimension N = 1, Dimension... Ns>
	requires(_IMPL gt_v<0, N, Ns...>)
	struct Tensor
	{
	public:
		static constexpr Dimension Rank = static_cast<Dimension>(N > 0 ? 1 : 0) + sizeof...(Ns);
		using SubTensor = _STD conditional_t<(Rank >= 2), Tensor<S, Ns...>, S>;

		constexpr Tensor() noexcept = default;

		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor(S2 scalar) noexcept
		{
			for (Dimension n{}; n < N; ++n)
				m_Scalars[n] = static_cast<S>(scalar);
		}

		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		requires(_STD is_convertible_v<S2, S> && N != N2 && (... && (Ns != N2s)))
		constexpr Tensor(const Tensor<S2, N2, N2s...>& tensor) noexcept
		{
			for (Dimension n{}; n < _GCEM min(N, N2); ++n)
				Set(n, tensor);
		}

		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator=(S2 scalar) noexcept
		{
			for (Dimension n{}; n < N; ++n)
				m_Scalars[n] = static_cast<S>(scalar);
			return *this;
		}

		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		requires(_STD is_convertible_v<S2, S> && N != N2 && (... && (Ns != N2s)))
		constexpr Tensor<S, N, Ns...>& operator=(const Tensor<S2, N2, N2s...>& tensor) noexcept
		{
			Dimension n{};
			for (; n < _GCEM min(N, N2); ++n)
				Set(n, tensor);
			if constexpr (N2 < N)
				for (; n < N; ++n)
					m_Scalars[n] = {};
			return *this;
		}
	public:
		// Scalar addition.
		template <Scalar S2>
		constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator+(S2 scalar) const noexcept
		{
			return static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(*this) += static_cast<_IMPL CT<S, S2>>(scalar);
		}

		// Scalar addition.
		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator+=(S2 scalar) noexcept
		{
			for (auto& s : m_Scalars)
				s += static_cast<S>(scalar);
			return *this;
		}

		// Element-wise addition.
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator+(const Tensor<S2, N2, N2s...>& tensor) const noexcept
		{
			return static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(*this) += static_cast<Tensor<_IMPL CT<S, S2>, N2, N2s...>>(tensor);
		}

		// Element-wise addition.
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator+=(const Tensor<S2, N2, N2s...>& tensor) noexcept
		{
			Dimension n{};
			for (; n < _GCEM min(N, N2); ++n)
				m_Scalars[n] += tensor[n];
			return *this;
		}
	public:
		// Scalar subtraction.
		template <Scalar S2>
		constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator-(S2 scalar) const noexcept
		{
			return static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(*this) -= static_cast<_IMPL CT<S, S2>>(scalar);
		}

		// Scalar subtraction.
		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator-=(S2 scalar) noexcept
		{
			for (auto& s : m_Scalars)
				s -= static_cast<S>(scalar);
			return *this;
		}

		// Element-wise subtraction.
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator-(const Tensor<S2, N2, N2s...>& tensor) const noexcept
		{
			return static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(*this) -= static_cast<Tensor<_IMPL CT<S, S2>, N2, N2s...>>(tensor);
		}

		// Element-wise subtraction.
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator-=(const Tensor<S2, N2, N2s...>& tensor) noexcept
		{
			Dimension n{};
			for (; n < _GCEM min(N, N2); ++n)
				m_Scalars[n] -= tensor[n];
			return *this;
		}
	public:
		// Scalar multiplication.
		template <Scalar S2>
		constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator*(S2 scalar) const noexcept
		{
			return static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(*this) *= static_cast<_IMPL CT<S, S2>>(scalar);
		}

		// Scalar multiplication.
		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator*=(S2 scalar) noexcept
		{
			for (auto& s : m_Scalars)
				s *= static_cast<S>(scalar);
			return *this;
		}

		// Element-wise multiplication.
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator*(const Tensor<S2, N2, N2s...>& tensor) const noexcept
		{
			return static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(*this) *= static_cast<Tensor<_IMPL CT<S, S2>, N2, N2s...>>(tensor);
		}

		// Element-wise multiplication.
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator*=(const Tensor<S2, N2, N2s...>& tensor) noexcept
		{
			Dimension n{};
			for (; n < _GCEM min(N, N2); ++n)
				m_Scalars[n] *= tensor[n];
			return *this;
		}
	public:
		// Scalar division.
		template <Scalar S2>
		constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator/(S2 scalar) const noexcept
		{
			return static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(*this) /= static_cast<_IMPL CT<S, S2>>(scalar);
		}

		// Scalar division.
		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator/=(S2 scalar) noexcept
		{
			for (auto& s : m_Scalars)
				s /= static_cast<S>(scalar);
			return *this;
		}

		// Element-wise division.
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator/(const Tensor<S2, N2, N2s...>& tensor) const noexcept
		{
			return static_cast<Tensor<_IMPL CT<S, S2>, N, Ns...>>(*this) /= static_cast<Tensor<_IMPL CT<S, S2>, N2, N2s...>>(tensor);
		}

		// Element-wise division.
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Tensor<S, N, Ns...>& operator/=(const Tensor<S2, N2, N2s...>& tensor) noexcept
		{
			Dimension n{};
			for (; n < _GCEM min(N, N2); ++n)
				m_Scalars[n] /= tensor[n];
			return *this;
		}
	public:
		constexpr Tensor<decltype(+S{}), N, Ns... > operator+() const noexcept
		{
			Tensor<decltype(+S{}), N, Ns... > result;
			for (Dimension n{}; n < N; ++n)
				result[n] = +m_Scalars[n];
			return result;
		}

		constexpr Tensor<decltype(-S{}), N, Ns...> operator-() const noexcept
		{
			Tensor<decltype(-S{}), N, Ns... > result;
			for (Dimension n{}; n < N; ++n)
				result[n] = -m_Scalars[n];
			return result;
		}
	public:
		template <Dimension N2>
		requires(N2 < N)
		constexpr SubTensor& at() noexcept
		{
			return m_Scalars[N2];
		}

		template <Dimension N2>
		requires(N2 < N)
		constexpr const SubTensor& at() const noexcept
		{
			return m_Scalars[N2];
		}

		constexpr SubTensor& operator[](Dimension n) noexcept
		{
			//__assume(n < N);
			return m_Scalars[n];
		}

		constexpr const SubTensor& operator[](Dimension n) const noexcept
		{
			//__assume(n < N);
			return m_Scalars[n];
		}
	public:
		template <Scalar S2>
		requires(!_STD is_same_v<S, S2> && _STD is_convertible_v<S, S2>)
		constexpr operator Tensor<S2, N, Ns...>() const noexcept
		{
			Tensor<S2, N, Ns...> result;
			for (Dimension n{}; n < N; ++n)
				result.Set(n, *this);
			return result;
		}

		template <Scalar S2>
		constexpr bool operator==(const Tensor<S2, N, Ns...>& tensor) const noexcept
		{
			if (this != (void*)&tensor)
				for (Dimension n{}; n < N; ++n)
					if (static_cast<Tensor<_IMPL CT<S, S2>, Ns...>>(m_Scalars[n]) != static_cast<Tensor<_IMPL CT<S, S2>, Ns...>>(tensor[n]))
						return false;
			return true;
		}

		template <Scalar S2>
		constexpr bool operator!=(const Tensor<S2, N, Ns...>& tensor) const noexcept
		{
			return !(*this == tensor);
		}
	public:
		template <typename SHasher = _STD hash<S>>
		struct Hasher
		{
			constexpr size_t operator()(const Tensor<S, N, Ns...>& tensor) const noexcept
			{
				size_t hash{};
				for (Dimension n{}; n < N; ++n)
				{
					size_t tensorHash;
					if constexpr (Rank == 1)
						tensorHash = SHasher{}(tensor[n]);
					else
						tensorHash = typename SubTensor::Hasher{}(tensor[n]);
					hash ^= (tensorHash << ((7 * n) & ((1ull << 6) - 1))) | ((tensorHash >> 56) & 0x7F);
				}
				return hash;
			}
		};

		constexpr size_t Hash() const noexcept
		{
			return Hasher{}(*this);
		}
	private:
		template <Scalar S2, Dimension N2 = 0, Dimension... N2s>
		requires(_STD is_convertible_v<S2, S>)
		constexpr auto Set(Dimension n, const Tensor<S2, N2, N2s...>& tensor) noexcept
		{
			if constexpr (Rank == 1)
				m_Scalars[n] = static_cast<S>(tensor[n]);
			else
				m_Scalars[n] = static_cast<Tensor<S, N2s...>>(tensor[n]);
		}
	protected:
		_STD array<SubTensor, N> m_Scalars{};

		template <Scalar S2, Dimension N2, Dimension... N2s>
		requires(_IMPL gt_v<0, N2, N2s...>)
		friend struct Tensor;
	};

	// External Operators

	template <Scalar S, Scalar S2, Dimension N = 0, Dimension... Ns>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator+(S scalar, const Tensor<S2, N, Ns...>& tensor) noexcept
	{
		return Tensor<_IMPL CT<S, S2>, N, Ns...>{scalar} += tensor;
	}

	template <Scalar S, Scalar S2, Dimension N = 0, Dimension... Ns>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator-(S scalar, const Tensor<S2, N, Ns...>& tensor) noexcept
	{
		return Tensor<_IMPL CT<S, S2>, N, Ns...>{scalar} -= tensor;
	}

	template <Scalar S, Scalar S2, Dimension N = 0, Dimension... Ns>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator*(S scalar, const Tensor<S2, N, Ns...>& tensor) noexcept
	{
		return Tensor<_IMPL CT<S, S2>, N, Ns...>{scalar} *= tensor;
	}

	template <Scalar S, Scalar S2, Dimension N = 0, Dimension... Ns>
	constexpr Tensor<_IMPL CT<S, S2>, N, Ns...> operator/(S scalar, const Tensor<S2, N, Ns...>& tensor) noexcept
	{
		return Tensor<_IMPL CT<S, S2>, N, Ns...>{scalar} /= tensor;
	}

	// Static Methods

	template <Scalar S1, Scalar S2, Dimension N1 = 0, Dimension N2 = 0, Dimension... N1s, Dimension... N2s>
	constexpr Tensor<_IMPL CT<S1, S2>, N1, N1s..., N2, N2s...> OuterProduct(const Tensor<S1, N1, N1s...>& tensor1, const Tensor<S2, N2, N2s...>& tensor2) noexcept
	{
		Tensor<_IMPL CT<S1, S2>, N1, N1s..., N2, N2s...> result;
		[]<Dimension NR, Dimension NL1, Dimension NL2, Dimension... NRs, Dimension... NL1s, Dimension... NL2s>(this auto&& self, Tensor< _IMPL CT<S1, S2>, NR, NRs...>& resultL, const Tensor<_IMPL CT<S1, S2>, NL1, NL1s...>& tensorL1, const Tensor<_IMPL CT<S1, S2>, NL2, NL2s...>& tensorL2) -> void
		{
			for (Dimension n{}; n < NL1; ++n)
			{
				if constexpr (_STD remove_cvref_t<decltype(tensorL1)>::Rank == 1)
					resultL[n] = static_cast<_IMPL CT<S1, S2>>(tensorL1[n]) * tensorL2;
				else 
					self(resultL[n], tensorL1[n], tensorL2);
			}
		}(result, static_cast<Tensor<_IMPL CT<S1, S2>, N1, N1s...>>(tensor1), static_cast<Tensor<_IMPL CT<S1, S2>, N2, N2s...>>(tensor2));
		return result;
	}

	template <Scalar S1, Scalar S2, Dimension N = 0, Dimension... Ns>
	constexpr _IMPL CT<S1, S2> InnerProduct(const Tensor<S1, N, Ns...>& tensor1, const Tensor<S2, N, Ns...>& tensor2) noexcept
	{
		_IMPL CT<S1, S2> sum{};
		[&sum]<Dimension NL, Dimension... NLs>(this auto&& self, const Tensor<_IMPL CT<S1, S2>, NL, NLs...>& tensor) -> void
		{
			for (Dimension n{}; n < NL; ++n)
			{
				if constexpr (_STD remove_cvref_t<decltype(tensor)>::Rank == 1)
					sum += static_cast<_IMPL CT<S1, S2>>(tensor[n]);
				else
					self(tensor[n]);
			}
		}(Tensor<_IMPL CT<S1, S2>, N, Ns...>{tensor1 * tensor2});
		return sum;
	}
}
