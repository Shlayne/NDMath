#pragma once

#include "ExtendedSTL.h"
#include <gcem.hpp>
#include <algorithm>
#include <ostream>
#include <utility>

namespace nd
{
	template<typename T, uint32_t N>
	struct Vector
	{
	public:
		static_assert(N > 0, "Vector must have at least one dimension.");
	public:
		constexpr Vector() noexcept = default;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Vector(const T2& crValue) noexcept;

		template<typename... Args, std::ext::disable_if_t<std::ext::are_convertible_v<T, Args...> && sizeof...(Args) == 1, int> = 0>
		constexpr Vector(Args&&... rrArgs) noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Vector<T, N>& operator=(const T2& crValue) noexcept;

		template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T> && N2 <= N, int> = 0>
		constexpr Vector<T, N>& operator=(const Vector<T2, N2>& crVector) noexcept;
	public:
		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Vector<T, N> operator+(const T2& crValue) const noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Vector<T, N>& operator+=(const T2& crValue) noexcept;

		template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T> && N2 <= N, int> = 0>
		constexpr Vector<T, N> operator+(const Vector<T2, N2>& crVector) const noexcept;

		template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T> && N2 <= N, int> = 0>
		constexpr Vector<T, N>& operator+=(const Vector<T2, N2>& crVector) noexcept;
	public:
		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Vector<T, N> operator-(const T2& crValue) const noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Vector<T, N>& operator-=(const T2& crValue) noexcept;

		template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T> && N2 <= N, int> = 0>
		constexpr Vector<T, N> operator-(const Vector<T2, N2>& crVector) const noexcept;

		template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T> && N2 <= N, int> = 0>
		constexpr Vector<T, N>& operator-=(const Vector<T2, N2>& crVector) noexcept;
	public:
		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Vector<T, N> operator*(const T2& crValue) const noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Vector<T, N>& operator*=(const T2& crValue) noexcept;

		template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T> && N2 <= N, int> = 0>
		constexpr Vector<T, N> operator*(const Vector<T2, N2>& crVector) const noexcept;

		template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T> && N2 <= N, int> = 0>
		constexpr Vector<T, N>& operator*=(const Vector<T2, N2>& crVector) noexcept;
	public:
		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Vector<T, N> operator/(const T2& crValue) const noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Vector<T, N>& operator/=(const T2& crValue) noexcept;

		template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T> && N2 <= N, int> = 0>
		constexpr Vector<T, N> operator/(const Vector<T2, N2>& crVector) const noexcept;

		template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T> && N2 <= N, int> = 0>
		constexpr Vector<T, N>& operator/=(const Vector<T2, N2>& crVector) noexcept;
	public:
		constexpr Vector<T, N> operator+() const noexcept;
		constexpr Vector<T, N> operator-() const noexcept;
	public:
		constexpr T& operator[](uint32_t index) noexcept;
		constexpr const T& operator[](uint32_t index) const noexcept;
	public:
		template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T, T2> && (N2 <= N), int> = 0>
		constexpr operator Vector<T2, N2>() const noexcept;
	private:
		template<uint32_t N2, typename T2, typename... Args, std::enable_if_t<std::is_convertible_v<T2, T> && (N2 >= 1), int> = 0>
		constexpr void Fill(const T2& crValue, Args&&... rrArgs) noexcept;

		template<uint32_t N2, typename T2, uint32_t N3, typename... Args, std::enable_if_t<std::is_convertible_v<T2, T> && (N2 >= N3), int> = 0>
		constexpr void Fill(const Vector<T2, N3>& crVector, Args&&... rrArgs) noexcept;
	private:
		T m_pElements[N]{ T(0) };
	};

	// External Operators

	template<typename T, uint32_t N, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Vector<T, N> operator+(const T2& crValue, const Vector<T, N>& crVector) noexcept;

	template<typename T, uint32_t N, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Vector<T, N> operator-(const T2& crValue, const Vector<T, N>& crVector) noexcept;

	template<typename T, uint32_t N, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Vector<T, N> operator*(const T2& crValue, const Vector<T, N>& crVector) noexcept;

	template<typename T, uint32_t N, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Vector<T, N> operator/(const T2& crValue, const Vector<T, N>& crVector) noexcept;

	// Static Methods

	template<typename T, uint32_t N, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Vector<T, N> Dot(const Vector<T, N>& crVector1, const Vector<T2, N>& crVector2) noexcept;

	template<typename T, uint32_t N>
	constexpr T Length2(const Vector<T, N>& crVector) noexcept;

	template<typename T, uint32_t N>
	constexpr T Length(const Vector<T, N>& crVector) noexcept;

	template<typename T, uint32_t N>
	constexpr Vector<T, N> Normalize(const Vector<T, N>& crVector) noexcept;

	template<typename T, uint32_t N, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Vector<T, N> Min(const Vector<T, N>& crVector1, const Vector<T2, N>& crVector2) noexcept;

	template<typename T, uint32_t N, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Vector<T, N> Max(const Vector<T, N>& crVector1, const Vector<T2, N>& crVector2) noexcept;

	template<typename T, uint32_t N>
	constexpr Vector<T, N> Abs(const Vector<T, N>& crVector1) noexcept;

	template<typename T, uint32_t N>
	constexpr Vector<T, N> Floor(const Vector<T, N>& crVector1) noexcept;

	template<typename T, uint32_t N>
	constexpr Vector<T, N> Ceil(const Vector<T, N>& crVector1) noexcept;

	template<typename T, uint32_t N, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Vector<T, N> Mod(const Vector<T, N>& crVector1, const Vector<T2, N>& crVector2) noexcept;

	template<typename T, uint32_t N, typename T2, typename T3, std::enable_if_t<std::ext::are_convertible_v<T, T2, T3>, int> = 0>
	constexpr Vector<T, N> Clamp(const Vector<T, N>& crVector, const Vector<T2, N>& crMin, const Vector<T3, N>& crMax) noexcept;

	template<typename T, uint32_t N>
	constexpr Vector<T, N> Round(const Vector<T, N>& crVector1) noexcept;

	template<typename T, uint32_t N, typename T2, typename T3, std::enable_if_t<std::ext::are_convertible_v<T, T2, T3>, int> = 0>
	constexpr Vector<T, N> Lerp(const Vector<T, N>& crVector1, const Vector<T2, N>& crVector2, const T3& crTValue) noexcept;

	template<typename T, uint32_t N>
	std::ostream& operator<<(std::ostream& rOstream, const Vector<T, N>& crVector);

	// Aliases

	template<uint32_t N>
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

#include "Vector.inl"
