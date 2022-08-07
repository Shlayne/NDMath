#pragma once

#include "ND/Vector.h"
#include <iomanip>

namespace nd
{
	template<typename T, uint32_t C, uint32_t R = C>
	struct Matrix
	{
	public:
		static_assert(C > 0 && R > 0, "Matrix must have at least one element.");
	public:
		constexpr Matrix() noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Matrix(const T2& crValue) noexcept;

		template<typename... Args, std::enable_if_t<std::ext::are_convertible_v<T, Args...> && sizeof...(Args) == C * R && C * R != 1, int> = 0>
		constexpr Matrix(Args&&... rrArgs) noexcept;

		template<uint32_t C2, uint32_t R2, typename T2, std::enable_if_t<std::is_convertible_v<T2, T> && (C2 <= C) && (R2 <= R), int> = 0>
		constexpr Matrix(const Matrix<T2, C2, R2>& crMatrix) noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Matrix<T, C, R>& operator=(const T2& crValue) noexcept;

		template<uint32_t C2, uint32_t R2, typename T2, std::enable_if_t<std::is_convertible_v<T2, T> && (C2 <= C) && (R2 <= R), int> = 0>
		constexpr Matrix& operator=(const Matrix<T2, C2, R2>& crMatrix) noexcept;
	public:
		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Matrix<T, C, R> operator+(const T2& crValue) const noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Matrix<T, C, R>& operator+=(const T2& crValue) noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Matrix<T, C, R> operator+(const Matrix<T2, C, R>& crMatrix) const noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Matrix<T, C, R>& operator+=(const Matrix<T2, C, R>& crMatrix) noexcept;
	public:
		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Matrix<T, C, R> operator-(const T2& crValue) const noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Matrix<T, C, R>& operator-=(const T2& crValue) noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Matrix<T, C, R> operator-(const Matrix<T2, C, R>& crMatrix) const noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Matrix<T, C, R>& operator-=(const Matrix<T2, C, R>& crMatrix) noexcept;
	public:
		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Matrix<T, C, R> operator*(const T2& crValue) const noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Matrix<T, C, R>& operator*=(const T2& crValue) noexcept;

		template<typename T2, uint32_t C2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Matrix<T, C2, R> operator*(const Matrix<T2, C2, C>& crMatrix) const noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T> && C == R, int> = 0>
		constexpr Matrix<T, C, R>& operator*=(const Matrix<T2, C, R>& crMatrix) noexcept;
	public:
		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Matrix<T, C, R> operator/(const T2& crValue) const noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Matrix<T, C, R>& operator/=(const T2& crValue) noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Matrix<T, C, R> operator/(const Matrix<T2, C, C>& crMatrix) const noexcept;

		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T> && C == R, int> = 0>
		constexpr Matrix<T, C, R>& operator/=(const Matrix<T2, C, C>& crMatrix) noexcept;
	public:
		template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
		constexpr Vector<T, R> operator*(const Vector<T2, C>& crVector) const noexcept;
	public:
		constexpr Matrix<T, C, R> operator+() const noexcept;
		constexpr Matrix<T, C, R> operator-() const noexcept;
	public:
		constexpr Vector<T, R>& operator[](uint32_t col) noexcept;
		constexpr const Vector<T, R>& operator[](uint32_t col) const noexcept;

		template<uint32_t C2, std::enable_if_t<(C2 < C), int> = 0>
		constexpr Vector<T, R> Col() const noexcept;

		template<uint32_t R2, std::enable_if_t<(R2 < R), int> = 0>
		constexpr Vector<T, C> Row() const noexcept;
	private:
		template<uint32_t C2, uint32_t R2, typename T2, typename... Args>
		constexpr void Fill(const T2& crParam, Args&&... rrArgs);
	private:
		Vector<T, R> m_pElements[C]{ Vector<T, R>() };
	private:
		template<typename T2, uint32_t C2, uint32_t R2>
		friend struct Matrix;
	};

	// External Operators

	template<typename T, uint32_t C, uint32_t R, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Matrix<T, C, R> operator+(const T2& crValue, const Matrix<T, C, R>& crMatrix) noexcept;

	template<typename T, uint32_t C, uint32_t R, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Matrix<T, C, R> operator-(const T2& crValue, const Matrix<T, C, R>& crMatrix) noexcept;

	template<typename T, uint32_t C, uint32_t R, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Matrix<T, C, R> operator*(const T2& crValue, const Matrix<T, C, R>& crMatrix) noexcept;

	template<typename T, uint32_t CR, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Matrix<T, CR> operator/(const T2& crValue, const Matrix<T, CR>& crMatrix) noexcept;

	template<typename T, uint32_t C, uint32_t R, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Vector<T, C> operator*(const Vector<T2, R>& crVector, const Matrix<T, C, R>& crMatrix) noexcept;

	// Static Methods

	template<uint32_t C, uint32_t R, typename T, uint32_t C2, uint32_t R2, std::enable_if_t<(C < C2) && (R < R2) && (C2 > 1) && (R2 > 1), int> = 0>
	constexpr Matrix<T, C2 - 1, R2 - 1> Submatrix(const Matrix<T, C2, R2>& crMatrix) noexcept;

	template<typename T, uint32_t C2, uint32_t R2, std::enable_if_t<(C2 > 1) && (R2 > 1), int> = 0>
	constexpr Matrix<T, C2 - 1, R2 - 1> Submatrix(uint32_t C, uint32_t R, const Matrix<T, C2, R2>& crMatrix) noexcept;

	template<typename T, uint32_t CR>
	constexpr T Trace(const Matrix<T, CR>& crMatrix) noexcept;

	template<typename T, uint32_t CR>
	constexpr T Determinant(const Matrix<T, CR>& crMatrix) noexcept;

	template<typename T, uint32_t C, uint32_t R>
	constexpr Matrix<T, R, C> Transpose(const Matrix<T, C, R>& crMatrix) noexcept;

	template<typename T, uint32_t CR>
	constexpr Matrix<T, CR> Minors(const Matrix<T, CR>& crMatrix) noexcept;

	template<typename T, uint32_t CR>
	constexpr Matrix<T, CR> Cofactors(const Matrix<T, CR>& crMatrix) noexcept;

	template<typename T, uint32_t CR>
	constexpr Matrix<T, CR> Adjugate(const Matrix<T, CR>& crMatrix) noexcept;

	template<typename T, uint32_t CR>
	constexpr Matrix<T, CR> Inverse(const Matrix<T, CR>& crMatrix) noexcept;

	template<typename T, uint32_t CR, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Matrix<T, CR> Inverse(const Matrix<T, CR>& crMatrix, const T2& crDeterminant) noexcept;

	template<typename T, uint32_t CR, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Matrix<T, CR> Translate(const Matrix<T, CR>& crMatrix, const Vector<T2, CR - 1>& crTranslation);

	template<uint32_t A1, uint32_t A2, typename T, uint32_t CR, typename T2, std::enable_if_t<std::is_convertible_v<T2, T> && (A1 < CR) && (A2 < CR) && A1 != A2, int> = 0>
	constexpr Matrix<T, CR> Rotate(const Matrix<T, CR, CR>& crMatrix, const T2& crRadians);

	template<typename T, uint32_t CR, typename T2, std::enable_if_t<std::is_convertible_v<T2, T> && (CR >= 2), int> = 0>
	constexpr Matrix<T, CR> Rotate(uint32_t A1, uint32_t A2, const Matrix<T, CR>& crMatrix, const T2& crRadians);

	template<typename T, uint32_t CR, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Matrix<T, CR> Scale(const Matrix<T, CR>& crMatrix, const Vector<T2, CR - 1>& crScale);

	template<typename T, uint32_t CR, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int> = 0>
	constexpr Matrix<T, CR> Scale(const Matrix<T, CR>& crMatrix, const T2& crScale);

	template<typename T, uint32_t C, uint32_t R>
	std::ostream& operator<<(std::ostream& rOstream, const Matrix<T, C, R>& crMatrix);

	// Aliases

	template<uint32_t C, uint32_t R = C>
	using MatrixCxRf = Matrix<float, C, R>;

	using Matrix1f = MatrixCxRf<1>;
	using Matrix2f = MatrixCxRf<2>;
	using Matrix3f = MatrixCxRf<3>;
	using Matrix4f = MatrixCxRf<4>;
	using Matrix5f = MatrixCxRf<5>;
	using Matrix6f = MatrixCxRf<6>;
	using Matrix7f = MatrixCxRf<7>;
	using Matrix8f = MatrixCxRf<8>;
}

#include "Matrix.inl"
