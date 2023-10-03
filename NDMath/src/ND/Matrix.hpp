#pragma once

#include "Impl.hpp"
#include "Tensor.hpp"
#include "Vector.hpp"
#include <gcem.hpp>

#define _ND ::nd::
#define _IMPL ::nd::impl::
#define _GCEM ::gcem::

namespace nd
{
	// Columns are the first index, rows are second.
	template <Scalar S, Dimension C, Dimension R = C>
	struct Matrix : public Tensor<S, 2, C, R>
	{
	public:
		// Constructs an identity matrix.
		constexpr Matrix() noexcept;

		// Constructs an identity matrix scaled by the given scalar.
		template <Scalar S2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Matrix(S2 scalar) noexcept;

		// Copies as many elements as fit into this matrix from the given matrix.
		template <Scalar S2, Dimension C2, Dimension R2>
		requires(_STD is_convertible_v<S2, S>)
		constexpr Matrix(const Matrix<S2, C2, R2>& matrix) noexcept;

		// Constructs a matrix with the given scalars.
		// NOTE: Elements are in order (col0, row0), (col1, row0), (col2, row0), ... , (col0, row1), ...
		template <Scalar... S2s>
		requires(sizeof...(S2s) == C * R && _IMPL all_true_v<_STD is_convertible_v<S2s, S>...>)
		constexpr Matrix(S2s&&... scalars) noexcept;

		using Tensor<S, 2, C, R>::operator=;
	public:
		using Tensor<S, 2, C, R>::operator+; // also handles unary plus.
		using Tensor<S, 2, C, R>::operator+=;
		using Tensor<S, 2, C, R>::operator-; // also handles unary plus.
		using Tensor<S, 2, C, R>::operator-=;
		using Tensor<S, 2, C, R>::operator*;
		using Tensor<S, 2, C, R>::operator*=;
		using Tensor<S, 2, C, R>::operator/;
		using Tensor<S, 2, C, R>::operator/=;
		using Tensor<S, 2, C, R>::at;
		using Tensor<S, 2, C, R>::operator[];
		using Tensor<S, 2, C, R>::operator==;
		using Tensor<S, 2, C, R>::operator!=;
	public:
		template <Scalar S2, Dimension C2>
		constexpr Matrix<_IMPL CT<S, S2>, C2, R> operator*(const Matrix<S2, C2, C>& matrix) const noexcept;
		template <Scalar S2>
		constexpr Vector<_IMPL CT<S, S2>, R> operator*(const Vector<S2, C>& vector) const noexcept;
	private:
		template <Dimension C2, Dimension R2, Scalar S2, Scalar... S2s>
		constexpr void Fill(S2 scalar, S2s&&... scalars);
	private:
		using Tensor<S, 2, C, R>::m_Scalars;
	};

	// External Operators

	//template <Scalar S1, Scalar S2, Dimension C, Dimension R>
	//constexpr Vector<_IMPL CT<S1, S2>, C> operator*(const Vector<S1, R>& vector, const Matrix<S2, C, R>& matrix) noexcept;

	// Static Methods

	template <Scalar S, Dimension C, Dimension R, Dimension C2, Dimension R2>
	requires(C < C2 && R < R2 && C2 > 1 && R2 > 1)
	constexpr Matrix<S, C2 - 1, R2 - 1> Submatrix(const Matrix<S, C2, R2>& matrix) noexcept;

	template <Scalar S, Dimension C2, Dimension R2>
	requires(/*C < C2 && R < R2 && */C2 > 1 && R2 > 1)
	constexpr Matrix<S, C2 - 1, R2 - 1> Submatrix(Dimension C, Dimension R, const Matrix<S, C2, R2>& matrix) noexcept;

	template <Scalar S, Dimension CR>
	constexpr S Trace(const Matrix<S, CR, CR>& matrix) noexcept;

	template <Scalar S, Dimension CR>
	constexpr S Determinant(const Matrix<S, CR, CR>& matrix) noexcept;

	template <Scalar S, Dimension C, Dimension R>
	constexpr Matrix<S, R, C> Transpose(const Matrix<S, C, R>& matrix) noexcept;

	template <Scalar S, Dimension CR>
	constexpr Matrix<S, CR, CR> Minors(const Matrix<S, CR, CR>& matrix) noexcept;

	template <Scalar S, Dimension CR>
	constexpr Matrix<S, CR, CR> Cofactors(const Matrix<S, CR, CR>& matrix) noexcept;

	template <Scalar S, Dimension CR>
	constexpr Matrix<S, CR, CR> Adjugate(const Matrix<S, CR, CR>& matrix) noexcept;

	template <Scalar S, Dimension CR>
	constexpr Matrix<S, CR, CR> Inverse(const Matrix<S, CR, CR>& matrix) noexcept;

	template <Scalar S1, Scalar S2, Dimension CR>
	constexpr Matrix<_IMPL CT<S1, S2>, CR, CR> Inverse(const Matrix<S1, CR, CR>& matrix, S2 determinant) noexcept;

	template <Scalar S1, Scalar S2, Dimension CR>
	constexpr Matrix<_IMPL CT<S1, S2>, CR, CR> Translate(const Matrix<S1, CR, CR>& matrix, const Vector<S2, CR - 1>& translation);

	template <Dimension A1, Dimension A2, Scalar S1, Scalar S2, Dimension CR>
	requires(A1 < CR - 1 && A2 < CR - 1 && A1 != A2)
	constexpr Matrix<_IMPL CT<S1, S2>, CR, CR> Rotate(const Matrix<S1, CR, CR>& matrix, S2 radians);

	template <Scalar S1, Scalar S2, Dimension CR>
	//requires(A1 < CR - 1 && A2 < CR - 1 && A1 != A2)
	constexpr Matrix<_IMPL CT<S1, S2>, CR, CR> Rotate(Dimension A1, Dimension A2, const Matrix<S1, CR, CR>& matrix, S2 radians);

	template <Scalar S1, Scalar S2, Dimension CR>
	constexpr Matrix<_IMPL CT<S1, S2>, CR, CR> Scale(const Matrix<S1, CR, CR>& matrix, const Vector<S2, CR - 1>& scale);

	template <Scalar S1, Scalar S2, Dimension CR>
	constexpr Matrix<_IMPL CT<S1, S2>, CR, CR> Scale(const Matrix<S1, CR, CR>& matrix, S2 scale);

	// Aliases

	template <Dimension C, Dimension R = C>
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

namespace nd
{
	template <Scalar S, Dimension C, Dimension R>
	constexpr Matrix<S, C, R>::Matrix() noexcept
		: Matrix<S, C, R>{1}
	{

	}

	template <Scalar S, Dimension C, Dimension R>
	template <Scalar S2>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Matrix<S, C, R>::Matrix(S2 scalar) noexcept
	{
		if (scalar != S{})
			for (Dimension cr{}; cr < _GCEM min(C, R); ++cr)
				m_Scalars[cr][cr] = static_cast<S>(scalar);
	}

	template <Scalar S, Dimension C, Dimension R>
	template <Scalar S2, Dimension C2, Dimension R2>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Matrix<S, C, R>::Matrix(const Matrix<S2, C2, R2>& matrix) noexcept
		: Tensor<S, 2, C, R>{matrix}
	{
		for (Dimension cr{_GCEM min(C2, R2)}; cr < _GCEM min(C, R); ++cr)
			m_Scalars[cr][cr] = S{1};
	}

	template <Scalar S, Dimension C, Dimension R>
	template <Scalar... S2s>
	requires(sizeof...(S2s) == C * R && _IMPL all_true_v<_STD is_convertible_v<S2s, S>...>)
	constexpr Matrix<S, C, R>::Matrix(S2s&&... scalars) noexcept
	{
		Fill<C, R>(_STD forward<S2s>(scalars)...);
	}

	template <Scalar S, Dimension C, Dimension R>
	template <Scalar S2, Dimension C2>
	constexpr Matrix<_IMPL CT<S, S2>, C2, R> Matrix<S, C, R>::operator*(const Matrix<S2, C2, C>& matrix) const noexcept
	{
		Matrix<_IMPL CT<S, S2>, C2, R> result{_IMPL CT<S, S2>{}};
		for (Dimension i{}; i < C2; ++i)
			for (Dimension j{}; j < C; ++j)
				result[i] += Vector<_IMPL CT<S, S2>, R>{m_Scalars[j]} * _IMPL CT<S, S2>{matrix[i][j]};
		return result;
	}

	template <Scalar S, Dimension C, Dimension R>
	template <Scalar S2>
	constexpr Vector<_IMPL CT<S, S2>, R> Matrix<S, C, R>::operator*(const Vector<S2, C>& vector) const noexcept
	{
		Vector<_IMPL CT<S, S2>, R> result;
		for (Dimension c{}; c < C; ++c)
			result += _IMPL CT<S, S2>{m_Scalars[c]} * _IMPL CT<S, S2>{vector[c]};
		return result;
	}

	template <Scalar S, Dimension C, Dimension R>
	template <Dimension C2, Dimension R2, Scalar S2, Scalar... S2s>
	constexpr void Matrix<S, C, R>::Fill(S2 scalar, S2s&&... scalars)
	{
		m_Scalars[C - C2][R - R2] = static_cast<S>(scalar);
		if constexpr (C2 > 1)
			Fill<C2 - 1, R2>(_STD forward<S2s>(scalars)...);
		else if constexpr (R2 > 1)
			Fill<C, R2 - 1>(_STD forward<S2s>(scalars)...);
	}

	//template <Scalar S1, Scalar S2, Dimension C, Dimension R>
	//constexpr Vector<_IMPL CT<S1, S2>, C> operator*(const Vector<S1, R>& vector, const Matrix<S2, C, R>& matrix) noexcept
	//{
	//	Vector<_IMPL CT<S1, S2>, C> result;
	//	Matrix<_IMPL CT<S1, S2>, R, C> m{Transpose(Matrix<_IMPL CT<S1, S2>, C, R>{matrix})};
	//	for (Dimension r{}; r < R; ++r)
	//		result += Vector<_IMPL CT<S1, S2>, R>{m[r]} * Vector<_IMPL CT<S1, S2>, R>{vector[r]};
	//	return result;
	//}
	
	template <Scalar S, Dimension C, Dimension R, Dimension C2, Dimension R2>
	requires(C < C2 && R < R2 && C2 > 1 && R2 > 1)
	constexpr Matrix<S, C2 - 1, R2 - 1> Submatrix(const Matrix<S, C2, R2>& matrix) noexcept
	{
		Matrix<S, C2 - 1, R2 - 1> result{};

		for (Dimension c{}; c < C; ++c)
			for (Dimension r{}; r < R; ++r)
				result[c][r] = matrix[c][r];

		for (Dimension c{C + 1}; c < C2; ++c)
			for (Dimension r{}; r < R; ++r)
				result[c - 1][r] = matrix[c][r];

		for (Dimension c{}; c < C; ++c)
			for (Dimension r{R + 1}; r < R2; ++r)
				result[c][r - 1] = matrix[c][r];

		for (Dimension c{C + 1}; c < C2; ++c)
			for (Dimension r{R + 1}; r < R2; ++r)
				result[c - 1][r - 1] = matrix[c][r];

		return result;
	}
	
	template <Scalar S, Dimension C2, Dimension R2>
	requires(/*C < C2 && R < R2 && */C2 > 1 && R2 > 1)
	constexpr Matrix<S, C2 - 1, R2 - 1> Submatrix(Dimension C, Dimension R, const Matrix<S, C2, R2>& matrix) noexcept
	{
		__assume(C < C2 && R < R2);

		Matrix<S, C2 - 1, R2 - 1> result{};

		for (Dimension c{}; c < C; ++c)
			for (Dimension r{}; r < R; ++r)
				result[c][r] = matrix[c][r];

		for (Dimension c{C + 1}; c < C2; ++c)
			for (Dimension r{}; r < R; ++r)
				result[c - 1][r] = matrix[c][r];

		for (Dimension c{}; c < C; ++c)
			for (Dimension r{R + 1}; r < R2; ++r)
				result[c][r - 1] = matrix[c][r];

		for (Dimension c{C + 1}; c < C2; ++c)
			for (Dimension r{R + 1}; r < R2; ++r)
				result[c - 1][r - 1] = matrix[c][r];

		return result;
	}

	template <Scalar S, Dimension CR>
	constexpr S Trace(const Matrix<S, CR, CR>& matrix) noexcept
	{
		S result{};
		for (Dimension cr{}; cr < CR; ++cr)
			result += matrix[cr][cr];
		return result;
	}

	template <Scalar S, Dimension CR>
	constexpr S Determinant(const Matrix<S, CR, CR>& matrix) noexcept
	{
		S result{};

		for (Dimension c{}; c < CR; ++c)
		{
			S subDeterminant{matrix[c][0] * Determinant(Submatrix(c, 0, matrix))};
			if (c % 2 == 0)
				result += subDeterminant;
			else
				result -= subDeterminant;
		};

		return result;
	}

	template <Scalar S>
	constexpr S Determinant(const Matrix<S, 1, 1>& matrix) noexcept
	{
		return matrix[0][0];
	}
	
	template <Scalar S, Dimension C, Dimension R>
	constexpr Matrix<S, R, C> Transpose(const Matrix<S, C, R>& matrix) noexcept
	{
		Matrix<S, R, C> result{};
		for (Dimension c{}; c < C; ++c)
			for (Dimension r{}; r < R; ++r)
				result[r][c] = matrix[c][r];
		return result;
	}

	template <Scalar S, Dimension CR>
	constexpr Matrix<S, CR, CR> Minors(const Matrix<S, CR, CR>& matrix) noexcept
	{
		Matrix<S, CR, CR> result{};
		for (Dimension c{}; c < CR; ++c)
			for (Dimension r{}; r < CR; ++r)
				result[c][r] = Determinant(Submatrix(c, r, matrix));
		return result;
	}

	template <Scalar S, Dimension CR>
	constexpr Matrix<S, CR, CR> Cofactors(const Matrix<S, CR, CR>& matrix) noexcept
	{
		Matrix<S, CR, CR> result{Minors(matrix)};
		for (Dimension c{}; c < CR; ++c)
			for (Dimension r{}; r < CR; ++r)
				if ((c + r) % 2 != 0)
					result[c][r] = -result[c][r];
		return result;
	}

	template <Scalar S, Dimension CR>
	constexpr Matrix<S, CR, CR> Adjugate(const Matrix<S, CR, CR>& matrix) noexcept
	{
		return Transpose(Cofactors(matrix));
	}

	template <Scalar S, Dimension CR>
	constexpr Matrix<S, CR, CR> Inverse(const Matrix<S, CR, CR>& matrix) noexcept
	{
		return Adjugate(matrix) / Determinant(matrix);
	}
	
	template <Scalar S, Dimension CR, Scalar S2>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Matrix<_IMPL CT<S, S2>, CR, CR> Inverse(const Matrix<S, CR, CR>& matrix, S2 determinant) noexcept
	{
		return Adjugate(Matrix<_IMPL CT<S, S2>, CR, CR>{matrix}) / _IMPL CT<S, S2>{determinant};
	}
	
	template <Scalar S, Dimension CR, Scalar S2>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Matrix<_IMPL CT<S, S2>, CR, CR> Translate(const Matrix<S, CR, CR>& matrix, const Vector<S2, CR - 1>& translation)
	{
		Matrix<_IMPL CT<S, S2>, CR, CR> t{};
		for (Dimension r{}; r < CR - 1; ++r)
			t[CR - 1][r] += _IMPL CT<S, S2>{translation[r]};
		return matrix * t;
	}
	
	template <Dimension A1, Dimension A2, Dimension CR, Scalar S, Scalar S2>
	requires(A1 < CR - 1 && A2 < CR - 1 && A1 != A2 && _STD is_convertible_v<S2, S>)
	constexpr Matrix<_IMPL CT<S, S2>, CR, CR> Rotate(const Matrix<S, CR, CR>& matrix, S2 radians)
	{
		Matrix<_IMPL CT<S, S2>, CR, CR> rotation{};

		auto sin{CT<S, S2>(_GCEM sin(CT<S, S2>(radians)))};
		auto cos{CT<S, S2>(_GCEM cos(CT<S, S2>(radians)))};

		rotation[A1][A1] = cos; rotation[A2][A1] = -sin;
		rotation[A1][A2] = sin; rotation[A2][A2] = cos;

		return Matrix<_IMPL CT<S, S2>, CR, CR>{matrix} * rotation;
	}
	
	template <Scalar S, Dimension CR, Scalar S2>
	requires(/*A1 < CR - 1 && A2 < CR - 1 && A1 != A2 && */_STD is_convertible_v<S2, S>)
	constexpr Matrix<_IMPL CT<S, S2>, CR, CR> Rotate(Dimension A1, Dimension A2, const Matrix<S, CR, CR>& matrix, S2 radians)
	{
		__assume(A1 < CR - 1 && A2 < CR - 1 && A1 != A2);

		Matrix<_IMPL CT<S, S2>, CR, CR> rotation{};

		auto sin{CT<S, S2>(_GCEM sin(CT<S, S2>(radians)))};
		auto cos{CT<S, S2>(_GCEM cos(CT<S, S2>(radians)))};

		rotation[A1][A1] = cos; rotation[A2][A1] = -sin;
		rotation[A1][A2] = sin; rotation[A2][A2] = cos;

		return Matrix<_IMPL CT<S, S2>, CR, CR>{matrix} * rotation;
	}
	
	template <Scalar S, Dimension CR, Scalar S2>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Matrix<_IMPL CT<S, S2>, CR, CR> Scale(const Matrix<S, CR, CR>& matrix, const Vector<S2, CR - 1>& scale)
	{
		Matrix<_IMPL CT<S, S2>, CR, CR> s{};
		for (Dimension cr{}; cr < CR - 1; ++cr)
			s[cr][cr] = CT<S, S2>(scale[cr]);
		return Matrix<_IMPL CT<S, S2>, CR, CR>{matrix} * s;
	}
	
	template <Scalar S, Dimension CR, Scalar S2>
	requires(_STD is_convertible_v<S2, S>)
	constexpr Matrix<_IMPL CT<S, S2>, CR, CR> Scale(const Matrix<S, CR, CR>& matrix, S2 scale)
	{
		Matrix<_IMPL CT<S, S2>, CR, CR> s{_IMPL CT<S, S2>{scale}};
		s[CR - 1][CR - 1] = _IMPL CT<S, S2>{1};
		return Matrix<_IMPL CT<S, S2>, CR, CR>{matrix} * s;
	}
}
