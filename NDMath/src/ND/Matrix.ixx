module;

#include <gcem.hpp>

export module nd.matrix;

import nd.dimension;
import nd.scalar;
import nd.vector;
import nd.impl;
import <iomanip>;
import <iostream>;

#define _ND ::nd::
#define _IMPL ::nd::impl::

namespace nd
{
	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	struct Matrix;
}

template<_ND Scalar... Scalars>
using CT = _STD common_type_t<Scalars...>;

// This is so small, I'm fine with duplicating it.
template<_ND Dimension N, _ND Scalar... Scalars>
using CVT = _ND Vector<N, CT<Scalars...>>;

template<_ND Dimension C, _ND Dimension R, _ND Scalar... Scalars>
using CMT = _ND Matrix<C, R, CT<Scalars...>>;

export namespace nd
{
	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	struct Matrix
	{
	public:
		constexpr Matrix() noexcept;

		template<Dimension C2, Dimension R2, Scalar S2> requires(C2 <= C && R2 <= R && _STD is_convertible_v<S2, S>)
		constexpr Matrix(const Matrix<C2, R2, S2>& matrix) noexcept;

		// Collides with copy constructor if requires is used instead of _STD enable_if_t.
		template<Scalar S2, _STD enable_if_t<_STD is_convertible_v<S2, S>, int> = 0>
		constexpr Matrix(const S2& scalar) noexcept;

		// NOTE: Elements are in order (col0, row0), (col1, row0), (col2, row0), ... , (col0, row1), ...
		// Collides with copy constructor if requires is used instead of _STD enable_if_t.
		template<Scalar... S2s, _STD enable_if_t<sizeof...(S2s) == R * C && _IMPL all_true_v<_STD is_convertible_v<S2s, S>...>, int> = 0>
		constexpr Matrix(S2s&&... scalars) noexcept;

		// Collides with copy operator= if requires is used instead of _STD enable_if_t.
		template<Scalar S2, _STD enable_if_t<_STD is_convertible_v<S2, S>, int> = 0>
		constexpr auto operator=(const S2& scalar) noexcept -> Matrix<C, R, S>&;

		template<Dimension C2, Dimension R2, Scalar S2> requires(C2 <= C && R2 <= R && _STD is_convertible_v<S2, S>)
		constexpr auto operator=(const Matrix<C2, R2, S2>& matrix) noexcept -> Matrix<C, R, S>&;
	public:
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator+(const S2& scalar) const noexcept -> CMT<C, R, S, S2>;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator+=(const S2& scalar) noexcept -> Matrix<C, R, S>&;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator+(const Matrix<C, R, S2>& matrix) const noexcept -> CMT<C, R, S, S2>;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator+=(const Matrix<C, R, S2>& matrix) noexcept -> Matrix<C, R, S>&;
	public:
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator-(const S2& scalar) const noexcept -> CMT<C, R, S, S2>;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator-=(const S2& scalar) noexcept -> Matrix<C, R, S>&;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator-(const Matrix<C, R, S2>& matrix) const noexcept -> CMT<C, R, S, S2>;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator-=(const Matrix<C, R, S2>& matrix) noexcept -> Matrix<C, R, S>&;
	public:
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator*(const S2& scalar) const noexcept -> CMT<C, R, S, S2>;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator*=(const S2& scalar) noexcept -> Matrix<C, R, S>&;

		template<Dimension C2, Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator*(const Matrix<C2, C, S2>& matrix) const noexcept -> CMT<C2, R, S, S2>;

		template<Scalar S2> requires(R == C && _STD is_convertible_v<S2, S>)
		constexpr auto operator*=(const Matrix<C, C, S2>& matrix) noexcept -> Matrix<C, R, S>&;
	public:
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator/(const S2& scalar) const noexcept -> CMT<C, R, S, S2>;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator/=(const S2& scalar) noexcept -> Matrix<C, R, S>&;

		template<Scalar S2> requires(R == C && _STD is_convertible_v<S2, S>)
		constexpr auto operator/(const Matrix<C, C, S2>& matrix) const noexcept -> CMT<C, R, S, S2>;

		template<Scalar S2> requires(R == C && _STD is_convertible_v<S2, S>)
		constexpr auto operator/=(const Matrix<C, C, S2>& matrix) noexcept -> Matrix<C, R, S>&;
	public:
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator*(const Vector<C, S2>& vector) const noexcept -> CVT<R, S, S2>;
	public:
		constexpr auto operator+() const noexcept -> Matrix<C, R, S>;
		constexpr auto operator-() const noexcept -> Matrix<C, R, S>;
	public:
		constexpr auto operator[](Dimension col) noexcept -> Vector<R, S>&;
		constexpr auto operator[](Dimension col) const noexcept -> const Vector<R, S>&;

		template<Dimension C2> requires(C2 < C)
		constexpr auto Col() const noexcept -> Vector<C, S>;

		template<Dimension R2> requires(R2 < R)
		constexpr auto Row() const noexcept -> Vector<C, S>;

		constexpr auto Col(Dimension C2) const noexcept -> Vector<C, S>;
		constexpr auto Row(Dimension R2) const noexcept -> Vector<C, S>;
	private:
		template<Dimension C2, Dimension R2, Scalar S2, Scalar... Scalars>
		constexpr auto Fill(const S2& scalar, Scalars&&... scalars) -> void;
	private:
		Vector<R, S> m_Scalars[C]{Vector<R, S>{}};
	};

	// External Operators

	template<Dimension C, Dimension R, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator+(const S2& scalar, const Matrix<C, R, S>& matrix) noexcept -> CMT<C, R, S, S2>;

	template<Dimension C, Dimension R, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator-(const S2& scalar, const Matrix<C, R, S>& matrix) noexcept -> CMT<C, R, S, S2>;

	template<Dimension C, Dimension R, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator*(const S2& scalar, const Matrix<C, R, S>& matrix) noexcept -> CMT<C, R, S, S2>;

	template<Dimension CR, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator/(const S2& scalar, const Matrix<CR, CR, S>& matrix) noexcept -> CMT<CR, CR, S, S2>;

	template<Dimension C, Dimension R, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator*(const Vector<R, S2>& vector, const Matrix<C, R, S>& matrix) noexcept -> CVT<C, S, S2>;

	// Static Methods

	template<Dimension C, Dimension R, Scalar S, Dimension C2, Dimension R2> requires(C < C2 && R < R2 && C2 > 1 && R2 > 1)
	constexpr auto Submatrix(const Matrix<C2, R2, S>& matrix) noexcept -> Matrix<C2 - 1, R2 - 1, S>;

	template<Dimension C2, Dimension R2, Scalar S> requires(/*C < C2 && R < R2 && */C2 > 1 && R2 > 1)
	constexpr auto Submatrix(Dimension C, Dimension R, const Matrix<C2, R2, S>& matrix) noexcept -> Matrix<C2 - 1, R2 - 1, S>;

	template<Dimension CR, Scalar S>
	constexpr auto Trace(const Matrix<CR, CR, S>& matrix) noexcept -> S;

	template<Dimension CR, Scalar S>
	constexpr auto Determinant(const Matrix<CR, CR, S>& matrix) noexcept -> S;

	template<Dimension C, Dimension R, Scalar S>
	constexpr auto Transpose(const Matrix<C, R, S>& matrix) noexcept -> Matrix<R, C, S>;

	template<Dimension CR, Scalar S>
	constexpr auto Minors(const Matrix<CR, CR, S>& matrix) noexcept -> Matrix<CR, CR, S>;

	template<Dimension CR, Scalar S>
	constexpr auto Cofactors(const Matrix<CR, CR, S>& matrix) noexcept -> Matrix<CR, CR, S>;

	template<Dimension CR, Scalar S>
	constexpr auto Adjugate(const Matrix<CR, CR, S>& matrix) noexcept -> Matrix<CR, CR, S>;

	template<Dimension CR, Scalar S>
	constexpr auto Inverse(const Matrix<CR, CR, S>& matrix) noexcept -> Matrix<CR, CR, S>;

	template<Dimension CR, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Inverse(const Matrix<CR, CR, S>& matrix, const S2& determinant) noexcept -> CMT<CR, CR, S, S2>;

	template<Dimension CR, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Translate(const Matrix<CR, CR, S>& matrix, const Vector<CR - 1, S2>& translation) -> CMT<CR, CR, S, S2>;

	template<Dimension A1, Dimension A2, Dimension CR, Scalar S, Scalar S2>
		requires(A1 < CR - 1 && A2 < CR - 1 && A1 != A2 && _STD is_convertible_v<S2, S>)
	constexpr auto Rotate(const Matrix<CR, CR, S>& matrix, const S2& radians) -> CMT<CR, CR, S, S2>;

	template<Dimension CR, Scalar S, Scalar S2>
		requires(/*A1 < CR - 1 && A2 < CR - 1 && A1 != A2 && */_STD is_convertible_v<S2, S>)
	constexpr auto Rotate(Dimension A1, Dimension A2, const Matrix<CR, CR, S>& matrix, const S2& radians) -> CMT<CR, CR, S, S2>;

	template<Dimension CR, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Scale(const Matrix<CR, CR, S>& matrix, const Vector<CR - 1, S2>& scale) -> CMT<CR, CR, S, S2>;

	template<Dimension CR, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Scale(const Matrix<CR, CR, S>& matrix, const S2& scale) -> CMT<CR, CR, S, S2>;

	template<Dimension C, Dimension R, Scalar S>
	auto operator<<(_STD ostream& ostream, const Matrix<C, R, S>& matrix) -> _STD ostream&;

	// Aliases

	template<Dimension C, Dimension R = C>
	using MatrixCRf = Matrix<C, R, float>;

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
	template<Dimension C, Dimension R, Scalar S>
	auto PrintRow(_STD ostream& ostream, const Matrix<C, R, S>& matrix, long long width, Dimension r) noexcept -> void
	{
		ostream << '[' << _STD setw(width) << +matrix[0][r];
		for (Dimension c{1}; c < C; c++)
			ostream << _STD setw(width + 1) << +matrix[c][r];
		ostream << ']';
	};
}

export namespace nd
{
	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	constexpr Matrix<C, R, S>::Matrix() noexcept
		: Matrix<C, R, S>{1}
	{

	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Dimension C2, Dimension R2, Scalar S2> requires(C2 <= C && R2 <= R && _STD is_convertible_v<S2, S>)
	constexpr Matrix<C, R, S>::Matrix(const Matrix<C2, R2, S2>& matrix) noexcept
		: Matrix<C, R, S>{1}
	{
		for (Dimension c{}; c < C2; c++)
			for (Dimension r{}; r < R2; r++)
				m_Scalars[c][r] = S{matrix[c][r]};
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2, _STD enable_if_t<_STD is_convertible_v<S2, S>, int>>
	constexpr Matrix<C, R, S>::Matrix(const S2& scalar) noexcept
	{
		if (scalar != S{})
			for (Dimension i{}; i < gcem::min(C, R); i++)
				m_Scalars[i][i] = S{scalar};
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar... S2s, _STD enable_if_t<sizeof...(S2s) == R * C && _IMPL all_true_v<_STD is_convertible_v<S2s, S>...>, int>>
	constexpr Matrix<C, R, S>::Matrix(S2s&&... scalars) noexcept
	{
		Fill<C, R>(_STD forward<S2s>(scalars)...);
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2, _STD enable_if_t<_STD is_convertible_v<S2, S>, int>>
	constexpr auto Matrix<C, R, S>::operator=(const S2& scalar) noexcept -> Matrix<C, R, S>&
	{
		return *this = Matrix<C, R, S>(scalar);
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Dimension C2, Dimension R2, Scalar S2> requires(C2 <= C && R2 <= R && _STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator=(const Matrix<C2, R2, S2>& matrix) noexcept -> Matrix<C, R, S>&
	{
		if (this != (void*)&matrix)
		{
			for (Dimension c{}; c < C2; c++)
				for (Dimension r{}; r < R2; r++)
					m_Scalars[c][r] = S{matrix[c][r]};

			for (Dimension c{C2}; c < C; c++)
				for (Dimension r{}; r < R2 + c - C2; r++)
					m_Scalars[c][r] = S{};

			for (Dimension r{R2}; r < R; r++)
				for (Dimension c{}; c < C2 + r - R2; c++)
					m_Scalars[c][r] = S{};

			for (Dimension cr{gcem::min(C2, R2)}; cr < gcem::min(C, R); cr++)
				m_Scalars[cr][cr] = S{1};
		}
		return *this;
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator+(const S2& scalar) const noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>{*this} += CT<S, S2>{scalar};
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator+=(const S2& scalar) noexcept -> Matrix<C, R, S>&
	{
		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				m_Scalars[c][r] += S{scalar};
		return *this;
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator+(const Matrix<C, R, S2>& matrix) const noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>{*this} += CMT<C, R, S, S2>{matrix};
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator+=(const Matrix<C, R, S2>& matrix) noexcept -> Matrix<C, R, S>&
	{
		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				m_Scalars[c][r] += S{matrix[c][r]};
		return *this;
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator-(const S2& scalar) const noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>{*this} -= CT<S, S2>{scalar};
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator-=(const S2& scalar) noexcept -> Matrix<C, R, S>&
	{
		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				m_Scalars[c][r] -= S{scalar};
		return *this;
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator-(const Matrix<C, R, S2>& matrix) const noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>{*this} -= CMT<C, R, S, S2>{matrix};
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator-=(const Matrix<C, R, S2>& matrix) noexcept -> Matrix<C, R, S>&
	{
		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				m_Scalars[c][r] -= S{matrix[c][r]};
		return *this;
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator*(const S2& scalar) const noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>{*this} *= CT<S, S2>{scalar};
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator*=(const S2& scalar) noexcept -> Matrix<C, R, S>&
	{
		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				m_Scalars[c][r] *= S{scalar};
		return *this;
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Dimension C2, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator*(const Matrix<C2, C, S2>& matrix) const noexcept -> CMT<C2, R, S, S2>
	{
		CMT<C2, R, S, S2> result{CT<S, S2>{}};
		for (Dimension i{}; i < C2; i++)
			for (Dimension j{}; j < C; j++)
				result[i] += CVT<R, S, S2>{m_Scalars[j]} * CT<S, S2>{matrix[i][j]};
		return result;
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(R == C && _STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator*=(const Matrix<C, C, S2>& matrix) noexcept -> Matrix<C, R, S>&
	{
		return *this = CMT<C, R, S, S2>{*this} * CT<S, S2>{matrix};
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator/(const S2& scalar) const noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>{*this} /= CT<S, S2>{scalar};
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator/=(const S2& scalar) noexcept -> Matrix<C, R, S>&
	{
		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				m_Scalars[c][r] /= S{scalar};
		return *this;
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(R == C && _STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator/(const Matrix<C, C, S2>& matrix) const noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>{*this} * Inverse(CMT<C, R, S, S2>{matrix});
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(R == C && _STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator/=(const Matrix<C, C, S2>& matrix) noexcept -> Matrix<C, R, S>&
	{
		return *this = CMT<C, R, S, S2>{*this} / CMT<C, C, S, S2>{matrix};
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator*(const Vector<C, S2>& vector) const noexcept -> CVT<R, S, S2>
	{
		CVT<R, S, S2> result;
		for (Dimension c{}; c < C; c++)
			result += CT<S, S2>{m_Scalars[c]} * CT<S, S2>{vector[c]};
		return result;
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	constexpr auto Matrix<C, R, S>::operator+() const noexcept -> Matrix<C, R, S>
	{
		return Matrix<C, R, S>{*this};
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	constexpr auto Matrix<C, R, S>::operator-() const noexcept -> Matrix<C, R, S>
	{
		Matrix<C, R, S> result{+*this};
		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				m_Scalars[c][r] = -m_Scalars[c][r];
		return result;
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	constexpr auto Matrix<C, R, S>::operator[](Dimension col) noexcept -> Vector<R, S>&
	{
		return m_Scalars[col];
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	constexpr auto Matrix<C, R, S>::operator[](Dimension col) const noexcept -> const Vector<R, S>&
	{
		return m_Scalars[col];
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Dimension C2> requires(C2 < C)
	constexpr auto Matrix<C, R, S>::Col() const noexcept -> Vector<C, S>
	{
		return m_Scalars[C2];
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Dimension R2> requires(R2 < R)
	constexpr auto Matrix<C, R, S>::Row() const noexcept -> Vector<C, S>
	{
		Vector<C, S> row;
		for (Dimension c{}; c < C; c++)
			row[c] = m_Scalars[c][R2];
		return row;
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	constexpr auto Matrix<C, R, S>::Col(Dimension C2) const noexcept -> Vector<C, S>
	{
		__assume(C2 < C);
		return m_Scalars[C2];
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	constexpr auto Matrix<C, R, S>::Row(Dimension R2) const noexcept -> Vector<C, S>
	{
		__assume(R2 < R);
		Vector<C, S> row;
		for (Dimension c{}; c < C; c++)
			row[c] = m_Scalars[c][R2];
		return row;
	}

	template<Dimension C, Dimension R, Scalar S> requires(C > 0 && R > 0)
	template<Dimension C2, Dimension R2, Scalar S2, Scalar... Scalars>
	constexpr auto Matrix<C, R, S>::Fill(const S2& scalar, Scalars&&... scalars) -> void
	{
		m_Scalars[C - C2][R - R2] = static_cast<S>(scalar);
		if constexpr (C2 > 1)
			Fill<C2 - 1, R2>(_STD forward<Scalars>(scalars)...);
		else if constexpr (R2 > 1)
			Fill<C, R2 - 1>(_STD forward<Scalars>(scalars)...);
	}
	
	template<Dimension C, Dimension R, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator+(const S2& scalar, const Matrix<C, R, S>& matrix) noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>{matrix} += CT<S, S2>{scalar};
	}
	
	template<Dimension C, Dimension R, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator-(const S2& scalar, const Matrix<C, R, S>& matrix) noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>{-matrix} += CT<S, S2>{scalar};
	}
	
	template<Dimension C, Dimension R, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator*(const S2& scalar, const Matrix<C, R, S>& matrix) noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>{matrix} *= CT<S, S2>{scalar};
	}
	
	template<Dimension CR, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator/(const S2& scalar, const Matrix<CR, CR, S>& matrix) noexcept -> CMT<CR, CR, S, S2>
	{
		return Inverse(CMT<CR, CR, S, S2>{matrix}) *= CT<S, S2>{scalar};
	}
	
	template<Dimension C, Dimension R, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator*(const Vector<R, S2>& vector, const Matrix<C, R, S>& matrix) noexcept -> CVT<C, S, S2>
	{
		CVT<C, S, S2> result;

		CMT<R, C, S, S2> m{Transpose(CMT<C, R, S, S2>{matrix})};
		for (Dimension r{}; r < R; r++)
			result += CVT<R, S, S2>{m[r]} * CVT<R, S, S2>{vector[r]};

		return result;
	}
	
	template<Dimension C, Dimension R, Scalar S, Dimension C2, Dimension R2> requires(C < C2 && R < R2 && C2 > 1 && R2 > 1)
	constexpr auto Submatrix(const Matrix<C2, R2, S>& matrix) noexcept -> Matrix<C2 - 1, R2 - 1, S>
	{
		Matrix<C2 - 1, R2 - 1, S> result{S{}};

		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				result[c][r] = matrix[c][r];

		for (Dimension c{C + 1}; c < C2; c++)
			for (Dimension r{}; r < R; r++)
				result[c - 1][r] = matrix[c][r];

		for (Dimension c{}; c < C; c++)
			for (Dimension r{R + 1}; r < R2; r++)
				result[c][r - 1] = matrix[c][r];

		for (Dimension c{C + 1}; c < C2; c++)
			for (Dimension r{R + 1}; r < R2; r++)
				result[c - 1][r - 1] = matrix[c][r];

		return result;
	}
	
	template<Dimension C2, Dimension R2, Scalar S>
		requires(/*C < C2 && R < R2 && */C2 > 1 && R2 > 1)
	constexpr auto Submatrix(Dimension C, Dimension R, const Matrix<C2, R2, S>& matrix) noexcept -> Matrix<C2 - 1, R2 - 1, S>
	{
		__assume(C < C2 && R < R2);

		Matrix<C2 - 1, R2 - 1, S> result{S{}};

		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				result[c][r] = matrix[c][r];

		for (Dimension c{C + 1}; c < C2; c++)
			for (Dimension r{}; r < R; r++)
				result[c - 1][r] = matrix[c][r];

		for (Dimension c{}; c < C; c++)
			for (Dimension r{R + 1}; r < R2; r++)
				result[c][r - 1] = matrix[c][r];

		for (Dimension c{C + 1}; c < C2; c++)
			for (Dimension r{R + 1}; r < R2; r++)
				result[c - 1][r - 1] = matrix[c][r];

		return result;
	}

	template<Dimension CR, Scalar S>
	constexpr auto Trace(const Matrix<CR, CR, S>& matrix) noexcept -> S
	{
		S result{S{}};
		for (Dimension cr{}; cr < CR; cr++)
			result += matrix[cr][cr];
		return result;
	}

	template<Dimension CR, Scalar S>
	constexpr auto Determinant(const Matrix<CR, CR, S>& matrix) noexcept -> S
	{
		S result{S{}};

		for (Dimension c{}; c < CR; c++)
		{
			S subDeterminant{matrix[c][0] * Determinant(Submatrix(c, 0, matrix))};
			if (c % 2 == 0)
				result += subDeterminant;
			else
				result -= subDeterminant;
		};

		return result;
	}

	template<Scalar S>
	constexpr auto Determinant(const Matrix<1, 1, S>& matrix) noexcept -> S
	{
		return matrix[0][0];
	}
	
	template<Dimension C, Dimension R, Scalar S>
	constexpr auto Transpose(const Matrix<C, R, S>& matrix) noexcept -> Matrix<R, C, S>
	{
		Matrix<R, C, S> result = S();
		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				result[r][c] = matrix[c][r];
		return result;
	}

	template<Dimension CR, Scalar S>
	constexpr auto Minors(const Matrix<CR, CR, S>& matrix) noexcept -> Matrix<CR, CR, S>
	{
		Matrix<CR, CR, S> result = S();
		for (Dimension c{}; c < CR; c++)
			for (Dimension r{}; r < CR; r++)
				result[c][r] = Determinant(Submatrix(c, r, matrix));
		return result;
	}

	template<Dimension CR, Scalar S>
	constexpr auto Cofactors(const Matrix<CR, CR, S>& matrix) noexcept -> Matrix<CR, CR, S>
	{
		Matrix<CR, CR, S> result = Minors(matrix);
		for (Dimension c{}; c < CR; c++)
			for (Dimension r{}; r < CR; r++)
				if ((c + r) % 2 != 0)
					result[c][r] = -result[c][r];
		return result;
	}

	template<Dimension CR, Scalar S>
	constexpr auto Adjugate(const Matrix<CR, CR, S>& matrix) noexcept -> Matrix<CR, CR, S>
	{
		return Transpose(Cofactors(matrix));
	}

	template<Dimension CR, Scalar S>
	constexpr auto Inverse(const Matrix<CR, CR, S>& matrix) noexcept -> Matrix<CR, CR, S>
	{
		return Adjugate(matrix) / Determinant(matrix);
	}
	
	template<Dimension CR, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Inverse(const Matrix<CR, CR, S>& matrix, const S2& determinant) noexcept -> CMT<CR, CR, S, S2>
	{
		return Adjugate(CMT<CR, CR, S, S2>{matrix}) / CT<S, S2>{determinant};
	}
	
	template<Dimension CR, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Translate(const Matrix<CR, CR, S>& matrix, const Vector<CR - 1, S2>& translation) -> CMT<CR, CR, S, S2>
	{
		CMT<CR, CR, S, S2> t;
		for (Dimension r{}; r < CR - 1; r++)
			t[CR - 1][r] += CT<S, S2>{translation[r]};
		return matrix * t;
	}
	
	template<Dimension A1, Dimension A2, Dimension CR, Scalar S, Scalar S2>
		requires(A1 < CR - 1 && A2 < CR - 1 && A1 != A2 && _STD is_convertible_v<S2, S>)
	constexpr auto Rotate(const Matrix<CR, CR, S>& matrix, const S2& radians) -> CMT<CR, CR, S, S2>
	{
		CMT<CR, CR, S, S2> rotation;

		auto sin{CT<S, S2>(gcem::sin(CT<S, S2>(radians)))};
		auto cos{CT<S, S2>(gcem::cos(CT<S, S2>(radians)))};

		rotation[A1][A1] = cos; rotation[A2][A1] = -sin;
		rotation[A1][A2] = sin; rotation[A2][A2] = cos;

		return CMT<CR, CR, S, S2>{matrix} * rotation;
	}
	
	template<Dimension CR, Scalar S, Scalar S2>
		requires(/*A1 < CR - 1 && A2 < CR - 1 && A1 != A2 && */_STD is_convertible_v<S2, S>)
	constexpr auto Rotate(Dimension A1, Dimension A2, const Matrix<CR, CR, S>& matrix, const S2& radians) -> CMT<CR, CR, S, S2>
	{
		__assume(A1 < CR - 1 && A2 < CR - 1 && A1 != A2);

		CMT<CR, CR, S, S2> rotation;

		auto sin{CT<S, S2>(gcem::sin(CT<S, S2>(radians)))};
		auto cos{CT<S, S2>(gcem::cos(CT<S, S2>(radians)))};

		rotation[A1][A1] = cos; rotation[A2][A1] = -sin;
		rotation[A1][A2] = sin; rotation[A2][A2] = cos;

		return CMT<CR, CR, S, S2>{matrix} * rotation;
	}
	
	template<Dimension CR, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Scale(const Matrix<CR, CR, S>& matrix, const Vector<CR - 1, S2>& scale) -> CMT<CR, CR, S, S2>
	{
		CMT<CR, CR, S, S2> s;
		for (Dimension cr{}; cr < CR - 1; cr++)
			s[cr][cr] = CT<S, S2>(scale[cr]);
		return CMT<CR, CR, S, S2>{matrix} * s;
	}
	
	template<Dimension CR, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Scale(const Matrix<CR, CR, S>& matrix, const S2& scale) -> CMT<CR, CR, S, S2>
	{
		CMT<CR, CR, S, S2> s{CT<S, S2>{scale}};
		s[CR - 1][CR - 1] = CT<S, S2>{1};
		return CMT<CR, CR, S, S2>{matrix} * s;
	}

	template<Dimension C, Dimension R, _STD floating_point S>
	auto operator<<(_STD ostream& ostream, const Matrix<C, R, S>& matrix) -> _STD ostream&
	{
		ostream << _STD showpos << _STD scientific;
		_IMPL PrintRow(ostream, matrix, 13, 0);
		for (Dimension r{1}; r < R; r++)
			_IMPL PrintRow(ostream << '\n', matrix, 13, r);
		return ostream << _STD defaultfloat << _STD noshowpos;
	}

	template<Dimension C, Dimension R, _STD integral S>
	auto operator<<(_STD ostream& ostream, const Matrix<C, R, S>& matrix) -> _STD ostream&
	{
		constexpr long long width{_STD is_signed_v<S> + static_cast<long long>(gcem::max(
			S{gcem::ceil(gcem::log10(_STD make_unsigned_t<S>(_STD numeric_limits<S>::min())))},
			S{gcem::ceil(gcem::log10(_STD make_unsigned_t<S>(_STD numeric_limits<S>::max())))}
		))};

		_IMPL PrintRow(ostream, matrix, width, 0);
		for (Dimension r{1}; r < R; r++)
			_IMPL PrintRow(ostream << '\n', matrix, width, r);
		return ostream;
	}
}
