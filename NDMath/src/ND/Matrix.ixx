module;

#include <gcem.hpp>

export module nd.matrix;

import nd.types;
import nd.impl;
import nd.tensor;
import nd.vector;
import <iomanip>;
import <iostream>;

#define _ND ::nd::
#define _IMPL ::nd::impl::

export namespace nd
{
	template<Scalar S, Dimension C, Dimension R>
	struct Matrix : public Tensor<S, 2, C, R>
	{
	public:
		constexpr Matrix() noexcept;

		template<Scalar S2, Dimension C2, Dimension R2> requires(C2 <= C && R2 <= R && _STD is_convertible_v<S2, S>)
		constexpr Matrix(const Matrix<S2, C2, R2>& matrix) noexcept;

		// Collides with copy constructor if requires is used instead of _STD enable_if_t.
		template<Scalar S2, _STD enable_if_t<_STD is_convertible_v<S2, S>, int> = 0>
		constexpr Matrix(S2 scalar) noexcept;

		// NOTE: Elements are in order (col0, row0), (col1, row0), (col2, row0), ... , (col0, row1), ...
		// Collides with copy constructor if requires is used instead of _STD enable_if_t.
		template<Scalar... S2s, _STD enable_if_t<sizeof...(S2s) == R * C && _IMPL all_true_v<_STD is_convertible_v<S2s, S>...>, int> = 0>
		constexpr Matrix(S2s&&... scalars) noexcept;

		// Collides with copy operator= if requires is used instead of _STD enable_if_t.
		template<Scalar S2, _STD enable_if_t<_STD is_convertible_v<S2, S>, int> = 0>
		constexpr auto operator=(S2 scalar) noexcept -> Matrix<S, C, R>&;

		template<Scalar S2, Dimension C2, Dimension R2> requires(C2 <= C && R2 <= R && _STD is_convertible_v<S2, S>)
		constexpr auto operator=(const Matrix<S2, C2, R2>& matrix) noexcept -> Matrix<S, C, R>&;
	public:
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator+(S2 scalar) const noexcept -> Matrix<_IMPL CT<S, S2>, C, R>;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator+=(S2 scalar) noexcept -> Matrix<S, C, R>&;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator+(const Matrix<S2, C, R>& matrix) const noexcept -> Matrix<_IMPL CT<S, S2>, C, R>;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator+=(const Matrix<S2, C, R>& matrix) noexcept -> Matrix<S, C, R>&;
	public:
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator-(S2 scalar) const noexcept -> Matrix<_IMPL CT<S, S2>, C, R>;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator-=(S2 scalar) noexcept -> Matrix<S, C, R>&;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator-(const Matrix<S2, C, R>& matrix) const noexcept -> Matrix<_IMPL CT<S, S2>, C, R>;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator-=(const Matrix<S2, C, R>& matrix) noexcept -> Matrix<S, C, R>&;
	public:
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator*(S2 scalar) const noexcept -> Matrix<_IMPL CT<S, S2>, C, R>;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator*=(S2 scalar) noexcept -> Matrix<S, C, R>&;

		template<Scalar S2, Dimension C2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator*(const Matrix<S2, C2, C>& matrix) const noexcept -> Matrix<_IMPL CT<S, S2>, C2, R>;

		template<Scalar S2> requires(R == C && _STD is_convertible_v<S2, S>)
		constexpr auto operator*=(const Matrix<S2, C, C>& matrix) noexcept -> Matrix<S, C, R>&;
	public:
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator/(S2 scalar) const noexcept -> Matrix<_IMPL CT<S, S2>, C, R>;

		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator/=(S2 scalar) noexcept -> Matrix<S, C, R>&;

		template<Scalar S2> requires(R == C && _STD is_convertible_v<S2, S>)
		constexpr auto operator/(const Matrix<S2, C, C>& matrix) const noexcept -> Matrix<_IMPL CT<S, S2>, C, R>;

		template<Scalar S2> requires(R == C && _STD is_convertible_v<S2, S>)
		constexpr auto operator/=(const Matrix<S2, C, C>& matrix) noexcept -> Matrix<S, C, R>&;
	public:
		template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
		constexpr auto operator*(const Vector<S2, C>& vector) const noexcept -> Vector<_IMPL CT<S, S2>, R>;
	public:
		constexpr auto operator+() const noexcept -> Matrix<S, C, R>;
		constexpr auto operator-() const noexcept -> Matrix<S, C, R>;
	public:
		//constexpr auto operator[](Dimension col) noexcept -> Vector<S, R>&;
		//constexpr auto operator[](Dimension col) const noexcept -> const Vector<S, R>&;
		using Tensor<S, 2, C, R>::operator[];

		template<Dimension C2> requires(C2 < C)
		constexpr auto Col() const noexcept -> Vector<S, R>;

		template<Dimension R2> requires(R2 < R)
		constexpr auto Row() const noexcept -> Vector<S, C>;

		constexpr auto Col(Dimension C2) const noexcept -> Vector<S, R>;
		constexpr auto Row(Dimension R2) const noexcept -> Vector<S, C>;
	private:
		template<Dimension C2, Dimension R2, Scalar S2, Scalar... Scalars>
		constexpr auto Fill(S2 scalar, Scalars&&... scalars) -> void;
	private:
		using Tensor<S, 2, C, R>::m_Scalars;
	};

	// External Operators

	template<Scalar S, Dimension C, Dimension R, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator+(S2 scalar, const Matrix<S, C, R>& matrix) noexcept -> Matrix<_IMPL CT<S, S2>, C, R>;

	template<Scalar S, Dimension C, Dimension R, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator-(S2 scalar, const Matrix<S, C, R>& matrix) noexcept -> Matrix<_IMPL CT<S, S2>, C, R>;

	template<Scalar S, Dimension C, Dimension R, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator*(S2 scalar, const Matrix<S, C, R>& matrix) noexcept -> Matrix<_IMPL CT<S, S2>, C, R>;

	template<Scalar S, Dimension CR, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator/(S2 scalar, const Matrix<S, CR, CR>& matrix) noexcept -> Matrix<_IMPL CT<S, S2>, CR, CR>;

	template<Scalar S, Dimension C, Dimension R, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator*(const Vector<S2, R>& vector, const Matrix<S, C, R>& matrix) noexcept -> Vector<_IMPL CT<S, S2>, C>;

	// Static Methods

	template<Scalar S, Dimension C, Dimension R, Dimension C2, Dimension R2> requires(C < C2 && R < R2 && C2 > 1 && R2 > 1)
	constexpr auto Submatrix(const Matrix<S, C2, R2>& matrix) noexcept -> Matrix<S, C2 - 1, R2 - 1>;

	template<Scalar S, Dimension C2, Dimension R2> requires(/*C < C2 && R < R2 && */C2 > 1 && R2 > 1)
	constexpr auto Submatrix(Dimension C, Dimension R, const Matrix<S, C2, R2>& matrix) noexcept -> Matrix<S, C2 - 1, R2 - 1>;

	template<Scalar S, Dimension CR>
	constexpr auto Trace(const Matrix<S, CR, CR>& matrix) noexcept -> S;

	template<Scalar S, Dimension CR>
	constexpr auto Determinant(const Matrix<S, CR, CR>& matrix) noexcept -> S;

	template<Scalar S, Dimension C, Dimension R>
	constexpr auto Transpose(const Matrix<S, C, R>& matrix) noexcept -> Matrix<S, R, C>;

	template<Scalar S, Dimension CR>
	constexpr auto Minors(const Matrix<S, CR, CR>& matrix) noexcept -> Matrix<S, CR, CR>;

	template<Scalar S, Dimension CR>
	constexpr auto Cofactors(const Matrix<S, CR, CR>& matrix) noexcept -> Matrix<S, CR, CR>;

	template<Scalar S, Dimension CR>
	constexpr auto Adjugate(const Matrix<S, CR, CR>& matrix) noexcept -> Matrix<S, CR, CR>;

	template<Scalar S, Dimension CR>
	constexpr auto Inverse(const Matrix<S, CR, CR>& matrix) noexcept -> Matrix<S, CR, CR>;

	template<Scalar S, Dimension CR, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Inverse(const Matrix<S, CR, CR>& matrix, S2 determinant) noexcept -> Matrix<_IMPL CT<S, S2>, CR, CR>;

	template<Scalar S, Dimension CR, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Translate(const Matrix<S, CR, CR>& matrix, const Vector<S2, CR - 1>& translation) -> Matrix<_IMPL CT<S, S2>, CR, CR>;

	template<Dimension A1, Dimension A2, Scalar S, Dimension CR, Scalar S2>
	requires(A1 < CR - 1 && A2 < CR - 1 && A1 != A2 && _STD is_convertible_v<S2, S>)
	constexpr auto Rotate(const Matrix<S, CR, CR>& matrix, S2 radians) -> Matrix<_IMPL CT<S, S2>, CR, CR>;

	template<Scalar S, Dimension CR, Scalar S2>
	requires(/*A1 < CR - 1 && A2 < CR - 1 && A1 != A2 && */_STD is_convertible_v<S2, S>)
	constexpr auto Rotate(Dimension A1, Dimension A2, const Matrix<S, CR, CR>& matrix, S2 radians) -> Matrix<_IMPL CT<S, S2>, CR, CR>;

	template<Scalar S, Dimension CR, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Scale(const Matrix<S, CR, CR>& matrix, const Vector<S2, CR - 1>& scale) -> Matrix<_IMPL CT<S, S2>, CR, CR>;

	template<Scalar S, Dimension CR, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Scale(const Matrix<S, CR, CR>& matrix, S2 scale) -> Matrix<_IMPL CT<S, S2>, CR, CR>;

	template<Scalar S, Dimension C, Dimension R>
	auto operator<<(_STD ostream& ostream, const Matrix<S, C, R>& matrix) -> _STD ostream&;

	// Aliases

	template<Dimension C, Dimension R = C>
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

// Implementation: Don't export.
namespace nd::impl
{
	template<Scalar S, Dimension C, Dimension R>
	auto PrintRow(_STD ostream& ostream, const Matrix<S, C, R>& matrix, _STD streamsize width, Dimension r) noexcept -> void
	{
		ostream << '[' << _STD setw(width) << +matrix[0][r];
		for (Dimension c{1}; c < C; c++)
			ostream << _STD setw(width + 1) << +matrix[c][r];
		ostream << ']';
	};
}

export namespace nd
{
	template<Scalar S, Dimension C, Dimension R>
	constexpr Matrix<S, C, R>::Matrix() noexcept
		: Matrix<S, C, R>{1}
	{

	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2, Dimension C2, Dimension R2> requires(C2 <= C && R2 <= R && _STD is_convertible_v<S2, S>)
	constexpr Matrix<S, C, R>::Matrix(const Matrix<S2, C2, R2>& matrix) noexcept
		: Matrix<S, C, R>{1}
	{
		for (Dimension c{}; c < C2; c++)
			for (Dimension r{}; r < R2; r++)
				m_Scalars[c][r] = S{matrix[c][r]};
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2, _STD enable_if_t<_STD is_convertible_v<S2, S>, int>>
	constexpr Matrix<S, C, R>::Matrix(S2 scalar) noexcept
	{
		if (scalar != S{})
			for (Dimension i{}; i < gcem::min(C, R); i++)
				m_Scalars[i][i] = static_cast<S>(scalar);
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar... S2s, _STD enable_if_t<sizeof...(S2s) == R * C && _IMPL all_true_v<_STD is_convertible_v<S2s, S>...>, int>>
	constexpr Matrix<S, C, R>::Matrix(S2s&&... scalars) noexcept
	{
		Fill<C, R>(_STD forward<S2s>(scalars)...);
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2, _STD enable_if_t<_STD is_convertible_v<S2, S>, int>>
	constexpr auto Matrix<S, C, R>::operator=(S2 scalar) noexcept -> Matrix<S, C, R>&
	{
		return *this = Matrix<S, C, R>(scalar);
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2, Dimension C2, Dimension R2> requires(C2 <= C && R2 <= R && _STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator=(const Matrix<S2, C2, R2>& matrix) noexcept -> Matrix<S, C, R>&
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

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator+(S2 scalar) const noexcept -> Matrix<_IMPL CT<S, S2>, C, R>
	{
		return Matrix<_IMPL CT<S, S2>, C, R>{*this} += _IMPL CT<S, S2>{scalar};
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator+=(S2 scalar) noexcept -> Matrix<S, C, R>&
	{
		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				m_Scalars[c][r] += S{scalar};
		return *this;
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator+(const Matrix<S2, C, R>& matrix) const noexcept -> Matrix<_IMPL CT<S, S2>, C, R>
	{
		return Matrix<_IMPL CT<S, S2>, C, R>{*this} += Matrix<_IMPL CT<S, S2>, C, R>{matrix};
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator+=(const Matrix<S2, C, R>& matrix) noexcept -> Matrix<S, C, R>&
	{
		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				m_Scalars[c][r] += S{matrix[c][r]};
		return *this;
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator-(S2 scalar) const noexcept -> Matrix<_IMPL CT<S, S2>, C, R>
	{
		return Matrix<_IMPL CT<S, S2>, C, R>{*this} -= _IMPL CT<S, S2>{scalar};
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator-=(S2 scalar) noexcept -> Matrix<S, C, R>&
	{
		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				m_Scalars[c][r] -= S{scalar};
		return *this;
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator-(const Matrix<S2, C, R>& matrix) const noexcept -> Matrix<_IMPL CT<S, S2>, C, R>
	{
		return Matrix<_IMPL CT<S, S2>, C, R>{*this} -= Matrix<_IMPL CT<S, S2>, C, R>{matrix};
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator-=(const Matrix<S2, C, R>& matrix) noexcept -> Matrix<S, C, R>&
	{
		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				m_Scalars[c][r] -= S{matrix[c][r]};
		return *this;
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator*(S2 scalar) const noexcept -> Matrix<_IMPL CT<S, S2>, C, R>
	{
		return Matrix<_IMPL CT<S, S2>, C, R>{*this} *= _IMPL CT<S, S2>{scalar};
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator*=(S2 scalar) noexcept -> Matrix<S, C, R>&
	{
		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				m_Scalars[c][r] *= S{scalar};
		return *this;
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2, Dimension C2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator*(const Matrix<S2, C2, C>& matrix) const noexcept -> Matrix<_IMPL CT<S, S2>, C2, R>
	{
		Matrix<_IMPL CT<S, S2>, C2, R> result{_IMPL CT<S, S2>{}};
		for (Dimension i{}; i < C2; i++)
			for (Dimension j{}; j < C; j++)
				result[i] += Vector<_IMPL CT<S, S2>, R>{m_Scalars[j]} * _IMPL CT<S, S2>{matrix[i][j]};
		return result;
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2> requires(R == C && _STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator*=(const Matrix<S2, C, C>& matrix) noexcept -> Matrix<S, C, R>&
	{
		return *this = Matrix<_IMPL CT<S, S2>, C, R>{*this} * _IMPL CT<S, S2>{matrix};
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator/(S2 scalar) const noexcept -> Matrix<_IMPL CT<S, S2>, C, R>
	{
		return Matrix<_IMPL CT<S, S2>, C, R>{*this} /= _IMPL CT<S, S2>{scalar};
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator/=(S2 scalar) noexcept -> Matrix<S, C, R>&
	{
		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				m_Scalars[c][r] /= S{scalar};
		return *this;
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2> requires(R == C && _STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator/(const Matrix<S2, C, C>& matrix) const noexcept -> Matrix<_IMPL CT<S, S2>, C, R>
	{
		return Matrix<_IMPL CT<S, S2>, C, R>{*this} * Inverse(Matrix<_IMPL CT<S, S2>, C, R>{matrix});
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2> requires(R == C && _STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator/=(const Matrix<S2, C, C>& matrix) noexcept -> Matrix<S, C, R>&
	{
		return *this = Matrix<_IMPL CT<S, S2>, C, R>{*this} / Matrix<_IMPL CT<S, S2>, C, C>{matrix};
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Matrix<S, C, R>::operator*(const Vector<S2, C>& vector) const noexcept -> Vector<_IMPL CT<S, S2>, R>
	{
		Vector<_IMPL CT<S, S2>, R> result;
		for (Dimension c{}; c < C; c++)
			result += _IMPL CT<S, S2>{m_Scalars[c]} * _IMPL CT<S, S2>{vector[c]};
		return result;
	}

	template<Scalar S, Dimension C, Dimension R>
	constexpr auto Matrix<S, C, R>::operator+() const noexcept -> Matrix<S, C, R>
	{
		return Matrix<S, C, R>{*this};
	}

	template<Scalar S, Dimension C, Dimension R>
	constexpr auto Matrix<S, C, R>::operator-() const noexcept -> Matrix<S, C, R>
	{
		Matrix<S, C, R> result{+*this};
		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				m_Scalars[c][r] = -m_Scalars[c][r];
		return result;
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Dimension C2> requires(C2 < C)
	constexpr auto Matrix<S, C, R>::Col() const noexcept -> Vector<S, R>
	{
		return m_Scalars[C2];
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Dimension R2> requires(R2 < R)
	constexpr auto Matrix<S, C, R>::Row() const noexcept -> Vector<S, C>
	{
		Vector<S, C> row;
		for (Dimension c{}; c < C; c++)
			row[c] = m_Scalars[c][R2];
		return row;
	}

	template<Scalar S, Dimension C, Dimension R>
	constexpr auto Matrix<S, C, R>::Col(Dimension C2) const noexcept -> Vector<S, R>
	{
		__assume(C2 < C);
		return m_Scalars[C2];
	}

	template<Scalar S, Dimension C, Dimension R>
	constexpr auto Matrix<S, C, R>::Row(Dimension R2) const noexcept -> Vector<S, C>
	{
		__assume(R2 < R);
		Vector<S, C> row;
		for (Dimension c{}; c < C; c++)
			row[c] = m_Scalars[c][R2];
		return row;
	}

	template<Scalar S, Dimension C, Dimension R>
	template<Dimension C2, Dimension R2, Scalar S2, Scalar... Scalars>
	constexpr auto Matrix<S, C, R>::Fill(S2 scalar, Scalars&&... scalars) -> void
	{
		m_Scalars[C - C2][R - R2] = static_cast<S>(scalar);
		if constexpr (C2 > 1)
			Fill<C2 - 1, R2>(_STD forward<Scalars>(scalars)...);
		else if constexpr (R2 > 1)
			Fill<C, R2 - 1>(_STD forward<Scalars>(scalars)...);
	}
	
	template<Scalar S, Dimension C, Dimension R, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator+(S2 scalar, const Matrix<S, C, R>& matrix) noexcept -> Matrix<_IMPL CT<S, S2>, C, R>
	{
		return Matrix<_IMPL CT<S, S2>, C, R>{matrix} += _IMPL CT<S, S2>{scalar};
	}
	
	template<Scalar S, Dimension C, Dimension R, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator-(S2 scalar, const Matrix<S, C, R>& matrix) noexcept -> Matrix<_IMPL CT<S, S2>, C, R>
	{
		return Matrix<_IMPL CT<S, S2>, C, R>{-matrix} += _IMPL CT<S, S2>{scalar};
	}
	
	template<Scalar S, Dimension C, Dimension R, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator*(S2 scalar, const Matrix<S, C, R>& matrix) noexcept -> Matrix<_IMPL CT<S, S2>, C, R>
	{
		return Matrix<_IMPL CT<S, S2>, C, R>{matrix} *= _IMPL CT<S, S2>{scalar};
	}
	
	template<Scalar S, Dimension CR, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator/(S2 scalar, const Matrix<S, CR, CR>& matrix) noexcept -> Matrix<_IMPL CT<S, S2>, CR, CR>
	{
		return Inverse(Matrix<_IMPL CT<S, S2>, CR, CR>{matrix}) *= _IMPL CT<S, S2>{scalar};
	}
	
	template<Dimension C, Dimension R, Scalar S, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto operator*(const Vector<S2, R>& vector, const Matrix<S, C, R>& matrix) noexcept -> Vector<_IMPL CT<S, S2>, C>
	{
		Vector<_IMPL CT<S, S2>, C> result;

		Matrix<_IMPL CT<S, S2>, R, C> m{Transpose(Matrix<_IMPL CT<S, S2>, C, R>{matrix})};
		for (Dimension r{}; r < R; r++)
			result += Vector<_IMPL CT<S, S2>, R>{m[r]} * Vector<_IMPL CT<S, S2>, R>{vector[r]};

		return result;
	}
	
	template<Scalar S, Dimension C, Dimension R, Dimension C2, Dimension R2> requires(C < C2 && R < R2 && C2 > 1 && R2 > 1)
	constexpr auto Submatrix(const Matrix<S, C2, R2>& matrix) noexcept -> Matrix<S, C2 - 1, R2 - 1>
	{
		Matrix<S, C2 - 1, R2 - 1> result{S{}};

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
	
	template<Scalar S, Dimension C2, Dimension R2>
	requires(/*C < C2 && R < R2 && */C2 > 1 && R2 > 1)
	constexpr auto Submatrix(Dimension C, Dimension R, const Matrix<S, C2, R2>& matrix) noexcept -> Matrix<S, C2 - 1, R2 - 1>
	{
		__assume(C < C2 && R < R2);

		Matrix<S, C2 - 1, R2 - 1> result{S{}};

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

	template<Scalar S, Dimension CR>
	constexpr auto Trace(const Matrix<S, CR, CR>& matrix) noexcept -> S
	{
		S result{S{}};
		for (Dimension cr{}; cr < CR; cr++)
			result += matrix[cr][cr];
		return result;
	}

	template<Scalar S, Dimension CR>
	constexpr auto Determinant(const Matrix<S, CR, CR>& matrix) noexcept -> S
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
	constexpr auto Determinant(const Matrix<S, 1, 1>& matrix) noexcept -> S
	{
		return matrix[0][0];
	}
	
	template<Scalar S, Dimension C, Dimension R>
	constexpr auto Transpose(const Matrix<S, C, R>& matrix) noexcept -> Matrix<S, R, C>
	{
		Matrix<S, R, C> result = S();
		for (Dimension c{}; c < C; c++)
			for (Dimension r{}; r < R; r++)
				result[r][c] = matrix[c][r];
		return result;
	}

	template<Scalar S, Dimension CR>
	constexpr auto Minors(const Matrix<S, CR, CR>& matrix) noexcept -> Matrix<S, CR, CR>
	{
		Matrix<S, CR, CR> result = S();
		for (Dimension c{}; c < CR; c++)
			for (Dimension r{}; r < CR; r++)
				result[c][r] = Determinant(Submatrix(c, r, matrix));
		return result;
	}

	template<Scalar S, Dimension CR>
	constexpr auto Cofactors(const Matrix<S, CR, CR>& matrix) noexcept -> Matrix<S, CR, CR>
	{
		Matrix<S, CR, CR> result = Minors(matrix);
		for (Dimension c{}; c < CR; c++)
			for (Dimension r{}; r < CR; r++)
				if ((c + r) % 2 != 0)
					result[c][r] = -result[c][r];
		return result;
	}

	template<Scalar S, Dimension CR>
	constexpr auto Adjugate(const Matrix<S, CR, CR>& matrix) noexcept -> Matrix<S, CR, CR>
	{
		return Transpose(Cofactors(matrix));
	}

	template<Scalar S, Dimension CR>
	constexpr auto Inverse(const Matrix<S, CR, CR>& matrix) noexcept -> Matrix<S, CR, CR>
	{
		return Adjugate(matrix) / Determinant(matrix);
	}
	
	template<Scalar S, Dimension CR, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Inverse(const Matrix<S, CR, CR>& matrix, S2 determinant) noexcept -> Matrix<_IMPL CT<S, S2>, CR, CR>
	{
		return Adjugate(Matrix<_IMPL CT<S, S2>, CR, CR>{matrix}) / _IMPL CT<S, S2>{determinant};
	}
	
	template<Scalar S, Dimension CR, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Translate(const Matrix<S, CR, CR>& matrix, const Vector<S2, CR - 1>& translation) -> Matrix<_IMPL CT<S, S2>, CR, CR>
	{
		Matrix<_IMPL CT<S, S2>, CR, CR> t;
		for (Dimension r{}; r < CR - 1; r++)
			t[CR - 1][r] += _IMPL CT<S, S2>{translation[r]};
		return matrix * t;
	}
	
	template<Dimension A1, Dimension A2, Dimension CR, Scalar S, Scalar S2>
	requires(A1 < CR - 1 && A2 < CR - 1 && A1 != A2 && _STD is_convertible_v<S2, S>)
	constexpr auto Rotate(const Matrix<S, CR, CR>& matrix, S2 radians) -> Matrix<_IMPL CT<S, S2>, CR, CR>
	{
		Matrix<_IMPL CT<S, S2>, CR, CR> rotation;

		auto sin{CT<S, S2>(gcem::sin(CT<S, S2>(radians)))};
		auto cos{CT<S, S2>(gcem::cos(CT<S, S2>(radians)))};

		rotation[A1][A1] = cos; rotation[A2][A1] = -sin;
		rotation[A1][A2] = sin; rotation[A2][A2] = cos;

		return Matrix<_IMPL CT<S, S2>, CR, CR>{matrix} * rotation;
	}
	
	template<Scalar S, Dimension CR, Scalar S2>
	requires(/*A1 < CR - 1 && A2 < CR - 1 && A1 != A2 && */_STD is_convertible_v<S2, S>)
	constexpr auto Rotate(Dimension A1, Dimension A2, const Matrix<S, CR, CR>& matrix, S2 radians) -> Matrix<_IMPL CT<S, S2>, CR, CR>
	{
		__assume(A1 < CR - 1 && A2 < CR - 1 && A1 != A2);

		Matrix<_IMPL CT<S, S2>, CR, CR> rotation;

		auto sin{CT<S, S2>(gcem::sin(CT<S, S2>(radians)))};
		auto cos{CT<S, S2>(gcem::cos(CT<S, S2>(radians)))};

		rotation[A1][A1] = cos; rotation[A2][A1] = -sin;
		rotation[A1][A2] = sin; rotation[A2][A2] = cos;

		return Matrix<_IMPL CT<S, S2>, CR, CR>{matrix} * rotation;
	}
	
	template<Scalar S, Dimension CR, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Scale(const Matrix<S, CR, CR>& matrix, const Vector<S2, CR - 1>& scale) -> Matrix<_IMPL CT<S, S2>, CR, CR>
	{
		Matrix<_IMPL CT<S, S2>, CR, CR> s;
		for (Dimension cr{}; cr < CR - 1; cr++)
			s[cr][cr] = CT<S, S2>(scale[cr]);
		return Matrix<_IMPL CT<S, S2>, CR, CR>{matrix} * s;
	}
	
	template<Scalar S, Dimension CR, Scalar S2> requires(_STD is_convertible_v<S2, S>)
	constexpr auto Scale(const Matrix<S, CR, CR>& matrix, S2 scale) -> Matrix<_IMPL CT<S, S2>, CR, CR>
	{
		Matrix<_IMPL CT<S, S2>, CR, CR> s{_IMPL CT<S, S2>{scale}};
		s[CR - 1][CR - 1] = _IMPL CT<S, S2>{1};
		return Matrix<_IMPL CT<S, S2>, CR, CR>{matrix} * s;
	}

	template<_STD floating_point S, Dimension C, Dimension R>
	auto operator<<(_STD ostream& ostream, const Matrix<S, C, R>& matrix) -> _STD ostream&
	{
		ostream << _STD showpos << _STD scientific;
		_IMPL PrintRow(ostream, matrix, 13, 0);
		for (Dimension r{1}; r < R; r++)
			_IMPL PrintRow(ostream << '\n', matrix, 13, r);
		return ostream << _STD defaultfloat << _STD noshowpos;
	}

	template<_STD integral S, Dimension C, Dimension R>
	auto operator<<(_STD ostream& ostream, const Matrix<S, C, R>& matrix) -> _STD ostream&
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
