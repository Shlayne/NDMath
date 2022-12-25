module;

#include <gcem.hpp>

export module nd.matrix;

import nd.scalar;
import nd.vector;
import std.core;

namespace nd
{
	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	struct Matrix;
}

template<nd::Scalar... Scalars>
using CT = std::common_type_t<Scalars...>;

// This is so small, I'm fine with duplicating it.
template<size_t N, nd::Scalar... Scalars>
using CVT = nd::Vector<N, CT<Scalars...>>;

template<size_t C, size_t R, nd::Scalar... Scalars>
using CMT = nd::Matrix<C, R, CT<Scalars...>>;

// Implementation: Don't export.
namespace impl
{
	// https://stackoverflow.com/questions/31767645/how-to-check-that-all-types-in-variadic-template-are-convertible-to-size-t

	template<bool...>
	struct bools {};

	template<bool B, bool... Bools>
	constexpr bool all_same_bool_v = std::is_same_v<bools<B, Bools...>, bools<Bools..., B>>;

	template<bool... Bools>
	constexpr bool all_true_v = all_same_bool_v<true, Bools...>;
}

export namespace nd
{
	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	struct Matrix
	{
	public:
		constexpr Matrix() noexcept;

		template<size_t C2, size_t R2, Scalar S2> requires(C2 <= C && R2 <= R && std::is_convertible_v<S2, S>)
		constexpr Matrix(const Matrix<C2, R2, S2>& matrix) noexcept;

		// Collides with copy constructor if requires is used instead of std::enable_if_t.
		template<Scalar S2, std::enable_if_t<std::is_convertible_v<S2, S>, int> = 0>
		constexpr Matrix(const S2& scalar) noexcept;

		// Collides with copy constructor if requires is used instead of std::enable_if_t.
		template<Scalar... S2s, std::enable_if_t<sizeof...(S2s) == R * C && impl::all_true_v<std::is_convertible_v<S2s, S>...>, int> = 0>
		constexpr Matrix(S2s&&... scalars) noexcept;

		// Collides with copy operator= if requires is used instead of std::enable_if_t.
		template<Scalar S2, std::enable_if_t<std::is_convertible_v<S2, S>, int> = 0>
		constexpr auto operator=(const S2& scalar) noexcept -> Matrix<C, R, S>&;

		template<size_t C2, size_t R2, Scalar S2> requires(C2 <= C && R2 <= R && std::is_convertible_v<S2, S>)
		constexpr auto operator=(const Matrix<C2, R2, S2>& matrix) noexcept -> Matrix<C, R, S>&;
	public:
		template<Scalar S2> requires(std::is_convertible_v<S2, S>)
		constexpr auto operator+(const S2& scalar) const noexcept -> CMT<C, R, S, S2>;

		template<Scalar S2> requires(std::is_convertible_v<S2, S>)
		constexpr auto operator+=(const S2& scalar) noexcept -> Matrix<C, R, S>&;

		template<Scalar S2> requires(std::is_convertible_v<S2, S>)
		constexpr auto operator+(const Matrix<C, R, S2>& matrix) const noexcept -> CMT<C, R, S, S2>;

		template<Scalar S2> requires(std::is_convertible_v<S2, S>)
		constexpr auto operator+=(const Matrix<C, R, S2>& matrix) noexcept -> Matrix<C, R, S>&;
	public:
		template<Scalar S2> requires(std::is_convertible_v<S2, S>)
		constexpr auto operator-(const S2& scalar) const noexcept -> CMT<C, R, S, S2>;

		template<Scalar S2> requires(std::is_convertible_v<S2, S>)
		constexpr auto operator-=(const S2& scalar) noexcept -> Matrix<C, R, S>&;

		template<Scalar S2> requires(std::is_convertible_v<S2, S>)
		constexpr auto operator-(const Matrix<C, R, S2>& matrix) const noexcept -> CMT<C, R, S, S2>;

		template<Scalar S2> requires(std::is_convertible_v<S2, S>)
		constexpr auto operator-=(const Matrix<C, R, S2>& matrix) noexcept -> Matrix<C, R, S>&;
	public:
		template<Scalar S2> requires(std::is_convertible_v<S2, S>)
		constexpr auto operator*(const S2& scalar) const noexcept -> CMT<C, R, S, S2>;

		template<Scalar S2> requires(std::is_convertible_v<S2, S>)
		constexpr auto operator*=(const S2& scalar) noexcept -> Matrix<C, R, S>&;

		template<size_t C2, Scalar S2> requires(std::is_convertible_v<S2, S>)
		constexpr auto operator*(const Matrix<C2, C, S2>& matrix) const noexcept -> CMT<C2, R, S, S2>;

		template<Scalar S2> requires(R == C && std::is_convertible_v<S2, S>)
		constexpr auto operator*=(const Matrix<C, C, S2>& matrix) noexcept -> Matrix<C, R, S>&;
	public:
		template<Scalar S2> requires(std::is_convertible_v<S2, S>)
		constexpr auto operator/(const S2& scalar) const noexcept -> CMT<C, R, S, S2>;

		template<Scalar S2> requires(std::is_convertible_v<S2, S>)
		constexpr auto operator/=(const S2& scalar) noexcept -> Matrix<C, R, S>&;

		template<Scalar S2> requires(R == C && std::is_convertible_v<S2, S>)
		constexpr auto operator/(const Matrix<C, C, S2>& matrix) const noexcept -> CMT<C, R, S, S2>;

		template<Scalar S2> requires(R == C && std::is_convertible_v<S2, S>)
		constexpr auto operator/=(const Matrix<C, C, S2>& matrix) noexcept -> Matrix<C, R, S>&;
	public:
		template<Scalar S2> requires(std::is_convertible_v<S2, S>)
		constexpr auto operator*(const Vector<C, S2>& vector) const noexcept -> CVT<R, S, S2>;
	public:
		constexpr auto operator+() const noexcept -> Matrix<C, R, S>;
		constexpr auto operator-() const noexcept -> Matrix<C, R, S>;
	public:
		constexpr auto operator[](size_t col) noexcept -> Vector<R, S>&;
		constexpr auto operator[](size_t col) const noexcept -> const Vector<R, S>&;

		template<size_t C2> requires(C2 < C)
		constexpr auto Col() const noexcept -> Vector<C, S>;

		template<size_t R2> requires(R2 < R)
		constexpr auto Row() const noexcept -> Vector<C, S>;
	private:
		template<size_t C2, size_t R2, Scalar S2, Scalar... Scalars>
		constexpr auto Fill(const S2& scalar, Scalars&&... scalars) -> void;
	private:
		Vector<R, S> m_Scalars[C]{ Vector<R, S>() };
	};

	// External Operators

	template<size_t C, size_t R, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto operator+(const S2& scalar, const Matrix<C, R, S>& matrix) noexcept -> CMT<C, R, S, S2>;

	template<size_t C, size_t R, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto operator-(const S2& scalar, const Matrix<C, R, S>& matrix) noexcept -> CMT<C, R, S, S2>;

	template<size_t C, size_t R, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto operator*(const S2& scalar, const Matrix<C, R, S>& matrix) noexcept -> CMT<C, R, S, S2>;

	template<size_t CR, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto operator/(const S2& scalar, const Matrix<CR, CR, S>& matrix) noexcept -> CMT<CR, CR, S, S2>;

	template<size_t C, size_t R, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto operator*(const Vector<R, S2>& vector, const Matrix<C, R, S>& matrix) noexcept -> CVT<C, S, S2>;

	// Static Methods

	template<size_t C, size_t R, Scalar S, size_t C2, size_t R2> requires(C < C2 && R < R2 && C2 > 1 && R2 > 1)
	constexpr auto Submatrix(const Matrix<C2, R2, S>& matrix) noexcept -> Matrix<C2 - 1, R2 - 1, S>;

	template<size_t C2, size_t R2, Scalar S>
	constexpr auto Submatrix(size_t C, size_t R, const Matrix<C2, R2, S>& matrix) noexcept -> Matrix<C2 - 1, R2 - 1, S>
		requires(C < C2 && R < R2 && C2 > 1 && R2 > 1);

	template<size_t CR, Scalar S>
	constexpr auto Trace(const Matrix<CR, CR, S>& matrix) noexcept -> S;

	template<size_t CR, Scalar S>
	constexpr auto Determinant(const Matrix<CR, CR, S>& matrix) noexcept -> S;

	template<size_t C, size_t R, Scalar S>
	constexpr auto Transpose(const Matrix<C, R, S>& matrix) noexcept -> Matrix<R, C, S>;

	template<size_t CR, Scalar S>
	constexpr auto Minors(const Matrix<CR, CR, S>& matrix) noexcept -> Matrix<CR, CR, S>;

	template<size_t CR, Scalar S>
	constexpr auto Cofactors(const Matrix<CR, CR, S>& matrix) noexcept -> Matrix<CR, CR, S>;

	template<size_t CR, Scalar S>
	constexpr auto Adjugate(const Matrix<CR, CR, S>& matrix) noexcept -> Matrix<CR, CR, S>;

	template<size_t CR, Scalar S>
	constexpr auto Inverse(const Matrix<CR, CR, S>& matrix) noexcept -> Matrix<CR, CR, S>;

	template<size_t CR, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Inverse(const Matrix<CR, CR, S>& matrix, const S2& determinant) noexcept -> CMT<CR, CR, S, S2>;

	template<size_t CR, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Translate(const Matrix<CR, CR, S>& matrix, const Vector<CR - 1, S2>& translation) -> CMT<CR, CR, S, S2>;

	template<size_t A1, size_t A2, size_t CR, Scalar S, Scalar S2>
		requires(A1 < CR - 1 && A2 < CR - 1 && A1 != A2 && std::is_convertible_v<S2, S>)
	constexpr auto Rotate(const Matrix<CR, CR, S>& matrix, const S2& radians) -> CMT<CR, CR, S, S2>;

	template<size_t CR, Scalar S, Scalar S2>
	constexpr auto Rotate(size_t A1, size_t A2, const Matrix<CR, CR, S>& matrix, const S2& radians) -> CMT<CR, CR, S, S2>
		requires(A1 < CR - 1 && A2 < CR - 1 && A1 != A2 && std::is_convertible_v<S2, S>);

	template<size_t CR, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Scale(const Matrix<CR, CR, S>& matrix, const Vector<CR - 1, S2>& scale) -> CMT<CR, CR, S, S2>;

	template<size_t CR, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Scale(const Matrix<CR, CR, S>& matrix, const S2& scale) -> CMT<CR, CR, S, S2>;

	template<size_t C, size_t R, Scalar S>
	auto operator<<(std::ostream& ostream, const Matrix<C, R, S>& matrix) -> std::ostream&;

	// Aliases

	template<size_t C, size_t R = C>
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
namespace impl
{
	using namespace ::nd;

	template<size_t C, size_t R, Scalar S>
	auto PrintRow(std::ostream& ostream, const Matrix<C, R, S>& matrix, long long width, size_t r) noexcept -> void
	{
		ostream << '[' << std::setw(width) << +matrix[0][r];
		for (size_t c = 1; c < C; c++)
			ostream << std::setw(width + 1) << +matrix[c][r];
		ostream << ']';
	};
}

export namespace nd
{
	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	constexpr Matrix<C, R, S>::Matrix() noexcept
		: Matrix<C, R, S>(1)
	{

	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<size_t C2, size_t R2, Scalar S2> requires(C2 <= C && R2 <= R && std::is_convertible_v<S2, S>)
	constexpr Matrix<C, R, S>::Matrix(const Matrix<C2, R2, S2>& matrix) noexcept
		: Matrix<C, R, S>(1)
	{
		for (size_t c = 0; c < C2; c++)
			for (size_t r = 0; r < R2; r++)
				m_Scalars[c][r] = S(matrix[c][r]);
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2, std::enable_if_t<std::is_convertible_v<S2, S>, int>>
	constexpr Matrix<C, R, S>::Matrix(const S2& scalar) noexcept
	{
		if (scalar != S())
			for (size_t i = 0; i < gcem::min(C, R); i++)
				m_Scalars[i][i] = S(scalar);
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar... S2s, std::enable_if_t<sizeof...(S2s) == R * C && impl::all_true_v<std::is_convertible_v<S2s, S>...>, int>>
	constexpr Matrix<C, R, S>::Matrix(S2s&&... scalars) noexcept
	{
		Fill<C, R>(std::forward<S2s>(scalars)...);
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2, std::enable_if_t<std::is_convertible_v<S2, S>, int>>
	constexpr auto Matrix<C, R, S>::operator=(const S2& scalar) noexcept -> Matrix<C, R, S>&
	{
		return *this = Matrix<C, R, S>(scalar);
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<size_t C2, size_t R2, Scalar S2> requires(C2 <= C && R2 <= R && std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator=(const Matrix<C2, R2, S2>& matrix) noexcept -> Matrix<C, R, S>&
	{
		if (this != (void*)&matrix)
		{
			for (size_t c = 0; c < C2; c++)
				for (size_t r = 0; r < R2; r++)
					m_Scalars[c][r] = S(matrix[c][r]);

			for (size_t c = C2; c < C; c++)
				for (size_t r = 0; r < R2 + c - C2; r++)
					m_Scalars[c][r] = S();

			for (size_t r = R2; r < R; r++)
				for (size_t c = 0; c < C2 + r - R2; c++)
					m_Scalars[c][r] = S();

			for (size_t cr = gcem::min(C2, R2); cr < gcem::min(C, R); cr++)
				m_Scalars[cr][cr] = S(1);
		}
		return *this;
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator+(const S2& scalar) const noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>(*this) += CT<S, S2>(scalar);
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator+=(const S2& scalar) noexcept -> Matrix<C, R, S>&
	{
		for (size_t c = 0; c < C; c++)
			for (size_t r = 0; r < R; r++)
				m_Scalars[c][r] += S(scalar);
		return *this;
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator+(const Matrix<C, R, S2>& matrix) const noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>(*this) += CMT<C, R, S, S2>(matrix);
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator+=(const Matrix<C, R, S2>& matrix) noexcept -> Matrix<C, R, S>&
	{
		for (size_t c = 0; c < C; c++)
			for (size_t r = 0; r < R; r++)
				m_Scalars[c][r] += S(matrix[c][r]);
		return *this;
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator-(const S2& scalar) const noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>(*this) -= CT<S, S2>(scalar);
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator-=(const S2& scalar) noexcept -> Matrix<C, R, S>&
	{
		for (size_t c = 0; c < C; c++)
			for (size_t r = 0; r < R; r++)
				m_Scalars[c][r] -= S(scalar);
		return *this;
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator-(const Matrix<C, R, S2>& matrix) const noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>(*this) -= CMT<C, R, S, S2>(matrix);
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator-=(const Matrix<C, R, S2>& matrix) noexcept -> Matrix<C, R, S>&
	{
		for (size_t c = 0; c < C; c++)
			for (size_t r = 0; r < R; r++)
				m_Scalars[c][r] -= S(matrix[c][r]);
		return *this;
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator*(const S2& scalar) const noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>(*this) *= CT<S, S2>(scalar);
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator*=(const S2& scalar) noexcept -> Matrix<C, R, S>&
	{
		for (size_t c = 0; c < C; c++)
			for (size_t r = 0; r < R; r++)
				m_Scalars[c][r] *= S(scalar);
		return *this;
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<size_t C2, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator*(const Matrix<C2, C, S2>& matrix) const noexcept -> CMT<C2, R, S, S2>
	{
		CMT<C2, R, S, S2> result = CT<S, S2>();
		for (size_t i = 0; i < C2; i++)
			for (size_t j = 0; j < C; j++)
				result[i] += CT<S, S2>(m_Scalars[j]) * CT<S, S2>(matrix[i][j]);
		return result;
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(R == C && std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator*=(const Matrix<C, C, S2>& matrix) noexcept -> Matrix<C, R, S>&
	{
		return *this = CMT<C, R, S, S2>(*this) * CT<S, S2>(matrix);
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator/(const S2& scalar) const noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>(*this) /= CT<S, S2>(scalar);
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator/=(const S2& scalar) noexcept -> Matrix<C, R, S>&
	{
		for (size_t c = 0; c < C; c++)
			for (size_t r = 0; r < R; r++)
				m_Scalars[c][r] /= S(scalar);
		return *this;
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(R == C && std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator/(const Matrix<C, C, S2>& matrix) const noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>(*this) * Inverse(CMT<C, R, S, S2>(matrix));
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(R == C && std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator/=(const Matrix<C, C, S2>& matrix) noexcept -> Matrix<C, R, S>&
	{
		return *this = CMT<C, R, S, S2>(*this) / CMT<C, C, S, S2>(matrix);
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Matrix<C, R, S>::operator*(const Vector<C, S2>& vector) const noexcept -> CVT<R, S, S2>
	{
		CVT<R, S, S2> result;
		for (size_t c = 0; c < C; c++)
			result += CT<S, S2>(m_Scalars[c]) * CT<S, S2>(vector[c]);
		return result;
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	constexpr auto Matrix<C, R, S>::operator+() const noexcept -> Matrix<C, R, S>
	{
		return Matrix<C, R, S>(*this);
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	constexpr auto Matrix<C, R, S>::operator-() const noexcept -> Matrix<C, R, S>
	{
		Matrix<C, R, S> result = +*this;
		for (size_t c = 0; c < C; c++)
			for (size_t r = 0; r < R; r++)
				m_Scalars[c][r] = -m_Scalars[c][r];
		return result;
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	constexpr auto Matrix<C, R, S>::operator[](size_t col) noexcept -> Vector<R, S>&
	{
		return m_Scalars[col];
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	constexpr auto Matrix<C, R, S>::operator[](size_t col) const noexcept -> const Vector<R, S>&
	{
		return m_Scalars[col];
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<size_t C2> requires(C2 < C)
	constexpr auto Matrix<C, R, S>::Col() const noexcept -> Vector<C, S>
	{
		return m_Scalars[C2];
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<size_t R2> requires(R2 < R)
	constexpr auto Matrix<C, R, S>::Row() const noexcept -> Vector<C, S>
	{
		Vector<C, S> row;
		for (size_t c = 0; c < C; c++)
			row[c] = m_Scalars[c][R2];
		return row;
	}

	template<size_t C, size_t R, Scalar S> requires(C > 0 && R > 0)
	template<size_t C2, size_t R2, Scalar S2, Scalar... Scalars>
	constexpr auto Matrix<C, R, S>::Fill(const S2& scalar, Scalars&&... scalars) -> void
	{
		m_Scalars[C - C2][R - R2] = S(scalar);
		if constexpr (C2 > 1)
			Fill<C2 - 1, R2>(std::forward<Scalars>(scalars)...);
		else if constexpr (R2 > 1)
			Fill<C, R2 - 1>(std::forward<Scalars>(scalars)...);
	}
	
	template<size_t C, size_t R, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto operator+(const S2& scalar, const Matrix<C, R, S>& matrix) noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>(matrix) += CT<S, S2>(scalar);
	}
	
	template<size_t C, size_t R, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto operator-(const S2& scalar, const Matrix<C, R, S>& matrix) noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>(-matrix) += CT<S, S2>(scalar);
	}
	
	template<size_t C, size_t R, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto operator*(const S2& scalar, const Matrix<C, R, S>& matrix) noexcept -> CMT<C, R, S, S2>
	{
		return CMT<C, R, S, S2>(matrix) *= CT<S, S2>(scalar);
	}
	
	template<size_t CR, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto operator/(const S2& scalar, const Matrix<CR, CR, S>& matrix) noexcept -> CMT<CR, CR, S, S2>
	{
		return Inverse(CMT<CR, CR, S, S2>(matrix)) *= CT<S, S2>(scalar);
	}
	
	template<size_t C, size_t R, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto operator*(const Vector<R, S2>& vector, const Matrix<C, R, S>& matrix) noexcept -> CVT<C, S, S2>
	{
		CVT<C, S, S2> result;

		CMT<R, C, S, S2> m = Transpose(CMT<C, R, S, S2>(matrix));
		for (size_t r = 0; r < R; r++)
			result += CVT<R, S, S2>(m[r]) * CVT<R, S, S2>(vector[r]);

		return result;
	}
	
	template<size_t C, size_t R, Scalar S, size_t C2, size_t R2> requires(C < C2 && R < R2 && C2 > 1 && R2 > 1)
	constexpr auto Submatrix(const Matrix<C2, R2, S>& matrix) noexcept -> Matrix<C2 - 1, R2 - 1, S>
	{
		Matrix<C2 - 1, R2 - 1, S> result = S();

		for (size_t c = 0; c < C; c++)
			for (size_t r = 0; r < R; r++)
				result[c][r] = matrix[c][r];

		for (size_t c = C + 1; c < C2; c++)
			for (size_t r = 0; r < R; r++)
				result[c - 1][r] = matrix[c][r];

		for (size_t c = 0; c < C; c++)
			for (size_t r = R + 1; r < R2; r++)
				result[c][r - 1] = matrix[c][r];

		for (size_t c = C + 1; c < C2; c++)
			for (size_t r = R + 1; r < R2; r++)
				result[c - 1][r - 1] = matrix[c][r];

		return result;
	}
	
	template<size_t C2, size_t R2, Scalar S>
	constexpr auto Submatrix(size_t C, size_t R, const Matrix<C2, R2, S>& matrix) noexcept -> Matrix<C2 - 1, R2 - 1, S>
		requires(C < C2 && R < R2 && C2 > 1 && R2 > 1)
	{
		//assert((C < C2) && (R < R2));

		Matrix<C2 - 1, R2 - 1, S> result = S();

		for (size_t c = 0; c < C; c++)
			for (size_t r = 0; r < R; r++)
				result[c][r] = matrix[c][r];

		for (size_t c = C + 1; c < C2; c++)
			for (size_t r = 0; r < R; r++)
				result[c - 1][r] = matrix[c][r];

		for (size_t c = 0; c < C; c++)
			for (size_t r = R + 1; r < R2; r++)
				result[c][r - 1] = matrix[c][r];

		for (size_t c = C + 1; c < C2; c++)
			for (size_t r = R + 1; r < R2; r++)
				result[c - 1][r - 1] = matrix[c][r];

		return result;
	}

	template<size_t CR, Scalar S>
	constexpr auto Trace(const Matrix<CR, CR, S>& matrix) noexcept -> S
	{
		S result = S();
		for (size_t cr = 0; cr < CR; cr++)
			result += matrix[cr][cr];
		return result;
	}

	template<size_t CR, Scalar S>
	constexpr auto Determinant(const Matrix<CR, CR, S>& matrix) noexcept -> S
	{
		S result = S();

		for (size_t c = 0; c < CR; c++)
		{
			S subDeterminant = matrix[c][0] * Determinant(Submatrix(c, 0, matrix));
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
	
	template<size_t C, size_t R, Scalar S>
	constexpr auto Transpose(const Matrix<C, R, S>& matrix) noexcept -> Matrix<R, C, S>
	{
		Matrix<R, C, S> result = S();
		for (size_t c = 0; c < C; c++)
			for (size_t r = 0; r < R; r++)
				result[r][c] = matrix[c][r];
		return result;
	}

	template<size_t CR, Scalar S>
	constexpr auto Minors(const Matrix<CR, CR, S>& matrix) noexcept -> Matrix<CR, CR, S>
	{
		Matrix<CR, CR, S> result = S();
		for (size_t c = 0; c < CR; c++)
			for (size_t r = 0; r < CR; r++)
				result[c][r] = Determinant(Submatrix(c, r, matrix));
		return result;
	}

	template<size_t CR, Scalar S>
	constexpr auto Cofactors(const Matrix<CR, CR, S>& matrix) noexcept -> Matrix<CR, CR, S>
	{
		Matrix<CR, CR, S> result = Minors(matrix);
		for (size_t c = 0; c < CR; c++)
			for (size_t r = 0; r < CR; r++)
				if ((c + r) % 2 != 0)
					result[c][r] = -result[c][r];
		return result;
	}

	template<size_t CR, Scalar S>
	constexpr auto Adjugate(const Matrix<CR, CR, S>& matrix) noexcept -> Matrix<CR, CR, S>
	{
		return Transpose(Cofactors(matrix));
	}

	template<size_t CR, Scalar S>
	constexpr auto Inverse(const Matrix<CR, CR, S>& matrix) noexcept -> Matrix<CR, CR, S>
	{
		return Adjugate(matrix) / Determinant(matrix);
	}
	
	template<size_t CR, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Inverse(const Matrix<CR, CR, S>& matrix, const S2& determinant) noexcept -> CMT<CR, CR, S, S2>
	{
		return Adjugate(CMT<CR, CR, S, S2>(matrix)) / CT<S, S2>(determinant);
	}
	
	template<size_t CR, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Translate(const Matrix<CR, CR, S>& matrix, const Vector<CR - 1, S2>& translation) -> CMT<CR, CR, S, S2>
	{
		CMT<CR, CR, S, S2> t;
		for (size_t r = 0; r < CR - 1; r++)
			t[CR - 1][r] += CT<S, S2>(translation[r]);
		return matrix * t;
	}
	
	template<size_t A1, size_t A2, size_t CR, Scalar S, Scalar S2>
		requires(A1 < CR - 1 && A2 < CR - 1 && A1 != A2 && std::is_convertible_v<S2, S>)
	constexpr auto Rotate(const Matrix<CR, CR, S>& matrix, const S2& radians) -> CMT<CR, CR, S, S2>
	{
		CMT<CR, CR, S, S2> rotation;

		auto sin = CT<S, S2>(gcem::sin(CT<S, S2>(radians)));
		auto cos = CT<S, S2>(gcem::cos(CT<S, S2>(radians)));

		rotation[A1][A1] = cos; rotation[A2][A1] = -sin;
		rotation[A1][A2] = sin; rotation[A2][A2] = cos;

		return CMT<CR, CR, S, S2>(matrix) * rotation;
	}
	
	template<size_t CR, Scalar S, Scalar S2>
	constexpr auto Rotate(size_t A1, size_t A2, const Matrix<CR, CR, S>& matrix, const S2& radians) -> CMT<CR, CR, S, S2>
		requires(A1 < CR - 1 && A2 < CR - 1 && A1 != A2 && std::is_convertible_v<S2, S>)
	{
		//assert(A1 < CR - 1 && A2 < CR - 1 && A1 != A2);

		CMT<CR, CR, S, S2> rotation;

		auto sin = CT<S, S2>(gcem::sin(CT<S, S2>(radians)));
		auto cos = CT<S, S2>(gcem::cos(CT<S, S2>(radians)));

		rotation[A1][A1] = cos; rotation[A2][A1] = -sin;
		rotation[A1][A2] = sin; rotation[A2][A2] = cos;

		return CMT<CR, CR, S, S2>(matrix) * rotation;
	}
	
	template<size_t CR, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Scale(const Matrix<CR, CR, S>& matrix, const Vector<CR - 1, S2>& scale) -> CMT<CR, CR, S, S2>
	{
		CMT<CR, CR, S, S2> s;
		for (size_t cr = 0; cr < CR - 1; cr++)
			s[cr][cr] = CT<S, S2>(scale[cr]);
		return CMT<CR, CR, S, S2>(matrix) * s;
	}
	
	template<size_t CR, Scalar S, Scalar S2> requires(std::is_convertible_v<S2, S>)
	constexpr auto Scale(const Matrix<CR, CR, S>& matrix, const S2& scale) -> CMT<CR, CR, S, S2>
	{
		CMT<CR, CR, S, S2> s = CT<S, S2>(scale);
		s[CR - 1][CR - 1] = CT<S, S2>(1);
		return CMT<CR, CR, S, S2>(matrix) * s;
	}

	template<size_t C, size_t R, std::floating_point S>
	auto operator<<(std::ostream& ostream, const Matrix<C, R, S>& matrix) -> std::ostream&
	{
		ostream << std::showpos << std::scientific;
		impl::PrintRow(ostream, matrix, 13, 0);
		for (size_t r = 1; r < R; r++)
			impl::PrintRow(ostream << '\n', matrix, 13, r);
		return ostream << std::defaultfloat << std::noshowpos;
	}

	template<size_t C, size_t R, std::integral S>
	auto operator<<(std::ostream& ostream, const Matrix<C, R, S>& matrix) -> std::ostream&
	{
		constexpr long long width = std::is_signed_v<S> + (long long)(gcem::max(
			S(gcem::ceil(gcem::log10(std::make_unsigned_t<S>(std::numeric_limits<S>::min())))),
			S(gcem::ceil(gcem::log10(std::make_unsigned_t<S>(std::numeric_limits<S>::max()))))
		));

		impl::PrintRow(ostream, matrix, width, 0);
		for (size_t r = 1; r < R; r++)
			impl::PrintRow(ostream << '\n', matrix, width, r);
		return ostream;
	}
}
