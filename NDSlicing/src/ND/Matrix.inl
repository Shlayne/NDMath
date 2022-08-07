namespace nd
{
	template<typename T, uint32_t C, uint32_t R>
	inline constexpr Matrix<T, C, R>::Matrix() noexcept
		: Matrix<T, C, R>(1)
	{

	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R>::Matrix(const T2& crValue) noexcept
	{
		if (crValue != T(0))
			for (uint32_t i = 0; i < gcem::min(C, R); i++)
				m_pElements[i][i] = T(crValue);
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename... Args, std::enable_if_t<std::ext::are_convertible_v<T, Args...> && sizeof...(Args) == C * R && C * R != 1, int>>
	inline constexpr Matrix<T, C, R>::Matrix(Args&&... rrArgs) noexcept
	{
		Fill<C, R>(std::forward<Args>(rrArgs)...);
	}

	template<typename T, uint32_t C, uint32_t R>
	template<uint32_t C2, uint32_t R2, typename T2, std::enable_if_t<std::is_convertible_v<T2, T> && (C2 <= C) && (R2 <= R), int>>
	inline constexpr Matrix<T, C, R>::Matrix(const Matrix<T2, C2, R2>& crMatrix) noexcept
		: Matrix<T, C, R>(1)
	{
		for (uint32_t c = 0; c < C2; c++)
			for (uint32_t r = 0; r < R2; r++)
				m_pElements[c][r] = T(crMatrix.m_pElements[c][r]);
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R>& Matrix<T, C, R>::operator=(const T2& crValue) noexcept
	{
		return *this = Matrix<T, C, R>(crValue);
	}

	template<typename T, uint32_t C, uint32_t R>
	template<uint32_t C2, uint32_t R2, typename T2, std::enable_if_t<std::is_convertible_v<T2, T> && (C2 <= C) && (R2 <= R), int>>
	inline constexpr Matrix<T, C, R>& Matrix<T, C, R>::operator=(const Matrix<T2, C2, R2>& crMatrix) noexcept
	{
		if (this != (void*)&crMatrix)
		{
			for (uint32_t c = 0; c < C2; c++)
				for (uint32_t r = 0; r < R2; r++)
					m_pElements[c][r] = T(crMatrix.m_pElements[c][r]);

			if constexpr (C2 < C)
				for (uint32_t c = C2; c < C; c++)
					for (uint32_t r = 0; r < R2 + c - C2; r++)
						m_pElements[c][r] = T(0);

			if constexpr (R2 < R)
				for (uint32_t r = R2; r < R; r++)
					for (uint32_t c = 0; c < C2 + r - R2; c++)
						m_pElements[c][r] = T(0);

			if constexpr (C2 < C && R2 < R)
				for (uint32_t cr = gcem::min(C2, R2); cr < gcem::min(C, R); cr++)
					m_pElements[cr][cr] = T(1);
		}
		return *this;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R> Matrix<T, C, R>::operator+(const T2& crValue) const noexcept
	{
		return +*this += crValue;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R>& Matrix<T, C, R>::operator+=(const T2& crValue) noexcept
	{
		for (uint32_t c = 0; c < C; c++)
			for (uint32_t r = 0; r < R; r++)
				m_pElements[c][r] += T(crValue);
		return *this;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R> Matrix<T, C, R>::operator+(const Matrix<T2, C, R>& crMatrix) const noexcept
	{
		return +*this += crMatrix;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R>& Matrix<T, C, R>::operator+=(const Matrix<T2, C, R>& crMatrix) noexcept
	{
		for (uint32_t c = 0; c < C; c++)
			for (uint32_t r = 0; r < R; r++)
				m_pElements[c][r] += T(crMatrix[c][r]);
		return *this;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R> Matrix<T, C, R>::operator-(const T2& crValue) const noexcept
	{
		return +*this -= crValue;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R>& Matrix<T, C, R>::operator-=(const T2& crValue) noexcept
	{
		for (uint32_t c = 0; c < C; c++)
			for (uint32_t r = 0; r < R; r++)
				m_pElements[c][r] -= T(crValue);
		return *this;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R> Matrix<T, C, R>::operator-(const Matrix<T2, C, R>& crMatrix) const noexcept
	{
		return +*this -= crMatrix;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R>& Matrix<T, C, R>::operator-=(const Matrix<T2, C, R>& crMatrix) noexcept
	{
		for (uint32_t c = 0; c < C; c++)
			for (uint32_t r = 0; r < R; r++)
				m_pElements[c][r] -= T(crMatrix[c][r]);
		return *this;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R> Matrix<T, C, R>::operator*(const T2& crValue) const noexcept
	{
		return +*this *= crValue;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R>& Matrix<T, C, R>::operator*=(const T2& crValue) noexcept
	{
		for (uint32_t c = 0; c < C; c++)
			for (uint32_t r = 0; r < R; r++)
				m_pElements[c][r] *= T(crValue);
		return *this;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, uint32_t C2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C2, R> Matrix<T, C, R>::operator*(const Matrix<T2, C2, C>& crMatrix) const noexcept
	{
		Matrix<T, C2, R> result = T(0);
		for (uint32_t i = 0; i < C2; i++)
			for (uint32_t j = 0; j < C; j++)
				result[i] += m_pElements[j] * crMatrix.m_pElements[i][j];
		return result;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T> && C == R, int>>
	inline constexpr Matrix<T, C, R>& Matrix<T, C, R>::operator*=(const Matrix<T2, C, R>& crMatrix) noexcept
	{
		return *this = *this * crMatrix;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R> Matrix<T, C, R>::operator/(const T2& crValue) const noexcept
	{
		return +*this /= crValue;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R>& Matrix<T, C, R>::operator/=(const T2& crValue) noexcept
	{
		for (uint32_t c = 0; c < C; c++)
			for (uint32_t r = 0; r < R; r++)
				m_pElements[c][r] /= T(crValue);
		return *this;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R> Matrix<T, C, R>::operator/(const Matrix<T2, C, C>& crMatrix) const noexcept
	{
		return +*this * Inverse(crMatrix);
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T> && C == R, int>>
	inline constexpr Matrix<T, C, R>& Matrix<T, C, R>::operator/=(const Matrix<T2, C, C>& crMatrix) noexcept
	{
		return *this = *this / crMatrix;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, R> Matrix<T, C, R>::operator*(const Vector<T2, C>& crVector) const noexcept
	{
		Vector<T, R> result;

		for (uint32_t c = 0; c < C; c++)
			result += m_pElements[c] * crVector[c];

		return result;
	}

	template<typename T, uint32_t C, uint32_t R>
	inline constexpr Matrix<T, C, R> Matrix<T, C, R>::operator+() const noexcept
	{
		return Matrix<T, C, R>(*this);
	}

	template<typename T, uint32_t C, uint32_t R>
	inline constexpr Matrix<T, C, R> Matrix<T, C, R>::operator-() const noexcept
	{
		Matrix<T, C, R> result = +*this;
		for (uint32_t c = 0; c < C; c++)
			for (uint32_t r = 0; r < R; r++)
				m_pElements[c][r] = -m_pElements[c][r];
		return result;
	}

	template<typename T, uint32_t C, uint32_t R>
	inline constexpr Vector<T, R>& Matrix<T, C, R>::operator[](uint32_t col) noexcept
	{
		return m_pElements[col];
	}

	template<typename T, uint32_t C, uint32_t R>
	inline constexpr const Vector<T, R>& Matrix<T, C, R>::operator[](uint32_t col) const noexcept
	{
		return m_pElements[col];
	}

	template<typename T, uint32_t C, uint32_t R>
	template<uint32_t C2, std::enable_if_t<(C2 < C), int>>
	inline constexpr Vector<T, R> Matrix<T, C, R>::Col() const noexcept
	{
		return m_pElements[C2];
	}

	template<typename T, uint32_t C, uint32_t R>
	template<uint32_t R2, std::enable_if_t<(R2 < R), int>>
	inline constexpr Vector<T, C> Matrix<T, C, R>::Row() const noexcept
	{
		Vector<T, C> row;
		for (uint32_t c = 0; c < C; c++)
			row[c] = m_pElements[c][R2];
		return row;
	}

	template<typename T, uint32_t C, uint32_t R>
	template<uint32_t C2, uint32_t R2, typename T2, typename... Args>
	inline constexpr void Matrix<T, C, R>::Fill(const T2& crParam, Args&&... rrArgs)
	{
		m_pElements[C - C2][R - R2] = T(crParam);
		if constexpr (C2 > 1)
			Fill<C2 - 1, R2>(std::forward<Args>(rrArgs)...);
		else if constexpr (R2 > 1)
			Fill<C, R2 - 1>(std::forward<Args>(rrArgs)...);
	}

	template<typename T, uint32_t C, uint32_t R, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R> operator+(const T2& crValue, const Matrix<T, C, R>& crMatrix) noexcept
	{
		return crMatrix + crValue;
	}

	template<typename T, uint32_t C, uint32_t R, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R> operator-(const T2& crValue, const Matrix<T, C, R>& crMatrix) noexcept
	{
		return -crMatrix + crValue;
	}

	template<typename T, uint32_t C, uint32_t R, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, C, R> operator*(const T2& crValue, const Matrix<T, C, R>& crMatrix) noexcept
	{
		return crMatrix * crValue;
	}

	template<typename T, uint32_t CR, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, CR> operator/(const T2& crValue, const Matrix<T, CR>& crMatrix) noexcept
	{
		return Inverse(crMatrix) * crValue;
	}

	template<typename T, uint32_t C, uint32_t R, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, C> operator*(const Vector<T2, R>& crVector, const Matrix<T, C, R>& crMatrix) noexcept
	{
		Vector<T, C> result;

		Matrix<T, R, C> matrix = Transpose(crMatrix);
		for (uint32_t r = 0; r < R; r++)
			result += matrix[r] * crVector[r];

		return result;
	}

	template<uint32_t C, uint32_t R, typename T, uint32_t C2, uint32_t R2, std::enable_if_t<(C < C2) && (R < R2) && (C2 > 1) && (R2 > 1), int>>
	inline constexpr Matrix<T, C2 - 1, R2 - 1> Submatrix(const Matrix<T, C2, R2>& crMatrix) noexcept
	{
		Matrix<T, C2 - 1, R2 - 1> result = T();

		if constexpr (C > 0 && R > 0)
			for (uint32_t c = 0; c < C; c++)
				for (uint32_t r = 0; r < R; r++)
					result[c][r] = crMatrix[c][r];

		if constexpr (C + 1 < C2 && R > 0)
			for (uint32_t c = C + 1; c < C2; c++)
				for (uint32_t r = 0; r < R; r++)
					result[c - 1][r] = crMatrix[c][r];

		if constexpr (C > 0 && R + 1 < R2)
			for (uint32_t c = 0; c < C; c++)
				for (uint32_t r = R + 1; r < R2; r++)
					result[c][r - 1] = crMatrix[c][r];

		if constexpr (C + 1 < C2 && R + 1 < R2)
			for (uint32_t c = C + 1; c < C2; c++)
				for (uint32_t r = R + 1; r < R2; r++)
					result[c - 1][r - 1] = crMatrix[c][r];

		return result;
	}

	template<typename T, uint32_t C2, uint32_t R2, std::enable_if_t<(C2 > 1) && (R2 > 1), int>>
	inline constexpr Matrix<T, C2 - 1, R2 - 1> Submatrix(uint32_t C, uint32_t R, const Matrix<T, C2, R2>& crMatrix) noexcept
	{
		//assert((C < C2) && (R < R2));

		Matrix<T, C2 - 1, R2 - 1> result = T();

		if (C > 0 && R > 0)
			for (uint32_t c = 0; c < C; c++)
				for (uint32_t r = 0; r < R; r++)
					result[c][r] = crMatrix[c][r];

		if (C + 1 < C2 && R > 0)
			for (uint32_t c = C + 1; c < C2; c++)
				for (uint32_t r = 0; r < R; r++)
					result[c - 1][r] = crMatrix[c][r];

		if (C > 0 && R + 1 < R2)
			for (uint32_t c = 0; c < C; c++)
				for (uint32_t r = R + 1; r < R2; r++)
					result[c][r - 1] = crMatrix[c][r];

		if (C + 1 < C2 && R + 1 < R2)
			for (uint32_t c = C + 1; c < C2; c++)
				for (uint32_t r = R + 1; r < R2; r++)
					result[c - 1][r - 1] = crMatrix[c][r];

		return result;
	}

	template<typename T, uint32_t CR>
	inline constexpr T Trace(const Matrix<T, CR>& crMatrix) noexcept
	{
		T result = T(0);

		for (uint32_t cr = 0; cr < CR; cr++)
			result += crMatrix[cr][cr];

		return result;
	}

	template<typename T, uint32_t CR>
	inline constexpr T Determinant(const Matrix<T, CR>& crMatrix) noexcept
	{
		T result = T(0);

		for (uint32_t c = 0; c < CR; c++)
		{
			T subDeterminant = crMatrix[c][0] * Determinant(Submatrix(c, 0, crMatrix));
			if (c % 2 == 0)
				result += subDeterminant;
			else
				result -= subDeterminant;
		};

		return result;
	}

	template<typename T>
	inline constexpr T Determinant(const Matrix<T, 1>& crMatrix) noexcept
	{
		return crMatrix[0][0];
	}

	template<typename T, uint32_t C, uint32_t R>
	inline constexpr Matrix<T, R, C> Transpose(const Matrix<T, C, R>& crMatrix) noexcept
	{
		Matrix<T, R, C> result = T(0);

		for (uint32_t c = 0; c < C; c++)
			for (uint32_t r = 0; r < R; r++)
				result[r][c] = crMatrix[c][r];

		return result;
	}

	template<typename T, uint32_t CR>
	inline constexpr Matrix<T, CR> Minors(const Matrix<T, CR>& crMatrix) noexcept
	{
		Matrix<T, CR> result = T(0);

		for (uint32_t c = 0; c < CR; c++)
			for (uint32_t r = 0; r < CR; r++)
				result[c][r] = Determinant(Submatrix(c, r, crMatrix));

		return result;
	}

	template<typename T, uint32_t CR>
	inline constexpr Matrix<T, CR> Cofactors(const Matrix<T, CR>& crMatrix) noexcept
	{
		Matrix<T, CR> result = Minors(crMatrix);

		for (uint32_t c = 0; c < CR; c++)
			for (uint32_t r = 0; r < CR; r++)
				if ((c + r) % 2 != 0)
					result[c][r] = -result[c][r];

		return result;
	}

	template<typename T, uint32_t CR>
	inline constexpr Matrix<T, CR> Adjugate(const Matrix<T, CR>& crMatrix) noexcept
	{
		return Transpose(Cofactors(crMatrix));
	}

	template<typename T, uint32_t CR>
	inline constexpr Matrix<T, CR> Inverse(const Matrix<T, CR>& crMatrix) noexcept
	{
		return Adjugate(crMatrix) / Determinant(crMatrix);
	}

	template<typename T, uint32_t CR, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, CR> Inverse(const Matrix<T, CR>& crMatrix, const T2& crDeterminant) noexcept
	{
		return Adjugate(crMatrix) / T(crDeterminant);
	}

	template<typename T, uint32_t CR, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, CR> Translate(const Matrix<T, CR>& crMatrix, const Vector<T2, CR - 1>& crTranslation)
	{
		Matrix<T, CR> translation;

		for (uint32_t r = 0; r < CR - 1; r++)
			translation[CR - 1][r] += T(crTranslation[r]);

		return crMatrix * translation;
	}

	template<uint32_t A1, uint32_t A2, typename T, uint32_t CR, typename T2, std::enable_if_t<std::is_convertible_v<T2, T> && (A1 < CR) && (A2 < CR) && A1 != A2, int>>
	inline constexpr Matrix<T, CR> Rotate(const Matrix<T, CR, CR>& crMatrix, const T2& crRadians)
	{
		Matrix<T, CR> rotation;

		T sin = T(gcem::sin(T(crRadians)));
		T cos = T(gcem::cos(T(crRadians)));

		rotation[A1][A1] = cos; rotation[A2][A1] = -sin;
		rotation[A1][A2] = sin; rotation[A2][A2] = cos;

		return crMatrix * rotation;
	}

	template<typename T, uint32_t CR, typename T2, std::enable_if_t<std::is_convertible_v<T2, T> && (CR >= 2), int>>
	inline constexpr Matrix<T, CR> Rotate(uint32_t A1, uint32_t A2, const Matrix<T, CR>& crMatrix, const T2& crRadians)
	{
		//assert((A1 < CR) && (A2 < CR) && A1 != A2);

		Matrix<T, CR> rotation;

		T sin = T(gcem::sin(T(crRadians)));
		T cos = T(gcem::cos(T(crRadians)));

		rotation[A1][A1] = cos; rotation[A2][A1] = -sin;
		rotation[A1][A2] = sin; rotation[A2][A2] = cos;

		return crMatrix * rotation;
	}

	template<typename T, uint32_t CR, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, CR> Scale(const Matrix<T, CR>& crMatrix, const Vector<T2, CR - 1>& crScale)
	{
		Matrix<T, CR> scale;

		for (uint32_t cr = 0; cr < CR - 1; cr++)
			scale[cr][cr] = T(crScale[cr]);

		return crMatrix * scale;
	}

	template<typename T, uint32_t CR, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Matrix<T, CR> Scale(const Matrix<T, CR>& crMatrix, const T2& crScale)
	{
		Matrix<T, CR> scale = T(crScale);
		scale[CR - 1][CR - 1] = T(1);
		return crMatrix * scale;
	}

	template<typename T, uint32_t C, uint32_t R>
	std::ostream& operator<<(std::ostream& rOstream, const Matrix<T, C, R>& crMatrix)
	{
		rOstream << std::scientific << '[' << std::setw(13) << crMatrix[0][0];
		for (uint32_t c = 1; c < C; c++)
			rOstream << std::setw(14) << crMatrix[c][0];
		rOstream << ']';

		for (uint32_t r = 1; r < R; r++)
		{
			rOstream << "\n[" << std::setw(13) << crMatrix[0][r];
			for (uint32_t c = 1; c < C; c++)
				rOstream << std::setw(14) << crMatrix[c][r];
			rOstream << ']';
		}

		return rOstream << std::defaultfloat;
	}
}
