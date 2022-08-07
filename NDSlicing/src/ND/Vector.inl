#include "Vector.h"
namespace nd
{
	template<typename T, uint32_t N>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N>::Vector(const T2& crValue) noexcept
	{
		for (T& rDim : m_pElements)
			rDim = T(crValue);
	}

	template<typename T, uint32_t N>
	template<typename... Args, std::ext::disable_if_t<std::ext::are_convertible_v<T, Args...> && sizeof...(Args) == 1, int>>
	inline constexpr Vector<T, N>::Vector(Args&&... rrArgs) noexcept
	{
		Fill<N>(std::forward<Args>(rrArgs)...);
	}

	template<typename T, uint32_t N>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N>& Vector<T, N>::operator=(const T2& crValue) noexcept
	{
		for (T& rDim : m_pElements)
			rDim = T(crValue);
		return *this;
	}

	template<typename T, uint32_t N>
	template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T> && N2 <= N, int>>
	inline constexpr Vector<T, N>& Vector<T, N>::operator=(const Vector<T2, N2>& crVector) noexcept
	{
		if (this != (void*)&crVector)
		{
			uint32_t n = 0;
			for (; n < N2; n++)
				m_pElements[n] = T(crVector[n]);
			for (; n < N; n++)
				m_pElements[n] = T(0);
		}
		return *this;
	}

	template<typename T, uint32_t N>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N> Vector<T, N>::operator+(const T2& crValue) const noexcept
	{
		return +*this += crValue;
	}

	template<typename T, uint32_t N>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N>& Vector<T, N>::operator+=(const T2& crValue) noexcept
	{
		for (T& rDim : m_pElements)
			rDim += T(crValue);
		return *this;
	}

	template<typename T, uint32_t N>
	template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T> && N2 <= N, int>>
	inline constexpr Vector<T, N> Vector<T, N>::operator+(const Vector<T2, N2>& crVector) const noexcept
	{
		return +*this += crVector;
	}

	template<typename T, uint32_t N>
	template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T> && N2 <= N, int>>
	inline constexpr Vector<T, N>& Vector<T, N>::operator+=(const Vector<T2, N2>& crVector) noexcept
	{
		for (uint32_t n = 0; n < N2; n++)
			m_pElements[n] += T(crVector[n]);
		return *this;
	}

	template<typename T, uint32_t N>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N> Vector<T, N>::operator-(const T2& crValue) const noexcept
	{
		return +*this -= crValue;
	}

	template<typename T, uint32_t N>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N>& Vector<T, N>::operator-=(const T2& crValue) noexcept
	{
		for (T& rDim : m_pElements)
			rDim -= T(crValue);
		return *this;
	}

	template<typename T, uint32_t N>
	template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T> && N2 <= N, int>>
	inline constexpr Vector<T, N> Vector<T, N>::operator-(const Vector<T2, N2>& crVector) const noexcept
	{
		return +*this -= crVector;
	}

	template<typename T, uint32_t N>
	template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T> && N2 <= N, int>>
	inline constexpr Vector<T, N>& Vector<T, N>::operator-=(const Vector<T2, N2>& crVector) noexcept
	{
		for (uint32_t n = 0; n < N2; n++)
			m_pElements[n] -= T(crVector[n]);
		return *this;
	}

	template<typename T, uint32_t N>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N> Vector<T, N>::operator*(const T2& crValue) const noexcept
	{
		return +*this *= crValue;
	}

	template<typename T, uint32_t N>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N>& Vector<T, N>::operator*=(const T2& crValue) noexcept
	{
		for (T& rDim : m_pElements)
			rDim *= T(crValue);
		return *this;
	}

	template<typename T, uint32_t N>
	template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T> && N2 <= N, int>>
	inline constexpr Vector<T, N> Vector<T, N>::operator*(const Vector<T2, N2>& crVector) const noexcept
	{
		return +*this *= crVector;
	}

	template<typename T, uint32_t N>
	template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T> && N2 <= N, int>>
	inline constexpr Vector<T, N>& Vector<T, N>::operator*=(const Vector<T2, N2>& crVector) noexcept
	{
		for (uint32_t n = 0; n < N2; n++)
			m_pElements[n] *= T(crVector[n]);
		return *this;
	}

	template<typename T, uint32_t N>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N> Vector<T, N>::operator/(const T2& crValue) const noexcept
	{
		return +*this /= crValue;
	}

	template<typename T, uint32_t N>
	template<typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N>& Vector<T, N>::operator/=(const T2& crValue) noexcept
	{
		for (T& rDim : m_pElements)
			rDim /= T(crValue);
		return *this;
	}

	template<typename T, uint32_t N>
	template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T>&& N2 <= N, int>>
	inline constexpr Vector<T, N> Vector<T, N>::operator/(const Vector<T2, N2>& crVector) const noexcept
	{
		return +*this /= crVector;
	}

	template<typename T, uint32_t N>
	template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T2, T>&& N2 <= N, int>>
	inline constexpr Vector<T, N>& Vector<T, N>::operator/=(const Vector<T2, N2>& crVector) noexcept
	{
		for (uint32_t n = 0; n < N2; n++)
			m_pElements[n] /= T(crVector[n]);
		return *this;
	}

	template<typename T, uint32_t N>
	inline constexpr Vector<T, N> Vector<T, N>::operator+() const noexcept
	{
		return Vector<T, N>(*this);
	}

	template<typename T, uint32_t N>
	inline constexpr Vector<T, N> Vector<T, N>::operator-() const noexcept
	{
		Vector<T, N> vector = +*this;
		for (T& rDim : m_pElements)
			rDim = -rDim;
		return vector;
	}

	template<typename T, uint32_t N>
	inline constexpr T& Vector<T, N>::operator[](uint32_t index) noexcept
	{
		return m_pElements[index];
	}

	template<typename T, uint32_t N>
	inline constexpr const T& Vector<T, N>::operator[](uint32_t index) const noexcept
	{
		return m_pElements[index];
	}

	template<typename T, uint32_t N>
	template<typename T2, uint32_t N2, std::enable_if_t<std::is_convertible_v<T, T2> && (N2 <= N), int>>
	inline constexpr Vector<T, N>::operator Vector<T2, N2>() const noexcept
	{
		Vector<T2, N2> result;

		for (uint32_t n = 0; n < N2; n++)
			result[n] = T2(m_pElements[n]);

		return result;
	}

	template<typename T, uint32_t N>
	template<uint32_t N2, typename T2, typename... Args, std::enable_if_t<std::is_convertible_v<T2, T> && (N2 >= 1), int>>
	inline constexpr void Vector<T, N>::Fill(const T2& crValue, Args&&... rrArgs) noexcept
	{
		m_pElements[N - N2] = T(crValue);
		if constexpr (sizeof...(Args) > 0)
			Fill<N2 - 1>(std::forward<Args>(rrArgs)...);
	}

	template<typename T, uint32_t N>
	template<uint32_t N2, typename T2, uint32_t N3, typename... Args, std::enable_if_t<std::is_convertible_v<T2, T> && (N2 >= N3), int>>
	inline constexpr void Vector<T, N>::Fill(const Vector<T2, N3>& crVector, Args&&... rrArgs) noexcept
	{
		for (uint32_t n = 0; n < N3; n++)
			m_pElements[N - N2 + n] = T(crVector[n]);
		if constexpr (sizeof...(Args) > 0)
			Fill<N2 - N3>(std::forward<Args>(rrArgs)...);
	}

	template<typename T, uint32_t N, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N> operator+(const T2& crValue, const Vector<T, N>& crVector) noexcept
	{
		return Vector<T, N>(crValue) += crVector;
	}

	template<typename T, uint32_t N, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N> operator-(const T2& crValue, const Vector<T, N>& crVector) noexcept
	{
		return Vector<T, N>(crValue) -= crVector;
	}

	template<typename T, uint32_t N, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N> operator*(const T2& crValue, const Vector<T, N>& crVector) noexcept
	{
		return Vector<T, N>(crValue) *= crVector;
	}

	template<typename T, uint32_t N, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N> operator/(const T2& crValue, const Vector<T, N>& crVector) noexcept
	{
		return Vector<T, N>(crValue) /= crVector;
	}

	template<typename T, uint32_t N, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N> Dot(const Vector<T, N>& crVector1, const Vector<T2, N>& crVector2) noexcept
	{
		T dot = T(0);
		Vector<T, N> product = crVector1 * crVector2;
		for (uint32_t n = 0; n < N; n++)
			dot += product[n];
		return dot;
	}
	
	template<typename T, uint32_t N>
	inline constexpr T Length2(const Vector<T, N>& crVector) noexcept
	{
		return Dot(crVector, crVector);
	}

	template<typename T, uint32_t N>
	inline constexpr T Length(const Vector<T, N>& crVector) noexcept
	{
		return T(gcem::sqrt(Length2(crVector)));
	}

	template<typename T, uint32_t N>
	inline constexpr Vector<T, N> Normalize(const Vector<T, N>& crVector) noexcept
	{
		return crVector / Length(crVector);
	}

	namespace impl
	{
		template<typename T, uint32_t N>
		inline constexpr Vector<T, N> Apply(T(*pFunc)(T), const Vector<T, N>& crVector)
		{
			Vector<T, N> result;

			for (uint32_t n = 0; n < N; n++)
				result[n] = pFunc(crVector[n]);

			return result;
		}

		template<typename T, uint32_t N>
		inline constexpr Vector<T, N> Apply(T(*pFunc)(T, T), const Vector<T, N>& crVector1, const Vector<T, N>& crVector2)
		{
			Vector<T, N> result;

			for (uint32_t n = 0; n < N; n++)
				result[n] = pFunc(crVector1[n], crVector2[n]);

			return result;
		}

		template<typename T, uint32_t N>
		inline constexpr Vector<T, N> Apply(T(*pFunc)(const T&, const T&, const T&), const Vector<T, N>& crVector1, const Vector<T, N>& crVector2, const Vector<T, N>& crVector3)
		{
			Vector<T, N> result;

			for (uint32_t n = 0; n < N; n++)
				result[n] = pFunc(crVector1[n], crVector2[n], crVector3[n]);

			return result;
		}

		template<typename T, uint32_t N>
		inline constexpr Vector<T, N> Apply(T(*pFunc)(const T&, const T&, const T&), const Vector<T, N>& crVector1, const Vector<T, N>& crVector2, const T& crValue)
		{
			Vector<T, N> result;

			for (uint32_t n = 0; n < N; n++)
				result[n] = pFunc(crVector1[n], crVector2[n], crValue);

			return result;
		}
	}

	template<typename T, uint32_t N, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N> Min(const Vector<T, N>& crVector1, const Vector<T2, N>& crVector2) noexcept
	{
		return impl::Apply(gcem::min, crVector1, Vector<T, N>(crVector2));
	}

	template<typename T, uint32_t N, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N> Max(const Vector<T, N>& crVector1, const Vector<T2, N>& crVector2) noexcept
	{
		return impl::Apply(gcem::max, crVector1, Vector<T, N>(crVector2));
	}

	template<typename T, uint32_t N>
	inline constexpr Vector<T, N> Abs(const Vector<T, N>& crVector1) noexcept
	{
		return impl::Apply(gcem::abs, crVector1);
	}

	template<typename T, uint32_t N>
	inline constexpr Vector<T, N> Floor(const Vector<T, N>& crVector1) noexcept
	{
		return impl::Apply(gcem::floor, crVector1);
	}

	template<typename T, uint32_t N>
	inline constexpr Vector<T, N> Ceil(const Vector<T, N>& crVector1) noexcept
	{
		return impl::Apply(gcem::ceil, crVector1);
	}

	template<typename T, uint32_t N, typename T2, std::enable_if_t<std::is_convertible_v<T2, T>, int>>
	inline constexpr Vector<T, N> Mod(const Vector<T, N>& crVector1, const Vector<T2, N>& crVector2) noexcept
	{
		return impl::Apply(gcem::fmod, crVector1, Vector<T, N>(crVector2));
	}

	template<typename T, uint32_t N, typename T2, typename T3, std::enable_if_t<std::ext::are_convertible_v<T, T2, T3>, int>>
	constexpr Vector<T, N> Clamp(const Vector<T, N>& crVector, const Vector<T2, N>& crMin, const Vector<T3, N>& crMax) noexcept
	{
		return impl::Apply(std::clamp, crVector, Vector<T, N>(crMin), Vector<T, N>(crMax));
	}

	template<typename T, uint32_t N>
	inline constexpr Vector<T, N> Trunc(const Vector<T, N>& crVector1) noexcept
	{
		return impl::Apply(gcem::trunc, crVector1);
	}

	template<typename T, uint32_t N>
	inline constexpr Vector<T, N> Round(const Vector<T, N>& crVector1) noexcept
	{
		return impl::Apply(gcem::round, crVector1);
	}

	namespace impl
	{
		template<typename T>
		inline constexpr T Lerp(const T& crA, const T& crB, const T& crT)
		{
			return (T(1) - crT) * crA + crT * crB;
		};
	}

	template<typename T, uint32_t N, typename T2, typename T3, std::enable_if_t<std::ext::are_convertible_v<T, T2, T3>, int>>
	inline constexpr Vector<T, N> Lerp(const Vector<T, N>& crVector1, const Vector<T2, N>& crVector2, const T3& crTValue) noexcept
	{
		return impl::Apply(impl::Lerp, crVector1, crVector2, crTValue);
	}

	template<typename T, uint32_t N>
	std::ostream& operator<<(std::ostream& rOstream, const Vector<T, N>& crVector)
	{
		rOstream << '<' << crVector[0];
		if constexpr (N > 1)
			for (uint32_t n = 1; n < N; n++)
				rOstream << ',' << crVector[n];
		return rOstream << '>';
	}
}
