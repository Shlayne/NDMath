namespace nd
{
	template<typename T, uint32_t N>
	inline constexpr Vertex<T, N>::Vertex(const Position& crPosition) noexcept
		: position(crPosition)
	{

	}
}
