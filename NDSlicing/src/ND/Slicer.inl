namespace nd
{
	template<typename T, uint32_t N>
	inline bool Slicer<T, N>::Slice(const MeshIn& crMeshIn, const T& crSlice) noexcept
	{
		Reset();

		if (!gcem::internal::is_finite(crSlice))
			return Fail(SlicerErrorCode_InvalidSlice);

		if (crMeshIn.vertices.empty() || crMeshIn.indices.empty())
			return Fail(SlicerErrorCode_EmptyInputMesh);

		if ((crMeshIn.indices.size() % N) != 0)
			return Fail(SlicerErrorCode_InvalidIndexCount);

		m_cpMeshIn = &crMeshIn;
		m_cpSlice = &crSlice;

		std::vector<Vertex<T, N - 1>> vertices;
		CalculateIntersections<N - 1>(vertices);

		return Succeed();
	}

	template<typename T, uint32_t N>
	template<uint32_t N2>
	inline void Slicer<T, N>::CalculateIntersections(std::vector<Vertex<T, N2>>& rOutVertices)
	{
		constexpr std::array<Index, N2 * (N2 + 1)> simplexIndices = GenerateSimplexIndices<N2>();

		const std::vector<Index>& crIndicesIn = m_cpMeshIn->indices;
		const std::vector<Vertex<T, N>>& crVerticesIn = m_cpMeshIn->vertices;

		for (size_t i = 0; i < crIndicesIn.size(); i += N2 + 1)
		{
			if constexpr (N2 == 1)
			{
				const Vector<T, N>& crVertexA = crVerticesIn[crIndicesIn[i + 0]].position;
				const Vector<T, N>& crVertexB = crVerticesIn[crIndicesIn[i + 1]].position;

				// If intersection...
				if (gcem::signbit(crVertexA[N2]) != gcem::signbit(crVertexB[N2]) || crVertexA[N2] == crVertexB[N2])
				{
					T t = crVertexA[N2] / (crVertexA[N2] * crVertexB[N2]);

					(Vector<T, N2>)crVertexA;
					//rOutVertices.emplace_back(nd::Lerp((Vector<T, N2>)crVertexA, (Vector<T, N2>)crVertexB, t));
				}
			}
			else
			{

			}
		}

		//std::vector<Vertex<T, N2 - 1>> vertices;
		//if constexpr (N2 > 1)
		//	CalculateIntersections<N2 - 1>(vertices);

		//	Get cross sections of (n - 1)-simplex
		//		...
		//			Get cross sections of 5-simplex
		//				Get cross sections of pentachora (4-simplex)
		//					Get cross sections of tetrahedra (3-simplex)
		//						Get cross sections of triangles (2-simplex)
		//							Get cross sections of lines (1-simplex)
	}

	template<typename T, uint32_t N>
	template<uint32_t N2, uint32_t I, uint32_t J>
	inline constexpr auto Slicer<T, N>::GenerateSimplexIndicesImpl()
	{
		if constexpr (J < N2)
			return std::ext::combine_integer_sequences(Indices<J + (J >= I)>(), GenerateSimplexIndicesImpl<N2, I, J + 1>());
		else if constexpr (I > 0)
			return GenerateSimplexIndicesImpl<N2, I - 1>();
		else
			return Indices<>();
	}

	template<typename T, uint32_t N>
	template<uint32_t N2>
	inline constexpr std::array<Index, N2 * (N2 + 1)> Slicer<T, N>::GenerateSimplexIndices()
	{
		return std::ext::integer_sequence_to_array(GenerateSimplexIndicesImpl<N2>());

		// Emulates the following, but at compile time:
		//std::array<Index, N2 * (N2 + 1)> indices;
		//Index index = 0;
		//for (uint32_t i = N2 + 1; i-- > 0; )
		//	for (uint32_t j = 0; j < N2; j++)
		//		indices[index++] = j + (j >= i ? 1 : 0);
		//return indices;
	}

	template<typename T, uint32_t N>
	inline SlicerErrorCode Slicer<T, N>::GetErrorCode() const noexcept
	{
		return m_ErrorCode;
	}

	template<typename T, uint32_t N>
	inline const char* Slicer<T, N>::GetErrorMessage() const noexcept
	{
		switch (m_ErrorCode)
		{
			case SlicerErrorCode_InvalidSlice:
				return "The provided hyperplane slice was nan or infinite.";
			case SlicerErrorCode_EmptyInputMesh:
				return "The provided input mesh had no vertices or had no indices.";
			case SlicerErrorCode_InvalidIndexCount:
				return "Too many or not enough indices in mesh.";
		}

		return nullptr;
	}

	template<typename T, uint32_t N>
	inline Slicer<T, N>::MeshOut Slicer<T, N>::GetOutputMesh() const noexcept
	{
		return m_MeshOut;
	}

	template<typename T, uint32_t N>
	inline void Slicer<T, N>::Reset() noexcept
	{
		m_cpMeshIn = nullptr;
		m_cpSlice = nullptr;
		m_MeshOut.vertices.clear();
		m_MeshOut.indices.clear();
		m_ErrorCode = SlicerErrorCode_None;
	}

	template<typename T, uint32_t N>
	inline bool Slicer<T, N>::Fail(SlicerErrorCode errorCode) noexcept
	{
		m_ErrorCode = errorCode;
		return false;
	}

	template<typename T, uint32_t N>
	inline bool Slicer<T, N>::Succeed() noexcept
	{
		m_ErrorCode = SlicerErrorCode_None;
		return true;
	}
}
