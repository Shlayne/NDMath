#pragma once

#include "ND/Matrix.h"
#include "ND/Mesh.h"
#include <array>
#include <set>

namespace nd
{
	using SlicerErrorCode = uint8_t;
	enum : SlicerErrorCode
	{
		SlicerErrorCode_None = 0,
		SlicerErrorCode_InvalidSlice = 1,
		SlicerErrorCode_EmptyInputMesh = 2,
		SlicerErrorCode_InvalidIndexCount = 3,
	};

	// Calculates intersections of meshes with hyperplanes.
	template<typename T, uint32_t N>
	class Slicer
	{
	public:
		static_assert(N > 1, "Slicer must have enough output dimensions.");
	public:
		using TransformIn = Matrix<T, N + 1>;
		using MeshIn = Mesh<T, N>;
		using MeshOut = Mesh<T, N - 1>;
	public:
		// Returns true if slicing resulted in a mesh, false if not or an error occurred.
		// You only need to check for errors if this returns false, which even then,
		// there may be no error. In this case, the input mesh did not intersect the
		// hyperplane described by the input transformation.
		// NOTE: Requires input mesh to contain only N - 1 simplices (N indices per).
		bool Slice(const MeshIn& crMeshIn, const T& crSlice = T(0)) noexcept;

		// Gets the error code of the most recent slice performed.
		SlicerErrorCode GetErrorCode() const noexcept;

		// Gets the error message of the most recent slice performed.
		const char* GetErrorMessage() const noexcept;

		// Gets the output mesh of the most recent slice performed.
		MeshOut GetOutputMesh() const noexcept;

		// Releases memory allocated by the most recent slice performed.
		// This is also done automatically at the beginning of Slice, so
		// there's no need to call this in between successive Slice calls.
		// But, after you're done using this Slicer, calling this is recommended.
		void Reset() noexcept;
	private:
		void CalculateIntersections(std::set<Vector<T, N - 1>>& rOutVertices);

		template<uint32_t N2>
		void CalculateIntersectionsImpl(std::set<Vector<T, N2 - 1>>& rOutVertices, const std::array<Index, N2>& crSimplexIndices);

		template<uint32_t N2>
		static constexpr std::array<Index, N2 * (N2 + 1)> GenerateSimplexIndices();

		template<uint32_t N2, uint32_t I = N2, uint32_t J = 0>
		static constexpr auto GenerateSimplexIndicesImpl();
	private:
		bool Fail(SlicerErrorCode errorCode) noexcept;
		bool Succeed() noexcept;
	private:
		const MeshIn* m_cpMeshIn = nullptr;
		const T* m_cpSlice = nullptr;
		MeshOut m_MeshOut;
		SlicerErrorCode m_ErrorCode = SlicerErrorCode_None;
	};
}

#include "Slicer.inl"
