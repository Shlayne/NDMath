#include "Test.h"

import nd;
import <iostream>;
// These can't be imported 'cause their header units aren't found by visual studio for some reason.
//#include <array>
//#include <set>
//#include <unordered_set>
//#include <functional>

using namespace nd;
// TODO: Loader.inl completion.

//template<size_t N, Scalar S>
//using LineSegment = std::pair<Vector<N, S>, Vector<N, S>>;
//
//template<size_t N, Scalar S>
//struct LineSegmentEqual
//{
//	constexpr auto operator()(const LineSegment<N, S>& lineSegment1, const LineSegment<N, S>& lineSegment2) const noexcept -> bool
//	{
//		std::cout << "pvv\n";
//		auto hash11 = Vector<N, S>::Hash()(lineSegment1.first);
//		auto hash12 = Vector<N, S>::Hash()(lineSegment1.second);
//		auto hash21 = Vector<N, S>::Hash()(lineSegment2.first);
//		auto hash22 = Vector<N, S>::Hash()(lineSegment2.second);
//		return (hash11 == hash21 && hash12 == hash22) ||
//			   (hash11 == hash22 && hash12 == hash21);
//	}
//};
//
//template< size_t N, size_t N2, Scalar S>
//struct LineSegmentArrayEqual
//{
//	constexpr auto operator()(const std::array<LineSegment<N, S>, N2>& lineSegmentArray1, const std::array<LineSegment<N, S>, N2>& lineSegmentArray2) const noexcept -> bool
//	{
//		std::cout << "apvv\n";
//		// TODO: too sensitive.
//		LineSegmentEqual<N, S> vpe;
//		return (vpe(lineSegmentArray1[0], lineSegmentArray2[0]) && vpe(lineSegmentArray1[1], lineSegmentArray2[1]) && vpe(lineSegmentArray1[2], lineSegmentArray2[2])) ||
//			   (vpe(lineSegmentArray1[0], lineSegmentArray2[1]) && vpe(lineSegmentArray1[1], lineSegmentArray2[0]) && vpe(lineSegmentArray1[2], lineSegmentArray2[2])) ||
//			   (vpe(lineSegmentArray1[0], lineSegmentArray2[0]) && vpe(lineSegmentArray1[1], lineSegmentArray2[2]) && vpe(lineSegmentArray1[2], lineSegmentArray2[1])) ||
//			   (vpe(lineSegmentArray1[0], lineSegmentArray2[2]) && vpe(lineSegmentArray1[1], lineSegmentArray2[0]) && vpe(lineSegmentArray1[2], lineSegmentArray2[1])) ||
//			   (vpe(lineSegmentArray1[0], lineSegmentArray2[1]) && vpe(lineSegmentArray1[1], lineSegmentArray2[2]) && vpe(lineSegmentArray1[2], lineSegmentArray2[0])) ||
//			   (vpe(lineSegmentArray1[0], lineSegmentArray2[2]) && vpe(lineSegmentArray1[1], lineSegmentArray2[1]) && vpe(lineSegmentArray1[2], lineSegmentArray2[0]));
//	}
//};
//
//template<size_t N, Scalar S>
//struct LineSegmentHash
//{
//	constexpr size_t operator()(const LineSegment<N, S>& lineSegment) const noexcept
//	{
//		return Vector<N, S>::Hash()(lineSegment.first) ^ Vector<N, S>::Hash()(lineSegment.second);
//	}
//};
//
//template<size_t N, size_t N2, Scalar S>
//struct LineSegmentArrayHash
//{
//	constexpr size_t operator()(const std::array<LineSegment<N, S>, N2>& lineSegmentArray) const noexcept
//	{
//		size_t hash = 0;
//		for (uint32_t n = 0; n < N2; n++)
//			hash ^= Vector<N, S>::Hash()(lineSegmentArray[n].first) | Vector<N, S>::Hash()(lineSegmentArray[n].second);
//		return hash;
//	}
//};
//
//using LineSegment3f = LineSegment<3, float>;
//using LineSegmentEqual3f = LineSegmentEqual<3, float>;
//using LineSegmentArrayEqual3f = LineSegmentArrayEqual<3, 3, float>;
//using LineSegmentHash3f = LineSegmentHash<3, float>;
//using LineSegmentArrayHash3f = LineSegmentArrayHash<3, 3, float>;

namespace test
{
	void Perform()
	{
		//constexpr std::array<LineSegment3f, 3> array1{ LineSegment3f{ Vector3f(2, 0, 0), Vector3f(0, 2, 0) }, LineSegment3f{ Vector3f(0, 2, 0), Vector3f(0, 0, 2) }, LineSegment3f{ Vector3f(0, 0, 2), Vector3f(2, 0, 0) } };
		//auto hash1 = LineSegmentArrayHash3f()(array1);
		//constexpr std::array<LineSegment3f, 3> array2{ LineSegment3f{ Vector3f(0, 2, 0), Vector3f(0, 0, 2) }, LineSegment3f{ Vector3f(1, 0, 0), Vector3f(0, 2, 0) }, LineSegment3f{ Vector3f(0, 0, 2), Vector3f(1, 0, 0) } };
		//auto hash2 = LineSegmentArrayHash3f()(array2);
		//auto hash3 = Vector3f::Hash<>()(Vector3f(2, 0, 0));
		//auto hash4 = Vector3f::Hash<>()(Vector3f(0, 2, 0));

		//constexpr Matrix5f transform = Rotate<0, 3>(Matrix5f(), 20.0f * 3.1415926535f / 180.0f);

		//Mesh<float, 4> inputMesh;
		//inputMesh.vertices.emplace_back(transform * Vector5f(1.0f, 1.0f / gcem::sqrt(3.0f), 1.0f / gcem::sqrt(6.0f), 1.0f / gcem::sqrt(10.0f), 1.0f));
		//inputMesh.vertices.emplace_back(transform * Vector5f(-1.0f, 1.0f / gcem::sqrt(3.0f), 1.0f / gcem::sqrt(6.0f), 1.0f / gcem::sqrt(10.0f), 1.0f));
		//inputMesh.vertices.emplace_back(transform * Vector5f(0.0f, -2.0f / gcem::sqrt(3.0f), 1.0f / gcem::sqrt(6.0f), 1.0f / gcem::sqrt(10.0f), 1.0f));
		//inputMesh.vertices.emplace_back(transform * Vector5f(0.0f, 0.0f, -gcem::sqrt(3.0f / 2.0f), 1.0f / gcem::sqrt(10.0f), 1.0f));
		//inputMesh.vertices.emplace_back(transform * Vector5f(0.0f, 0.0f, 0.0f, -2.0f * gcem::sqrt(2.0f / 5.0f), 1.0f));
		//inputMesh.indices = { 0, 1, 2, 3,  0, 1, 2, 4,  0, 1, 3, 4,  0, 2, 3, 4,  1, 2, 3, 4 };

		//std::unordered_set<std::array<LineSegment3f, 3>, LineSegmentArrayHash3f, LineSegmentArrayEqual3f> intersectionTriangles;

		//// For each tetrahedra...
		//for (size_t i = 0; i < inputMesh.indices.size(); i += 4)
		//{
		//	std::unordered_set<LineSegment3f, LineSegmentHash3f, LineSegmentEqual3f> intersectionLines;

		//	constexpr uint8_t tetrahedronPointCount = 4;
		//	Index tetrahedronIndices[tetrahedronPointCount]
		//	{
		//		inputMesh.indices[i + 0],
		//		inputMesh.indices[i + 1],
		//		inputMesh.indices[i + 2],
		//		inputMesh.indices[i + 3],
		//	};

		//	// For each triangle in the tetrahedron...
		//	for (uint8_t j = 0; j < tetrahedronPointCount; j++)
		//	{
		//		std::unordered_set<LineSegment3f, LineSegmentHash3f, LineSegmentEqual3f> tempIntersectionLines;
		//		std::unordered_set<Vector3f, Vector3f::Hash<>> intersectionPoints;

		//		constexpr uint8_t trianglePointCount = 3;
		//		Index triangleIndices[trianglePointCount]
		//		{
		//			tetrahedronIndices[0 + (j <= 0)],
		//			tetrahedronIndices[1 + (j <= 1)],
		//			tetrahedronIndices[2 + (j <= 2)],
		//		};

		//		// For each line segment in the triangle...
		//		for (uint8_t k = 0; k < trianglePointCount; k++)
		//		{
		//			Vector4f pointA = inputMesh.vertices[triangleIndices[k]].position;
		//			Vector4f pointB = inputMesh.vertices[triangleIndices[(k + 1) % trianglePointCount]].position;

		//			constexpr float slice = 0.0f;
		//			float pointASlice = pointA[3] - slice;
		//			float pointBSlice = pointB[3] - slice;

		//			// If hyperplane intersects entire line segment...
		//			if (pointASlice == pointBSlice)
		//				tempIntersectionLines.emplace(static_cast<Vector3f>(pointA), static_cast<Vector3f>(pointB));
		//			// If hyperplane intersects with only point A...
		//			else if (pointASlice == 0.0f)
		//				intersectionPoints.insert(pointA);
		//			// If hyperplane intersects with only point B...
		//			else if (pointBSlice == 0.0f)
		//				intersectionPoints.insert(pointB);
		//			// If hyperplane intersects point along line segment...
		//			else if (gcem::signbit(pointASlice) != gcem::signbit(pointBSlice))
		//			{
		//				float t = pointASlice / (pointASlice - pointBSlice);
		//				if (0.0f <= t && t <= 1.0f)
		//					intersectionPoints.insert(Lerp(static_cast<Vector3f>(pointA), static_cast<Vector3f>(pointB), t));
		//				else
		//					__debugbreak(); // TODO: I don't think this is possible?
		//			}
		//		}

		//		// If hyperplane intersects entire triangle...
		//		if (tempIntersectionLines.size() == 3)
		//		{
		//			auto it = tempIntersectionLines.cbegin();
		//			const auto& crLine1 = *it;
		//			const auto& crLine2 = *++it;
		//			const auto& crLine3 = *++it;
		//			intersectionTriangles.insert({ crLine1, crLine2, crLine3 });
		//			tempIntersectionLines.clear();
		//		}
		//		else if (!tempIntersectionLines.empty())
		//		{
		//			for (const auto& crLine : tempIntersectionLines)
		//				intersectionLines.insert(crLine);
		//			tempIntersectionLines.clear();
		//		}

		//		if (intersectionPoints.size() == 2)
		//		{
		//			auto it = intersectionPoints.cbegin();
		//			const auto& crPointA = *it;
		//			const auto& crPointB = *++it;
		//			intersectionLines.emplace(crPointA, crPointB);
		//			intersectionPoints.clear();
		//		}
		//		else if (!intersectionPoints.empty())
		//		{
		//			//__debugbreak();
		//			intersectionPoints.clear();
		//		}
		//	}

		//	// If intersection is a triangle...
		//	if (intersectionLines.size() == 3)
		//	{
		//		auto it = intersectionLines.cbegin();
		//		const auto& crLine1 = *it;
		//		const auto& crLine2 = *++it;
		//		const auto& crLine3 = *++it;
		//		intersectionTriangles.insert({ crLine1, crLine2, crLine3 });
		//		intersectionLines.clear();
		//	}
		//	// If intersection is a quad or tetrahedron...
		//	else if (intersectionLines.size() == 4)
		//	{

		//	}
		//	else if (!intersectionLines.empty())
		//	{
		//		//__debugbreak();
		//		intersectionLines.clear();
		//	}
		//}



		//__debugbreak();

		//Slicer<float, 4> slicer;

		//if (slicer.Slice(inputMesh))
		//{
		//	std::cout << "Mesh to be rendered.\n";
		//}
		//else if (SlicerErrorCode errorCode = slicer.GetErrorCode(); errorCode != SlicerErrorCode_None)
		//{
		//	std::cout << "Error while slicing mesh: code = " << +errorCode << '\n' << slicer.GetErrorMessage() << '\n';
		//}
		//else
		//{
		//	std::cout << "Mesh did not intersect hyperplane.\n";
		//}

		Matrix<4, 4, signed char> a
		(
			0, true, 2.0, 3.0f,
			4, 5, 6, 7,
			8, 9, 10, CHAR_MIN,
			12, 13, 14, UCHAR_MAX
		);

		Matrix4f b = a;

		std::cout << a << '\n' << '\n';
		std::cout << b << '\n' << '\n';

		std::cin.get();
	}
}
