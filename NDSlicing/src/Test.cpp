#include "Test.h"
#include "ND/Loader.h"
#include "ND/Slicer.h"
#include <iostream>

using namespace nd;

// TODO: Loader.inl completion.

namespace test
{
	void Perform()
	{
		constexpr auto aaaaaa = Vector4f(1, 2, 3, 4);
		constexpr Vector3f bbbbbb = aaaaaa.operator Vector3f();
		constexpr Vector3f cccccc = (Vector3f)aaaaaa;

		Mesh<float, 2> inputMesh;
		Matrix<float, 3> transformation;

		Slicer<float, 2> slicer;

		if (slicer.Slice(inputMesh))
		{
			std::cout << "Mesh to be rendered.\n";
		}
		else if (SlicerErrorCode errorCode = slicer.GetErrorCode(); errorCode != SlicerErrorCode_None)
		{
			std::cout << "Error while slicing mesh: code = " << +errorCode << '\n' << slicer.GetErrorMessage() << '\n';
		}
		else
		{
			std::cout << "Mesh did not intersect hyperplane.\n";
		}

		std::cin.get();
	}
}
