export module test;

import nd;
import <iostream>;

using namespace nd;

export namespace test
{
	void Perform()
	{
		Matrix<4, 4, signed char> a
		(
			0, true, 2.0, 3.0f,
			4, 5, 6, 7,
			8, 9, 10, -128,
			12, 13, 14, -1
		);

		Matrix4f b = a;

		std::cout << a << '\n' << '\n';
		std::cout << b << '\n' << '\n';

		std::cin.get();
	}
}
