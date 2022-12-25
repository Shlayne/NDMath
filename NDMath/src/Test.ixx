export module test;

import nd;
import <iostream>;

using namespace nd;

export namespace test
{
	auto Perform() -> void
	{
		Matrix4f a
		(
			0, true, 2.0, 3.0f,
			4, 5, 6, 7,
			8, 9, 10, -128,
			0, 13, 14, -1
		);

		auto b = Determinant(a);
		auto c = Inverse(a);

		std::cout << a << '\n' << '\n';
		std::cout << b << '\n' << '\n';
		std::cout << c << '\n' << '\n';

		std::cin.get();
	}
}
