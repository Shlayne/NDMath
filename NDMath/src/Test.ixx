export module test;

import nd;
import <iostream>;

using namespace nd;

export namespace test
{
	auto Perform() -> void
	{
		constexpr Tensor<float, 0> a;
		constexpr Tensor<float, 1, 5> b;
		constexpr Tensor<float, 2, 5, 3> c;
		constexpr Tensor<float, 3, 5, 3, 6> d;
		constexpr Tensor<float, 4, 5, 3, 6, 2> e;

		constexpr Tensor<float, 3, 3, 6, 2> e0 = e[4];
		constexpr Tensor<float, 2, 6, 2> e1 = e0[2];
		constexpr Tensor<float, 1, 2> e2 = e1[5];
		constexpr Tensor<float, 0> e3 = e2[1];
		constexpr float e3a = e3;

		//Matrix4f a
		//{
		//	0, true, 2.0, 3.0f,
		//	4, 5, 6, 7,
		//	8, 9, 10, -128,
		//	0, 13, 14, -1
		//};

		//auto b{Determinant(a)};
		//auto c{Inverse(a)};

		//_STD cout << a << '\n' << '\n';
		//_STD cout << b << '\n' << '\n';
		//_STD cout << c << '\n' << '\n';

		_STD cin.get();
	}
}
