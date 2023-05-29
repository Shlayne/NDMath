export module test;

import nd;
import <iostream>;

using namespace nd;

export namespace test
{
	auto Perform() -> void
	{
		//constexpr Tensor<float, 0> a;
		//constexpr Tensor<float, 1, 5> b;
		//constexpr Tensor<float, 2, 5, 3> c;
		//constexpr Tensor<float, 3, 5, 3, 6> d;
		//Tensor<float, 4, 5, 3, 6, 2> e;

		//constexpr Tensor<float, 3, 3, 6, 2> e0 = e[4];
		//constexpr Tensor<float, 2, 6, 2> e1 = e0[2];
		//constexpr Tensor<float, 1, 2> e2 = e1[5];
		//constexpr Tensor<float, 0> e3 = e2[1];
		//constexpr float e3a = e2[1];
		//constexpr float e3b = e3;

		Tensor<float, 1, 4> a;
		a[0] = 1;
		a[1] = 3;
		a[2] = -4;
		a[3] = 2.0;
		//a[4] = 5.7; // aborts (good)
		_STD cout << (a != a) << '\n';

		_STD cout << a.at<3>() << '\n'; // compiles (good)
		//a.at<4>(); // doesn't compile (good)

		//Vector4f b{a * 2}; // for some reason, this crashes intellisense... :/

		Tensor<float, 1, 5> b;
		b = 6;
		Vector3f e{Vector4f{5.9, 2.6f}};
		Vector1f f = b;

		_STD cout << a.Hash() << '\n';

		//_STD cout << b << '\n';

		//Tensor<long, 1, 4> c = static_cast<Tensor<int, 1, 4>>(a); // TODO: doesn't compile.
		//auto c = static_cast<Tensor<float, 1, 4>>(a);

		e = f;
		f = e;

		a = -a;
		auto d = +a;

		a += 5;
		a = a + 5;
		a += b;
		a = a + b;

		a -= 5;
		a = a - 5;
		a -= b;
		a = a - b;

		a *= 5;
		a = a * 5;
		a *= b;
		a = a * b;

		a /= 5;
		a = a / 5;
		a /= b;
		a = a / b;

		a = 2 + a;
		a = 2 - a;
		a = 2 * a;
		a = 2 / a;

		//Matrix4f x
		//{
		//	0, true, 2.0, 3.0f,
		//	4, 5, 6, 7,
		//	8, 9, 10, -128,
		//	0, 13, 14, -1
		//};

		//auto y{Determinant(x)};
		//auto z{Inverse(x)};

		//_STD cout << x << '\n' << '\n';
		//_STD cout << y << '\n' << '\n';
		//_STD cout << z << '\n' << '\n';

		_STD cin.get();
	}
}
