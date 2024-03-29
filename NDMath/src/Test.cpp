#include "ND/ND.hpp"
#include <iostream>

using namespace nd;

namespace test
{
	template <Scalar S>
	S RandomRange(S min, S max)
	{
		return (rand() / static_cast<std::common_type_t<S, float>>(RAND_MAX)) * (max - min) + min;
	}

	void Perform()
	{
		{
			time_t time_ = time(nullptr);
			srand(static_cast<unsigned int>(time_) ^ static_cast<unsigned int>(time_ >> 32));
			static_cast<void>(rand());
		}

		//constexpr Tensor<float> a;
		//constexpr Tensor<float, 5> b;
		//constexpr Tensor<float, 5, 3> c;
		//constexpr Tensor<float, 5, 3, 6> d;
		//Tensor<float, 4, 5, 3, 6, 2> e;

		//constexpr Tensor<float, 3, 6, 2> e0 = e[4];
		//constexpr Tensor<float, 6, 2> e1 = e0[2];
		//constexpr Tensor<float, 2> e2 = e1[5];
		//constexpr Tensor<float> e3 = e2[1];
		//constexpr float e3a = e2[1];
		//constexpr float e3b = e3;

		Tensor<float, 4> a;

		a[0] = 1;
		a[1] = 3;
		a[2] = -4;
		a[3] = 2.0;
		//a[4] = 5.7; // aborts (good)
		std::cout << (a != a) << '\n';

		std::cout << a.at<3>() << '\n'; // compiles (good)
		//a.at<4>(); // doesn't compile (good)

		Vector4f _ab{a * 2};

		Tensor<float, 5> b;
		b = 6;
		Vector3f e{Vector4f{5.9, 2.6f}};
		Vector2f f;
		Vector2f _af = b;

		std::cout << a.Hash() << '\n';

		//std::cout << b << '\n';

		Tensor<long, 4> _ac = a;
		auto c = static_cast<Tensor<float, 4>>(a);

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

		Matrix4f m;
		Matrix6f n = 3;
		Matrix8f o = n;
		Matrix2f p = o;

		m = 2;
		n = 2;
		m = n;
		n = m;

		m += 5;
		m = m + 5;
		m += n;
		m = m + n;

		m -= 5;
		m = m - 5;
		m -= n;
		m = m - n;

		m *= 5;
		m = m * 5;
		m *= n;
		m = m * n;

		m /= 5;
		m = m / 5;
		m /= n;
		m = m / n;

		m = 2 + m;
		m = 2 - m;
		m = 2 * m;
		m = 2 / m;

		std::cout << m.at<3>() << '\n'; // compiles (good)
		//m.at<4>(); // doesn't compile (good)
		std::cout << n.at<5>() << '\n'; // compiles (good)
		//n.at<6>(); // doesn't compile (good)
		std::cout << (m != m) << '\n';

		Matrix4f x
		{
			0, true, 2.0, 3.0f,
			4, 5, 6, 7,
			8, 9, 10, -128,
			0, 13, 14, -1
		};

		//auto y{Determinant(x)};
		//auto z{Inverse(x)};

		//std::cout << x << '\n' << '\n';
		//std::cout << y << '\n' << '\n';
		//std::cout << z << '\n' << '\n';

		Vector4f _a{1, 2, 3, 4};
		Vector3f _b{5, 6, 7};

		auto _c{OuterProduct(_a, x)};
		std::cout << InnerProduct(_a, _a) << ' ' << InnerProduct(_b, _b) << '\n';


		// example hyperplane basis vectors
		Matrix4f hyperplaneBases
		{
			1, 0, 3, 2,
			0, 0, 8, 0,
			2, 1, 0, 0,
			0, 4, 0, 7,
		};

		for (Dimension n{}; n < 4; ++n)
			// TODO: remove <float, 4> requirement.
			// It's only there cause it can't tell how to convert a Tensor<S, N> to a Vector<S, N>
			// without the Vector<S, N> being predefined. In this case, it's not predefined, it's
			// templated which takes the template args of the parameter, but can't convert.
			// The conversion operator works, but if it doesn't know what to convert to, it's stuck.
			hyperplaneBases[n] = Normalize(hyperplaneBases[n]);
		std::cout << "Hyperplane Bases:\n" << hyperplaneBases << '\n';

		// example bounds where
		// x is bounded to [-3, 10],
		// y is bounded to [1, 3],
		// z is bounded to [0, 0],
		// w is bounded to [0, 0],
		Vector4f minBounds = {0.0f, 0.0f, 10.0f, -10.0f};
		Vector4f maxBounds = {1.0f, -1.0f, 100.0f, -100.0f};

		const auto nml = Normalize(maxBounds);
		std::cout << "Normal.x: " << nml[0] << '\n';
		//nml = maxBounds;
		std::cout << "Normal: " << nml << '\n';
		std::cout << "Normal Length: " << Length(nml) << '\n';
		auto nml_1 = nml + 2; // compiles (good)
		//auto nml_2 = nml += 2; // doesn't compile (good)
		auto nml_3 = nml - 2; // compiles (good)
		//auto nml_4 = nml -= 2; // doesn't compile (good)
		auto nml_5 = nml * 2; // compiles (good)
		//auto nml_6 = nml *= 2; // doesn't compile (good)
		auto nml_7 = nml / 2; // compiles (good)
		//auto nml_8 = nml /= 2; // doesn't compile (good)

		Vector4f randomPoint = RandomPointInBoundedHyperplane(hyperplaneBases, minBounds, maxBounds, test::RandomRange<float>);
		std::cout << "Random Point: " << randomPoint << '\n';

		std::cin.get();
	}
}
