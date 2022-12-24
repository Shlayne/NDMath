export module nd.scalar;

import <concepts>;

export namespace nd
{
	template<typename T>
	concept Scalar = std::integral<T> || std::floating_point<T>;
}
