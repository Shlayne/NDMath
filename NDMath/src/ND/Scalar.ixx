export module nd.scalar;

import std.core;

export namespace nd
{
	template<typename T>
	concept Scalar = std::integral<T> || std::floating_point<T>;
}
