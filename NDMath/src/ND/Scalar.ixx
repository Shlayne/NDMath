export module nd.scalar;

import <concepts>;

export namespace nd
{
	template<typename T>
	concept Scalar = _STD integral<T> || _STD floating_point<T>;
}
