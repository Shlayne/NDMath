#pragma once

#include "Impl.hpp"
#include "Matrix.hpp"
#include <gcem.hpp>
#include <random>

#define _ND ::nd::
#define _IMPL ::nd::impl::
#define _GCEM ::gcem::

namespace nd
{
	// The basis vectors for the hyperplane are each one column in the given bases matrix.
	// They are expected to be normalized. If not, the random values will be scaled by their magnitudes.
	template <Scalar S, Dimension N, _IMPL Function<S, S, S> RandomFunction>
	constexpr Vector<S, N> RandomPointInBoundedHyperplane(const Matrix<S, N, N>& bases, const Vector<S, N>& minBounds, const Vector<S, N>& maxBounds, const RandomFunction& randomFunction) noexcept
	{
		Vector<S, N> randomPoint;
		for (Dimension n{}; n < N; ++n)
			randomPoint += bases[n] * randomFunction(minBounds[n], maxBounds[n]);
		return randomPoint;
	}
}
