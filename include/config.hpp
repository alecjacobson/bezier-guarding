#ifndef BZMSH_CONFIG_HPP
#define BZMSH_CONFIG_HPP


// - in ALL_EXACT mode, exact predicates and constructions are used.
// - otherwise, exact predicates are used (for intersection tests),
//   but control points and envelopes are constructed in double precision.
//   This comes with the risk of numerical degeneracies in complicated configurations


// Note: ALL_EXACT is defined via CMake
#ifdef ALL_EXACT
typedef mpq_class Scalar;
#else
typedef double Scalar;
#endif

#endif
