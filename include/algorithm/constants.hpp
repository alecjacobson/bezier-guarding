#ifndef BZMSH_CONSTANTS_HPP
#define BZMSH_CONSTANTS_HPP

#include "gmpxx.h"

namespace bzmsh
{
/**
 * @brief Maximum number of bisections per curve in intersection tests.
 */
extern const int MAX_INTERSECTION_DEPTH;
/**
 * @brief Maximum number of bisections per curve during envelope bisection.
 */
extern const int MAX_GUARDING_DEPTH;
/**
 * @brief Minimum bounding box circumference for intersection pre-tests.
 */
extern const mpq_class MIN_INTERSECTION_CIRCUMFERENCE;
/**
 * @brief Cache depth of bisected curves in curve-curve-intersection test.
 */
extern const int INTERSECTION_CACHE_DEPTH;
/**
 * @brief t-Parameter margin around previously found curve intersection-t inside of which any new
 *        intersection will be snapped to the same point as the first intersection.
 */
extern const mpq_class INTERSECTION_T_EPSILON;
/**
 * @brief Minimum allowed angle between the tangents of any two curves with common endpoints (in
 * degrees).
 */
extern const double MIN_TANGENT_ANGLE_DEGREES;
/**
 * @brief Amount by which to shift guarding points of curves (non-curve quad corners) outwards
 */
extern mpq_class GUARDING_NORMAL_SHIFT; // mu
/**
 * @brief Maximum allowed turning angle of a curve towards any direction (CW/CCW)
 */
extern const double GUARDING_MAX_CURVATURE_DEGREES;
/**
 * @brief How much padding will be between the (bounding) box created for CDT and the contained curves.
 *
 * Will be applied as fraction of total height/width.
 *
 */
extern const double CDT_BOX_PADDING;
} // namespace bzmsh

#endif
