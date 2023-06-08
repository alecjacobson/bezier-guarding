#include "algorithm/constants.hpp"

namespace bzmsh
{
const int MAX_INTERSECTION_DEPTH = 30;
const int MAX_GUARDING_DEPTH = 200;
const int INTERSECTION_CACHE_DEPTH = 0;
const mpq_class MIN_INTERSECTION_CIRCUMFERENCE = mpq_class(1, (int)1e8);
const mpq_class INTERSECTION_T_EPSILON = mpq_class(1, (int)1e6);
const double MIN_TANGENT_ANGLE_DEGREES = 0.001;
const double GUARDING_MAX_CURVATURE_DEGREES = 179;
mpq_class GUARDING_NORMAL_SHIFT = mpq_class(1, 100); //mu
const double CDT_BOX_PADDING = 0.05;

} // namespace BZMSH
