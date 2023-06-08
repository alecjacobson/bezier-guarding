#ifndef BZMSH_INTERSECTIONTESTS_HPP
#define BZMSH_INTERSECTIONTESTS_HPP

#include "datastructure/bezier_curve.hpp"
#include "datastructure/line_segment.hpp"
#include "util/helper_functions.hpp"

namespace bzmsh
{

/**
 * @brief Find intersection point between two infinite lines if it exists.
 *
 * @tparam T float-type (used for double and mpq_class)
 * @param l1Pt1 First point on line 1
 * @param l1Pt2 Second point on line 1
 * @param l2Pt1 First point on line 2
 * @param l2Pt2 Second point on line 2
 * @param intersection resulting intersection (output)
 * @return false if none exists
 * @return true else
 */
template <typename T>
bool getLineIntersection(const Vec2<T>& l1Pt1,
                         const Vec2<T>& l1Pt2,
                         const Vec2<T>& l2Pt1,
                         const Vec2<T>& l2Pt2,
                         Vec2<T>& intersection)
{
    T a1 = l1Pt2.y - l1Pt1.y;
    T b1 = l1Pt1.x - l1Pt2.x;
    T c1 = a1 * l1Pt1.x + b1 * l1Pt1.y;

    T a2 = l2Pt2.y - l2Pt1.y;
    T b2 = l2Pt1.x - l2Pt2.x;
    T c2 = a2 * l2Pt1.x + b2 * l2Pt1.y;

    T determinant = a1 * b2 - a2 * b1;

    if (determinant == T(0))
    {
        return false;
    }
    intersection.x = (b2 * c1 - b1 * c2) / determinant;
    intersection.y = (a1 * c2 - a2 * c1) / determinant;
    return true;
}

/**
 * @brief Simple test if point lies in convex polygon.
 *
 * @param poly polygon given as ordered vector of lines (CW or CCW)
 * @param pt point to query
 * @param clockWise whether poly is CW or CCW
 * @param includeBorder whether polygon border should be handled as inside
 * @return true if inside
 * @return false if outside
 */
bool inConvexPoly(const vector<LineSegment<mpq_class>>& poly,
                  const Vec2<mpq_class>& pt,
                  bool clockWise = false,
                  bool includeBorder = false);

bool halfOpenSegmentsIntersect(const LineSegment<mpq_class>& seg1,
                       const LineSegment<mpq_class>& seg2,
                       bool start1open,
                       bool start2open);

bool guardingQuadsIntersect(const vector<LineSegment<mpq_class>>& poly1,
                          const vector<LineSegment<mpq_class>>& poly2);

/**
 * @brief Auxiliary method to register a new intersection between two
 *        curve segments.
 *
 * @param segmentI first segment
 * @param segmentJ second segment
 * @param epsilon to determine whether an intersection is too close to an endpoint to be a true
 *                intersection
 * @param curveIdToSplits the map the curve id, intersection t and curve(t) will be saved in
 * @return true if true intersection and map has been updated
 * @return false else
 */
bool addNewIntersection(const BezierCurve<mpq_class>& segmentI,
                        const BezierCurve<mpq_class>& segmentJ,
                        std::map<uint, SetOfSplits>& curveIdToSplits);

/**
 * @brief Gather all the intersection t-values and points in a map,
 *        that maps from curve id to a set of intersection points
 *
 * @param curves curves to check for intersections
 * @param curveIdToSplits map from curve id to a set of intersection points
 * @return false if there was curve overlap (i.e. possibly infinite intersections)
 * @return true else
 */
bool getIntersections(const vector<BezierCurve<mpq_class>> curves,
                      std::map<uint, SetOfSplits>& curveIdToSplits);

} // namespace bzmsh
#endif
