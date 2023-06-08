#ifndef BZMSH_HELPERFUNCTIONS_HPP
#define BZMSH_HELPERFUNCTIONS_HPP

#include "algorithm/constants.hpp"
#include "util/logger.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <gmpxx.h>
#include <list>
#include <map>
#include <set>
#include "algorithm/timer.hpp"

namespace bzmsh
{
using std::list;
using std::map;
using std::pair;
using std::set;
using std::vector;

// Forward declarations
template <typename T> struct BezierMesh;
template <typename T> struct BezierCurve;
template <typename T> struct GuardedBezierCurve;
template <typename T> struct Vec2;

/**
 * @brief Templated absolute function for use with both double and mpq_class
 *
 * @tparam T float type (double or mpq_class)
 * @param t
 * @return T abs of @p t
 */
template <typename T> T utilAbs(T t)
{
    if (t < 0)
        return -t;
    return t;
}

/**
 * @brief Utility function because mpq_class does not support typecasting to double
 *
 * @tparam T_IN type to cast from
 * @tparam T_OUT type to cast to
 * @param t
 * @return T_OUT t cast to type T_OUT
 */
template <typename T_IN, typename T_OUT> inline T_OUT utilCast(T_IN t)
{
    return static_cast<T_OUT>(t);
}

template <> inline double utilCast<mpq_class, double>(mpq_class t)
{
    return t.get_d();
}

/**
 * @brief Check equality within epsilon margin.
 *
 * @tparam T float type (to be used with double and mpq_class)
 * @param t1 value 1
 * @param t2 value 2
 * @param epsilon margin
 * @return true if equal within margin
 * @return false else
 */
template <typename T> bool eqWithEpsilon(T t1, T t2, T epsilon)
{
    return utilAbs<T>(t2 - t1) < epsilon;
}

/**
 * @brief Factorial
 *
 * @param n number
 * @return int n!
 */
double fac(int n);

/**
 * @brief Angle between 2 vectors
 *
 * @tparam T float-type (to be used with double and mpq_class)
 * @param v1 vector 1
 * @param v2 vector 2
 * @return double angle
 */
template <typename T> double angleRad(const Vec2<T>& v1, const Vec2<T>& v2)
{
    double dot = utilCast<T, double>(v1 * v2);
    double norm = v1.length() * v2.length();
    if (norm == 0.0)
        dot = 0.0;
    else
        dot /= norm;
    return acos(dot);
}

/**
 * @brief Conversion from degrees to radian
 *
 * @param t angle in degrees
 * @return double t in radian
 */
double toRadian(double t);

/**
 * @brief Conversion from radian to degrees
 *
 * @param t angle in radian
 * @return double t in degrees
 */
double toDegree(double t);

/**
 * @brief Convert a vector from one float type to another
 *
 * @tparam T_IN type to convert from
 * @tparam T_OUT type to convert to
 * @param v vector
 * @return Vec2<T_OUT> v converted to Vec2<T_OUT>
 */
template <typename T_IN, typename T_OUT> Vec2<T_OUT> convertVec(const Vec2<T_IN>& v)
{
    return Vec2<T_OUT>(utilCast<T_IN, T_OUT>(v.x), utilCast<T_IN, T_OUT>(v.y));
}

template <typename T_IN, typename T_OUT> std::vector<Vec2<T_OUT>> convertVecVec(const std::vector<Vec2<T_IN>>& vv)
{
    vector<Vec2<T_OUT>> converted;
    for (const Vec2<T_IN>& v : vv)
    {
        converted.emplace_back(convertVec<T_IN, T_OUT>(v));
    }
    return converted;
}

/**
 * @brief Convert a curve from one float type to another
 *
 * @tparam T_IN type to convert from
 * @tparam T_OUT type to convert to
 * @param c curve
 * @return BezierCurve<T_OUT> c converted to BezierCurve<T_OUT>
 */
template <typename T_IN, typename T_OUT> BezierCurve<T_OUT> convertCurve(const BezierCurve<T_IN>& c)
{
    vector<Vec2<T_OUT>> ctrlptsConverted;
    for (const Vec2<T_IN>& pt : c.m_ctrlpts)
    {
        ctrlptsConverted.emplace_back(convertVec<T_IN, T_OUT>(pt));
    }
    BezierCurve<T_OUT> convC(ctrlptsConverted,
                             utilCast<T_IN, T_OUT>(c.m_tStart),
                             utilCast<T_IN, T_OUT>(c.m_tEnd),
                             c.m_id);
    convC.m_origTStart = utilCast<T_IN, T_OUT>(c.m_origTStart);
    convC.m_origTEnd = utilCast<T_IN, T_OUT>(c.m_origTEnd);
    convC.m_origId = c.m_origId;
    convC.m_origLen = utilCast<T_IN, T_OUT>(c.m_origLen);
    return convC;
}

/**
 * @brief Convert a guarded curve from one float type to another
 *
 * @tparam T_IN type to convert from
 * @tparam T_OUT type to convert to
 * @param c guarded curve
 * @return GuardedBezierCurve<T_OUT> c converted to GuardedBezierCurve<T_OUT>
 */
template <typename T_IN, typename T_OUT>
GuardedBezierCurve<T_OUT> convertGuardedCurve(const GuardedBezierCurve<T_IN>& c)
{
    BezierCurve<T_OUT> convC = convertCurve<T_IN, T_OUT>(c.m_curve);
    GuardedBezierCurve<T_OUT> convGuardedC(convC, false);
    convGuardedC.m_valid = c.m_valid;
    convGuardedC.m_slope = convertVec<T_IN, T_OUT>(c.m_slope);
    convGuardedC.m_slope1 = convertVec<T_IN, T_OUT>(c.m_slope1);
    convGuardedC.m_slope2 = convertVec<T_IN, T_OUT>(c.m_slope2);
    convGuardedC.m_guardPtLeft = convertVec<T_IN, T_OUT>(c.m_guardPtLeft);
    convGuardedC.m_guardPtRight = convertVec<T_IN, T_OUT>(c.m_guardPtRight);
    convGuardedC.m_hullPtLeft = convertVec<T_IN, T_OUT>(c.m_hullPtLeft);
    convGuardedC.m_hullPtRight = convertVec<T_IN, T_OUT>(c.m_hullPtRight);
    convGuardedC.m_priority = utilCast<T_IN, T_OUT>(c.m_priority);
    convGuardedC.m_timesCut = c.m_timesCut;
    return convGuardedC;
}

/**
 * @brief Guard all curves in a given vector
 *
 * @tparam T float type (to be used with double and mpq_class)
 * @param curves curves to guard
 * @param guardedCurves resulting guarded curves
 * @return true if successful
 * @return false else
 */
template <typename T, typename T2>
bool guardAll(const vector<BezierCurve<T>>& curves, vector<GuardedBezierCurve<T2>>& guardedCurves)
{
    for (uint i = 0; i < curves.size(); i++)
    {
        list<BezierCurve<T2>> segments = {convertCurve<T,T2>(curves[i])};
        
        uint numsplits = 0;
        while (!segments.empty())
        {
           
            BezierCurve<T2> segment = segments.front();
            segments.pop_front();

            guardedCurves.emplace_back(segment);
            if(!guardedCurves.back().m_valid)
            {
                guardedCurves.pop_back();
                numsplits++;
                vector<BezierCurve<T2>> bisectedSegment;
                if (numsplits > MAX_GUARDING_DEPTH || !bisect(segment, T2(0.5), bisectedSegment))
                {
                    Logger::lout(Logger::DEBUG)
                        << "Guarding failed for curve " << curves[i] << endl;
                    return false;
                }
                
                segments.insert(segments.begin(), bisectedSegment.begin(), bisectedSegment.end());
            }
        }
    }
    return true;
}

/**
 * @brief Get the shortest distance of a point from an infinite line
 *
 * @tparam T float-type (to be used with double and mpq_class)
 * @param linePt1 first point on line
 * @param linePt2 second point on line
 * @param pt point
 * @return T distance of pt frome line
 */
template <typename T>
T distanceToLine(const Vec2<T>& linePt1, const Vec2<T>& linePt2, const Vec2<T>& pt)
{
    Vec2<T> lineDir = linePt2 - linePt1;
    T norm = lineDir.length();
    return utilAbs<T>(pt.cross(lineDir) + linePt2.cross(linePt1)) / norm;
}

/**
 * @brief Get the perpendicular projection of a point on a line
 *
 * @tparam T float-type (to be used with double and mpq_class)
 * @param linePt1 first point on line
 * @param linePt2 second point on line
 * @param pt point
 * @return Vec2<T> pt projected onto line
 */
template <typename T>
Vec2<T> projectOntoLine(const Vec2<T>& linePt1, const Vec2<T>& linePt2, const Vec2<T>& pt)
{
    Vec2<T> lineDir = linePt2 - linePt1;
    return linePt1 + lineDir * (lineDir * (pt - linePt1)) / (lineDir * lineDir);
}

/**
 * @brief Relative orientation of a point to a directed line. AKA counterclockwise.
 *
 */
const int ORI_LEFT = 1;
/**
 * @brief Relative orientation of a point to a directed line. AKA clockwise.
 *
 */
const int ORI_RIGHT = -1;
/**
 * @brief Relative orientation of a point to a directed line.
 *
 */
const int ORI_COLINEAR = 0;
/**
 * @brief Get the relative orientation of a point to a directed line
 *
 * @tparam T float-type (to be used with double and mpq_class)
 * @param linePt1 first point on line
 * @param linePt2 second point on line
 * @param pt point
 * @return int orientation (-1, 0 or 1)
 */
template <typename T> int orientation(Vec2<T> linePt1, Vec2<T> linePt2, Vec2<T> pt)
{
    T val = (linePt2.x - linePt1.x) * (pt.y - linePt1.y)
            - (linePt2.y - linePt1.y) * (pt.x - linePt1.x);

    if (val == T(0))
        return ORI_COLINEAR; // colinear

    return (val > 0) ? ORI_LEFT : ORI_RIGHT; // left or right
}

/**
 * @brief Utility function used to get subtriangle control points from original triangle
 *
 */
template <typename T> Vec2<T> interpolate(const std::vector<Vec2<T>>& C,
                         std::vector<int>& indC,
                         int i,
                         int j,
                         int n1,
                         int n2,
                         int n3);

/**
 * FOR ANYTHING ORDER RELATED REFER TO THIS
 *
 *  RECURSIVE            ROWWISE
 *
 *  2                    14
 *  | \                  | \
 *  9   8                12  13
 *  |     \              |     \
 * 10  14   7            9  10  11
 *  |         \          |         \
 * 11  12  13   6        5   6   7   8
 *  |             \      |             \
 *  0---3---4---5---1    0---1---2---3---4
 *
 * REVERSE rowwise means going from right to left in every row
 * REVERSE recursive means going from clockwise instead of counterclockwise
 *
 * The recursive ordering of inner control points for higher degrees (as the name implies) uses the
 * same order (vertices -> edges -> inner)
 */

/**
 * @brief Retrieve sequential recursive ordering of control points from per-row control points
 *
 * @tparam T any control point descriptor (indices or positions etc.)
 * @param rows rows of control point descriptors
 * @param recursive resulting sequential control point descriptors in recursive ordering
 * @param reverse whether ordering should be reversed
 */
template <typename T>
void recursiveFromRows(std::vector<std::vector<T>>& rows,
                       std::vector<T>& recursive,
                       bool reverse = false)
{
    if (reverse)
    {
        for (std::vector<T>& row : rows)
        {
            std::reverse(row.begin(), row.end());
        }
    }
    recursive.clear();
    while (!rows.empty())
    {
        uint nRows = rows.size();
        if (nRows == 1)
        {
            recursive.emplace_back(rows[0][0]);
            break;
        }
        if (nRows == 2)
        {
            recursive.insert(recursive.end(), rows[0].begin(), rows[0].end());
            recursive.emplace_back(rows[1][0]);
            break;
        }

        recursive.emplace_back(rows[0][0]);
        recursive.emplace_back(rows[0][nRows - 1]);
        recursive.emplace_back(rows[nRows - 1][0]);
        for (uint i = 1; i < nRows - 1; i++)
        {
            recursive.emplace_back(rows[0][i]);
        }
        for (uint i = 1; i < nRows - 1; i++)
        {
            recursive.emplace_back(rows[i].back());
        }
        for (uint i = nRows - 2; i >= 1; i--)
        {
            recursive.emplace_back(rows[i].front());
        }

        rows.erase(rows.end() - 1);
        rows.erase(rows.end() - 1);
        rows.erase(rows.begin());
        for (uint i = 0; i < rows.size(); i++)
        {
            rows[i].erase(rows[i].end() - 1);
            rows[i].erase(rows[i].begin());
        }
    }
}

/**
 * @brief Retrieve sequential recursive ordering of control points from sequential rowwise control
 * points
 *
 * @tparam T any control point descriptor (indices or positions etc.)
 * @param rowwise control points in sequential rowwise ordering
 * @param recursive resulting control points in sequential recursive ordering
 * @param degree curve degree
 * @param reverse whether ordering should be reversed
 */
template <typename T>
void reorderToRecursive(const std::vector<T>& rowwise,
                        std::vector<T>& recursive,
                        int degree,
                        bool reverse = false)
{
    recursive.clear();
    uint nRows = degree + 1;
    std::vector<std::vector<T>> rows(nRows);
    uint pos = 0;
    uint row = 0;
    while (pos < rowwise.size())
    {
        uint elementsPerRow = nRows - row;
        for (uint element = 0; element < elementsPerRow; element++)
        {
            rows[row].emplace_back(rowwise[pos++]);
        }
        row++;
    }

    recursiveFromRows(rows, recursive, reverse);
    assert(rowwise.size() == recursive.size());
}

/**
 * @brief Retrieve sequential rowwise ordering of control points from sequential recursive control
 * points
 *
 * @tparam T any control point descriptor (indices or positions etc.)
 * @param recursive resulting control points in sequential recursive ordering
 * @param rowwise control points in sequential rowwise ordering
 * @param degree curve degree
 * @param reverse whether ordering should be reversed
 */
template <typename T>
void reorderToRowwise(const std::vector<T> recursive,
                      std::vector<T>& rowwise,
                      int degree,
                      bool reverse = false)
{
    rowwise = std::vector<T>(recursive.size());
    std::vector<int> indices(recursive.size());
    for (uint i = 0; i < recursive.size(); i++)
        indices[i] = i;

    std::vector<int> indicesRecursive;
    reorderToRecursive(indices, indicesRecursive, degree, reverse);

    int i = 0;
    for (int index : indicesRecursive)
    {
        rowwise[index] = recursive[i++];
    }
}

/**
 * @brief Check validity of bezier triangles for each triangle of a mesh.
 *
 * @param mesh Bezier mesh
 * @return true if all triangles are valid
 * @return false if any triangle is invalid (flipped)
 */
template <typename T>
bool checkMesh(BezierMesh<T>& mesh, bool& uncertain);

using Split = std::pair<const mpq_class, Vec2<mpq_class>>;

/**
 * @brief Comparator for mpq_class with epsilon margin
 *
 */
struct cmpMPQWithIntersectionEpsilon
{
    bool operator()(const mpq_class& lhs, const mpq_class& rhs) const
    {
        return !eqWithEpsilon(lhs, rhs, INTERSECTION_T_EPSILON) && lhs < rhs;
    }
};

using SetOfSplits = std::map<mpq_class, Vec2<mpq_class>, cmpMPQWithIntersectionEpsilon>;

} // namespace bzmsh

#endif
