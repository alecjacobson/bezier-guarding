#ifndef BZMSH_CURVE_ENVELOPE_HPP
#define BZMSH_CURVE_ENVELOPE_HPP

#include "algorithm/constants.hpp"
#include "algorithm/intersection_tests.hpp"
#include "datastructure/bezier_curve.hpp"
#include "util/helper_functions.hpp"
#include "util/logger.hpp"
#include "config.hpp"

namespace bzmsh
{
/**
 * @brief Represents a bezier curve that is completely enclosed in a quad made from two triangles.
 *        The center-facing edge of each of those two (higher-order) triangles is formed by this
 *        curve
 *
 * @tparam T float-type
 */
template <typename T> struct GuardedBezierCurve
{
    /**
     * @brief Constructs a new guard curve around a given bezier curve respecting some epsilon
     *        values
     *
     * @param c curve to guard
     */
    GuardedBezierCurve(const BezierCurve<T>& c, bool constr = true)
        : m_valid(false), m_curve(c), m_timesCut(0)
    {
        if (constr)
            construct();
    }

    GuardedBezierCurve() : m_valid(false)
    {
    }

    /**
     * @brief Construct the guarding quad for this curve
     *
     * @param epsilon amount by which to enlarge the guarding quad in direction of the curves
     *                average normal
     * @param maxDeviation maximum tolerated average curve turning angle IN DEGREES
     * @return true if successful
     * @return false else
     */
    bool construct()
    {
        m_valid = false;
        int d = m_curve.degree();
        assert(d > 0);

        if (!m_curve.isGuardable())
        {
            Logger::lout(Logger::DEBUG)
                << "Curve " << m_curve << " was not actually guardable" << endl;
            return false;
        }
        
        m_slope = m_curve.curveAxis();

        T len = (m_curve.endPt() - m_curve.startPt()).norm1();
        if(m_curve.m_origLen <= 0)
            m_curve.m_origLen = len;
        T shift = utilCast<mpq_class, T>(GUARDING_NORMAL_SHIFT) * (len*len / m_curve.m_origLen);

        Vec2<T> normal = m_slope.rot90ccw().normalizedNorm1();
        
        m_slope1 = m_curve.splus();
        m_slope2 = m_curve.sminus();

        T maxDistLeft = 0, maxDistRight = 0;
        m_hullPtLeft = m_curve[0];
        m_hullPtRight = m_curve[0];

        Vec2<T> nn = m_slope.rot90ccw();
        maxDistLeft = nn * m_curve[0];
        maxDistRight = nn * m_curve[0];
        for (int i = 1; i <= d; i++)
        {
          T dist = nn * m_curve[i];
          if (dist > maxDistLeft)
          {
              maxDistLeft = dist;
              m_hullPtLeft = m_curve[i];
          }
          if (dist < maxDistRight)
          {
              maxDistRight = dist;
              m_hullPtRight = m_curve[i];
          }
        }
        
        if(!getLineIntersection(m_curve[0], m_curve[0]+m_slope1, m_curve[d], m_curve[d]+m_slope2, m_guardPtLeft)
             || !getLineIntersection(m_curve[0], m_curve[0]+m_slope2, m_curve[d], m_curve[d]+m_slope1, m_guardPtRight))
        {
            //special case for straight lines
            m_guardPtLeft = (m_curve[0] + m_curve[d])/2.0;
            m_guardPtRight = m_guardPtLeft;
        }
        
        // add mu-shift:
        m_guardPtLeft = m_guardPtLeft + (normal * shift);
        m_guardPtRight = m_guardPtRight - (normal * shift);
          

        setPriority();

        m_valid = true;
        return true;
    }

    /**
     * @brief Set this curves priority based on the diameter of the guarding quad
     *
     */
    void setPriority()
    {
        m_priority = m_curve.curveAxis().rot90ccw().normalizedNorm1() * (m_guardPtLeft - m_guardPtRight);
    }

    /**
     * @brief Get an axis aligned bounding box that wraps this curves guarding quad
     *
     * @return BBox2 axis aligned bounding box that wraps this curves guarding quad
     */
    BBox<T> bbox() const
    {
        BBox<T> box(m_curve[0]);
        box.widen(m_curve[m_curve.degree()]);
        box.widen(m_guardPtLeft);
        box.widen(m_guardPtRight);
        return box;
    }

    /**
     * @brief Get the edges of this curves guarding quad IN COUNTERCLOCKWISE ORDER
     *
     * @return vector<LineSegment<T>> the edges of this curves guarding quad IN COUNTERCLOCKWISE
     * ORDER
     */
    vector<LineSegment<T>> quadLines() const
    {
        return {LineSegment<T>(m_curve.startPt(), m_guardPtRight),
                LineSegment<T>(m_guardPtRight, m_curve.endPt()),
                LineSegment<T>(m_curve.endPt(), m_guardPtLeft),
                LineSegment<T>(m_guardPtLeft, m_curve.startPt())};
    }

    /**
     * @brief Get the vertices of this curves guarding quad IN COUNTERCLOCKWISE ORDER
     *
     * @return vector<Vec2<T>> the vertices of this curves guarding quad IN COUNTERCLOCKWISE ORDER
     */
    vector<Vec2<T>> quadPoints() const
    {
        return {Vec2<T>(m_curve.startPt()),
                Vec2<T>(m_guardPtRight),
                Vec2<T>(m_curve.endPt()),
                Vec2<T>(m_guardPtLeft)};
    }

    /**
     * @brief Get the edges of the guarding triangles forming this curves guarding quad IN
     * COUNTERCLOCKWISE ORDER
     *
     * @return vector<vector<LineSegment<T>>> the edges of the guarding triangles forming this
     *         curves guarding quad IN COUNTERCLOCKWISE ORDER
     */
    vector<vector<LineSegment<T>>> triLines() const
    {
        return {{LineSegment<T>(m_curve.startPt(), m_guardPtRight),
                 LineSegment<T>(m_guardPtRight, m_curve.endPt()),
                 LineSegment<T>(m_curve.endPt(), m_curve.startPt())},
                {LineSegment<T>(m_curve.startPt(), m_curve.endPt()),
                 LineSegment<T>(m_curve.endPt(), m_guardPtLeft),
                 LineSegment<T>(m_guardPtLeft, m_curve.startPt())}};
    }

    /**
     * @brief Must be true, after constructing the guard and hull points (if successful)
     *        and false else.
     *
     */
    bool m_valid;
    /**
     * @brief Bezier curve that is getting guarded
     *
     */
    BezierCurve<T> m_curve;
    /**
     * @brief Slope vector used for construction
     *
     */
    Vec2<T> m_slope;
    Vec2<T> m_slope1, m_slope2;
    /**
     * @brief The two guard points oriented on left and right side of the ray defined by first and
     *        last curve controlpoints.
     *
     */
    Vec2<T> m_guardPtLeft, m_guardPtRight;
    /**
     * @brief Intermediate hulll point (outermost control points on left/right side of the ray
     *        defined by first and last curve controlpoints).
     *
     */
    Vec2<T> m_hullPtLeft, m_hullPtRight;
    /**
     * @brief This curves priority for bisections (bigger quad diameter -> higher priority)
     *
     */
    T m_priority;
    /**
     * @brief How many times this guarded curve was cut. Used for detecting unsresolvable quad
     *        intersections.
     *
     */
    int m_timesCut;
    
};

} // namespace bzmsh

#endif
