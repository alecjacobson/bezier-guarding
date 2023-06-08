#ifndef BZMSH_CURVE_HPP
#define BZMSH_CURVE_HPP

#include <cassert>
#include <vector>

#include "datastructure/bbox.hpp"
#include "datastructure/vec2.hpp"
#include "util/logger.hpp"

namespace bzmsh
{
using std::vector;

/**
 * @brief For representing 2D curves (especially bezier curves and line segments)
 *
 * @tparam T float-type
 */
template <typename T> struct Curve
{
    /**
     * @brief Construct a new curve
     *
     * @param ctrlpts The curves control points
     * @param tStart Starting t
     * @param tEnd Ending t
     * @param id curve id
     */
    Curve(const vector<Vec2<T>>& ctrlpts, T tStart = 0, T tEnd = 1, int id = 0)
        : m_tStart(tStart), m_tEnd(tEnd), m_origTStart(tStart), m_origTEnd(tEnd), m_id(id),
          m_origId(id), m_ctrlpts(ctrlpts), m_origLen(-1.0)
    {
        if (ctrlpts.size() < 2)
        {
            Logger::lout(Logger::WARN) << "Curve with <2 control points created!" << endl;
        }
    }

    /**
     * @brief Construct a new default curve
     *
     */
    Curve() : m_tStart(0), m_tEnd(1), m_origTStart(0), m_origTEnd(1), m_id(-1), m_origId(-1), m_origLen(-1.0)
    {
    }
    
    virtual ~Curve() {}

    /**
     * @brief Return the curves degree (line = 1 etc.)
     *
     * @return int the curves degree
     */
    inline int degree() const
    {
        return m_ctrlpts.size() - 1;
    }

    /**
     * @brief Return the first controlpoint
     *
     * @return Vec2<T> this curves first controlpoint
     */
    Vec2<T> startPt() const
    {
        assert(degree() >= 0);
        return m_ctrlpts[0];
    }

    /**
     * @brief Return the last controlpoint
     *
     * @return Vec2<T> this curves last controlpoint
     */
    Vec2<T> endPt() const
    {
        assert(degree() >= 0);
        return m_ctrlpts[degree()];
    }

    /**
     * @brief Indexed access to controlpoints
     *
     * @param i control point index
     * @return Vec2<T> this curves @p i th control point
     */
    Vec2<T> operator[](uint i) const
    {
        assert(i >= 0 && i < m_ctrlpts.size());
        return m_ctrlpts[i];
    }

    /**
     * @brief Indexed access to controlpoints
     *
     * @param i control point index
     * @return Vec2<T> this curves @p i th control point
     */
    Vec2<T>& operator[](uint i)
    {
        assert(i >= 0 && i < m_ctrlpts.size());
        return m_ctrlpts[i];
    }

    /**
     * @brief Return a bounding box that exactly wraps all control points
     *
     * @return BBox<T> bounding box
     */
    virtual BBox<T> bbox() const
    {
        BBox<T> box(m_ctrlpts[0]);
        for (uint i = 1; i < m_ctrlpts.size(); i++)
            box.widen(m_ctrlpts[i]);
        return box;
    }

    /**
     * @brief Evaluate this curve for parameter @p t
     * @pre t must be in [0, 1]
     *
     * @param t parametrization value
     * @return Vec2<T> this curves value at @p t
     */
    virtual Vec2<T> operator()(const T& t) const = 0;

    /**
     * @brief These define the subinterval of [0, 1] this curve is actually defined in.
     *        After cropping, these should always be 0, 1.
     *
     */
    T m_tStart, m_tEnd;
    /**
     * @brief If this curve was derived from another curve by cutting, these will refer
     *        to the part of the original curve this curve represents.
     *        E.g. m_origTStart = 0, m_origTEnd = 0.25 will represent the first quarter
     *        of some original curve.
     *
     */
    T m_origTStart, m_origTEnd;
    /**
     * @brief This curves id (this should be equal to the curve's index in a vector in this
     *        algorithm library!)
     *
     */
    int m_id;
    /**
     * @brief The curve id of the original curve this curve was cut from at some point (if a read
     *        json file contained curve ids, this is it and will not change during algorithm, check
     * this for easier debugging)
     *
     */
    int m_origId;
    /**
     * @brief This curves controlpoints
     *
     */
    vector<Vec2<T>> m_ctrlpts;

    T m_origLen;
};

} // namespace bzmsh

#endif
