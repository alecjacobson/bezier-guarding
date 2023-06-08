#ifndef BZMSH_LINE_CPP
#define BZMSH_LINE_CPP

#include "datastructure/curve.hpp"

namespace bzmsh
{

/**
 * @brief Represents a line segment, on the basis of linear parametrization
 *
 * @tparam T float-type (used for double and mpq_class)
 */
template <typename T> struct LineSegment : public Curve<T>
{
    using Curve<T>::m_ctrlpts;
    using Curve<T>::m_tStart;
    using Curve<T>::m_tEnd;
    using Curve<T>::m_id;

    /**
     * @brief Construct a new bezier curve
     *
     * @param ctrlpts The bezier curves control points
     * @param tStart Starting t
     * @param tEnd Ending t
     * @param id curve id
     */
    LineSegment(const vector<Vec2<T>>& ctrlpts, T tStart = 0, T tEnd = 1, int id = 0)
        : Curve<T>(ctrlpts, tStart, tEnd, id)
    {
        assert(ctrlpts.size() > 1);
        if (ctrlpts.size() > 2)
        {
            Logger::lout(Logger::WARN)
                << "LineSegment with more than 2 control points created, only first and last "
                   "will be used!"
                << endl;
            m_ctrlpts = {ctrlpts[0], ctrlpts[ctrlpts.size() - 1]};
        }
    }

    /**
     * @brief Constructs a new line segment from its two endpoints
     *
     * @param startPt
     * @param endPt
     * @param id
     */
    LineSegment(Vec2<T> startPt, Vec2<T> endPt, int id = 0)
        : Curve<T>(vector<Vec2<T>>({startPt, endPt}), 0, 1, id)
    {
    }

    LineSegment() : Curve<T>()
    {
    }

    /**
     * @brief Return the value of this curve for parameter @p t
     * @pre t must be in [0, 1]
     *
     * @param t parametrization value
     * @return Vec2<T> this curves value at @p t
     */
    virtual Vec2<T> operator()(const T& t) const
    {
        assert(t >= 0 && t <= 1);
        if (t == T(0))
            return m_ctrlpts[0];
        if (t == T(1))
            return m_ctrlpts[this->degree()];
        return this->endPt() * t + this->startPt() * (T(1) - t);
    }
};

} // namespace bzmsh

#endif
