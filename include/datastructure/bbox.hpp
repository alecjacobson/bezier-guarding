#ifndef BZMSH_BBOX_HPP
#define BZMSH_BBOX_HPP

#include "datastructure/vec2.hpp"
#include <algorithm>

namespace bzmsh
{
using std::max;
using std::min;
/**
 * @brief Axis aligned bounding box (mainly used for preliminary intersection checks)
 *
 * @tparam T float-type (used for double and mpq_class)
 */
template <typename T> struct BBox
{
    T minX;
    T minY;
    T maxX;
    T maxY;

    /**
     * @brief Initialize a bounding box with a 2D point
     *
     * @param pt point
     */
    BBox(const Vec2<T>& pt) : minX(pt.x), minY(pt.y), maxX(pt.x), maxY(pt.y)
    {
    }

    /**
     * @brief Widen the bbox to enclose a new point
     *
     * @param pt point
     */
    void widen(const Vec2<T>& pt)
    {
        minX = min(pt.x, minX);
        minY = min(pt.y, minY);
        maxX = max(pt.x, maxX);
        maxY = max(pt.y, maxY);
    }

    /**
     * @brief Widen the bbox to enclose another bounding box
     *
     * @param bb bounding box
     */
    void widen(const BBox& bb)
    {
        minX = min(bb.minX, minX);
        minY = min(bb.minY, minY);
        maxX = max(bb.maxX, maxX);
        maxY = max(bb.maxY, maxY);
    }

    /**
     * @brief Area of this bounding box
     *
     * @return T area
     */
    T area()
    {
        return (maxX - minX) * (maxY - minY);
    }

    /**
     * @brief circumference of this bounding box
     *
     * @return T circumference
     */
    T circumference()
    {
        return T(2) * ((maxX - minX) + (maxY - minY));
    }

    /**
     * @brief Quick intersection test with another bounding box
     *
     * @param other other bounding box to test intersection with
     * @param includeBorder whether common border point means intersecting
     * @return true if intersecting
     * @return false if not intersecting
     */
    bool intersectsWith(const BBox<T>& other, bool includeBorder = true)
    {
        if (includeBorder)
            return !(other.minX > maxX || other.maxX < minX || other.maxY < minY
                     || other.minY > maxY);
        else
            return !(other.minX >= maxX || other.maxX <= minX || other.maxY <= minY
                     || other.minY >= maxY);
    }
};
} // namespace bzmsh
#endif
