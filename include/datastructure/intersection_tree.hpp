#ifndef BZMSH_INTERSECTION_TREE_HPP
#define BZMSH_INTERSECTION_TREE_HPP

#include "aabbcc/AABB.h"
#include "datastructure/bbox.hpp"
//#include <gmpxx.h>

namespace bzmsh
{
/**
 * @brief Represents a tree for quick preliminary intersection checks based on axis aligned bounding
 * boxes. Boxes only touching each other at their border will be handled as intersecting.
 *
 */
class IntersectionTree
{
  public:
    /**
     * @brief Construct a tree with an expected number of elements
     * @pre @p elements > 0
     *
     * @param elements expected number of elements
     */
    IntersectionTree(uint elements);

    /**
     * @brief Insert a new bounding box into the tree
     *
     * @param bbox bounding box
     * @param id id to assign to the newly added bounding box
     * @return uint @p id assigned id for later querying and removing
     */
    template <typename T> uint insert(const BBox<T>& bbox, uint id)
    {
        vector<double> lower({utilCast<T, double>(bbox.minX), utilCast<T, double>(bbox.minY)});
        vector<double> upper({utilCast<T, double>(bbox.maxX), utilCast<T, double>(bbox.maxY)});
        m_tree.insertParticle(id, lower, upper);
        m_maxId = std::max(id + 1, m_maxId + 1);
        return id;
    }

    /**
     * @brief Insert a new bounding box into the tree
     *
     * @param bbox bounding box
     * @return uint assigned id for later querying and removing
     */
    template <typename T> uint insert(const BBox<T>& bbox)
    {
        uint id = m_maxId;
        return insert(bbox, id);
    }

    /**
     * @brief Remove bounding box with id == @p id
     *
     * @param id bounding box id
     */
    void remove(uint id);

    /**
     * @brief Query the ids of bounding boxes intersecting with a given bounding box
     *
     * @param id id of bounding box to be queried
     * @return vector<uint> ids of bounding boxes the given one intersects with
     */
    vector<uint> query(uint id);

  private:
    uint m_maxId;
    aabb::Tree m_tree;
};

} // namespace bzmsh
#endif
