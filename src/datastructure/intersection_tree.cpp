#include "datastructure/intersection_tree.hpp"

namespace bzmsh
{
IntersectionTree::IntersectionTree(uint elements) : m_maxId(0), m_tree(2, 0.0, elements, true)
{
    assert(elements > 0);
}

void IntersectionTree::remove(uint id)
{
    m_tree.removeParticle(id);
}

vector<uint> IntersectionTree::query(uint id)
{
    return m_tree.query(id);
}

} // namespace bzmsh
