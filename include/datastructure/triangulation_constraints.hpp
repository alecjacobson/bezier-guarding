#ifndef BZMSH_TRIANGULATION_CONSTRAINTS_HPP
#define BZMSH_TRIANGULATION_CONSTRAINTS_HPP

#include "datastructure/bezier_mesh.hpp"
#include "datastructure/curve.hpp"
#include "datastructure/guarded_bezier_curve.hpp"
#include <map>

namespace bzmsh
{

/**
 * @brief Constraints to pass into a constrained delaunay triangulation (CDT)
 *
 */
template <typename T> struct TriangulationConstraints
{
    static const int MARK_OFFSET = 2;

    /**
     * @brief Add vertex (called point in cdt library)
     *
     * @param p point
     * @return uint index in sequential m_seqPoints vector
     */
    uint addPoint(const Vec2<T>& p);

    /**
     * @brief Add a segment constraint
     *
     * @param c curve (start and endpoint will form the segment)
     * @param markCurve segment marker
     * @return true if successful
     * @return false if invalid curve
     */
    bool addSegment(const Curve<T>& c, int markCurve);

    /**
     * @brief Add a whole guarding quad (and its segments) as constraints
     *
     * @param gc guarded curve
     * @return true if successful
     * @return false if invalid quad
     */
    bool addQuad(const GuardedBezierCurve<T>& gc);

    /**
     * @brief Add 4 segment constraints enclosing all vertices in an axis-aligned bbox
     *
     * @return true if successfull
     * @return false if no vertices to box in
     */
    bool addBbox();


    /**
     * @brief Helper functions to keep track of point id
     * @param id
     * @return
     */
    Vec2<T> get_point(int id)
    {
        if(id < m_point2index.size())
            return Vec2<T>( m_seqPoints[2*id], m_seqPoints[2*id+1] );
        return Vec2<T>( T(0.0), T(0.0));
    }

    int get_point_index(Vec2<T> p)
    {
        if(m_point2index.find(p) == m_point2index.end())
            return -1;
        return m_point2index[p];
    }
    
    vector<T> m_seqPoints;
    vector<int> m_seqSegments;
    vector<int> m_seqSegmentMarkers;
    std::map<struct Vec2<T>, uint, struct Vec2<T>::compareLex> m_point2index;
    std::map<int, vector<int>> m_curveIdToVertices;
};
} // namespace bzmsh

#endif
