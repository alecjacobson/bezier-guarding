#include "datastructure/triangulation_constraints.hpp"

#include "algorithm/constants.hpp"
#include "datastructure/bbox.hpp"
#include "datastructure/bezier_mesh.hpp"
#include "datastructure/curve.hpp"
#include "datastructure/guarded_bezier_curve.hpp"
#include "util/helper_functions.hpp"
#include <map>

namespace bzmsh
{
using std::map;

template <typename T>
uint TriangulationConstraints<T>::addPoint(const Vec2<T>& p)
{
    if (m_point2index.find(p) == m_point2index.end())
    {
        m_point2index[p] = m_seqPoints.size() / 2;
        m_seqPoints.emplace_back(p.x);
        m_seqPoints.emplace_back(p.y);
    }
    return m_point2index[p];
}

template uint TriangulationConstraints<double>::addPoint(const Vec2<double>& p);
template uint TriangulationConstraints<mpq_class>::addPoint(const Vec2<mpq_class>& p);

template <typename T>
bool TriangulationConstraints<T>::addSegment(const Curve<T>& c, int markCurve)
{
    if (c.degree() < 1)
    {
        Logger::lout(Logger::ERROR) << "Curve of degree <1 can't form a segment for CDT" << endl;
        return false;
    }
    if (c.startPt() == c.endPt())
    {
        Logger::lout(Logger::ERROR)
            << "Curve with identical start and endpoint can't form a segment for CDT" << endl;
        return false;
    }

    m_seqSegments.emplace_back(addPoint(c.startPt()));
    m_seqSegments.emplace_back(addPoint(c.endPt()));
    m_seqSegmentMarkers.emplace_back(markCurve + MARK_OFFSET);
    return true;
}


template bool TriangulationConstraints<double>::addSegment(const Curve<double>& c, int markCurve);
template bool TriangulationConstraints<mpq_class>::addSegment(const Curve<mpq_class>& c, int markCurve);

template <typename T>
bool TriangulationConstraints<T>::addQuad(const GuardedBezierCurve<T>& gc)
{
    for (const Vec2<T>& pt : gc.quadPoints())
        m_curveIdToVertices[gc.m_curve.m_id].emplace_back(addPoint(pt));

    for (const LineSegment<T>& l : gc.quadLines())
        if (!addSegment(l, Edge::MARK_DEFAULT))
            return false;

    return addSegment(gc.m_curve, gc.m_curve.m_id);
}

template bool TriangulationConstraints<double>::addQuad(const GuardedBezierCurve<double>& gc);
template bool TriangulationConstraints<mpq_class>::addQuad(const GuardedBezierCurve<mpq_class>& gc);

template <typename T>
bool TriangulationConstraints<T>::addBbox()
{
    if (m_point2index.empty())
        return false;

    BBox<T> bbox((*m_point2index.begin()).first);
    for (auto kv : m_point2index)
    {
        bbox.widen(kv.first);
    }

    T xrange = bbox.maxX - bbox.minX;
    T yrange = bbox.maxY - bbox.minY;
    bbox.minX -= CDT_BOX_PADDING * xrange;
    bbox.maxX += CDT_BOX_PADDING * xrange;
    bbox.minY -= CDT_BOX_PADDING * yrange;
    bbox.maxY += CDT_BOX_PADDING * yrange;

    addSegment(
        LineSegment<T>(Vec2<T>(bbox.minX, bbox.minY), Vec2<T>(bbox.maxX, bbox.minY)),
        Edge::MARK_BBOX);
    addSegment(
        LineSegment<T>(Vec2<T>(bbox.maxX, bbox.minY), Vec2<T>(bbox.maxX, bbox.maxY)),
        Edge::MARK_BBOX);
    addSegment(
        LineSegment<T>(Vec2<T>(bbox.maxX, bbox.maxY), Vec2<T>(bbox.minX, bbox.maxY)),
        Edge::MARK_BBOX);
    addSegment(
        LineSegment<T>(Vec2<T>(bbox.minX, bbox.maxY), Vec2<T>(bbox.minX, bbox.minY)),
        Edge::MARK_BBOX);

    return true;
}

template bool TriangulationConstraints<double>::addBbox();
template bool TriangulationConstraints<mpq_class>::addBbox();

} // namespace bzmsh
