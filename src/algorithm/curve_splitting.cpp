#include "algorithm/curve_splitting.hpp"

#include "algorithm/intersection_tests.hpp"
#include "datastructure/bezier_curve.hpp"
#include "datastructure/intersection_tree.hpp"
#include "datastructure/line_segment.hpp"
#include "fileio/eps_writer.hpp"

#include <list>
#include <map>
#include <set>
#include <tuple>

namespace bzmsh
{
using std::list;

template <typename T>
bool split(const BezierCurve<T>& curve,
           const SetOfSplits& splitSet,
           vector<BezierCurve<T>>& curvesOut)
{
    curvesOut.clear();
    if (splitSet.size() == 0)
    {
        curvesOut.emplace_back(curve);
        return false;
    }

    int d = curve.degree();

    vector<T> splitParams;
    vector<Vec2<T>> splitPoints;
    if (splitSet.find(T(0)) == splitSet.end())
    {
        splitParams.emplace_back(0);
        splitPoints.emplace_back(curve.startPt());
    }
    for (const Split& s : splitSet)
    {
        splitParams.emplace_back(s.first.get_d());
        splitPoints.emplace_back(convertVec<mpq_class, T>(s.second));
    }
    if (splitSet.find(T(1)) == splitSet.end())
    {
        splitParams.emplace_back(1);
        splitPoints.emplace_back(curve.endPt());
    }

    for (uint i = 0; i < splitParams.size() - 1; i++)
    {
        BezierCurve<T> curveSegment = curve;
        curveSegment.m_tStart = splitParams[i];
        curveSegment.m_tEnd = splitParams[i + 1];
        if(crop(curveSegment))
        {
            curveSegment.m_origTStart
                    = splitParams[i] * (curve.m_origTEnd - curve.m_origTStart) + curve.m_origTStart;
            curveSegment.m_origTEnd
                    = splitParams[i + 1] * (curve.m_origTEnd - curve.m_origTStart) + curve.m_origTStart;
            curveSegment[0] = splitPoints[i];
            curveSegment[d] = splitPoints[i + 1];
            curvesOut.emplace_back(curveSegment);
        }
    }

    return true;
}

template bool split<double>(const BezierCurve<double>& curve,
           const SetOfSplits& splitSet,
           vector<BezierCurve<double>>& curvesOut);

template bool split<mpq_class>(const BezierCurve<mpq_class>& curve,
           const SetOfSplits& splitSet,
           vector<BezierCurve<mpq_class>>& curvesOut);
           

template <typename T>
void splitNonWellGuardables(vector<BezierCurve<T>>& curves)
{
    list<BezierCurve<T>> curveSegments;
    curveSegments.insert(curveSegments.begin(), curves.begin(), curves.end());
    curves.clear();

    while (!curveSegments.empty())
    {
        BezierCurve<T> front = curveSegments.front();
        curveSegments.pop_front();

        if (!front.isWellGuardable())
        {
            Logger::lout(Logger::DEBUG) << "Bisecting curve " << front << " with angleRange "
                                          << toDegree(front.angleRange()) << "°" << endl;
                                          
            vector<BezierCurve<T>> bisectedC;
            bisect(front, T(0.5), bisectedC);

            Logger::lout(Logger::DEBUG)
                << "Subcurve 1: " << bisectedC[0] << " with angleRange "
                << toDegree(bisectedC[0].angleRange()) << "°" << endl;
            Logger::lout(Logger::DEBUG)
                << "Subcurve 2: " << bisectedC[1] << " with angleRange "
                << toDegree(bisectedC[1].angleRange()) << "°" << endl;
            curveSegments.emplace_front(bisectedC[0]);
            curveSegments.emplace_front(bisectedC[1]);
        }
        else
        {
            curves.emplace_back(front);
        }
    }
}

template void splitNonWellGuardables<double>(vector<BezierCurve<double>>& curves);
template void splitNonWellGuardables<mpq_class>(vector<BezierCurve<mpq_class>>& curves);

template <typename T>
bool splitIntersecting(vector<GuardedBezierCurve<T>>& guardedCurves)
{
    map<uint, GuardedBezierCurve<mpq_class>> idToSegment;
    list<GuardedBezierCurve<mpq_class>> curveSegments;

    for (const GuardedBezierCurve<T>& c : guardedCurves)
    {
        curveSegments.emplace_back(convertGuardedCurve<T, mpq_class>(c));
    }

    IntersectionTree tree(guardedCurves.size() * 2);

    while (!curveSegments.empty())
    {
        GuardedBezierCurve<mpq_class> segmentI = curveSegments.front();
        curveSegments.pop_front();

        if (segmentI.m_timesCut > MAX_GUARDING_DEPTH)
        {
            Logger::lout(Logger::DEBUG)
                << "Bisection of guarded curves reached set max refinement limit. Aborting." << endl;
            return false;
        }

        uint i = tree.insert(segmentI.bbox());
        idToSegment.emplace(i, segmentI);

        for (uint j : tree.query(i))
        {
            GuardedBezierCurve<mpq_class> segmentJ = idToSegment[j];

            vector<LineSegment<mpq_class>> quadLineI = segmentI.quadLines();
            vector<LineSegment<mpq_class>> quadLineJ = segmentJ.quadLines();

#ifndef ALL_EXACT
            for(auto s : quadLineI)
            {
                s.m_ctrlpts[0] = convertVec<T, mpq_class>( convertVec<mpq_class, T> (s.m_ctrlpts[0]));
                s.m_ctrlpts[1] = convertVec<T, mpq_class>( convertVec<mpq_class, T> (s.m_ctrlpts[1]));
            }
            for(auto s : quadLineJ)
            {
                s.m_ctrlpts[0] = convertVec<T, mpq_class>( convertVec<mpq_class, T> (s.m_ctrlpts[0]));
                s.m_ctrlpts[1] = convertVec<T, mpq_class>( convertVec<mpq_class, T> (s.m_ctrlpts[1]));
            }
#endif
            if(guardingQuadsIntersect(quadLineI, quadLineJ))
            {
                GuardedBezierCurve<mpq_class>* toBisect;
                if (segmentI.m_priority > segmentJ.m_priority)
                {
                    toBisect = &segmentI;
                    tree.remove(i);
                    idToSegment.erase(i);
                }
                else
                {
                    toBisect = &segmentJ;
                    tree.remove(j);
                    idToSegment.erase(j);
                }

                vector<BezierCurve<mpq_class>> bisectedCurve;
                vector<GuardedBezierCurve<mpq_class>> bisectedGuardedCurve;

                if (!bisect(toBisect->m_curve, mpq_class(1, 2), bisectedCurve))
                {
                    Logger::lout(Logger::DEBUG)
                        << "Bisecting " << (*toBisect) << " failed " << endl;
                    return false;
                }
                
                if (!guardAll(bisectedCurve, bisectedGuardedCurve))
                {
                    Logger::lout(Logger::DEBUG)
                        << "Guarding " << (*toBisect) << " failed " << endl;
                    return false;
                }
                for (GuardedBezierCurve<mpq_class>& bisectedCurve : bisectedGuardedCurve)
                {
                    bisectedCurve.m_timesCut = toBisect->m_timesCut + 1;
                }
                curveSegments.insert(curveSegments.begin(),
                                     bisectedGuardedCurve.begin(),
                                     bisectedGuardedCurve.end());

                // If we bisected the segment for which we are querying intersections right now, we
                // need to restart with its subsegments
                if (segmentI.m_priority > segmentJ.m_priority)
                {
                    break;
                }
            }
        }
    }

    guardedCurves.clear();
    for (const auto& kv : idToSegment)
    {
        guardedCurves.emplace_back(convertGuardedCurve<mpq_class, T>(kv.second));
    }
    return true;
}


template bool splitIntersecting<double>(vector<GuardedBezierCurve<double>>& guardedCurves);
template bool splitIntersecting<mpq_class>(vector<GuardedBezierCurve<mpq_class>>& guardedCurves);

} // namespace bzmsh
