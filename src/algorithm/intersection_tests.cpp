#include "algorithm/intersection_tests.hpp"

#include "algorithm/curve_splitting.hpp"
#include "algorithm/triangulating.hpp"
#include "datastructure/intersection_tree.hpp"
#include "datastructure/line_segment.hpp"
#include "fileio/eps_writer.hpp"
#include "util/logger.hpp"

#include "algorithm/timer.hpp"

#include <list>
#include <map>
#include <set>
#include <tuple>

namespace bzmsh
{
using std::list;
using std::map;
using std::set;
using std::tuple;

bool inConvexPoly(const vector<LineSegment<mpq_class>>& poly,
                  const Vec2<mpq_class>& pt,
                  bool clockWise,
                  bool includeBorder)
{
    assert(poly.size() >= 3);
    int outside = clockWise ? ORI_LEFT : ORI_RIGHT;
    for (const LineSegment<mpq_class>& l : poly)
    {
        int ori = orientation(l.startPt(), l.endPt(), pt);
        if ((ori == ORI_COLINEAR && !includeBorder) || ori == outside)
            return false;
    }
    return true;
}

bool halfOpenSegmentsIntersect(const LineSegment<mpq_class>& seg1,
                       const LineSegment<mpq_class>& seg2,
                       bool start1open,
                       bool start2open)
{
    Vec2<mpq_class> seg1Pt1 = seg1.startPt();
    Vec2<mpq_class> seg1Pt2 = seg1.endPt();
    Vec2<mpq_class> seg2Pt1 = seg2.startPt();
    Vec2<mpq_class> seg2Pt2 = seg2.endPt();

    int o1 = orientation(seg1Pt1, seg1Pt2, seg2Pt1);
    int o2 = orientation(seg1Pt1, seg1Pt2, seg2Pt2);
    if(o1*o2 > 0) return false; // early rejection; one segment lies strictly on one side of the other segment
    if(!start2open) std::swap(o1, o2);
    int o3 = orientation(seg2Pt1, seg2Pt2, seg1Pt1);
    int o4 = orientation(seg2Pt1, seg2Pt2, seg1Pt2);
    if(o3*o4 > 0) return false; // early rejection; one segment lies strictly on one side of the other segment
    if(!start1open) std::swap(o3, o4);
    // orientation is such that (o1=0 AND o3=0) indicates coincident "open ends"
    
    if(o1 == 0 && o3 == 0 && o2 != 0) { assert(o4 != 0); return false; } // only intersecting at open ends
    
    if(o1 != o2 && o3 != o4) return true; // unique intersection point (interior or endpoint)
    
    // remains to check overlaps
    if(o1 == 0 && o2 == 0)
    {
      mpq_class p1,p2,p3,p4;
      mpq_class xdiff = seg1Pt2.x-seg1Pt1.x;
      if(xdiff < 0) xdiff = -xdiff;
      mpq_class ydiff = seg1Pt2.y-seg1Pt1.y;
      if(ydiff < 0) ydiff = -ydiff;
      bool use_x = (ydiff < xdiff);
      if(use_x)
      {
        p1 = seg2Pt1.x - seg1Pt1.x;
        p2 = seg2Pt1.x - seg1Pt2.x;
        p3 = seg2Pt2.x - seg1Pt1.x;
        p4 = seg2Pt2.x - seg1Pt2.x;
      }
      else
      {
        p1 = seg2Pt1.y - seg1Pt1.y;
        p2 = seg2Pt1.y - seg1Pt2.y;
        p3 = seg2Pt2.y - seg1Pt1.y;
        p4 = seg2Pt2.y - seg1Pt2.y;
      }
        
      if(start1open && start2open)
        if(p1 == 0 && p2*p3 > 0) return false; // only overlapping at open ends
      if(!start1open && !start2open)
        if(p4 == 0 && p2*p3 > 0) return false; // only overlapping at open ends
      if(start1open && !start2open)
        if(p3 == 0 && p1*p4 > 0) return false; // only overlapping at open ends
      if(!start1open && start2open)
        if(p2 == 0 && p1*p4 > 0) return false; // only overlapping at open ends
      
      if(p1*p2 <= 0) return true; //seg2Pt1 on seg1
      if(p3*p4 <= 0) return true; //seg2Pt2 on seg1
      if(p1*p4 <= 0) return true; //seg2 contained in or including seg1
    }

    return false;
}

bool guardingQuadsIntersect(const vector<LineSegment<mpq_class>>& poly1,
                          const vector<LineSegment<mpq_class>>& poly2)
{
  assert(poly1.size() == 4 && poly2.size() == 4);
  
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
    {
      if(halfOpenSegmentsIntersect(poly1[i], poly2[j], (i%2==0), (j%2==0)))
        return true;
    }
  
  // q1 is enclosed within q2
  if (inConvexPoly(poly1, poly2[1].startPt(), false, false))
      return true;

  // q2 is enclosed within q1
  if (inConvexPoly(poly2, poly1[1].startPt(), false, false))
      return true;
      
  return false;
}

bool addNewIntersection(const BezierCurve<mpq_class>& segmentI,
                        const BezierCurve<mpq_class>& segmentJ,
                        map<uint, SetOfSplits>& curveIdToSplits)
{
    // This accounts for common endpoints of curve segments that were split beforehand
    if (segmentI.m_id == segmentJ.m_id
        && (segmentI.m_origTEnd == segmentJ.m_origTStart
            || segmentJ.m_origTEnd == segmentI.m_origTStart))
    {
        return false;
    }

    mpq_class tI = (segmentI.m_origTStart + segmentI.m_origTEnd) / mpq_class(2);
    mpq_class tJ = (segmentJ.m_origTStart + segmentJ.m_origTEnd) / mpq_class(2);

    bool isFirstpointI = segmentI.m_origTStart == mpq_class(0);
    bool isLastpointI = segmentI.m_origTEnd == mpq_class(1);
    bool isFirstpointJ = segmentJ.m_origTStart == mpq_class(0);
    bool isLastpointJ = segmentJ.m_origTEnd == mpq_class(1);
    bool isEndpointI = isFirstpointI || isLastpointI;
    bool isEndpointJ = isFirstpointJ || isLastpointJ;
    // This prevents intersections too close to an endpoint
    // -> register the intersection at the endpoint instead
    // (this should also handle tiny curve overlap)
    if ((!isFirstpointI && eqWithEpsilon(tI, mpq_class(0), INTERSECTION_T_EPSILON))
        || (!isLastpointI && eqWithEpsilon(tI, mpq_class(1), INTERSECTION_T_EPSILON))
        || (!isFirstpointJ && eqWithEpsilon(tJ, mpq_class(0), INTERSECTION_T_EPSILON))
        || (!isLastpointJ && eqWithEpsilon(tJ, mpq_class(1), INTERSECTION_T_EPSILON)))
        return false;

    Vec2<mpq_class> endpointI, endpointJ;
    if (isFirstpointI)
    {
        tI = mpq_class(0);
        endpointI = segmentI.startPt();
    }
    else if (isLastpointI)
    {
        tI = mpq_class(1);
        endpointI = segmentI.endPt();
    }
    if (isFirstpointJ)
    {
        tJ = mpq_class(0);
        endpointJ = segmentJ.startPt();
    }
    else if (isLastpointJ)
    {
        tJ = mpq_class(1);
        endpointJ = segmentJ.endPt();
    }
    // Prevent common endpoints being registered as intersection
    if (isEndpointI && isEndpointJ && endpointI == endpointJ)
    {
        return false;
    }

    uint idI = segmentI.m_id;
    uint idJ = segmentJ.m_id;
    bool splitExistsI = curveIdToSplits.find(idI) != curveIdToSplits.end()
                        && curveIdToSplits[idI].find(tI) != curveIdToSplits[idI].end();
    bool splitExistsJ = curveIdToSplits.find(idJ) != curveIdToSplits.end()
                        && curveIdToSplits[idJ].find(tJ) != curveIdToSplits[idJ].end();
    // Prevent new intersections getting registered too close to already registered ones
    if (splitExistsI && splitExistsJ)
    {
        return false;
    }

    // Snap the intersection to a common point
    Vec2<mpq_class> commonPoint;
    if (splitExistsI)
    {
        commonPoint = curveIdToSplits[idI][tI];
    }
    else if (splitExistsJ)
    {
        commonPoint = curveIdToSplits[idJ][tJ];
    }
    else if (isEndpointI)
    {
        commonPoint = endpointI;
    }
    else if (isEndpointJ)
    {
        commonPoint = endpointJ;
    }
    else
    {
        commonPoint = (segmentI(mpq_class(1, 2)) + segmentJ(mpq_class(1, 2))) / mpq_class(2);
    }

    curveIdToSplits[idI][tI] = commonPoint;
    curveIdToSplits[idJ][tJ] = commonPoint;

    return true;
}

bool getIntersections(const vector<BezierCurve<mpq_class>> curves,
                      map<uint, SetOfSplits>& curveIdToSplits)
{
    curveIdToSplits.clear();

    const uint maxDepth = MAX_INTERSECTION_DEPTH;
    const uint cacheDepth = INTERSECTION_CACHE_DEPTH;

    map<uint, pair<BezierCurve<mpq_class>, uint>> idToSegmentAndDepth;
    uint cacheCounter = curves.size();

    IntersectionTree tree(curves.size() * 2);

    // Fill tree with curvesegments one by one and query intersections
    // with segments already in the tree.
    // On POSSIBLE intersection, bisect recursively until no more intersection gets reported
    // or we have an intersection for sure, which we will register.
    for (uint i = 0; i < curves.size(); i++)
    {
        tree.insert(curves[i].bbox(), i);
        idToSegmentAndDepth[i] = std::make_pair(curves[i], 0);
        bool cachedI = false;
        // Ask for possible intersections based on bbox
        for (uint j : tree.query(i))
        {
            bool cachedJ = false;
            int intersectionCount = 0;
            uint depthJ = idToSegmentAndDepth[j].second;
            uint depth = 0;
            BezierCurve<mpq_class> segmentI = curves[i];
            BezierCurve<mpq_class> segmentJ = idToSegmentAndDepth[j].first;
            list<tuple<uint, BezierCurve<mpq_class>, BezierCurve<mpq_class>>> intersectionStack;
            intersectionStack.emplace_front(std::make_tuple(depth, segmentI, segmentJ));
            while (!intersectionStack.empty())
            {
                depth = std::get<0>(intersectionStack.front());
                segmentI = std::get<1>(intersectionStack.front());
                segmentJ = std::get<2>(intersectionStack.front());
                intersectionStack.pop_front();
                if (depth < maxDepth
                    && (segmentI.bbox().circumference() + segmentJ.bbox().circumference())
                           > MIN_INTERSECTION_CIRCUMFERENCE)
                {
                    // Preliminary checks whether intersection would be redundant within epsilon
                    // margin
                    mpq_class tI = (segmentI.m_origTStart + segmentI.m_origTEnd) / mpq_class(2);
                    mpq_class tJ = (segmentJ.m_origTStart + segmentJ.m_origTEnd) / mpq_class(2);
                    uint idI = segmentI.m_id;
                    uint idJ = segmentJ.m_id;
                    bool isEndpointI = eqWithEpsilon(tI, mpq_class(0), INTERSECTION_T_EPSILON)
                                       || eqWithEpsilon(tI, mpq_class(1), INTERSECTION_T_EPSILON);
                    bool isEndpointJ = eqWithEpsilon(tJ, mpq_class(0), INTERSECTION_T_EPSILON)
                                       || eqWithEpsilon(tJ, mpq_class(1), INTERSECTION_T_EPSILON);

                    if (isEndpointI && isEndpointJ
                        && (segmentI.startPt() == segmentJ.startPt()
                            || segmentI.endPt() == segmentJ.startPt()
                            || segmentI.startPt() == segmentJ.endPt()
                            || segmentI.endPt() == segmentJ.startPt()))
                    {
                        continue;
                    }

                    bool splitExistsI
                        = curveIdToSplits.find(idI) != curveIdToSplits.end()
                          && curveIdToSplits[idI].find(tI) != curveIdToSplits[idI].end();
                    bool splitExistsJ
                        = curveIdToSplits.find(idJ) != curveIdToSplits.end()
                          && curveIdToSplits[idJ].find(tJ) != curveIdToSplits[idJ].end();
                    if (splitExistsI && splitExistsJ)
                    {
                        continue;
                    }

                    // Bisect
                    vector<BezierCurve<mpq_class>> bisectedSegmentI, bisectedSegmentJ;
                    bisect(segmentI, mpq_class(1, 2), bisectedSegmentI);
                    bisect(segmentJ, mpq_class(1, 2), bisectedSegmentJ);

                    // Replace curves by 2 bisected parts if caching is enabled
                    if (!cachedI && cacheDepth != 0 && !isEndpointI)
                    {
                        cachedI = true;
                        tree.remove(i);
                        idToSegmentAndDepth.erase(i);
                        idToSegmentAndDepth[tree.insert(bisectedSegmentI[0].bbox(), cacheCounter++)]
                            = std::make_pair(bisectedSegmentI[0], 1);
                        idToSegmentAndDepth[tree.insert(bisectedSegmentI[1].bbox(), cacheCounter++)]
                            = std::make_pair(bisectedSegmentI[1], 1);
                    }

                    if (depthJ < cacheDepth && !cachedJ && !isEndpointJ)
                    {
                        cachedJ = true;
                        tree.remove(j);
                        idToSegmentAndDepth.erase(j);
                        idToSegmentAndDepth[tree.insert(bisectedSegmentJ[0].bbox(), cacheCounter++)]
                            = std::make_pair(bisectedSegmentJ[0], depthJ + 1);
                        idToSegmentAndDepth[tree.insert(bisectedSegmentJ[1].bbox(), cacheCounter++)]
                            = std::make_pair(bisectedSegmentJ[1], depthJ + 1);
                    }

                    // Check bisected parts for intersections
                    BBox<mpq_class> boxI[2]
                        = {bisectedSegmentI[0].bbox(), bisectedSegmentI[1].bbox()};
                    BBox<mpq_class> boxJ[2]
                        = {bisectedSegmentJ[0].bbox(), bisectedSegmentJ[1].bbox()};

                    for (uint idI : {0, 1})
                    {
                        for (uint idJ : {0, 1})
                        {
                            if (boxI[idI].intersectsWith(boxJ[idJ]))
                            {
                                intersectionStack.emplace_front(std::make_tuple(
                                    depth + 1, bisectedSegmentI[idI], bisectedSegmentJ[idJ]));
                            }
                        }
                    }
                }
                // register intersection only if not redundant
                else if (addNewIntersection(segmentI, segmentJ, curveIdToSplits))
                {
                    // Count the intersections for a given segment pair, if higher than possible
                    // for two non-overlapping Bezier curves, report overlap
                    // (actually 2 * degree, but to be safe we choose 3 * degree; intersections for
                    // overlapping curves will be infinitely many anyways)
                    intersectionCount++;
                    if (intersectionCount > 3 * segmentI.degree())
                    {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

} // namespace bzmsh
