#include "util/helper_functions.hpp"

#include "algorithm/curve_splitting.hpp"
#include "datastructure/triangulation_constraints.hpp"
#include "datastructure/vec2.hpp"
#include "util/bezier_injectivity_checker.hpp"
#include "util/logger.hpp"
#include <algorithm>
#include <cassert>
#include <list>
#include <map>
#include <set>

namespace bzmsh
{

double fac(int n)
{
    if (n == 0)
        return 1.0;
    if (n == 1)
        return 1.0;
    if (n == 2)
        return 2.0;
    if (n == 3)
        return 6.0;
    if (n == 4)
        return 24.0;
    if (n == 5)
        return 120.0;
    return n * fac(n - 1);
}

double toRadian(double t)
{
    return M_PI * t / 180.0;
}

double toDegree(double t)
{
    return 180.0 * t / M_PI;
}

template <typename T>
Vec2<T> interpolate(const std::vector<Vec2<T>>& C,
                         std::vector<int>& indC,
                         int i,
                         int j,
                         int n1,
                         int n2,
                         int n3)
{
    if (n1 > 0)
        return (interpolate(C, indC, i, j, n1 - 1, n2, n3)
                + interpolate(C, indC, i - 1, j + 1, n1 - 1, n2, n3))
               * T(0.5);
    if (n2 > 0)
        return (interpolate(C, indC, i, j, n1, n2 - 1, n3)
                + interpolate(C, indC, i, j - 1, n1, n2 - 1, n3))
               * T(0.5);
    if (n3 > 0)
        return (interpolate(C, indC, i, j, n1, n2, n3 - 1)
                + interpolate(C, indC, i + 1, j, n1, n2, n3 - 1))
               * T(0.5);
    return C[indC[i] + j];
}

template Vec2<double> interpolate<double>(const std::vector<Vec2<double>>& C,
                         std::vector<int>& indC,
                         int i,
                         int j,
                         int n1,
                         int n2,
                         int n3);
template Vec2<mpq_class> interpolate<mpq_class>(const std::vector<Vec2<mpq_class>>& C,
                        std::vector<int>& indC,
                        int i,
                        int j,
                        int n1,
                        int n2,
                        int n3);

template <typename T>
bool checkMesh(BezierMesh<T>& mesh, bool& uncertain)
{
    assert(mesh.degree >= 1);
    int degree = mesh.degree;
    int badTriangles = 0;
    int uncertainTriangles = 0;
    
    mesh.ctrlpointFlipped.resize(mesh.allCtrlpts.size());

    // Check whether all the triangles of the mesh are injective
    for (const Triangle& t : mesh.allTriangles)
    {
        assert(t.edges.size() == 3);
        assert(t.vertices.size() == 3);
        assert((int)t.ctrlpts.size() == (degree + 1) * (degree + 2) / 2 - 3 * degree);
        vector<Vec2<T>> ctrlptsRecursive;
        for (int vindex : t.vertices)
        {
            ctrlptsRecursive.emplace_back(mesh.allVertices[vindex]);
        }
        assert(orientation(ctrlptsRecursive[0], ctrlptsRecursive[1], ctrlptsRecursive[2])
               == ORI_LEFT);
        for (int i = 0; i < 3; i++)
        {
            int eindex = t.edges[i];
            int vindexfrom = t.vertices[i];
            vector<int> edgeCtrlptIndices = mesh.allEdges[eindex].ctrlpts;
            if (vindexfrom != mesh.allEdges[eindex].from)
            {
                std::reverse(edgeCtrlptIndices.begin(), edgeCtrlptIndices.end());
            }
            for (int edgeCtrlptIndex : edgeCtrlptIndices)
            {
                ctrlptsRecursive.emplace_back(mesh.allCtrlpts[edgeCtrlptIndex]);
            }
        }
        for (int innerCtrlptIndex : t.ctrlpts)
        {
            ctrlptsRecursive.emplace_back(mesh.allCtrlpts[innerCtrlptIndex]);
        }

        assert((int)ctrlptsRecursive.size() == (degree + 1) * (degree + 2) / 2);
        vector<Vec2<T>> ctrlptsRowwise;
        reorderToRowwise(ctrlptsRecursive, ctrlptsRowwise, degree, true);

        BezierTriangleInjectivityChecker<T> ck(degree);
        std::vector<T> xx, yy;
        for (const Vec2<T>& ctrlpt : ctrlptsRowwise)
        {
            xx.emplace_back(ctrlpt.x);
            yy.emplace_back(ctrlpt.y);
        }
        int result = ck.is_definitely_injective_by_determinant(xx, yy, true);
        if(result == 0) badTriangles++;
        if(result == 2) uncertainTriangles++;
        
        
        for(int ctrlPoint : t.ctrlpts)
        {
           mesh.ctrlpointFlipped[ctrlPoint] = (result == 0);
        }
    }

    if (badTriangles > 0)
    {
        Logger::lout(Logger::INFO) << badTriangles << " out of " << mesh.allTriangles.size()
        << " triangles were IRREGULAR" << endl;
        return false;
    }
    if (uncertainTriangles > 0)
    {
        uncertain = true;
        Logger::lout(Logger::INFO) << "Regularity of " << uncertainTriangles << " out of " << mesh.allTriangles.size()
                                 << " triangles could not be certified." << endl;
    }

    return true;
}

template bool checkMesh<double>(BezierMesh<double>& mesh, bool& uncertain);
template bool checkMesh<mpq_class>(BezierMesh<mpq_class>& mesh, bool& uncertain);

} // namespace bzmsh
