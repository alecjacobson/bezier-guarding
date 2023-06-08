#include "fileio/msh_writer.hpp"

#include <fstream>
#include <iostream>

namespace bzmsh
{
using std::ofstream;
using std::ostream;
using std::string;

template <typename T>
bool MSHWriter<T>::write(const string& filename, BezierMesh<T>& mesh)
{
    assert(mesh.degree >= 1);
    ofstream out(filename);
    if (!out.good() || mesh.allEdges.empty())
    {
        return false;
    }
    // precision for double output
    out.precision(10);

    int triType = 0;
    int degree = mesh.degree;
    std::vector<int> reorder;
    switch (degree)
    {
    case 1:
        triType = 2;
        break;
    case 2:
        triType = 9;
        break;
    case 3:
        triType = 21;
        break;
    case 4:
        triType = 23;
        break;
    case 5:
        triType = 25;
        break;
    case 6:
        triType = 42;
        break;
    case 7:
        triType = 43;
        break;
    case 8:
        triType = 44;
        break;
    case 9:
        triType = 45;
        break;
    case 10:
        triType = 46;
        break;
    default:
        Logger::lout(Logger::ERROR)
            << "No support for .msh output of degree " << degree << " yet" << endl;
        return false;
    }

    /// Header
    out << "$MeshFormat"
        << "\n"
        << "4.1 0 " << std::to_string(sizeof(double)) << "\n"
        << "$EndMeshFormat"
        << "\n"
        << endl;

    /// Nodes
    int nNodes = mesh.allVertices.size() + mesh.allCtrlpts.size();
    int ctrlptOffset = mesh.allVertices.size();
    out << "$Nodes"
        << "\n"
        << indent(1) << "1 " << nNodes << " 1 " << nNodes << "\n"
        << indent(1) << "2 1 0 " << nNodes << "\n";
    for (int id = offset(0, baseOffset); id < offset(nNodes, baseOffset); id++)
        out << indent(1) << id << "\n";
    for (Vec2<T> vT : mesh.allVertices)
    {
        Vec2<double> v = convertVec<T,double>(vT);
        out << indent(1) << v.x << " " << v.y << " " << 0 << "\n";
    }
    for (Vec2<T> ctrlptT : mesh.allCtrlpts)
    {
        Vec2<double> ctrlpt = convertVec<T,double>(ctrlptT);
        out << indent(1) << ctrlpt.x << " " << ctrlpt.y << " " << 0 << "\n";
    }
    out << "$EndNodes"
        << "\n"
        << endl;

    /// Elements
    int nTriangles = mesh.allTriangles.size();
    out << "$Elements"
        << "\n"
        << indent(1) << "1 " << nTriangles << " 1 " << nTriangles << "\n"
        << indent(1) << "2 1 " << triType << " " << nTriangles << "\n";

    for (int i = 0; i < nTriangles; i++)
    {
        const Triangle& t = mesh.allTriangles[i];
        assert(t.edges.size() == 3);
        assert(t.vertices.size() == 3);
        out << indent(1) << offset(i, baseOffset);
        for (int vindex : t.vertices)
        {
            out << " " << offset(vindex, baseOffset);
        }
        for (int j = 0; j < 3; j++)
        {
            vector<int> edgeCtrlptIndices = mesh.allEdges[t.edges[j]].ctrlpts;
            int vfrom = t.vertices[j];
            if (mesh.allEdges[t.edges[j]].from != vfrom)
            {
                std::reverse(edgeCtrlptIndices.begin(), edgeCtrlptIndices.end());
            }
            assert((int)edgeCtrlptIndices.size() == degree - 1);
            for (int edgeCtrlptIndex : edgeCtrlptIndices)
            {
                out << " " << offset(edgeCtrlptIndex, baseOffset + ctrlptOffset);
            }
        }
        assert((int)t.ctrlpts.size() == (degree + 1) * (degree + 2) / 2 - 3 * degree);
        for (int innerCtrlptIndex : t.ctrlpts)
        {
            out << " " << offset(innerCtrlptIndex, baseOffset + ctrlptOffset);
        }

        out << "\n";
    }
    out << "$EndElements"
        << "\n"
        << endl;

    return true;
}

template bool MSHWriter<double>::write(const string& filename, BezierMesh<double>& mesh);
template bool MSHWriter<mpq_class>::write(const string& filename, BezierMesh<mpq_class>& mesh);

} // namespace bzmsh
