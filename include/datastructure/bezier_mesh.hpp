#ifndef BZMSH_BEZIERMESH_HPP
#define BZMSH_BEZIERMESH_HPP

#include "datastructure/bbox.hpp"
#include "datastructure/curve.hpp"
#include "datastructure/guarded_bezier_curve.hpp"
#include "datastructure/triangulation_constraints.hpp"
#include "util/helper_functions.hpp"
#include <map>

namespace bzmsh
{
using std::map;

/**
 * @brief Undirected edge between 2 vertices. Can be marked.
 */
struct Edge
{
    /**
     * @brief Mark for bbox edge
     */
    const static int MARK_BBOX = -5;
    /**
     * @brief Mark for edge from subdivision (for eps visualization)
     */
    const static int MARK_SUBDIV = -4;
    /**
     * @brief Mark for center triangle edges in subdivision (for eps visualization)
     */
    const static int MARK_SUBDIV_CENTER = -3;
    /**
     * @brief Default marker. Standard edge. Non-boundary. Non-curved.
     */
    const static int MARK_DEFAULT = -2;
    /**
     * @brief Boundary edge. Non-curved.
     */
    const static int MARK_BOUNDARY = -1;
    /**
     * @brief IDs >= this will refer to the curve that this edge represents.
     */
    const static int MARK_CURVEIDLIMIT = 0;

    /**
     * @brief Construct a new edge
     *
     * @param _from vertex index 1
     * @param _to vertex inex 2
     * @param _marker marker to assign
     */
    Edge(int _from, int _to, int _marker = MARK_DEFAULT) : from(_from), to(_to), marker(_marker)
    {
        if (from == to)
        {
            Logger::lout(Logger::ERROR) << "invalid edge created" << endl;
            from = -1;
            to = -1;
        }
        if (from > to)
        {
            std::swap(from, to);
        }
    }

    int from = -1;
    int to = -1;
    int marker = MARK_DEFAULT;

    /**
     * @brief Edge ctrlpt indices. Do not include the endpoints!
     * Therefore number of controlpoints = degree -1
     */
    vector<int> ctrlpts;

    /**
     * @brief Comparator for ordering edges based on first vertex
     */
    struct compare
    {
        bool operator()(const Edge& lhs, const Edge& rhs) const
        {
            return (lhs.from < rhs.from) || ((lhs.from == rhs.from) && (lhs.to < rhs.to));
        }
    };
};

struct Triangle
{
    vector<int> vertices; // counterclockwise
    vector<int> edges;    // counterclockwise
    vector<int> ctrlpts;  // only inner controlpoints!

    std::vector<int> allctrlpts;  // All controlpoints (includes offset for non vertex)
};

struct Jacobian
{
    double a;
    double b;
    double c;
    double d;
};

template<typename T>
struct BezierMesh
{
    BezierMesh(int _degree = -1);

    /**
     * @brief Add a vertex into the mesh
     *
     * @param v position vector
     * @return int index of vertex in allVertices vector
     */
    int addVertex(const Vec2<T>& v);

    /**
     * @brief Add a control point into the mesh
     *
     * @param v position vector
     * @return int index of control point in allCtrlpts vector
     */
    int addCtrlpt(const Vec2<T>& v);

    /**
     * @brief Add an edge into the mesh
     *
     * @param from index of first vertex
     * @param to index of second vertex
     * @param marker marker to assign to the edge
     * @return int index of edge in allEdges vector
     */
    int addEdge(int from, int to, int marker);

    /**
     * @brief Add an edge into the mesh
     *
     * @param vfrom first vertex
     * @param vto second vertex
     * @param marker marker to assign to the edge
     * @return int index of edge in allEdges vector
     */
    int addEdge(const Vec2<T>& vfrom, const Vec2<T>& vto, int marker);

    /**
     * @brief Add a triangle into the mesh
     *
     * @param vertexIndices the three vertex indices in counterclockwise order
     * @return int index of the triangle in allTriangles vector
     */
    int addTriangle(const vector<int>& vertexIndices);

    /**
     * @brief Add a triangle into the mesh
     *
     * @param vertices the three vertices in ccw order
     * @return int index of the triangle in allTriangles vector
     */
    int addTriangle(const vector<Vec2<T>>& vertices);

    int triangle_count() { return  allTriangles.size(); };

    int degree;

    map<struct Vec2<T>, int, struct Vec2<T>::compareLex> vertex2index;
    map<struct Vec2<T>, int, struct Vec2<T>::compareLex> ctrlpt2index;
    map<Edge, int, Edge::compare> edge2index;
    vector<Vec2<T>> allVertices;
    vector<Vec2<T>> allCtrlpts;
    vector<Edge> allEdges;
    vector<Triangle> allTriangles;
    
    vector<bool> ctrlpointFlipped; //remove


    //Vertex points followed by rest of control points
    vector< Vec2<T> > allPts;

    //Pre-computation: Local triangle neighborhood information
    std::vector< std::vector<Triangle> > neighbourFaces;
    std::vector< std::set<int> > neighbourFaces_fixedV;
    std::vector< std::map<int, int> > localVertexMap;

    //Pre-computation: Bezier coefficients
    std::vector< std::pair<double,double> > uv_sampling;
    std::vector< std::pair<int, int> > C_IJK, C1_IJK;
    std::vector< std::vector<double> > UV_IJK;

    /**
     * @brief Total number of control points:  (degree+1)*(degree+2)/2
     */
    int m_control_point_count;

    /**
     * @brief Control points (wrt global vector allPts) in row-wise ordered manner
     *
     * @param Triangle t
     * @return control points of a triangle in global sense
     */
    std::vector<int> contrlPts( Triangle t);

    /**
     * @brief i, j, k to local id conversion
     */
    int bz_control_point_local_id(int i, int j, int k);

    /**
     * @brief Bezier coefficients
     *
     * @param B^degree_{i, j, k} and (u,v) as barycentric coordinate
     * @return coefficients value in double
     */
    double bz_value(int degree, int i, int j, int k, double u, double v);

    /**
     * @brief Returns previously stored Bezier coefficients to avoid multiple computations over quadrature
     */
    double bz_get_uv_ijk(int ui, int i, int j, int k);

    /**
     * @brief Compute Jacobian of a triangle
     *
     * @param triangle t
     * @return vector of Jacobians at all quadrature points
     */
    std::vector<Jacobian> bz_jacobian(Triangle& t);

    /**
     * @brief Compute Jacobian of a triangle
     *
     * @param vector of all control points and control point id's of a triangle
     * @return vector of Jacobians at all quadrature points
     */
    std::vector<Jacobian> bz_jacobian(std::vector< Vec2<T> >& ctrPt, std::vector<int>& t_ctrl_pt);

    /**
     * @brief Computes gradient of our conformal energy
     *
     * @param control points, jacobian, and local map(control pt id -> local id)
     * @return Grad vector as Vec2 per control point
     */
    void bz_gradient(std::vector<int>& t_ctrl_pt, std::vector<Jacobian>& Jac, std::vector< Vec2<double> >& Grad, std::map<int, int>& map );

    /**
     * @brief Compute Conformal Energy of a triangle
     *
     * @param Jacobian vector over all quadrature points
     * @return Conformal energy
     */
    double bz_energy(std::vector<Jacobian>& Jac);

    /**
     * @brief Compute Conformal Energy of 1-ring triangles around a vertex
     *
     * @param vector of all control points and base vertex id
     * @return Conformal energy
     */
    double bz_energy(vector< Vec2<T> >& ctrPt, int v_id);

    /**
     * @brief Check orientation of all trianges around v_id
     *
     * @param control points of a triangle and the vector of all control points
     * @return True if flipped
     */
    bool check_flip(vector< Vec2<T> >& ctrPt, std::vector<int>& f_ctrPt);

    /**
     * @brief Check orientation of 1-ring trianges around v_id
     *
     * @param id of base vertex and the vector of all control points
     * @return True if flipped
     */
    bool check_flip(vector< Vec2<T> >& ctrPt, int v_id);


    /**
     * @brief Fill local triangulation information for optimization
     */
    void fillNeighbourInfo();

    /**
     * @brief Gradient step for scale invariant conformal distortion measure [Hormann and Greiner 2000]
     */
    void gradient_step( int v_id);

    /**
     * @brief Main function to optimize the mesh via gradient descent in a local manner (1 ring triangles)
     */
    void optimization_conformal();
};

} // namespace bzmsh

#endif
