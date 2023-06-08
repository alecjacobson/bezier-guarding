#include "datastructure/bezier_mesh.hpp"
#include "util/bezier_injectivity_checker.hpp"
#include "util/helper_functions.hpp"

#define SQRT3_INV 0.57735026919

namespace bzmsh
{

template <typename T>
BezierMesh<T>::BezierMesh(int _degree) : degree(_degree)
{
}

template <typename T>
int BezierMesh<T>::addVertex(const Vec2<T>& v)
{
    if (vertex2index.find(v) == vertex2index.end())
    {
        vertex2index[v] = allVertices.size();
        allVertices.emplace_back(v);
    }
    return vertex2index[v];
}

template <typename T>
int BezierMesh<T>::addCtrlpt(const Vec2<T>& v)
{
    if (ctrlpt2index.find(v) == ctrlpt2index.end())
    {
        ctrlpt2index[v] = allCtrlpts.size();
        allCtrlpts.emplace_back(v);
    }
    return ctrlpt2index[v];
}

template <typename T>
int BezierMesh<T>::addEdge(int from, int to, int marker)
{
    Edge e(from, to, marker);
    if (edge2index.find(e) == edge2index.end())
    {
        edge2index[e] = allEdges.size();
        allEdges.emplace_back(e);
    }
    return edge2index[e];
}

template <typename T>
int BezierMesh<T>::addEdge(const Vec2<T>& vfrom, const Vec2<T>& vto, int marker)
{
    int from = addVertex(vfrom);
    int to = addVertex(vto);
    Edge e(from, to, marker);
    if (edge2index.find(e) == edge2index.end())
    {
        edge2index[e] = allEdges.size();
        allEdges.emplace_back(e);
    }
    return edge2index[e];
}

template <typename T>
int BezierMesh<T>::addTriangle(const vector<int>& vertexIndices)
{
    assert(vertexIndices.size() == 3);
    Triangle tri;
    tri.vertices = vertexIndices;

    for (int i = 0; i < 3; i++)
    {
        tri.edges.emplace_back(addEdge(tri.vertices[i], tri.vertices[(i + 1) % 3], -2));
    }

    allTriangles.emplace_back(tri);

    return allTriangles.size() - 1;
}

template <typename T>
int BezierMesh<T>::addTriangle(const vector<Vec2<T>>& vertices)
{
    assert(vertices.size() == 3);
    vector<int> vertexIndices;
    for (int i = 0; i < 3; i++)
    {
        vertexIndices.emplace_back(addVertex(vertices[i]));
    }
    return addTriangle(vertexIndices);
}


template <typename T>
std::vector<int> BezierMesh<T>::contrlPts( Triangle t)
{
    int offset = allVertices.size();
    vector<int> ctrlpts;
    for (int vindex : t.vertices)
    {
        ctrlpts.emplace_back(vindex);
    }
    for (int i = 0; i < 3; i++)
    {
        int eindex = t.edges[i];
        int vindexfrom = t.vertices[i];
        vector<int> edgeCtrlptIndices = allEdges[eindex].ctrlpts;
        if (vindexfrom != allEdges[eindex].from)
        {
            std::reverse(edgeCtrlptIndices.begin(), edgeCtrlptIndices.end());
        }
        for (int edgeCtrlptIndex : edgeCtrlptIndices)
        {
            ctrlpts.emplace_back(edgeCtrlptIndex + offset);
        }
    }
    for (int innerCtrlptIndex : t.ctrlpts)
    {
        ctrlpts.emplace_back(innerCtrlptIndex + offset);
    }

    assert((int)ctrlpts.size() == (degree + 1) * (degree + 2) / 2);
    vector<int> ctrlptsRowwise;
    reorderToRowwise(ctrlpts, ctrlptsRowwise, degree, true);

    return ctrlptsRowwise;
}

template <typename T>
void BezierMesh<T>::fillNeighbourInfo()
{
    int offset = allVertices.size();
    neighbourFaces.clear();
    neighbourFaces_fixedV.clear();

    for(int i=0; i<vertex2index.size(); i++)
    {
        std::vector< Triangle > a;
        std::set<int> b;

        neighbourFaces.push_back(a);
        neighbourFaces_fixedV.push_back(b);
    }

    //Go over all Triangles and add local information
    for (Triangle& t : allTriangles)
    {
        t.allctrlpts = contrlPts(t);
        //for (int vindex : t.vertices)
        for (int i=0; i<3; i++)
        {
            int v_id = t.vertices[i];
            neighbourFaces[v_id].push_back(t);

            //mark fixed for optimization
            for (int j=0; j<3; j++)
            {
                Edge& e = allEdges[t.edges[j]];
                bool fixed = true;
                if(e.from == v_id || e.to == v_id)
                {
                    fixed = false;
                    if(e.marker == Edge::MARK_BBOX || e.marker >= Edge::MARK_CURVEIDLIMIT)
                        fixed = true;
                }
                if(fixed)
                {
                    neighbourFaces_fixedV[v_id].insert( e.from );
                    neighbourFaces_fixedV[v_id].insert( e.to );
                    for(int k : e.ctrlpts)
                        neighbourFaces_fixedV[v_id].insert( k+offset );
                }
            }
        }
    }

    allPts.clear();
    for( Vec2<T> p: allVertices )
        allPts.push_back(p);
    for( Vec2<T> p: allCtrlpts )
        allPts.push_back(p);

    localVertexMap.clear();
    for(int i=0; i<neighbourFaces.size(); i++)
    {
        int id = 0;
        std::map<int,int> m;
        for (const Triangle& t : neighbourFaces[i])
        {
            for(int c : t.allctrlpts)
            {
                if(m.end() == m.find(c))
                    m[c] = id++;
            }
        }
        localVertexMap.push_back(m);
    }
}

template <typename T>
int BezierMesh<T>::bz_control_point_local_id(int i, int j, int k)
{
    int d = i + j + k;
    return (d+1)*i - ((i*(i-1))/2) + j;
}

template <typename T>
double BezierMesh<T>::bz_value(int degree, int i, int j, int k, double u, double v)
{
    if(degree != i + j + k)
        printf("issue with bezier poly degree: %d = %d + %d + %d\n", degree, i, j, k);
    if(i<0 || j<0 || k<0)
        return 0.0;
    double b = fac(degree);
    b *= pow(u, i);
    b *= pow(v, j);
    b *= pow(1.0-u-v, k);
    b /= fac(i);
    b /= fac(j);
    b /= fac(k);
    return b;
}

template <typename T>
double BezierMesh<T>::bz_get_uv_ijk(int ui, int i, int j, int k)
{
    if(i<0 || j<0 || k<0)
        return 0.0;

    return UV_IJK[ui][ bz_control_point_local_id(i, j, k) ];
}

template <typename T>
std::vector<Jacobian> BezierMesh<T>::bz_jacobian(Triangle& t)
{
    return bz_jacobian(allPts, t.allctrlpts);
}

template <typename T>
std::vector<Jacobian> BezierMesh<T>::bz_jacobian(std::vector< Vec2<T> >& ctrPt, std::vector<int>& t_ctrl_pt)
{
    //       | j0*x_u  j1*x_u + j2*x_v |
    // J  =  |                         |
    //       | j0*y_u  j1*y_u + j2*y_v |
    std::vector<Jacobian> J;
    for(unsigned int ui=0; ui<uv_sampling.size(); ui++)
    {
        double x_u = 0.0, x_v = 0.0;
        double y_u = 0.0, y_v = 0.0;

        for(int cp=0; cp<m_control_point_count; cp++)
        {
            int i = C_IJK[cp].first;
            int j = C_IJK[cp].second;
            int k = degree - i - j;
            Vec2<double> bz = convertVec<T, double>(ctrPt[ t_ctrl_pt[cp] ]);
            double bz_vali = bz_get_uv_ijk(ui, i-1, j, k);
            double bz_valj = bz_get_uv_ijk(ui, i, j-1, k);
            double bz_valk = bz_get_uv_ijk(ui, i, j, k-1);

            x_u += bz.x * (bz_vali-bz_valk);
            y_u += bz.y * (bz_vali-bz_valk);
            x_v += bz.x * (bz_valj-bz_valk);
            y_v += bz.y * (bz_valj-bz_valk);
        }

        Jacobian j_local;
        j_local.a = x_u;
        j_local.b = SQRT3_INV * (-x_u + 2.0*x_v);
        j_local.c = y_u;
        j_local.d = SQRT3_INV * (-y_u + 2.0*y_v);
        J.push_back(j_local);
    }
    return J;
}

template <typename T>
void BezierMesh<T>::bz_gradient(std::vector<int>& t_ctrl_pt, std::vector<Jacobian>& Jac, std::vector< Vec2<double> >& Grad, std::map<int, int>& map )
{
    double wt = 1.0 / uv_sampling.size();
    for(unsigned int ui=0; ui<uv_sampling.size(); ui++)
    {
        double& a = Jac[ui].a;
        double& b = Jac[ui].b;
        double& c = Jac[ui].c;
        double& d = Jac[ui].d;
        double det_Jinv = 1.0 / ((a*d) - (b*c));
        double norm_J = (a*a) + (b*b) + (c*c) + (d*d);

        for(int i=0; i<=degree; i++)
        {
            for(int j=0; i+j<=degree; j++)
            {
                int k = degree - i - j;
                int id = t_ctrl_pt[ bz_control_point_local_id(i, j, k) ];

                double bz_vali = bz_get_uv_ijk(ui, i-1, j, k);
                double bz_valj = bz_get_uv_ijk(ui, i, j-1, k);
                double bz_valk = bz_get_uv_ijk(ui, i, j, k-1);

                double d_a = (bz_vali - bz_valk);
                double d_b = -SQRT3_INV*(bz_vali - bz_valk) + 2.0*SQRT3_INV*(bz_valj - bz_valk);

                Vec2<double> g = Grad[ map[id] ];
                g.x += wt * (det_Jinv*(2*a*d_a + 2*b*d_b - (d*d_a - c*d_b)*norm_J*det_Jinv));
                g.y += wt * (det_Jinv*(2*d*d_b + 2*c*d_a - (a*d_b - b*d_a)*norm_J*det_Jinv));
                Grad[ map[id] ] = g;
            }
        }
    }
}

template <typename T>
double BezierMesh<T>::bz_energy(std::vector<Jacobian>& Jac)
{
    double e = 0.0;
    for(Jacobian j : Jac)
    {
        double& a = j.a;
        double& b = j.b;
        double& c = j.c;
        double& d = j.d;

        double det = (a*d) - (b*c);
        double nume = (a*a) + (b*b) + (c*c) + (d*d);
        e += nume / det;
    }
    e /= uv_sampling.size(); //weight of the sample
    return e;
}

template <typename T>
double BezierMesh<T>::bz_energy(vector< Vec2<T> >& ctrPt, int v_id)
{
    double e = 0.0;
    for(Triangle t : neighbourFaces[v_id])
    {
        std::vector<Jacobian> jac = bz_jacobian(ctrPt, t.allctrlpts);
        e += bz_energy(jac);
    }
    return e;
}

template <typename T>
bool BezierMesh<T>::check_flip(vector< Vec2<T> >& ctrPt, std::vector<int>& f_ctrPt)
{
    BezierTriangleInjectivityChecker<T> ck(degree);
    std::vector<T> xx, yy;
    for (int p : f_ctrPt)
    {
        xx.emplace_back(ctrPt[p].x);
        yy.emplace_back(ctrPt[p].y);
    }
    int result = ck.is_definitely_injective_by_determinant(xx, yy, true);
    if(1 == result)
        return false;
    return true;
}

template <typename T>
bool BezierMesh<T>::check_flip(vector< Vec2<T> >& ctrPt, int v_id)
{
    for(Triangle& t : neighbourFaces[v_id])
    {
        if(check_flip(ctrPt, t.allctrlpts))
            return true;
    }
    return false;
}

template <typename T>
void BezierMesh<T>::gradient_step(int v_id)
{
    if(neighbourFaces[v_id].size() <= 0)
        return;

    int max_itr = 30;
    std::vector< Vec2<double> > Grad;

    for(unsigned int i=0; i<localVertexMap[v_id].size(); i++)
        Grad.push_back(Vec2<double>(0.0, 0.0));

   for(Triangle t : neighbourFaces[v_id])
   {
       std::vector<Jacobian> Jac = bz_jacobian(t);
       bz_gradient(t.allctrlpts, Jac, Grad, localVertexMap[v_id] );
   }

   //Mark fixed points
   for(int i : neighbourFaces_fixedV[v_id])
       Grad[ localVertexMap[v_id][i] ] = Vec2<double>(0.0, 0.0);

   vector< Vec2<T> > ctrPt_temp = allPts;

   double new_energy;
   double old_energy = bz_energy(ctrPt_temp, v_id);

   int itr = 0;
   double m_lambda = 2.0;
   for(int k=0; k<2; k++)
   {
       m_lambda *= 0.5;
       for(std::pair<int, int> m : localVertexMap[v_id])
       {
           ctrPt_temp[m.first].x = allPts[m.first].x - utilCast<double, T>(m_lambda*Grad[m.second].x);
           ctrPt_temp[m.first].y = allPts[m.first].y - utilCast<double, T>(m_lambda*Grad[m.second].y);
       }
       if(0 == k)
       {
           if(check_flip(ctrPt_temp, v_id))
               k = -1;
           else //flips checking done
               k = 1;
       }
       if(1 == k)
       {
           new_energy = bz_energy(ctrPt_temp, v_id);
           if(new_energy>old_energy)
               k = 0;
           else if(check_flip(ctrPt_temp, v_id)) // Unexpected flip with smaller step size
               k = -1;
       }

       itr++;
       if(itr>max_itr)
           break;
   }
   if(itr<=max_itr)
   {
       //printf("lambda(%d) = %0.10f\n", itr, m_lambda);
       //double percentage = 100.0*(old_energy-new_energy) / old_energy;
       //printf("E reduction = %0.10f %%  (%f/%f)\n", percentage, new_energy, old_energy);

       allPts = ctrPt_temp;
       old_energy = new_energy;
   }
}

template <typename T>
void BezierMesh<T>::optimization_conformal()
{
    int optItr = 10;
    m_control_point_count = (degree+1)*(degree+2)/2;

    //Set Quadrature and Local ids

    //Reverse map: id -> i,j,k
    C_IJK.resize(m_control_point_count);
    C1_IJK.resize(degree*(degree+1)/2);
    for(int i=0; i<=degree; i++)
    {
        for(int j=0; i+j<=degree; j++)
        {
            int k = degree - i - j;
            int id = bz_control_point_local_id(i, j, k);
            C_IJK[id] = std::make_pair(i, j);
        }
    }
    for(int i=0; i<degree; i++)
    {
        for(int j=0; i+j<degree; j++)
        {
            int k = degree - 1 - i - j;
            int id = bz_control_point_local_id(i, j, k);
            C1_IJK[id] = std::make_pair(i, j);
        }
    }

    //Quadrature for Energy computation
    int quadrature = 36*(degree-1);
    quadrature = (std::sqrt(8*quadrature+1) - 1.0)/2.0;
    if(quadrature<1)
        quadrature = 1;
    for(int i=0; i<quadrature; i++)
    {
        for(int j=0; j<quadrature-i; j++)
        {
            if(1 == degree)
            {
                uv_sampling.emplace_back( std::make_pair(0.0, 0.0) );
                i = quadrature+1;
                break;
            }
            uv_sampling.emplace_back( std::make_pair((1.0*i)/(quadrature-1.0), (1.0*j)/(quadrature-1.0)) );
        }
    }

    //Pre-compute Bezier coefficients
    for(unsigned int i=0; i<uv_sampling.size(); i++)
    {
        double& u = uv_sampling[i].first;
        double& v = uv_sampling[i].second;

        std::vector<double> uv_ijk;
        for(unsigned int cp=0; cp<C1_IJK.size(); cp++)
        {
            int i = C1_IJK[cp].first;
            int j = C1_IJK[cp].second;
            int k = degree - i - j - 1;
            uv_ijk.push_back( degree*bz_value(degree-1, i, j, k, u, v) );
        }
        UV_IJK.emplace_back(uv_ijk);
    }

    //Fill Neighborhood information for all the vertices
    fillNeighbourInfo();


    for(int i=0; i<optItr; i++)
    {
        //Logger::lout(Logger::INFO) << "Itr: " << i+1 << " / " << optItr << endl;
        for(int j=0; j<neighbourFaces.size(); j++)
        {
            gradient_step(j);
            //break;
        }
    }

    //Update triangulation variables
    int offset = allVertices.size();
    for( unsigned int i=0; i<allPts.size(); i++)
    {
        if(i < offset)
            allVertices[i] = allPts[i];
        else
            allCtrlpts[i-offset] = allPts[i];
    }
}

template struct BezierMesh<double>;
template struct BezierMesh<mpq_class>;

} // namespace bzmsh
