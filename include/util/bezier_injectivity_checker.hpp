#ifndef BZMSH_BEZIER_INJECTIVITY_CHECKER_HPP
#define BZMSH_BEZIER_INJECTIVITY_CHECKER_HPP

#include <queue>
#include <vector>

namespace bzmsh
{

/**
 * @brief Utility class to check validity of bezier triangles. 
 *
 */
template <typename T>
class BezierTriangleInjectivityChecker
{
  private:
    T min_area = 0.0;
    #ifdef ALL_EXACT
      int max_subdiv = 200; // number of subdivisions for conservative injectivity test
    #else
      int max_subdiv = 15;  // more does not make much sense in double precision
    #endif

    int current_level = 0;
    std::vector<T> v1, v2;
    int degreeB;
    int degreeC;

    //(degree+1)*i - ((i*(i-1))/2) + j;
    std::vector<int> indB;
    std::vector<int> indC;

    T b1(int i, int j)
    {
        return v1[indB[i] + j];
    }
    T b2(int i, int j)
    {
        return v2[indB[i] + j];
    }

    T xu(int i, int j)
    {
        return (b1(i + 1, j) - b1(i, j));
    } // constant factor does not matter: *3.0
    T xv(int i, int j)
    {
        return (b1(i, j + 1) - b1(i, j));
    } // constant factor does not matter: *3.0
    T yu(int i, int j)
    {
        return (b2(i + 1, j) - b2(i, j));
    } // constant factor does not matter: *3.0
    T yv(int i, int j)
    {
        return (b2(i, j + 1) - b2(i, j));
    } // constant factor does not matter: *3.0

    T P(const std::vector<T>& C, int i, int j, int n1, int n2, int n3);

    void bisect(const std::vector<T>& C,
                int degree,
                const std::vector<int>& ind,
                std::vector<T>& C0,
                std::vector<T>& C1);
    
    T facT(int n)
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
        return fac(n - 1) * n;
    }

  public:
    BezierTriangleInjectivityChecker()
    {
    }

    // control points are expected in lexicographic index order, i.e b003, b012, b021, b030, b102,
    // b111, b120, b201, b210, b300.
    BezierTriangleInjectivityChecker(int _degree) : degreeB(_degree), degreeC(2 * (_degree - 1))
    {
        set_degree(_degree);
    }

    void set_degree(int degree);
    void set_max_subdiv(int subdiv)
    {
        max_subdiv = subdiv;
    }

    int is_definitely_injective_by_determinant(const std::vector<T>& _b1,
                                                const std::vector<T>& _b2,
                                                bool via_bisect = true);

    int level()
    {
        return current_level;
    }
};

} // namespace bzmsh

#include "util/bezier_injectivity_checker.cpp"

#endif
