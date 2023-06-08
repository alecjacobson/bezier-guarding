#include "util/bezier_injectivity_checker.hpp"

#include <algorithm>

namespace bzmsh
{

template <typename T>
void BezierTriangleInjectivityChecker<T>::set_degree(int degree)
{
    degreeB = degree;
    degreeC = 2 * (degree - 1);

    indB = std::vector<int>(degreeB + 1);
    int val = 0;
    int add = degreeB + 1;
    for (int i = 0; i <= degreeB; i++)
    {
        indB[i] = val;
        val += add;
        add--;
    }
    indC = std::vector<int>(degreeC + 1);
    val = 0;
    add = degreeC + 1;
    for (int i = 0; i <= degreeC; i++)
    {
        indC[i] = val;
        val += add;
        add--;
    }
}

template <typename T>
T BezierTriangleInjectivityChecker<T>::P(
    const std::vector<T>& C, int i, int j, int n1, int n2, int n3)
{
    if (n1 > 0)
        return 0.5 * (P(C, i, j, n1 - 1, n2, n3) + P(C, i - 1, j + 1, n1 - 1, n2, n3));
    if (n2 > 0)
        return 0.5 * (P(C, i, j, n1, n2 - 1, n3) + P(C, i, j - 1, n1, n2 - 1, n3));
    if (n3 > 0)
        return 0.5 * (P(C, i, j, n1, n2, n3 - 1) + P(C, i + 1, j, n1, n2, n3 - 1));
    return C[indC[i] + j];
}

template <typename T>
void BezierTriangleInjectivityChecker<T>::bisect(const std::vector<T>& C,
                                                   int degree,
                                                   const std::vector<int>& ind,
                                                   std::vector<T>& C0,
                                                   std::vector<T>& C1)
{
    std::vector<T> dCrow(degree + 1);
    for (int i = 0; i <= degree; i++)
    {
        for (int j = 0; j <= degree - i; j++)
            dCrow[j] = C[ind[i] + j];
        C0[ind[0] + i] = dCrow[degree - i];
        C1[ind[0] + degree - i] = dCrow[0];
        for (int l = 1; l <= degree - i; l++)
        {
            for (int j = 0; j <= degree - i - l; j++)
                dCrow[j] = 0.5 * (dCrow[j] + dCrow[j + 1]);
            C0[ind[l] + i] = dCrow[degree - i - l];
            C1[ind[l] + degree - i - l] = dCrow[0];
        }
    }
}


template <typename T>
int BezierTriangleInjectivityChecker<T>::is_definitely_injective_by_determinant(
    const std::vector<T>& _b1, const std::vector<T>& _b2, bool via_bisect)
{
    v1 = _b1;
    v2 = _b2;

    // compute quartic control values
    int derivDegree = degreeB - 1;
    std::vector<T> C(((degreeC + 1) * (degreeC + 2)) / 2, 0.0);
    for (int i = 0; i <= degreeC; i++)
        for (int j = 0; j <= degreeC - i; j++)
        {
            for (int ri = 0; ri <= derivDegree; ri++)
                for (int rj = 0; rj <= derivDegree - ri; rj++)
                {
                    int si = i - ri;
                    int sj = j - rj;
                    int rk = derivDegree - ri - rj;
                    int sk = derivDegree - si - sj;
                    if (si < 0)
                        continue;
                    if (sj < 0)
                        continue;
                    if (si + sj > derivDegree)
                        continue;
                    if (degreeC - ri - rj - si - sj != degreeC - i - j)
                        continue;
                    T lijk = (facT(i) * facT(j) * facT(degreeC - i - j))
                                  / (facT(ri) * facT(rj) * facT(rk) * facT(si) * facT(sj)
                                     * facT(sk)); // constant factor does not matter:  * 1.0/6.0;
                                                 // (general: 1.0/facT(degreeB))
                    C[indC[i] + j] += lijk * (xu(ri, rj) * yv(si, sj) - xv(ri, rj) * yu(si, sj));
                }
        }

    std::queue<std::pair<std::vector<T>, int>> q;
    q.push(std::make_pair(C, 0));
    while (!q.empty())
    {
        std::vector<T> C = q.front().first;
        current_level = q.front().second;
        q.pop();

        T mi = *std::min_element(C.begin(), C.end());
        if (mi > min_area)
            continue; // determinant positive everywhere in this triangle

        T mic
            = std::min(C[indC[degreeC] + 0],
                       std::min(C[indC[0] + degreeC], C[indC[0] + 0])); // minimum of corner values
        if (mic <= min_area)
        {
            return 0; // determinant non-positive at some (corner) point
        }

        if (current_level >= max_subdiv)
        {
            return 2; // injectivity could not be certified
        }

        if (via_bisect)
        {
            // two subtriangles
            std::vector<T> C0(C.size());
            std::vector<T> C1(C.size());
            bisect(C, degreeC, indC, C0, C1);
            q.push(std::make_pair(C0, current_level + 1));
            q.push(std::make_pair(C1, current_level + 1));
        }
        else
        {
            // four subtriangles
            std::vector<T> C0, C1, C2, C3;
            for (int i = 0; i <= degreeC; i++)
                for (int m = 0; m <= i; m++)
                    C0.push_back(P(C, i, 0, m, 0, degreeC - i));
            for (int i = 0; i <= degreeC; i++)
                for (int m = 0; m <= degreeC - i; m++)
                    C1.push_back(P(C, i, degreeC - i, i, m, 0));
            for (int i = 0; i <= degreeC; i++)
                for (int m = 0; m <= degreeC - i; m++)
                    C2.push_back(P(C, 0, i, 0, i, m));
            for (int i = 0; i <= degreeC; i++)
                for (int j = 0; j <= degreeC - i; j++)
                    C3.push_back(P(C, i, j, i, j, degreeC - i - j));

            q.push(std::make_pair(C0, current_level + 1));
            q.push(std::make_pair(C1, current_level + 1));
            q.push(std::make_pair(C2, current_level + 1));
            q.push(std::make_pair(C3, current_level + 1));
        }
    }
    return 1;
}

} // namespace bzmsh
