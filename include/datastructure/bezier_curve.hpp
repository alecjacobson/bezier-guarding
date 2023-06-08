#ifndef BZMSH_BEZIER_CURVE_HPP
#define BZMSH_BEZIER_CURVE_HPP

#include "datastructure/curve.hpp"
#include "datastructure/line_segment.hpp"
#include "util/helper_functions.hpp"
#include <gmpxx.h>

namespace bzmsh
{
using std::max;
using std::pair;

/**
 * @brief Represents a bezier curve
 *
 * @tparam T float-type (used for double and mpq_class)
 */
template <typename T> struct BezierCurve : public Curve<T>
{
    using Curve<T>::m_ctrlpts;
    using Curve<T>::m_tStart;
    using Curve<T>::m_tEnd;
    using Curve<T>::m_id;
    using Curve<T>::degree;

    /**
     * @brief Construct a new bezier curve
     *
     * @param ctrlpts The bezier curves control points
     * @param tStart Starting t
     * @param tEnd Ending t
     * @param id curve id
     */
    BezierCurve(const vector<Vec2<T>>& ctrlpts, T tStart = 0, T tEnd = 1, int id = 0)
        : Curve<T>(ctrlpts, tStart, tEnd, id)
    {
    }

    /**
     * @brief Convert line segment to bezier curve
     *
     * @param l line
     */
    BezierCurve(const LineSegment<T>& l) : Curve<T>(l)
    {
    }

    BezierCurve() : Curve<T>()
    {
    }

    /**
     * @brief Check whether this curve is irregular
     *
     * @return true if first and second or last and second to last control points are equal
     * @return false else
     */
    bool irregular() const
    {
        size_t d = this->degree();
        assert(d > 0);
        return (m_ctrlpts[1] == m_ctrlpts[0]) || (m_ctrlpts[d] == m_ctrlpts[d - 1]);
    }

    bool irregular2() const
    {
        size_t d = this->degree();
        assert(d > 0);
        return (m_ctrlpts[1] == m_ctrlpts[0]) || (m_ctrlpts[d] == m_ctrlpts[d - 1]) || (m_ctrlpts[0] == m_ctrlpts[d]);
    }

    /**
     * @brief Evaluate this curve for parameter @p t
     * @pre t must be in [0, 1]
     *
     * @param t parametrization value
     * @return Vec2<T> this curves value at @p t
     */
    virtual Vec2<T> operator()(const T& t) const
    {
        assert(t >= T(0) && t <= T(1));
        if (t == T(0))
            return m_ctrlpts[0];
        if (t == T(1))
            return m_ctrlpts[this->degree()];

        Vec2<T> b(0, 0);
        uint d = this->degree();
        T fac_n = T(fac(d));
        vector<T> powT = {T(1)};
        vector<T> powOneMinusT = {T(1)};
        for (uint i = 1; i <= d; i++)
        {
            powT.emplace_back(powT.back() * t);
            powOneMinusT.emplace_back(powOneMinusT.back() * (T(1) - t));
        }
        for (uint i = 0; i <= d; i++)
        {
            b += m_ctrlpts[i]
                 * (fac_n * powT[i] * powOneMinusT[d - i] / (T(fac(i)) * T(fac(d - i))));
        }

        return b;
    }

    /**
     * @brief Returns true if the curve is guardable
     *
     * @return is the curve guardable?
     */
    bool isGuardable() const
    {
        return (curveAxis() != Vec2<T>(0,0));
    }
    
    bool isWellGuardable() const
    {
        Vec2<T> axis = curveAxis();
              
        if(axis.norm1() == 0) return false;
        
        Vec2<T> sp = splus();
        Vec2<T> sm = sminus();
        
        if(sp.AbsAngleWith(sm) > toRadian(GUARDING_MAX_CURVATURE_DEGREES)) return false;
        if(sp.AbsAngleWith(axis) > toRadian(GUARDING_MAX_CURVATURE_DEGREES)*0.5) return false;
        if(sm.AbsAngleWith(axis) > toRadian(GUARDING_MAX_CURVATURE_DEGREES)*0.5) return false;
        
        return true;
    }
    
    /**
    * @brief Returns an axis if the curve is guardable, (0,0) otherwise
    *
    * @return Vec2<T> curve axis vector
    */
    Vec2<T> curveAxis() const
       {
         size_t d = this->degree();
         Vec2<T> zero(0,0);
         if(d <= 0) return zero;
         if(d == 1) return m_ctrlpts[1] - m_ctrlpts[0];
          
          Vec2<T> splus = m_ctrlpts[1] - m_ctrlpts[0];
          Vec2<T> sminus = m_ctrlpts[2] - m_ctrlpts[1];
          if(sminus.ccwOf(splus)) std::swap(splus, sminus);
          
          for(uint i = 2; i < d; i++)
          {
            Vec2<T> si = m_ctrlpts[i+1] - m_ctrlpts[i];
                        
            if(si.ccwOf(splus))
            {
              if(si.ccwOf(sminus)) splus = si;
              else return zero;
            }
            else
            {
              if(si.ccwOf(sminus)) {} //within cone
              else sminus = si;
            }
          }
          
          // try nice axis (using non-rational 2-norm)
          Vec2<T> axis = (splus * sminus.length()) + (sminus * splus.length());

          if(axis*splus <= 0 || axis*sminus <= 0) // use safe axis (using rational 1-norm)
            axis = (splus * sminus.norm1()) + (sminus * splus.norm1());
            
          if(axis*splus <= 0 || axis*sminus <= 0)
            return zero;
          
          return axis;
       }
       
    Vec2<T> splus() const
    {
       size_t d = this->degree();
       Vec2<T> zero(0,0);
       if(d <= 0) return zero;
       Vec2<T> splus = m_ctrlpts[1] - m_ctrlpts[0];
      
       for(uint i = 1; i < d; i++)
       {
         Vec2<T> si = m_ctrlpts[i+1] - m_ctrlpts[i];
         if(si.ccwOf(splus)) splus = si;
       }
       
       return splus;
    }
    
    Vec2<T> sminus() const
    {
       size_t d = this->degree();
       Vec2<T> zero(0,0);
       if(d <= 0) return zero;
       Vec2<T> sminus = m_ctrlpts[1] - m_ctrlpts[0];
      
       for(uint i = 1; i < d; i++)
       {
         Vec2<T> si = m_ctrlpts[i+1] - m_ctrlpts[i];
         if(sminus.ccwOf(si)) sminus = si;
       }
       
       return sminus;
    }
    
    double angleRange() const
    {
      size_t d = this->degree();
      double angleMin = 0;
      double angleMax = 0;
      double angleCur = 0;
      for(uint i = 1; i < d; i++)
      {
        Vec2<T> s0 = m_ctrlpts[i] - m_ctrlpts[i-1];
        Vec2<T> s1 = m_ctrlpts[i+1] - m_ctrlpts[i];
        double angle = std::acos(utilCast<T,double>((s0.normalized())*(s1.normalized())));
        if(s0.ccwOf(s1)) angleCur -= angle;
        else angleCur += angle;
        if(angleCur < angleMin) angleMin = angleCur;
        if(angleCur > angleMax) angleMax = angleCur;
      }
      return angleMax - angleMin;
    }
      

    /**
     * @brief Increase this curves degree by creating new controlpoints. Curve itself remains
     *        unchanged
     *
     * @param delta by how much to increase degree
     * @return true if curve has controlpoints
     * @return false else
     */
    bool increaseDegree(uint delta)
    {
        if (degree() < 0)
        {
            return false;
        }
        while (delta-- > 0)
        {
            vector<Vec2<T>> ctrlpts = {m_ctrlpts[0]};
            for (uint i = 1; i < m_ctrlpts.size(); i++)
            {
                Vec2<T> pt = (m_ctrlpts[i - 1] * T(i) + m_ctrlpts[i] * T(degree() + 1 - i))
                             / T(degree() + 1);
                ctrlpts.emplace_back(pt);
            }
            ctrlpts.emplace_back(m_ctrlpts[degree()]);
            m_ctrlpts = ctrlpts;
        }
        return true;
    }

    /**
     * @brief check whether this and another (or the same) curve have equal endpoints and an
     *        (almost) equal tangent vector at those
     *
     * @param other other curve (may be this one)
     * @param angleEpsilon
     * @return true if common endpoint and similar tangent
     * @return false else
     */
    bool commonTangent(const BezierCurve<T>& other, double angleEpsilon) const
    {
        uint d = degree();
        uint otherd = other.degree();
        assert(d > 0 && otherd > 0);
        if (&other == this)
        {
            return (this->startPt() == this->endPt()
                    && angleRad((*this)[1] - (*this)[0], (*this)[d - 1] - (*this)[d])
                           < angleEpsilon);
        }

        return (this->startPt() == other.startPt()
                && angleRad((*this)[1] - (*this)[0], other[1] - other[0]) < angleEpsilon)
               || (this->startPt() == other.endPt()
                   && angleRad((*this)[1] - (*this)[0], other[otherd - 1] - other[otherd])
                          < angleEpsilon)
               || (this->endPt() == other.startPt()
                   && angleRad((*this)[d - 1] - (*this)[d], other[1] - other[0]) < angleEpsilon)
               || (this->endPt() == other.endPt()
                   && angleRad((*this)[d - 1] - (*this)[d], other[otherd - 1] - other[otherd])
                          < angleEpsilon);
    }
};

} // namespace bzmsh

#endif
