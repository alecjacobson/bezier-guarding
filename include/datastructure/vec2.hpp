#ifndef BZMSH_VEC2_HPP
#define BZMSH_VEC2_HPP

#include "util/helper_functions.hpp"
#include <cmath>

namespace bzmsh
{

/**
 * @brief For representing 2D points and vectors
 *
 * @tparam T coordinate datatype
 */
template <typename T> struct Vec2
{
    /**
     * @brief Construct a 2D Vector
     */
    Vec2(const T& _x = 0, const T& _y = 0) : x(_x), y(_y)
    {
    }

    /**
     * @brief Get the vectors length
     *
     * @return float length
     */
    double length() const
    {
        return sqrt(utilCast<T, double>(lengthSquared()));
    }
    
    Vec2<T> normalized() const
    {
        return *this/length();
    }
    
    T norm1() const
    {
        if(x > 0)
        {
            if(y > 0) return x+y;
            else return x-y;
        }
        else
        {
            if(y > 0) return y-x;
            else return -x-y;
        }
    }
    
    Vec2<T> normalizedNorm1() const
    {
        return *this/norm1();
    }

    /**
     * @brief Get the vectors length
     *
     * @return float length
     */
    T lengthSquared() const
    {
        return x * x + y * y;
    }

    /**
     * @brief Comparison of two Vec2s
     * @param vec vector
     * @return vector
     */
    bool operator==(const Vec2<T>& vec) const
    {
        return (x == vec.x) && (y == vec.y);
    }

    bool operator!=(const Vec2<T>& vec) const
    {
        return !(*this == vec);
    }

    /**
     * @brief Vector addition
     * @param vec vector
     * @return vector
     */
    Vec2<T> operator+(const Vec2<T>& vec) const
    {
        return Vec2<T>(x + vec.x, y + vec.y);
    }

    Vec2<T> operator+=(const Vec2<T>& vec)
    {
        return *this = *this + vec;
    }

    /**
     * @brief Vector subtraction
     * @param vec vector
     * @return vector
     */
    Vec2<T> operator-(const Vec2<T>& vec) const
    {
        return Vec2<T>(x - vec.x, y - vec.y);
    }

    Vec2<T> operator-=(const Vec2<T>& vec)
    {
        return *this = *this - vec;
    }

    /**
     * @brief Scale vector
     * @param scale scalar
     * @return vector
     */
    Vec2<T> operator*(const T& scale) const
    {
        return Vec2<T>(x * scale, y * scale);
    }

    Vec2<T> operator*=(const T& scale)
    {
        return *this = *this * scale;
    }

    /**
     * @brief Scale vector
     * @param scale scalar
     * @return vector
     */
    Vec2<T> operator/(const T& scale) const
    {
        return Vec2<T>(x / scale, y / scale);
    }

    Vec2<T> operator/=(const T& scale)
    {
        return *this = *this / scale;
    }

    /**
     * @brief Dot product
     * @param vec vector
     * @return result (as a float)vec2
     */
    T operator*(const Vec2<T>& vec) const
    {
        return x * vec.x + y * vec.y;
    }

    /**
     * @brief 2D analogon of cross product
     *
     * @param vec vector to form cross product with
     * @return T 2D cross product
     */
    T cross(const Vec2<T>& vec) const
    {
        return x * vec.y - y * vec.x;
    }
    
    /**
     * @brief ccw orientation test
     *
     * @param vec vector to check orientation with
     * @return true if this ccw of vec
     */
    bool ccwOf(const Vec2<T>& vec) const
    {
        return this->cross(vec) < 0;
    }
    
    double AbsAngleWith(const Vec2<T>& vec) const
    {
        return std::acos(utilCast<T, double>((*this*vec)/this->length()/vec.length()));
    }
    
    Vec2<T> rot90ccw()
    {
        return Vec2<T>(-y, x);
    }

    /**
     * @brief Comparator for ordering vectors lexicographically (x more important than y)
     */
    struct compareLex
    {
        bool operator()(const Vec2<T>& lhs, const Vec2<T>& rhs) const
        {
            return (lhs.x < rhs.x) || ((lhs.x == rhs.x) && (lhs.y < rhs.y));
        }
    };
    
    friend std::ostream &operator<<( std::ostream &outp, const Vec2<T> &v ) {
       outp << "(" << v.x << ", " << v.y << ")";
       return outp;
    }

    /**
     * @brief The internal values of the vector
     */
    T x, y;
};

} // namespace bzmsh

#endif
