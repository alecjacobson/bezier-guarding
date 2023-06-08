#ifndef BZMSH_LOGGER_HPP
#define BZMSH_LOGGER_HPP

#include <ctime>
#include <functional>
#include <gmpxx.h>
#include <iomanip>
#include <iostream>
#include <vector>

namespace bzmsh
{
using std::endl;
using std::vector;
template <typename T> struct Vec2;
template <typename T> struct BezierCurve;
template <typename T> struct LineSegment;
template <typename T> struct GuardedBezierCurve;

/**
 * @brief Use to log events of any severity level
 *
 */
class Logger
{
  public:
    Logger() = delete;
    Logger(const Logger&) = delete;
    Logger(Logger&&) = delete;
    Logger& operator=(const Logger&) = delete;
    Logger& operator=(Logger&&) = delete;

    enum Level : uint8_t
    {
        NONE = 0b000000000,
        DEBUG = 0b00000001,
        INFO = 0b00000010,
        WARN = 0b00000100,
        ERROR = 0b00001000,
        FATAL = 0b00010000,
        ALL = 0b00111111
    };

    /**
     * @brief To configure messages of which severity should appear in logs.
     *        If NDEBUG is defined this is initialized to only Warnings, Error and Fatal,
     *        else everything gets logged by default.
     */
#ifdef NDEBUG
    static const Level showLvl = static_cast<Level>(Level::ALL & ~Level::DEBUG);
#else
    static const Level showLvl = Level::ALL;
#endif
    ;

    /**
     * @brief Get the LogLevel of a message as string
     *
     */
    static std::string toString(Level msgLvl)
    {
        switch (msgLvl)
        {
        case Level::NONE:
            return "None";
        case Level::DEBUG:
            return "Debug";
        case Level::INFO:
            return "Info";
        case Level::WARN:
            return "Warning";
        case Level::ERROR:
            return "Error";
        case Level::FATAL:
            return "Fatal";
        case Level::ALL:
            return "All";
        default:
            return "Unknown";
        }
    }

    static Logger& lout(Level msgLvl = Level::INFO)
    {
        static Logger l(std::cout);
        l.m_lvl = msgLvl;
        if ((l.m_lvl & showLvl) != 0)
        {
            auto t = std::time(nullptr);
            auto tm = *std::localtime(&t);
            l << "[" << toString(msgLvl) << "] " << std::put_time(&tm, "%H:%M:%S") << ": ";
        }
        return l;
    }

    Level m_lvl = Level::DEBUG;

    Logger& operator<<(const mpq_class& t)
    {
        if ((m_lvl & showLvl) != 0)
        {
            if (t.get_den() > 1e6)
                m_out << t.get_d();
            else
                m_out << t;
        }
        return *this;
    }

    template <typename T> Logger& operator<<(const T& t)
    {
        if ((m_lvl & showLvl) != 0)
        {
            m_out << t;
        }
        return *this;
    }

    template <typename T> Logger& operator<<(const Vec2<T>& v)
    {
        if ((m_lvl & showLvl) != 0)
        {
            *this << "(" << v.x << ", " << v.y << ")";
        }
        return *this;
    }

    template <typename T> Logger& operator<<(const vector<Vec2<T>>& v)
    {
        if ((m_lvl & showLvl) != 0)
        {
            m_out << "[";
            for (uint i = 0; i < v.size() - 1; i++)
            {
                *this << v[i] << ", ";
            }
            *this << v[v.size() - 1] << "]";
        }
        return *this;
    }

    template <typename T> Logger& operator<<(const BezierCurve<T>& c)
    {
        if ((m_lvl & showLvl) != 0)
        {
            *this << "(v)id: (" << c.m_origId << ")" << c.m_id << ", (v)tStart: (" << c.m_origTStart
                  << ")" << c.m_tStart << ", (v)tEnd: (" << c.m_origTEnd << ")" << c.m_tEnd
                  << ", ctrlpts: " << c.m_ctrlpts;
        }
        return *this;
    }

    template <typename T> Logger& operator<<(const LineSegment<T>& c)
    {
        if ((m_lvl & showLvl) != 0)
        {
            *this << "(v)id: (" << c.m_origId << ")" << c.m_id << ", (v)tStart: (" << c.m_origTStart
                  << ")" << c.m_tStart << ", (v)tEnd: (" << c.m_origTEnd << ")" << c.m_tEnd
                  << ", start: " << c[0] << ", end: " << c[1];
        }
        return *this;
    }

    template <typename T> Logger& operator<<(const GuardedBezierCurve<T>& gc)
    {
        return *this << gc.m_curve;
    }

    Logger& operator<<(std::ostream& (*manip)(std::ostream&))
    {
        if ((m_lvl & showLvl) != 0)
        {
            m_out << manip;
        }
        return *this;
    }

  private:
    Logger(std::ostream& out) : m_out(out)
    {
    }

    /**
     * @brief Which stream to log to (default is cout)
     *
     */
    std::ostream& m_out;
}; // namespace bzmsh

} // namespace bzmsh

#endif
