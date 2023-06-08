#include "fileio/json_reader.hpp"

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <vector>

#include "algorithm/curve_splitting.hpp"
#include "datastructure/bezier_curve.hpp"
#include "datastructure/line_segment.hpp"
#include "util/logger.hpp"

bool ignore_irregular_curve = false;

namespace bzmsh
{
using std::ifstream;

bool JSONReader::read(const std::string& filename, vector<BezierCurve<double>>& curvesOut)
{
    using nlohmann::json;
    using namespace std;

    ifstream is(filename);
    Logger::lout(Logger::INFO) << "Trying to read curves from file " + filename << endl;

    if (!is.good())
    {
        Logger::lout(Logger::ERROR) << "Bad input stream/file!" << endl;
        return false;
    }

    curvesOut.clear();
    json feature_info = json({});
    is >> feature_info;

    if (feature_info.empty())
    {
        Logger::lout(Logger::ERROR) << "The json data you tried to read was empty" << endl;
        return false;
    }

    int ignored = 0;
    int curvesRead = 0;
    int linesRead = 0;
    for (size_t i = 0; i < feature_info.size(); i++)
    {
        if (feature_info[i]["type"] == "Line")
        {
            LineSegment<double> l = extractJsonLine(feature_info[i]);
            if (l.startPt() == l.endPt() || l.m_tStart == l.m_tEnd || !crop(l))
            {
                Logger::lout(Logger::INFO)
                    << "Input Line " << l << " was ill-defined and was ignored!" << endl;
                ignored++;
            }
            else
            {
                linesRead++;
                curvesOut.emplace_back(BezierCurve<double>(l));
            }
        }
        else if (feature_info[i]["type"] == "BezierCurve")
        {
            BezierCurve<double> c = extractJsonBezier(feature_info[i]);
            if (c.degree() <= 0 || c.m_tStart == c.m_tEnd || !crop(c))
            {
                Logger::lout(Logger::INFO)
                    << "Input Curve " << c << " was ill-defined and was ignored!"
                    << endl;
                ignored++;
            }
            else if(ignore_irregular_curve && (c[0] == c[1] || c[c.degree()] == c[c.degree()-1]))
            {
                Logger::lout(Logger::INFO)
                    << "Input Curve " << c << " was irregular and was ignored!"
                    << endl;
                ignored++;
            }
            else if(c.m_tStart != c.m_tEnd)
            {
                c = BezierCurve<double>(c.m_ctrlpts, 0.0, 1.0, c.m_id);
                curvesRead++;
                curvesOut.emplace_back(c);
            }
        }
    }

    Logger::lout(Logger::INFO) << "Line segments read: " << linesRead << ", "
                               << "Bezier curves read: " << curvesRead << ", "
                               << "Bad curves ignored: " << ignored << endl;

    if (curvesOut.size() == 0)
        return false;

    return true;
}

LineSegment<double> JSONReader::extractJsonLine(const nlohmann::json& attr)
{
    double t1, t2;
    int id;
    vector<Vec2<double>> ctrlpts;

    ctrlpts.emplace_back(extractJsonPoint(attr["start"]));
    ctrlpts.emplace_back(extractJsonPoint(attr["end"]));

    extractJsonCurve(attr, t1, t2, id);

    return LineSegment<double>(ctrlpts, t1, t2, id);
}

BezierCurve<double> JSONReader::extractJsonBezier(const nlohmann::json& attr)
{
    double t1 = 0;
    double t2 = 1;
    int id = -1;
    vector<Vec2<double>> ctrlpts;

    for (size_t i = 0; i < attr["poles"].size(); i++)
    {
        ctrlpts.emplace_back(extractJsonPoint(attr["poles"][i]));
    }

    extractJsonCurve(attr, t1, t2, id);

    size_t d = ctrlpts.size()-1;
    if(d<2)
        return BezierCurve<double>(ctrlpts, t1, t2, id);

    return BezierCurve<double>(ctrlpts, t1, t2, id);
}

void JSONReader::extractJsonCurve(const nlohmann::json& attr, double& t1, double& t2, int& id)
{
    static uint uid = 0;

    if (attr.contains("paras"))
    {
        t1 = static_cast<double>(attr["paras"][0]);
        t2 = static_cast<double>(attr["paras"][1]);
        if(t1>t2)
            std::swap(t1, t2);
    }
    else
    {
        t1 = 0.0;
        t2 = 1.0;
    }

    if (attr.contains("curve_id"))
    {
        id = attr["curve_id"];
        uid = std::max(id, (int)uid) + 1;
    }
    else
    {
        id = uid++;
    }
}

Vec2<double> JSONReader::extractJsonPoint(const nlohmann::json& array)
{
    return Vec2<double>(array[0], array[1]);
}

} // namespace bzmsh
