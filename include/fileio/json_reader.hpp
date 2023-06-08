#ifndef BZMSH_JSONIO_HPP
#define BZMSH_JSONIO_HPP

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include <vector>

#include "datastructure/bezier_curve.hpp"
#include "datastructure/line_segment.hpp"
#include "util/logger.hpp"

namespace bzmsh
{
using std::ifstream;

/**
 * @brief In/Output for .json files
 *
 */
class JSONReader
{
  public:
    /**
     * @brief Read curves from a .json file.
     *
     * @param filename file to read from
     * @param curves_out curves will be written here
     * @return true if reading successful
     * @return false else
     */
    bool read(const std::string& filename,
              vector<BezierCurve<double>>& curves_out);

  private:
    LineSegment<double> extractJsonLine(const nlohmann::json& attr);

    BezierCurve<double> extractJsonBezier(const nlohmann::json& attr);

    void extractJsonCurve(const nlohmann::json& attr, double& t1, double& t2, int& id);

    Vec2<double> extractJsonPoint(const nlohmann::json& array);
};

} // namespace bzmsh

#endif
