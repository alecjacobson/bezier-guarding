#ifndef BZMSH_EPSWRITER_HPP
#define BZMSH_EPSWRITER_HPP

#include "datastructure/bbox.hpp"
#include "datastructure/bezier_curve.hpp"
#include "datastructure/guarded_bezier_curve.hpp"
#include "datastructure/triangulation_constraints.hpp"
#include <fstream>
#include <ostream>

namespace bzmsh
{

template <typename T>
class EPSWriter
{
  public:
    struct Color
    {
        Color(double _r, double _g, double _b) : r(_r), g(_g), b(_b)
        {
        }
        double r = 0.0;
        double g = 0.0;
        double b = 0.0;
    };

    enum LineStyle
    {
        NONE,
        SOLID,
        DASHED,
        DOTTED,
        DASHDOTTED
    };

    Color colorCurves = {0.1, 0.5, 0.0};
    Color colorEdges = {0.6, 0.6, 0.6};
    Color colorSubcurves = {0.85, 0.85, 0.85};
    Color colorVertices = {0.0, 0.0, 1.0};
    Color colorCtrlpts = {1.0, 0.0, 0.0};
    LineStyle lineStyleCurves = LineStyle::SOLID;
    LineStyle lineStyleEdges = LineStyle::SOLID;
    LineStyle lineStyleSubcurves = LineStyle::DASHED;
    double lineWidthCurves = 0.5;
    double lineWidthEdges = 0.2;
    double lineWidthSubcurves = 0.1;
    double radiusVertices = 0.35;
    double radiusCtrlpts = 0.25;
    int samplesPerCurve = 20;
    int samplesPerSubcurve = 10;
    int subLevels = 1;

    double scalingFactor = 1.0;
    double shiftX = 0.0;
    double shiftY = 0.0;

    /**
     * @brief Write a visualization of a mesh to an .eps postscript file (vector format).
     *        Bezier curves of degree 3 are natively supported, other degrees will be sampled.
     *
     * @param filename file to write output to
     * @param mesh mesh to print
     * @param resetScaling if this writer was used to print another mesh, should scaling be reset?
     * @return true if writing was successful
     * @return false else
     */
    bool write(const std::string& filename, const BezierMesh<T>& mesh, bool resetScaling = true);

  private:
    /**
     * @brief Write the eps header including the bbox to a file stream.
     *
     * @param out output file stream
     * @param bbox bbox information to include in header
     */
    void writeHeader(std::ostream& out, const BBox<double>& bbox) const;

    /**
     * @brief Write all edges (curves and normal edges and subdivision edges) contained in the mesh
     * as curves/lines in eps output. Line width, style and color configurable
     *
     * @param out output file stream
     * @param allEdges edges to output
     * @param allVertices vertex position information
     * @param allCtrlpts control point position information
     */
    void writeEdges(std::ostream& out,
                    const vector<Edge>& allEdges,
                    const vector<Vec2<double>>& allVertices,
                    const vector<Vec2<double>>& allCtrlpts);

    /**
     * @brief Write all vertices contained in the mesh as circles. Radius and color configurable.
     *
     * @param out output file stream
     * @param allVertices vertex position information
     */
    void writeVertices(std::ostream& out, const vector<Vec2<double>>& allVertices) const;

    /**
     * @brief Write all control points contained in the mesh as circles. Radius and color
     * configurable.
     *
     * @param out output file stream
     * @param allCtrlpts control point position information
     */
    void writeCtrlpts(std::ostream& out, const vector<Vec2<double>>& allCtrlpts, const vector<bool>& flipped) const;

    /**
     * @brief Write a circle with given center and radius to file stream.
     *
     * @param out output file stream
     * @param center center position
     * @param radius circle radius
     */
    void writeCircle(std::ostream& out, const Vec2<double>& center, double radius) const;

    /**
     * @brief Write a certain linestyle to the output file stream (affects any following drawn
     * features)
     *
     * @param out output file stream
     * @param c color
     * @param lineStyle line style
     * @param lineWidth line width
     */
    void writeLineStyle(std::ostream& out, Color c, LineStyle lineStyle, double lineWidth) const;

    /**
     * @brief Write cubic bezier curves using its control points to output file stream
     *
     * @param out output file stream
     * @param c control points
     */
    void writeCurve(std::ostream& out, const vector<Vec2<double>>& ctrlpts) const;

    /**
     * @brief Write cubic bezier curves using its control points to output file stream
     *
     * @param out output file stream
     * @param c control points
     */
    void writeCubicBezier(std::ostream& out, const vector<Vec2<double>>& c) const;

    /**
     * @brief Write a bezier curve of any degree to ouput file stream by sampling it.
     *
     * @param out output file stream
     * @param c curve
     * @param sampleCount number of samples equally distributed
     */
    void writeSampledBezier(std::ostream& out, const BezierCurve<double>& c, int sampleCount) const;

    /**
     * @brief Write a straight line to output file stream.
     *
     * @param out output file stream
     * @param l line to draw
     */
    void writeLine(std::ostream& out, const LineSegment<double>& l) const;
};

template class EPSWriter<double>;
template class EPSWriter<mpq_class>;

} // namespace bzmsh


#endif
