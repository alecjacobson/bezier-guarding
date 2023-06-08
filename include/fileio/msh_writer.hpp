#ifndef BZMSH_MSH_WRITER_HPP
#define BZMSH_MSH_WRITER_HPP

#include "datastructure/triangulation_constraints.hpp"

namespace bzmsh
{

template <typename T>
class MSHWriter
{
  public:
    uint indentWidth = 2;
    char indentChar = ' ';
    const int baseOffset = 1; // 0 index is forbidden in MSH format

    /**
     * @brief Write all triangle elements of a given mesh (including all control points)
     *        to a mesh file in gmsh file format.
     *
     * @param filename file to write output to
     * @param mesh mesh to print
     * @return true if writing was successful
     * @return false else
     */
    bool write(const std::string& filename, BezierMesh<T>& mesh);

  private:
    inline std::string indent(int depth)
    {
        return std::string(indentWidth * depth, indentChar);
    }

    inline int offset(int i, int j)
    {
        return i + j;
    }
};
} // namespace bzmsh
#endif
