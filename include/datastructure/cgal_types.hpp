#ifndef CGAL_TYPES_H
#define CGAL_TYPES_H

#include <CGAL/gmpxx.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel_exact;

typedef Kernel_exact::FT FT_exact;
typedef Kernel_exact::Point_2 Point_2_exact;
typedef Kernel_exact::Point_3 Point_3_exact;
typedef Kernel_exact::Triangle_2 Triangle_2_exact;
typedef Kernel_exact::Triangle_3 Triangle_3_exact;
typedef Kernel_exact::Segment_2 Segment_2_exact;
typedef Kernel_exact::Ray_2 Ray_2_exact;
typedef Kernel_exact::Direction_2 Direction_2_exact;
typedef Kernel_exact::Vector_2 Vector_2_exact;
typedef Kernel_exact::Line_2 Line_2_exact;

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
typedef CGAL::Exact_predicates_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel_exact, CGAL::Default, Itag> CDT_2D_EXACT;

typedef CDT_2D_EXACT::Finite_faces_iterator CDT_Faces_iterator;
typedef CDT_2D_EXACT::Face_handle CDT_Face;
typedef CDT_2D_EXACT::Vertex_handle CDT_Vertex;
typedef CDT_2D_EXACT::Edge CDT_Edge;

#endif // CGAL_TYPES_H


