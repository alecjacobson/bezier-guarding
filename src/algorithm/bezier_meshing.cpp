#include "algorithm/bezier_meshing.hpp"

#include <gmpxx.h>
#include <iostream>
#include <set>

#include "algorithm/constants.hpp"
#include "algorithm/curve_splitting.hpp"
#include "algorithm/intersection_tests.hpp"
#include "algorithm/triangulating.hpp"
#include "datastructure/guarded_bezier_curve.hpp"
#include "datastructure/intersection_tree.hpp"
#include "fileio/eps_writer.hpp"
#include "fileio/json_reader.hpp"
#include "fileio/msh_writer.hpp"

#include "algorithm/timer.hpp"

#include "config.hpp"


namespace bzmsh
{
using std::max;
using std::multimap;
using std::set;
using std::string;

void save_log_file(int t, std::string f_name, int t_count=-1)
{
    std::string file_name;
    if(-1 == t) //CURVE_OVERLAPPING_NODES
        file_name = "common_nodes2.txt";
    else if(0 == t) //CURVE_OVERLAPPING_NODES
        file_name = "common_nodes.txt";
    else if(1 == t) //CURVE_OVERLAP
        file_name = "overlap.txt";
    else if(2 == t) //CURVE_CO_TANGENTS
        file_name = "cotangents.txt";
    else if(3 == t) //Flip Uncertainity
        file_name = "uncertain_flipped.txt";
    else if(4 == t) //CURVE_STRAIGHT_FLIPPED
        file_name = "others_flipped.txt";
    else if(5 == t) //CURVE_GOOD
        file_name = "good_curve.txt";
    else if(6 == t) //more intersections
        file_name = "more.txt";
    else if(7 == t) //guard still intersecting
        file_name = "guard_fail.txt";
    else if(8 == t) //guard still intersecting
        file_name = "guard_intrfail.txt";
    else if(10 == t) //cdt failed
        file_name = "cdt.txt";
    else if(11 == t) //control points filling failed
        file_name = "bz.txt";
    else if(19 == t) //total timings
        file_name = "timings.txt";
    else if(20 == t) //all curves for logging
        file_name = "process.txt";
    else if(21 == t || 22 == t) //triangle count
        file_name = "triangles.txt";
    else
        file_name = "xxx.txt";

    FILE *error_file = fopen(file_name.c_str(), "a");
    if(19 == t)
    {
        double time_spent = subtractTimes(t_start)/1e-3;
        fprintf( error_file, "%s %f\n", f_name.c_str(), time_spent );
    }
    else if(20 == t)
        fprintf( error_file, "%s %d\n", f_name.c_str(), t_count );
    else if(21 == t)
        fprintf( error_file, "%d", t_count );
    else if(22 == t)
        fprintf( error_file, "\n" );
    else
        fprintf( error_file, "%s\n", f_name.c_str() );
    fflush(error_file);
    fclose(error_file);
}

bool bezierMeshing(std::string filename, bool addBbox, bool epsOutput, bool gmshOutput, bool optimizeMesh, int exp_type)
{
#ifdef ALL_EXACT
    Logger::lout(Logger::INFO) << "Running in EXACT mode." << endl;
#else
    Logger::lout(Logger::INFO) << "Running in FLOAT mode." << endl;
#endif

    if(exp_type>=0)
    {
        int mu = 1;
        while(exp_type > 0)
        {
            mu *= 10;
            exp_type--;
        }
        printf("setting GUARDING_NORMAL_SHIFT = 1 / %d\n", mu);
        GUARDING_NORMAL_SHIFT = mpq_class(1, mu);
    }

    vector<BezierCurve<double>> curves;
    if (!JSONReader().read(filename, curves))
    {
        Logger::lout(Logger::FATAL) << "Could not read json file, aborting..." << endl;
        return false;
    }

    save_log_file(20, filename, curves.size());
    save_log_file(22, filename);

    Logger::lout(Logger::INFO) << "Total curves =  " << curves.size() << endl;
    
    t_start = std::chrono::high_resolution_clock::now();
    original_file_name = filename;

    int degree = 1;
    for (uint i = 0; i < curves.size(); i++)
    {
        // Algorithm relies on curve id matching vector index (m_origId remains unchanged)
        curves[i].m_id = i;
        degree = max(curves[i].degree(), degree);
        if (curves[i].irregular())
        {
            Logger::lout(Logger::FATAL) << "Irregular curve encountered, aborting..." << endl;
            save_log_file(0, filename);
            return false;
        }
    }

    for (BezierCurve<double>& curve : curves)
    {
        curve.increaseDegree(degree - curve.degree());
    }

    Logger::lout(Logger::INFO) << "Partitioning curves into guardable subcurves..."
                               << endl;

    vector<BezierCurve<mpq_class>> virtualCurveSegments;
    for (const BezierCurve<double>& c : curves)
    {
        virtualCurveSegments.emplace_back(convertCurve<double, mpq_class>(c));
    }

    splitNonWellGuardables(virtualCurveSegments);
    for (uint i = 0; i < virtualCurveSegments.size(); i++)
    {
        // Splitting might produce irregular curves in limited precision mode, so check again
        if (virtualCurveSegments[i].irregular())
        {
            Logger::lout(Logger::FATAL)
                << "Partitioning led to numerically irregular curves, aborting..." << endl;
            save_log_file(0, filename);
            return false;
        }
    }

    Logger::lout(Logger::INFO)
        << "Looking for intersections between exact curve segments, number of actual curves: "
        << curves.size() << "; number of curve segments: " << virtualCurveSegments.size() << endl;

    map<uint, SetOfSplits> curveIdToSplits;
    if (!getIntersections(virtualCurveSegments, curveIdToSplits))
    {
        Logger::lout(Logger::FATAL) << "Some curves were overlapping, aborting..." << endl;
        save_log_file(1, filename);
        return false;
    }
    if (!curveIdToSplits.empty())
    {
        Logger::lout(Logger::WARN) << curveIdToSplits.size()
                                   << " curves were intersecting, trying to split those..." << endl;

        // try to remove all intersections by splitting curves at intersection points
        vector<BezierCurve<double>> curvesCopy = curves;
        curves.clear();
        for (uint i = 0; i < curvesCopy.size(); i++)
        {
            vector<BezierCurve<double>> splitCurve;
            split(curvesCopy[i], curveIdToSplits[i], splitCurve);
            curves.insert(curves.end(), splitCurve.begin(), splitCurve.end());
        }
        curvesCopy.clear();

        // Reindex and check if still regular after splitting
        for (uint i = 0; i < curves.size(); i++)
        {
            curves[i].m_id = i;
            curves[i].m_origTStart = 0.0;
            curves[i].m_origTEnd = 1.0;

            if (curves[i].irregular())
            {
                Logger::lout(Logger::FATAL)
                    << "Splitting at intersections led to irregular curves, aborting..." << endl;
                save_log_file(0, filename);
                return false;
            }
        }
        Logger::lout(Logger::INFO)
            << "Number of actual curves after splitting: " << curves.size() << endl;
          

        // Repeat intersection test
        Logger::lout(Logger::INFO)
            << "Checking whether split curves are properly representable in "
               "the desired precision..."
            << endl;

        virtualCurveSegments.clear();
        for (const BezierCurve<double>& c : curves)
        {
            virtualCurveSegments.emplace_back(convertCurve<double, mpq_class>(c));
        }
        splitNonWellGuardables(virtualCurveSegments);
        Logger::lout(Logger::INFO)
            << "Looking for intersections between segments, "
            << "number of actual curves: " << curves.size()
            << "; number of curve segments: " << virtualCurveSegments.size() << endl;
        if (!getIntersections(virtualCurveSegments, curveIdToSplits))
        {
            Logger::lout(Logger::FATAL)
                << "There are overlaps, aborting..." << endl;
            save_log_file(1, filename);
            return false;
        }

        if (!curveIdToSplits.empty())
        {
            save_log_file(6, filename);
            Logger::lout(Logger::FATAL) << curveIdToSplits.size()
                                        << " curves are intersecting in the desired precision after attempting to "
                                           "resolve all intersections, aborting..."
                                        << endl;
            return false;
        }
        else
        {
            Logger::lout(Logger::INFO)
                << "Managed to resolve all intersections, number of curves now: " << curves.size()
                << endl;
        }
    }
    else // from "if (!curveIdToSplits.empty())" -> NO initial intersections
    {
        Logger::lout(Logger::INFO) << "No intersecting curves were found, continuing..." << endl;
    }
    

    // Splitting for guardability is already done, so we can use these results
    curves.clear();
    for (uint i = 0; i < virtualCurveSegments.size(); i++)
    {
        curves.emplace_back(convertCurve<mpq_class, double>(virtualCurveSegments[i]));
        curves[i].m_id = i;
        curves[i].m_origTStart = 0.0;
        curves[i].m_origTEnd = 1.0;
    }
    virtualCurveSegments.clear();

    Logger::lout(Logger::INFO) << "Checking for degenerate corners..." << endl;
    multimap<uint, uint> similarTangent;
    double angleEpsilon = toRadian(MIN_TANGENT_ANGLE_DEGREES);
    IntersectionTree tree(curves.size());
    for (uint i = 0; i < curves.size() - 1; i++)
    {
        tree.insert(curves[i].bbox());
        BezierCurve<double>& curveI = curves[i];

        for (uint j : tree.query(i))
        {
            BezierCurve<double>& curveJ = curves[j];
            if (curveI.commonTangent(curveJ, angleEpsilon))
            {
                similarTangent.emplace(std::make_pair(i, j));
            }
        }
    }
    if (!similarTangent.empty())
    {
        save_log_file(2, filename);
        Logger::lout(Logger::FATAL) << "Near-zero-angle corners encountered in " << similarTangent.size()
                                    << " cases; (perhaps try lower MIN_TANGENT_ANGLE_DEGREES). Aborting..." << endl;
        return false;
    }

    Logger::lout(Logger::INFO) << "Wrapping curves in guarding quad envelopes..." << endl;
    vector<GuardedBezierCurve<Scalar>> guardedCurves;
    if (!guardAll(curves, guardedCurves))
    {
        save_log_file(7, filename);
        Logger::lout(Logger::FATAL) << "Guarding of curves failed, aborting..." << endl;
        return false;
    }
    curves.clear();
    virtualCurveSegments.clear();

    Logger::lout(Logger::INFO) << "Splitting envelopes until no longer intersecting..."
                               << endl;
    if (!splitIntersecting(guardedCurves))
    {
        save_log_file(8, filename);
        Logger::lout(Logger::FATAL) << "Unresolvable intersections between guarded Bezier curves, "
                                       "perhaps numerical issue, or MAX_GUARDING_DEPTH set too low. Aborting..."
                                    << endl;
        return false;
    }

    Logger::lout(Logger::INFO)
        << "Successfully split all intersecting guarded curves, number of guarded curves now: "
        << guardedCurves.size() << endl;

    // Reindex once more
    for (uint i = 0; i < guardedCurves.size(); i++)
    {
        guardedCurves[i].m_curve.m_id = i;

        if (guardedCurves[i].m_curve.irregular2())
        {
            Logger::lout(Logger::FATAL)
                << "Guard splitting led to irregular curves, aborting..." << endl;
            save_log_file(-1, filename);
            return false;
        }
    }

    Logger::lout(Logger::INFO) << "Starting polygon triangulation" << endl;
    TriangulationConstraints<Scalar> constraints;
    for (uint i = 0; i < guardedCurves.size(); i++)
    {
        guardedCurves[i].m_curve.m_id = i;
        constraints.addQuad(guardedCurves[i]);
    }
    BezierMesh<Scalar> mesh(degree);
    if (!performCDT(constraints, mesh, addBbox))
    {
        save_log_file(10, filename);
        Logger::lout(Logger::FATAL) << "Polygon triangulation failed. There probably were envelopes that still "
                                       "intersected. Aborting..."
                                    << endl;
        return false;
    }
    Logger::lout(Logger::INFO) << "Finished polygon triangulation successfully."
                               << endl;


    //Save number of triangles to log
    save_log_file(21, filename, mesh.triangle_count());


    Logger::lout(Logger::INFO) << "Computing regular Bezier control points..." << endl;
    if (!insertMeshCtrlpts(guardedCurves, mesh))
    {
        save_log_file(11, filename);
        Logger::lout(Logger::FATAL)
            << "Computing Bezier control points failed, aborting..." << endl;
        return false;
    }
    Logger::lout(Logger::INFO) << "Computing Bezier control points done."
                               << endl;

    save_log_file(19, filename);

    Logger::lout(Logger::INFO) << "Checking resulting mesh validity..." << endl;
    bool uncertain = false;
    if (!checkMesh(mesh, uncertain))
    {
        save_log_file(4, filename);
        Logger::lout(Logger::FATAL) << "There were numerically irregular Bezier triangles in the mesh. "
                                       "Aborting..."
                                    << endl;
        return false;
    }
    Logger::lout(Logger::INFO) << "Mesh is valid."
                               << endl;
    if(uncertain)
        save_log_file(3, filename);
    else
        save_log_file(5, filename);

    string filenameBase = filename.substr(0, filename.find_last_of('.'));
    if (addBbox)
        filenameBase += "_boxed";
    if (epsOutput)
    {
        Logger::lout(Logger::INFO) << "Writing visualization output to .eps-file..." << endl;
        EPSWriter<Scalar>().write(filenameBase + ".eps", mesh);
    }
    if (gmshOutput)
    {
        Logger::lout(Logger::INFO) << "Writing final output to .msh-file..." << endl;
        MSHWriter<Scalar>().write(filenameBase + ".msh", mesh);
    }

    if(optimizeMesh)
    {
        filenameBase += "_opt";
        Logger::lout(Logger::INFO) << "Relaxing the triangulation (may take a while)..." << endl;
        mesh.optimization_conformal();

        if (epsOutput)
        {
            Logger::lout(Logger::INFO) << "Writing visualization output to .eps-file..." << endl;
            EPSWriter<Scalar>().write(filenameBase + ".eps", mesh);
        }
        if (gmshOutput)
        {
            Logger::lout(Logger::INFO) << "Writing final output to .msh-file..." << endl;
            MSHWriter<Scalar>().write(filenameBase + ".msh", mesh);
        }
    }
    Logger::lout(Logger::INFO) << "Finished successfully." << endl;
    return true;
}

} // namespace bzmsh
