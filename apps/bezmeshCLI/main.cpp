#include <iostream>

#include "algorithm/bezier_meshing.hpp"
#include "datastructure/bezier_curve.hpp"
#include "fileio/json_reader.hpp"

using namespace bzmsh;
using namespace std;

void print_usage(string executable)
{
    cout << "Usage: specify path to json file, whether to add a bounding box (b), whether to "
            "output eps "
            "visualization (e) and whether to output mesh in gmsh file format (g) (as string of "
            "charflags)"
         << endl;
    cout << "Use this pattern: " + executable + " <path_to_json_curve_file> [optional flagstring]"
         << endl;
    cout << "Example: " + executable + " ../../resource/my_curves.json \"beg\"" << endl;
}

int main(int argc, char** argv)
{
    if (argc != 2 && argc != 3 && argc != 4)
    {
        print_usage(argv[0]);
        return -1;
    }

    string filename = argv[1];
    bool addBbox = false;
    bool epsOutput = false;
    bool gmshOutput = false;
    bool optimizeMesh = false;
    int exp_type = -1;
    if (argc >= 3)
    {
        string flagstring(argv[2]);
        if (flagstring.find('b') != flagstring.npos)
        {
            addBbox = true;
        }
        if (flagstring.find('e') != flagstring.npos)
        {
            epsOutput = true;
        }
        if (flagstring.find('g') != flagstring.npos)
        {
            gmshOutput = true;
        }
        if (flagstring.find('o') != flagstring.npos)
        {
            optimizeMesh = true;
            addBbox = true; //Default
        }
        
        if (argc >= 4)
        {
            sscanf(argv[3], "%d", &exp_type);
            printf("exp type = %d\n", exp_type);
        }
    }
    
    return bezierMeshing(filename, addBbox, epsOutput, gmshOutput, optimizeMesh, exp_type);
}
