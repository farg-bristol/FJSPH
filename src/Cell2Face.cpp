#include <algorithm>
#include <omp.h>
#include <chrono>
#include <set>
#include <cmath>
#include <assert.h>
#include <TECIO.h>
#include "Convert.h"

using namespace std::chrono;
using namespace netCDF;

/*Define Simulation Dimension*/
#ifndef SIMDIM
#define SIMDIM 3
#endif

/****** Eigen vector definitions ************/
typedef std::array<real,3> StateVecD;
typedef std::array<int,3> StateVecI;

const std::string WHITESPACE = " \n\r\t\f\v";

struct FACE
{
    FACE():left(-1),right(-1),zone(-1),nodes({-1,-1,-1,-1}) {};

    inline void to_array(std::array<int,7>& face_array)
    {
        face_array[0] = left;
        face_array[1] = right;
        face_array[2] = zone;
        face_array[3] = nodes[0];
        face_array[4] = nodes[1];
        face_array[5] = nodes[2];
        face_array[6] = nodes[3];
    }

    inline void from_array(std::array<int,7> const& face_array)
    {
        left     = face_array[0];
        right    = face_array[1];
        zone     = face_array[2];
        nodes[0] = face_array[3];
        nodes[1] = face_array[4];
        nodes[2] = face_array[5];
        nodes[3] = face_array[6];
    }

    inline int get_node(std::array<int,7> const& face_array, int index)
    {
        return face_array[index+3];
    }

    int left, right;
    int zone;
    std::array<int,4> nodes;
};


struct ProgressBar
{
    ProgressBar():percent(0),bar_length(40)
    {   // Start with an empty bar
        string spaces = string((bar_length),' ');
        cout << "\rPercent: [" << spaces << "] " << int(round(percent * 100)) 
                << "%" << std::flush;};

    void reset()
    {
        percent = 0;
        string spaces = string((bar_length),' ');
        cout << "\rPercent: [" << spaces << "] " << int(round(percent * 100)) 
                << "%" << std::flush;
    }

    void update(float percent_)
    {
        int new_percent = round(percent_ * 100);
        if (new_percent > round(percent*100))
        {
            percent = percent_;
            int bar_fill = round(percent * bar_length);
            string hashes = string(bar_fill, '#');
            string spaces = string((bar_length - hashes.size()),' ');
            cout << "\rPercent: [" << hashes << spaces << "] " << int(round(percent * 100)) 
                << "%" << std::flush;
        }
    }

    float percent;
    int bar_length;
};


namespace TAUtoFace
{
    inline void reorder(std::array<int,4>& nodes)
    {
        std::array<int,4> reordered{{-1,-1,-1,-1}};

        int num_nodes = 4;
        if (nodes[3] == -1) 
            num_nodes = 3;

        int min_index = 0;
        int min_val = nodes[0];

        for (int ii = 1; ii < num_nodes; ii++)
        {
            if (nodes[ii] < min_val)
            {
               min_index = ii;
               min_val = nodes[ii];
            }
        }
          
        assert (min_val >= 0);

        int jj = 0;
        for (int ii = min_index; ii < num_nodes; ii++)
        {
            reordered[jj] = nodes[ii];
            jj++;
        }

        for (int ii = 0; ii < min_index; ii++)
        {
            reordered[jj] = nodes[ii];
            jj++;
        }

        for (int ii = 0; ii < num_nodes; ii++)
            nodes[ii] = reordered[ii];
    }

    inline void update_face_array(FACE& f, int& face_count, std::vector<std::array<int,7>>& faces)
    {
        f.to_array(faces[face_count]);
        face_count++;
    }


    size_t tetrahedra(int& meshID, size_t& num_cells, std::vector<std::array<int,7>>& faces)
    {
        /* 
         * Tetrahedron
         *
         *          3
         *         /|\
         *        / | \
         *       /  |  \
         *      0---|---2
         *        \ | /
         *          1
         *
         */

        cout << "  Tetrahedras" << endl; 
        // ProgressBar progress;
        int retval;
        /*Retrieve the cdata*/
        int num_tets_dim, num_points_per_tet_dim;
        size_t num_tets, num_points_per_tet;
        vector<vector<size_t>> points_of_tets; 
        if ((retval = nc_inq_dimid(meshID, "no_of_tetraeders", &num_tets_dim)))
        {
            cout << "\rNo tetraeders in mesh file." << std::flush;
            num_tets = 0;
        }
        else
        {
            Get_Dimension(meshID, "no_of_tetraeders", num_tets_dim, num_tets);
            Get_Dimension(meshID, "points_per_tetraeder", num_points_per_tet_dim, num_points_per_tet);
            points_of_tets = Get_Element(meshID, "points_of_tetraeders",num_tets,num_points_per_tet);
        }
        
        vector<std::array<int,7>> tets(num_tets*4);

        #pragma omp parallel for
        for(size_t ii = 0; ii < num_tets; ii++)
        {
            FACE f;
            f.left = num_cells + ii;

            /* Face 1 */
            f.nodes[0] = points_of_tets[ii][0];
            f.nodes[1] = points_of_tets[ii][2];
            f.nodes[2] = points_of_tets[ii][1];
            reorder(f.nodes);
            f.to_array(tets[ii*4]);

            /* Face 2 */            
            f.nodes[0] = points_of_tets[ii][0];
            f.nodes[1] = points_of_tets[ii][1];
            f.nodes[2] = points_of_tets[ii][3];
            reorder(f.nodes);
            f.to_array(tets[ii*4+1]);

            /* Face 3 */            
            f.nodes[0] = points_of_tets[ii][1];
            f.nodes[1] = points_of_tets[ii][2];
            f.nodes[2] = points_of_tets[ii][3];
            reorder(f.nodes);
            f.to_array(tets[ii*4+2]);

            /* Face 4 */            
            f.nodes[0] = points_of_tets[ii][2];
            f.nodes[1] = points_of_tets[ii][0];
            f.nodes[2] = points_of_tets[ii][3];
            reorder(f.nodes);
            f.to_array(tets[ii*4+3]);

            // progress.update(float(ii)/num_tets);
        }
        num_cells += num_tets;
        faces.insert(faces.end(), tets.begin(), tets.end());
        return num_tets;
    }

    size_t hexahedra(int& meshID, size_t& num_cells, std::vector<std::array<int,7>>& faces)
    {
        /* 
         * Hexahedron
         *                             faces:
         *    4-------7         0----------------3
         *    |\      |\        |\              /|
         *    | \     | \       | \      3     / |
         *    |  5-------6      |  \4--------7/  |
         *    |  |    |  |      |   |        |   |
         *    0--|----3  |      | 5 |   4    | 2 | 0
         *     \ |     \ |      |   |        |   |
         *      \|      \|      |  /5--------6\  |
         *       1-------2      | /     1      \ |
         *                      |/              \|
         *                      1----------------2
         *
         */

        cout << "  Hexahedras" << endl;
        // ProgressBar progress;
        int retval;
        /*Retrieve the cdata*/
        int num_hexs_dim, num_points_per_hex_dim;
        size_t num_hexs, num_points_per_hex;
        vector<vector<size_t>> points_of_hexs; 
        if ((retval = nc_inq_dimid(meshID, "no_of_hexaeders", &num_hexs_dim)))
        {
            cout << "\rNo hexaeders in mesh file." << std::flush;
            num_hexs = 0;
        }
        else
        {
            Get_Dimension(meshID, "no_of_hexaeders", num_hexs_dim, num_hexs);
            Get_Dimension(meshID, "points_per_hexaeder", num_points_per_hex_dim, num_points_per_hex);
            points_of_hexs = Get_Element(meshID, "points_of_hexaeders",num_hexs,num_points_per_hex);
        }
        
        vector<std::array<int,7>> hexs(num_hexs*6);

        #pragma omp parallel for
        for(size_t ii = 0; ii < num_hexs; ii++)
        {
            FACE f;
            f.left = num_cells + ii;

            /* Face 1 */
            f.nodes[0] = points_of_hexs[ii][0];
            f.nodes[1] = points_of_hexs[ii][3];
            f.nodes[2] = points_of_hexs[ii][2];
            f.nodes[3] = points_of_hexs[ii][1];
            reorder(f.nodes);
            f.to_array(hexs[ii*6]);

            /* Face 2 */            
            f.nodes[0] = points_of_hexs[ii][0];
            f.nodes[1] = points_of_hexs[ii][1];
            f.nodes[2] = points_of_hexs[ii][5];
            f.nodes[3] = points_of_hexs[ii][4];
            reorder(f.nodes);
            f.to_array(hexs[ii*6+1]);

            /* Face 3 */            
            f.nodes[0] = points_of_hexs[ii][1];
            f.nodes[1] = points_of_hexs[ii][2];
            f.nodes[2] = points_of_hexs[ii][6];
            f.nodes[3] = points_of_hexs[ii][5];
            reorder(f.nodes);
            f.to_array(hexs[ii*6+2]);

            /* Face 4 */            
            f.nodes[0] = points_of_hexs[ii][2];
            f.nodes[1] = points_of_hexs[ii][3];
            f.nodes[2] = points_of_hexs[ii][7];
            f.nodes[3] = points_of_hexs[ii][6];
            reorder(f.nodes);
            f.to_array(hexs[ii*6+3]);

            /* Face 5 */            
            f.nodes[0] = points_of_hexs[ii][0];
            f.nodes[1] = points_of_hexs[ii][4];
            f.nodes[2] = points_of_hexs[ii][7];
            f.nodes[3] = points_of_hexs[ii][3];
            reorder(f.nodes);
            f.to_array(hexs[ii*6+4]);

            /* Face 6 */            
            f.nodes[0] = points_of_hexs[ii][4];
            f.nodes[1] = points_of_hexs[ii][5];
            f.nodes[2] = points_of_hexs[ii][6];
            f.nodes[3] = points_of_hexs[ii][7];
            reorder(f.nodes);
            f.to_array(hexs[ii*6+5]);
            
            // progress.update(float(ii)/num_hexs);
        }
        num_cells += num_hexs;
        
        faces.insert(faces.end(), hexs.begin(), hexs.end());
        return num_hexs;
    }
    
    size_t pyramid(int& meshID, size_t& num_cells, std::vector<std::array<int,7>>& faces)
    {
        /* 
         * Pyramid                   faces
         *            4
         *          /|            0-------3
         *         / |            |\     /|
         *        / |.|           | \ 3 / |
         *       /  |.|           |  \ /  |
         *      /  |  .|          | 4 4 2 | 0
         *     /   |  .|          |  / \  |
         *    0---|---3 |         | /   \ |
         *     \  |    \|         |/  1  \|
         *      \|      \|        1-------2
         *       1-------2
         *
         */

        cout << "  Pyramids" << endl;
        // ProgressBar progress;
        int retval;
        /*Retrieve the cdata*/
        int num_pyrs_dim, num_points_per_pyr_dim;
        size_t num_pyrs, num_points_per_pyr;
        vector<vector<size_t>> points_of_pyrs; 
        if ((retval = nc_inq_dimid(meshID, "no_of_pyramids", &num_pyrs_dim)))
        {
            cout << "\rNo pyramids in mesh file." << std::flush;
            num_pyrs = 0;
        }
        else
        {
            Get_Dimension(meshID, "no_of_pyramids", num_pyrs_dim, num_pyrs);
            Get_Dimension(meshID, "points_per_pyramid", num_points_per_pyr_dim, num_points_per_pyr);
            points_of_pyrs = Get_Element(meshID, "points_of_pyramids",num_pyrs,num_points_per_pyr);
        }
        
        vector<std::array<int,7>> pyrs(num_pyrs*5);

        #pragma omp parallel for
        for(size_t ii = 0; ii < num_pyrs; ii++)
        {
            FACE f;
            f.left = num_cells + ii;

            /* Face 1 */
            f.nodes[0] = points_of_pyrs[ii][0];
            f.nodes[1] = points_of_pyrs[ii][4];
            f.nodes[2] = points_of_pyrs[ii][3];
            reorder(f.nodes);
            f.to_array(pyrs[ii*5]);

            /* Face 2 */            
            f.nodes[0] = points_of_pyrs[ii][0];
            f.nodes[1] = points_of_pyrs[ii][1];
            f.nodes[2] = points_of_pyrs[ii][4];
            reorder(f.nodes);
            f.to_array(pyrs[ii*5+1]);

            /* Face 3 */            
            f.nodes[0] = points_of_pyrs[ii][1];
            f.nodes[1] = points_of_pyrs[ii][2];
            f.nodes[2] = points_of_pyrs[ii][4];
            reorder(f.nodes);
            f.to_array(pyrs[ii*5+2]);

            /* Face 4 */            
            f.nodes[0] = points_of_pyrs[ii][2];
            f.nodes[1] = points_of_pyrs[ii][3];
            f.nodes[2] = points_of_pyrs[ii][4];
            reorder(f.nodes);
            f.to_array(pyrs[ii*5+3]);

            /* Face 5 */            
            f.nodes[0] = points_of_pyrs[ii][0];
            f.nodes[1] = points_of_pyrs[ii][3];
            f.nodes[2] = points_of_pyrs[ii][2];
            f.nodes[3] = points_of_pyrs[ii][1];
            reorder(f.nodes);
            f.to_array(pyrs[ii*5+4]);

            
            // progress.update(float(ii)/num_pyrs);
        }
        num_cells += num_pyrs;
        faces.insert(faces.end(), pyrs.begin(), pyrs.end());
        return num_pyrs;
    }

    
    size_t prism(int& meshID, size_t& num_cells, std::vector<std::array<int,7>>& faces)
    {
        /* 
         * Prism                  face schematic
         *                               0
         *     3-------5          0-------------2
         *     |\     /|           \     4     /
         *     | \  /  |             3-------5
         *     |  4    |              \  1  /
         *     |  |    |            2  \  /   3
         *     0--|----2                4
         *      \ |   /                 |
         *       \| /                   1
         *        1
         *
         */
        
        cout << "  Prisms" << endl;
        // ProgressBar progress;
        int retval;
        /*Retrieve the cdata*/
        int num_pris_dim, num_points_per_pri_dim;
        size_t num_pris, num_points_per_pri;
        vector<vector<size_t>> points_of_pris; 
        if ((retval = nc_inq_dimid(meshID, "no_of_prisms", &num_pris_dim)))
        {
            cout << "\rNo prisms in mesh file." << std::flush;
            num_pris = 0;
        }
        else
        {
            Get_Dimension(meshID, "no_of_prisms", num_pris_dim, num_pris);
            Get_Dimension(meshID, "points_per_prism", num_points_per_pri_dim, num_points_per_pri);
            points_of_pris = Get_Element(meshID, "points_of_prisms",num_pris,num_points_per_pri);
        }
        
        vector<std::array<int,7>> pris(num_pris*5);
        // size_t cellcount = 0;
        #pragma omp parallel for
        for(size_t ii = 0; ii < num_pris; ii++)
        {
            FACE f;
            f.left = num_cells + ii;

            /* Face 1 */
            f.nodes[0] = points_of_pris[ii][0];
            f.nodes[1] = points_of_pris[ii][2];
            f.nodes[2] = points_of_pris[ii][1];
            reorder(f.nodes);
            f.to_array(pris[ii*5]);

            /* Face 2 */            
            f.nodes[0] = points_of_pris[ii][3];
            f.nodes[1] = points_of_pris[ii][4];
            f.nodes[2] = points_of_pris[ii][5];
            reorder(f.nodes);
            f.to_array(pris[ii*5+1]);

            /* Face 3 */            
            f.nodes[0] = points_of_pris[ii][0];
            f.nodes[1] = points_of_pris[ii][1];
            f.nodes[2] = points_of_pris[ii][4];
            f.nodes[3] = points_of_pris[ii][3];
            reorder(f.nodes);
            f.to_array(pris[ii*5+2]);

            /* Face 4 */            
            f.nodes[0] = points_of_pris[ii][1];
            f.nodes[1] = points_of_pris[ii][2];
            f.nodes[2] = points_of_pris[ii][5];
            f.nodes[3] = points_of_pris[ii][4];
            reorder(f.nodes);
            f.to_array(pris[ii*5+3]);

            /* Face 5 */            
            f.nodes[0] = points_of_pris[ii][2];
            f.nodes[1] = points_of_pris[ii][0];
            f.nodes[2] = points_of_pris[ii][3];
            f.nodes[3] = points_of_pris[ii][5];
            reorder(f.nodes);
            f.to_array(pris[ii*5+4]);

            // #pragma omp atomic
            //     cellcount++;

            
                // progress.update(float(cellcount)/num_pris);
        }
        num_cells += num_pris;
        
        faces.insert(faces.end(), pris.begin(), pris.end());
        return num_pris;
    }

    /* Surface faces */
    size_t surface_quads(int& meshID, std::vector<std::array<int,7>>& faces)
    {
        
        cout << "  Surface quads" << endl;
        // ProgressBar progress;
        int retval;
        /*Retrieve the cdata*/
        int num_quad_dim, num_points_per_quad_dim;
        size_t num_quad, num_points_per_quad;
        vector<vector<size_t>> points_of_quad; 
        if ((retval = nc_inq_dimid(meshID, "no_of_surfacequadrilaterals", &num_quad_dim)))
        {
            cout << "No quads in mesh file." << endl;
            num_quad = 0;
        }
        else
        {
            Get_Dimension(meshID, "no_of_surfacequadrilaterals", num_quad_dim, num_quad);
            Get_Dimension(meshID, "points_per_surfacequadrilateral", num_points_per_quad_dim, num_points_per_quad);
            points_of_quad = Get_Element(meshID, "points_of_surfacequadrilaterals",num_quad,num_points_per_quad);
        }

        vector<std::array<int,7>> quads(num_quad);

        #pragma omp parallel for
        for(size_t ii = 0; ii < num_quad; ii++)
        {
            FACE f;

            /* Face 1 */
            f.nodes[0] = points_of_quad[ii][3];
            f.nodes[1] = points_of_quad[ii][2];
            f.nodes[2] = points_of_quad[ii][1];
            f.nodes[3] = points_of_quad[ii][0];
            reorder(f.nodes);
            f.to_array(quads[ii]);

            // progress.update(float(ii)/num_quad);
        }
        faces.insert(faces.end(), quads.begin(), quads.end());
        return num_quad;
    }

    size_t surface_trigs(int& meshID, std::vector<std::array<int,7>>& faces)
    {
        
        cout << "  Surface trigs" << endl;
        // ProgressBar progress;
        int retval;
        /*Retrieve the cdata*/
        int num_trig_dim, num_points_per_trig_dim;
        size_t num_trig, num_points_per_trig;
        vector<vector<size_t>> points_of_trig; 
        if ((retval = nc_inq_dimid(meshID, "no_of_prisms", &num_trig_dim)))
        {
            cout << "No trigs in mesh file." << endl;
            num_trig = 0;
        }
        else
        {
            Get_Dimension(meshID, "no_of_surfacetriangles", num_trig_dim, num_trig);
            Get_Dimension(meshID, "points_per_surfacetriangle", num_points_per_trig_dim, num_points_per_trig);
            points_of_trig = Get_Element(meshID, "points_of_surfacetriangles",num_trig,num_points_per_trig);
        }

        vector<std::array<int,7>> trigs(num_trig);

        #pragma omp parallel for
        for(size_t ii = 0; ii < num_trig; ii++)
        {
            FACE f;

            /* Face 1 */
            f.nodes[0] = points_of_trig[ii][2];
            f.nodes[1] = points_of_trig[ii][1];
            f.nodes[2] = points_of_trig[ii][0];
            reorder(f.nodes);
            f.to_array(trigs[ii]);

            // progress.update(float(ii)/num_trig);
        }
        faces.insert(faces.end(), trigs.begin(), trigs.end());
        return num_trig;
    }

    void write_points(int& fin, int& fout, size_t const& nPnts, 
        int& ptsxID, int& ptsyID, int& ptszID)
    {
        // # Write points
        cout << "Writing Points" << endl;

        /*Put them in the file*/
        double *coord = Get_Real_Scalar(fin, "points_xc", nPnts);
        Write_Variable_Scalar(fout, ptsxID, coord, "x coordinate");
            
        coord = Get_Real_Scalar(fin, "points_yc", nPnts);
        Write_Variable_Scalar(fout, ptsyID, coord, "y coordinate");

        coord = Get_Real_Scalar(fin, "points_zc", nPnts);
        Write_Variable_Scalar(fout, ptszID, coord, "z coordinate");

        cout << "Completed writing points" << endl;
    }

    // Face conversion for my own cell containment queries. Written in NetCDF as 
    // basically a face based TAU format.
    void Write_Face_Data(int& fin, string const& meshIn, string const& meshOut, 
                        std::vector<std::array<int,7>> const& faces,
                        size_t const& nPnts, size_t const& nElem, size_t const& nFace,
                        size_t const& nSurf, int* faceMarkers)
    { 
        #ifdef DEBUG
        dbout << "Entering Write_Face_Data..." << endl;
        dbout << "Output file: " << meshOut << endl;
        #endif

        cout << "Attempting to write output file." << endl;
        cout << "File: " << meshOut << endl;
        
        /* Create the file. */
        int retval;
        int fout;
        if ((retval = nc_create(meshOut.c_str(), NC_CLOBBER | NC_64BIT_OFFSET, &fout)))
        {	
            cout << "Error: Failed to open file \"" << meshOut << "\" when outputting." << endl;
            ERR(retval);
            exit(-1);
        }

        if ((retval = nc_put_att_text(fout, NC_GLOBAL, "converted_mesh_filename", meshIn.length(), meshIn.c_str())))
        {
            cout << "Error: Failed to attach filename attribute" << endl;
            ERR(retval);
            exit(-1);
        }

        // Separate faces to write
        size_t ntFaces = 0;
        size_t nqFaces = 0;
        int* tfaces;
        int* qfaces;
        int* left;
        int* right;
        // int* zone;

        // Find out how many of each face type exist
        #pragma omp parallel for
        for(size_t ii = 0; ii < nFace; ++ii)
        {
            if(faces[ii][6] == -1)
            {   
                #pragma omp atomic
                ntFaces++;
            }
            else
            {
                #pragma omp atomic
                nqFaces++;            
            }
        }

        // Now allocate face arrays to have data.
        tfaces = new int[ntFaces*3];
        qfaces = new int[nqFaces*4];
        left = new int[nFace];
        right = new int[nFace];

        // Go back through faces, placing them in the arrays
        size_t tFace = 0;
        size_t qFace = 0;
        for(size_t ii = 0; ii < nFace; ++ii)
        {
            if(faces[ii][6] == -1)
            {   
                left[tFace]       = faces[ii][0];
                right[tFace]      = faces[ii][1];
                // zone[tFace]       = faces[ii][2];
                tfaces[tFace*3]   = faces[ii][3];
                tfaces[tFace*3+1] = faces[ii][4];
                tfaces[tFace*3+2] = faces[ii][5];
                tFace++;
            }
            else
            { 
                left[ntFaces + qFace]  = faces[ii][0];
                right[ntFaces + qFace] = faces[ii][1];
                // zone[ntFace + qFace]  = faces[ii][2];
                qfaces[qFace*4]       = faces[ii][3];
                qfaces[qFace*4+1]     = faces[ii][4];
                qfaces[qFace*4+2]     = faces[ii][5];
                qfaces[qFace*4+3]     = faces[ii][6];
                qFace++;            
            }
        }

        /* Define the dimensions. */
        int nElemID, nFaceID, nTFaceID, nPpTFcID, nQFaceID, nPpQFcID, nMarkID, nPntsID;

        Define_Dimension(fout, "no_of_elements", nElem, &nElemID);
        Define_Dimension(fout, "no_of_faces", nFace, &nFaceID);

        if(ntFaces != 0)
        {
            Define_Dimension(fout, "no_of_triangles", ntFaces, &nTFaceID);
            Define_Dimension(fout, "points_per_triangle", 3, &nPpTFcID);
        }

        if(nqFaces != 0)
        {
            Define_Dimension(fout, "no_of_quadrilaterals", nqFaces, &nQFaceID);
            Define_Dimension(fout, "points_per_quadrilateral", 4, &nPpQFcID);
        }

        Define_Dimension(fout, "no_of_surfaceelements", nSurf, &nMarkID);
        Define_Dimension(fout, "no_of_points", nPnts, &nPntsID);

        /* Define the variables */
        int faceTVarID, faceQVarID, leftVarID, rightVarID, markID, ptsxID, ptsyID, ptszID;
        
        /*Define the faces*/
        int dimIDs[] = {nTFaceID,nPpTFcID};

        if(ntFaces != 0)
        {
            Define_Variable(fout, "points_of_triangles", NC_INT, 2,
                                    dimIDs, &faceTVarID);
        }

        if(nqFaces != 0)
        {
            dimIDs[0] = nQFaceID;
            dimIDs[1] = nPpQFcID;
            Define_Variable(fout, "points_of_quadrilaterals", NC_INT, 2,
                                    dimIDs, &faceQVarID);
        }

        Define_Variable(fout, "left_element_of_faces", NC_INT, 1,
                                    &nFaceID, &leftVarID);
        Define_Variable(fout, "right_element_of_faces", NC_INT, 1,
                                &nFaceID, &rightVarID);
        Define_Variable(fout, "boundarymarker_of_surfaces", NC_INT, 1,
                                &nMarkID, &markID);


        /*Define the points*/
        Define_Variable(fout, "points_xc", NC_DOUBLE, 1, &nPntsID, &ptsxID);
        Define_Variable(fout, "points_yc", NC_DOUBLE, 1, &nPntsID, &ptsyID);
        Define_Variable(fout, "points_zc", NC_DOUBLE, 1, &nPntsID, &ptszID);

        /* End define mode. */
        if ((retval = nc_enddef(fout)))
            ERR(retval);

        /*Create the C array for the faces*/
        size_t start[] = {0,0};
        size_t end[] = {ntFaces,3};

        if(ntFaces != 0)
        {
            /*Put faces into the file*/
            cout << "Writing triangle faces" << endl;
            Write_Variable_Array(fout, faceTVarID, start, end, &tfaces[0], "Triangle faces");
        }

        if(nqFaces != 0)
        {
            /*Put faces into the file*/
            end[0] = nqFaces;
            end[1] = 4;
            cout << "Writing quadrilateral faces" << endl;
            Write_Variable_Array(fout, faceQVarID, start, end, &qfaces[0], "Quadrilateral faces");
        }

        /*Put face left and right into the file*/
        cout << "Writing left elements of faces" << endl;
        Write_Variable_Scalar(fout, leftVarID, &left[0], "Left elements");
        cout << "Writing right elements of faces" << endl;
        Write_Variable_Scalar(fout, rightVarID, &right[0], "Right elements");
        cout << "Writing surface markers" << endl;
        Write_Variable_Scalar(fout, markID, &faceMarkers[0], "Surface markers");

        write_points(fin,fout,nPnts,ptsxID,ptsyID,ptszID);
        
        
        if ((retval = nc_close(fout)))
        {	
            cout << "Error: Failed to close file \"" << meshOut << "\" when outputting." << endl;
            ERR(retval);
            exit(-1);
        }
    }

    void Write_Tecplot_Binary(int& fin, string const& meshOut, 
                        std::vector<std::array<int,7>> const& faces,
                        size_t const& nPnts, size_t const& nElem, size_t const& nFace)
    {
        string meshOut_ = meshOut + ".plt";
        int32_t debug = 0;
        #ifdef DEBUG
        dbout << "Entering Write_ASCII_Face_Data..." << endl;
        cout << "Attempting write output file." << endl;
        cout << "File: " << meshOut_ << endl;
        debug = 1;
        #endif

        string varNames = "X,Y,Z";
        string title = "3D face based grid conversion output";
        string group = "Mesh data";
        int32_t fileType = 0;
        int32_t fileFormat = 0;
        int32_t fileIsDouble = 1;

        if(TECINI142(title.c_str(),varNames.c_str(),meshOut_.c_str(),".",
                    &fileFormat,&fileType,&debug,&fileIsDouble))
        {
            cout << "Failed to open " << meshOut_ << endl;
            exit(-1);
        }
        
        // Separate faces to write
        size_t ntFaces = 0;
        size_t nqFaces = 0;
        
        // Find out how many of each face type exist
        #pragma omp parallel for reduction(+:ntFaces,nqFaces)
        for(size_t ii = 0; ii < nFace; ++ii)
        {
            if(faces[ii][6] == -1)
            {   
                ntFaces++;
            }
            else
            {
                nqFaces++;            
            }
        }

        size_t TotalNumFaceNodes = ntFaces*3 + nqFaces*4;
        // Now allocate face arrays to have data.        
        int32_t* left = new int32_t[nFace];
        int32_t* right = new int32_t[nFace];
        int32_t* faceNodes = new int32_t[TotalNumFaceNodes];
        int32_t* faceNodeCount = new int32_t[nFace];
        // Go back through faces, placing them in the arrays
        size_t tFace = 0;
        size_t qFace = 0;
        for(size_t ii = 0; ii < nFace; ++ii)
        {
            if(faces[ii][6] == -1)
            {   
                left[tFace]            = faces[ii][0]+1;
                if(faces[ii][1] < 0)
                    right[tFace] = 0;
                else
                    right[tFace]       = faces[ii][1]+1;

                faceNodeCount[tFace] = 3;
                faceNodes[tFace*3]   = faces[ii][3]+1;
                faceNodes[tFace*3+1] = faces[ii][4]+1;
                faceNodes[tFace*3+2] = faces[ii][5]+1;
                tFace++;
            }
            else
            { 
                left[ntFaces + qFace]  = faces[ii][0]+1;
                if(faces[ii][1] < 0)
                    right[ntFaces + qFace] = 0;
                else
                    right[ntFaces + qFace] = faces[ii][1]+1;

                faceNodeCount[ntFaces + qFace] = 4;
                faceNodes[ntFaces*3 + qFace*4]   = faces[ii][3]+1;
                faceNodes[ntFaces*3 + qFace*4+1] = faces[ii][4]+1;
                faceNodes[ntFaces*3 + qFace*4+2] = faces[ii][5]+1;
                faceNodes[ntFaces*3 + qFace*4+3] = faces[ii][6]+1;
                qFace++;            
            }
        }


        int32_t zoneType = 7; /* FE Polyhedron */	
        int32_t numNodes = nPnts;
        int32_t numElems = nElem;
        int64_t numFaces = nFace;

        /* Need to find how many face nodes there are */
        int64_t numFaceNodes = TotalNumFaceNodes;
        double  solTime = 0.0;
        int32_t strandID = 0;     // Static Zone
        int32_t unused = 0;       // ParentZone is no longer used
        int32_t numBoundaryFaces = 0;
        int32_t numBoundaryConnections = 0;
        int32_t shareConnectivity = 0;

        if(TECPOLYZNE142(
            "Polyhedral Zone",
            &zoneType,
            &numNodes,
            &numElems,
            &numFaces,
            &numFaceNodes,
            &solTime,
            &strandID,
            &unused,
            &numBoundaryFaces,
            &numBoundaryConnections,
            NULL,
            NULL,  // All nodal variables
            NULL,
            &shareConnectivity))
        {
            printf("Polyhedron: error calling TECPOLYZNE\n");
            exit(-1);
        }

        int32_t nPnts_ = nPnts;
        double *coord = Get_Real_Scalar(fin, "points_xc", nPnts);
        TECDAT142(&nPnts_, coord, &fileIsDouble);
        
        coord = Get_Real_Scalar(fin, "points_yc", nPnts);
        TECDAT142(&nPnts_, coord, &fileIsDouble);
        
        coord = Get_Real_Scalar(fin, "points_zc", nPnts);
        TECDAT142(&nPnts_, coord, &fileIsDouble);

        delete[] coord;
        
        int32_t faceOffset = nFace;
        if(TECPOLYFACE142(
                    &faceOffset,
                    faceNodeCount,
                    &faceNodes[0],
                    &left[0],
                    &right[0]))
        {
            printf("Error calling TECPOLYFACE\n");
            exit(-1);
        }

        delete [] faceNodes;
        delete [] faceNodeCount;
        delete [] left;
        delete [] right;

        if(TECEND142())
        {
            printf("Polyquads: error calling TECEND\n");
            exit(-1);
        }

    }

    void Write_ASCII_Face_Data(int& fin, string const& meshOut, 
                        std::vector<std::array<int,7>> const& faces,
                        size_t const& nPnts, size_t const& nElem, size_t const& nFace)
    {
        string meshOut_ = meshOut + ".dat";
        #ifdef DEBUG
        dbout << "Entering Write_ASCII_Face_Data..." << endl;
        cout << "Attempting write output file." << endl;
        cout << "File: " << meshOut_ << endl;
        #endif
        std::ofstream fout(meshOut_,std::ios::out);
        if(!fout.is_open())
        {
            cout << "Couldn't open the output file." << endl;
            exit(-1);
        }

        // Separate faces to write
        size_t ntFaces = 0;
        size_t nqFaces = 0;
        int* tfaces;
        int* qfaces;
        int* left;
        int* right;
        
        // Find out how many of each face type exist
        #pragma omp parallel for reduction(+:ntFaces,nqFaces)
        for(size_t ii = 0; ii < nFace; ++ii)
        {
            if(faces[ii][6] == -1)
            {   
                ntFaces++;
            }
            else
            {
                nqFaces++;            
            }
        }

        // Now allocate face arrays to have data.
        tfaces = new int[ntFaces*3];
        qfaces = new int[nqFaces*4];
        left = new int[nFace];
        right = new int[nFace];

        // Go back through faces, placing them in the arrays
        size_t tFace = 0;
        size_t qFace = 0;
        for(size_t ii = 0; ii < nFace; ++ii)
        {
            if(faces[ii][6] == -1)
            {   
                left[tFace]       = faces[ii][0];
                right[tFace]      = faces[ii][1];
                // zone[tFace]       = faces[ii][2];
                tfaces[tFace*3]   = faces[ii][3];
                tfaces[tFace*3+1] = faces[ii][4];
                tfaces[tFace*3+2] = faces[ii][5];
                tFace++;
            }
            else
            { 
                left[ntFaces + qFace]  = faces[ii][0];
                right[ntFaces + qFace] = faces[ii][1];
                // zone[ntFace + qFace]  = faces[ii][2];
                qfaces[qFace*4]       = faces[ii][3];
                qfaces[qFace*4+1]     = faces[ii][4];
                qfaces[qFace*4+2]     = faces[ii][5];
                qfaces[qFace*4+3]     = faces[ii][6];
                qFace++;            
            }
        }

        size_t TotalNumFaceNodes = ntFaces*3 + nqFaces*4;
        

        fout << "VARIABLES= \"X\" \"Y\" \"Z\" " << endl;
        fout << "ZONE T=\"FEPOLYHEDRON Test\"" << endl;
        fout << "ZONETYPE=FEPOLYHEDRON" << endl;
        fout << "NODES=" << nPnts << " ELEMENTS=" << nElem << " FACES=" << nFace << endl;
        fout << "TotalNumFaceNodes=" << TotalNumFaceNodes << endl;
        fout << "NumConnectedBoundaryFaces=0 TotalNumBoundaryConnections=0" << endl;

        size_t w = 15;
        size_t preci = 6;
        fout << std::left << std::scientific << std::setprecision(preci);
        // fout << fdata.numElem << " " << fdata.numFaces << "  " <<  fdata.numPoint << endl;
        
        /*Write vertices in block format (Each dimension in turn)*/
        /*Get the coordinate data*/
        double *coord = Get_Real_Scalar(fin, "points_xc", nPnts);
        
        size_t newl = 0;
        fout << std::setw(1);
        for(size_t ii = 0; ii < nPnts; ++ii)
        {
            fout << std::setw(w) << coord[ii];
            newl++;

            if(newl>4)
            {
                fout << endl;
                fout << " ";
                newl=0;
            }
        }
        
        coord = Get_Real_Scalar(fin, "points_yc", nPnts);

        newl = 0;
        for(size_t ii = 0; ii < nPnts; ++ii)
        {
            fout << std::setw(w) << coord[ii];
            newl++;

            if(newl>4)
            {
                fout << endl;
                fout << " ";
                newl=0;
            }
        }
                
        coord = Get_Real_Scalar(fin, "points_zc", nPnts);

        newl = 0;
        for(size_t ii = 0; ii < nPnts; ++ii)
        {
            fout << std::setw(w) << coord[ii];
            newl++;

            if(newl>4)
            {
                fout << endl;
                fout << " ";
                newl=0;
            }
        }
        
        fout << endl;
        

        fout << std::left << std::fixed;
        w = 9;
        /*Inform of how many vertices in each face*/
        fout << "#node count per face" << endl;
        newl = 0;
        for (size_t ii = 0; ii < ntFaces; ++ii)
        {
            fout << std::setw(w) << 3;
            newl++;

            if(newl>4)
            {
                fout << endl;
                newl=0;
            }
        }

        for (size_t ii = 0; ii < nqFaces; ++ii)
        {
            fout << std::setw(w) << 4;
            newl++;

            if(newl>4)
            {
                fout << endl;
                newl=0;
            }
        }
        fout << endl;

        /*Write the face data*/
        fout << "#face nodes" << endl;
        for (size_t ii = 0; ii < ntFaces; ++ii)
        {
            for(size_t jj = 0; jj < 3; ++jj)
            {	/*Write face vertex indexes*/
                fout << std::setw(w) << tfaces[ii*3+jj]+1;
                // if (vertex > fdata.nPnts)
                // {
                // 	cout << "Trying to write a vertex outside of the number of points." << endl;
                // }
            }
            fout << endl;
        }

        for (size_t ii = 0; ii < nqFaces; ++ii)
        {
            for(size_t jj = 0; jj < 4; ++jj)
            {	/*Write face vertex indexes*/
                fout << std::setw(w) << qfaces[ii*4+jj]+1;
                // if (vertex > fdata.nPnts)
                // {
                // 	cout << "Trying to write a vertex outside of the number of points." << endl;
                // }
            }
            fout << endl;
        }

        /*Write face left and right*/
        newl = 0;
        fout << "#left elements" << endl;
        for (size_t ii = 0; ii < nFace; ++ii)
        {
            fout << std::setw(w) << left[ii]+1 ;
            newl++;

            if(newl>4)
            {
                fout << endl;
                newl=0;
            }
        }
        fout << endl;

        fout << "#right elements" << endl;
        newl = 0;
        for (size_t ii = 0; ii < nFace; ++ii)
        {
            if(right[ii] < 0)
                fout<< std::setw(w) << 0;
            else
                fout << std::setw(w) << right[ii]+1;

            newl++;

            if(newl>4)
            {
                fout << endl;
                newl=0;
            }
        }
        fout << endl;
        fout.close();

        #ifdef DEBUG
        dbout << "Exiting Write_Face_Data..." << endl;
        #endif
    }
}

int main(int argc, char *argv[])
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
	high_resolution_clock::time_point t2;

    if (argc > 3) 
	{	/*Check number of input arguments*/
		cout << "\tWARNING: only a maximum of two input arguments accepted,\n";
		cout << "1: Input mesh file\n";
		cout << "2: Output mesh file\n";
		cout << "Other inputs will be ignored." << endl << endl;
	}

	if (argc == 1)
    {	/*Check if input has been provided*/
    	cout << "\tERROR: No inputs provided. Stopping... \n";
    	exit(-1);    	
    }

    string meshIn = argv[1];
    string meshOut;
    if(argc == 2)
    {
        meshOut = meshIn + ".faces";
    }
    else if (argc == 3)
    {
        meshOut = argv[2];
    }

	#ifdef DEBUG 
		cout << "Attempting read of NetCDF file." << endl;
		cout << "Mesh file: " << meshIn << endl;
		// cout << "Solution file: " << solIn << endl;
	#endif

	int meshID;
	int retval;

	if ((retval = nc_open(meshIn.c_str(), NC_NOWRITE, &meshID)))
	{	
		cout << "Failed to open mesh file \"" << meshIn << "\"" << endl;
		ERR(retval);
		exit(-1);
	}

	cout << "Mesh file open. Reading cell data..." << endl;

    int ptDimID, elemDimID, surfElemDimID;
	size_t nPnts, nElem, nSurfElem;

    // Retrieve how many elements there are.
	Get_Dimension(meshID, "no_of_elements", elemDimID, nElem);
	Get_Dimension(meshID, "no_of_points", ptDimID, nPnts);
    
	Get_Dimension(meshID, "no_of_surfaceelements", surfElemDimID, nSurfElem);

	cout << "nElem : " << nElem << " nPts: " << nPnts << endl;

    /* Get the faces from each cell type */
    std::vector<std::array<int,7>> faces;
    size_t num_cells = 0;

    TAUtoFace::surface_trigs(meshID,faces);
    TAUtoFace::surface_quads(meshID,faces);

    // get boundary markers
    int* faceMarkers = new int[nSurfElem];
	faceMarkers = Get_Int_Scalar(meshID, "boundarymarker_of_surfaces", nSurfElem);

    // Assign boundary markers
    std::set<int> zones;
    for(size_t ii = 0; ii < nSurfElem; ++ii)
    {
        faces[ii][2] = faceMarkers[ii];
        zones.insert(faceMarkers[ii]);
    }
    
    ProgressBar progress;
    progress.reset();
    TAUtoFace::tetrahedra(meshID, num_cells, faces);
    progress.update(float(num_cells) / nElem);
    TAUtoFace::prism(meshID, num_cells, faces);
    progress.update(float(num_cells) / nElem);
    TAUtoFace::pyramid(meshID, num_cells, faces);
    progress.update(float(num_cells) / nElem);
    TAUtoFace::hexahedra(meshID, num_cells, faces);
    progress.update(float(num_cells) / nElem);

    assert(num_cells == nElem);

    size_t nFace = faces.size()/2;
    
    cout << std::scientific << std::left << std::setprecision(4);
    t2 = high_resolution_clock::now();
	float duration = duration_cast<microseconds>(t2-t1).count()/1e6;
    cout << "\nCompleted ingesting TAU mesh. Time taken: " << duration << endl;

    /* Match faces */
    cout << "Matching faces" << endl;

    /* Build node to face pointer */
    vector<vector<int>> node_to_face(nPnts,vector<int>());
    for(size_t ii = 0; ii < faces.size(); ++ii)
        node_to_face[faces[ii][3]].emplace_back(ii);

    std::vector<std::array<int,7>> f_list, bf_list;

    int face_count = faces.size();
    progress.reset();
    while (face_count > 0)
    {
        bool found = false;

        FACE f;

        face_count--;
        f.from_array(faces[face_count]);
        // pop face

        int n = f.nodes[0];
        if (n != -1)
        {        
            for (int& nn:node_to_face[n])
            {
                if (nn != -1 && nn < face_count)
                {
                    // Potential neighbour face
                    FACE fn;
                    fn.from_array(faces[nn]);
                
                    if (f.nodes[3] != -1)// Quad
                    {
                        if (f.nodes[0] == fn.nodes[0])
                        {
                            if (f.nodes[1] == fn.nodes[3] && f.nodes[2] == fn.nodes[2] && f.nodes[3] == fn.nodes[1])
                            {    
                                found = true;
                                f.right = fn.left;
                                f.zone  = fn.zone;
                                fn.nodes[0] = -1;
                                fn.to_array(faces[nn]);
                                nn = -1;
                                break;
                            }
                        }
                    }
                    else if(fn.nodes[3] == -1) // Triangle /* (f.nodes[3] == -1)*/
                    {
                        if (f.nodes[0] == fn.nodes[0])
                        {
                            if (f.nodes[1] == fn.nodes[2] && f.nodes[2] == fn.nodes[1])
                            {       
                                found = true;
                                f.right = fn.left;
                                f.zone  = fn.zone;
                                fn.nodes[0] = -1;
                                fn.to_array(faces[nn]);
                                nn = -1;
                                break;
                            }
                        }
                    }
                }
            }

            if(!found)
            {
                cout << "Face count: " << face_count << endl;
                printf("Unmatched face: %d %d %d %d  zone: %d\n",f.nodes[0],f.nodes[1],f.nodes[2],f.nodes[3], f.zone);
                cout << "Neighbour faces: " << endl;
                for (int& nn:node_to_face[n])
                {
                    if(nn != -1)
                    {
                        FACE fn;
                        fn.from_array(faces[nn]);
                        printf("%d %d %d %d zone: %d\n",fn.nodes[0],fn.nodes[1],fn.nodes[2],fn.nodes[3],fn.zone);
                    }
                }
            }
            assert(found);

            if (f.zone == -1)
            {
                assert(f.right != -1);
                std::array<int,7> temp;
                f.to_array(temp);
                f_list.emplace_back(temp);
            }
            else
            {
                assert (f.right == -1);
                std::array<int,7> temp;
                f.to_array(temp);
                bf_list.emplace_back(temp);
            }
            // #if len(f_list)%100000 == 0:
            // #    print len(f_list)
        }

        progress.update(float(nFace-face_count/2) / nFace);
    }
    // Add boundary faces
        
    f_list.insert(f_list.end(),bf_list.begin(),bf_list.end());
        
    // Check we have extracted correct number of faces
    assert (f_list.size() == nFace);
    // num_cells = 0;
    // Extra checks and set halo cells
    for (size_t ii = 0; ii < nFace; ++ii)
    {
        assert(f_list[ii][0] != -1);
        if (f_list[ii][2] == -1) // if zone is a surface
            assert (f_list[ii][1] != -1);
        // if (f_list[ii][1] == -1) // right cell is not defined yet
        // {
        //     f_list[ii][1] = num_cells;
        //     num_cells += 1;
        // }
    } 
    // Check number of cells
    // assert (num_cells == nElem + nSurfElem);

    size_t num_faces = f_list.size();
    nFace = num_faces;
    // num_cells = nElem;

    high_resolution_clock::time_point t3 = high_resolution_clock::now();
	duration = duration_cast<microseconds>(t3-t2).count()/1e6;
    cout << "\nCompleted face matching. Time taken: " << duration << endl;
    
    cout << "Writing face based mesh" << endl;
    TAUtoFace::Write_Face_Data(meshID,meshIn,meshOut,f_list,
                nPnts,nElem,nFace,nSurfElem,faceMarkers);

    // TAUtoFace::Write_Tecplot_Binary(meshID,meshOut,f_list,
    //             nPnts,nElem,nFace);

    // TAUtoFace::Write_ASCII_Face_Data(meshID,meshOut,f_list,
    //             nPnts,nElem,nFace);

    high_resolution_clock::time_point t4 = high_resolution_clock::now();
	duration = duration_cast<microseconds>(t4-t3).count()/1e6;
    cout << "Completed writing of face based mesh. Time taken: " << duration << endl;

    cout << "All complete! Total time taken: " << duration_cast<microseconds>(t4-t1).count()/1e6 << endl;

    return 0;
}
