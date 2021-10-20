/**
 * @brief Program for testing the containment of a point, rather than having to run the entire code.
 * 
 */

#include "Var.h"
#include "CDFIO.h"
#include "FOAMIO.h"
#include "Containment.h"
#include "Neighbours.h"

#include <chrono>
using namespace std::chrono;
using namespace nanoflann;

void read_params(int const& argc, char **argv, SIM &svar, StateVecD &testp)
{
    ifstream fin(argv[1]);

    if(!fin.is_open())
    {
        cout << "Could not open file \"" << argv[1] << "\"\n";
        exit(-1);
    }

    string line;
    while(getline(fin,line))
    {
        /* trim any comments away */
		size_t end = line.find_first_of('#');
		if(end != std::string::npos)
			line.substr(0,end+1);

        Get_String(line, "Primary grid face filename", svar.taumesh);
        Get_String(line, "Boundary mapping filename", svar.taubmap);
        Get_String(line, "OpenFOAM folder", svar.foamdir);
        Get_Number(line, "OpenFOAM binary (0/1)", svar.isBinary);
        Get_Number(line, "Label size (32/64)", svar.labelSize);
        Get_Number(line, "Scalar size (32/64)", svar.scalarSize);
        Get_Vector(line, "Test point", testp);
    }
    fin.close();

    if(svar.taumesh.empty())
    {
		
		if(svar.foamdir.empty())
		{
			cout << "No mesh file defined. Cannot continue." << endl;
            exit(-1);
		}
    }
    else if(svar.taubmap.empty())
    {
        cout << "Input TAU bmap file not defined." << endl;
        exit(-1);
    }
}

void test_containtment(MESH const& cells, Vec_Tree const& tree, StateVecD const& testp)
{
    cout << "Starting test at 5 cells..." << endl;
    size_t num_results = 5;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
	
    vector<size_t> ret_indexes(num_results);
    vector<real> out_dists_sqr(num_results);

    nanoflann::KNNResultSet<real> resultSet(num_results);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
    
    tree.index->findNeighbors(resultSet, &testp[0], nanoflann::SearchParams(10,0,true));

	high_resolution_clock::time_point t2 = high_resolution_clock::now();

    cout << "Time taken to get neighbours: " << 
        duration_cast<microseconds>(t2-t1).count() << " microseconds." << endl;

    uint pert = 0;
    size_t ii = 0;
    bool found = false;
    for( auto const& cell:ret_indexes)
    {
        if(CheckCell(cell,cells,testp,pert))
        {
            cout << "Successfully found the cell: " << cell << endl;
            cout << "Index in the search: " << ii << endl;
            found = true;
            break;
        }
        ++ii;
    }    

    high_resolution_clock::time_point t3 = high_resolution_clock::now();
    cout << "Time taken to search: " << duration_cast<microseconds>(t3-t2).count() << " microseconds." << endl;

    if(found == false)
    {
        num_results = 100;
        cout << "Expanding search to " << num_results << endl;

        t1 = high_resolution_clock::now();

        ret_indexes.resize(num_results);
        out_dists_sqr.resize(num_results);

        resultSet = nanoflann::KNNResultSet<real>(num_results);
        resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
        
        tree.index->findNeighbors(resultSet, &testp[0], nanoflann::SearchParams(10,0,true));

        t2 = high_resolution_clock::now();

        cout << "Time taken to get neighbours: " << 
        duration_cast<microseconds>(t2-t1).count() << " microseconds." << endl;

        ii = 0;
        for( auto const& cell:ret_indexes)
        {
            if(CheckCell(cell,cells,testp,pert))
            {
                cout << "Successfully found the cell: " << cell << endl;
                cout << "Index in the search: " << ii << endl;
                found = true;
                break;
            }
            ++ii;
        }    

        t3 = high_resolution_clock::now();
        cout << "Time taken to search: " << duration_cast<microseconds>(t3-t2).count() << " microseconds." << endl;

    }

    if(found == false)
    {
        num_results = 200;
        cout << "Expanding search to " << num_results << endl;

        t1 = high_resolution_clock::now();

        ret_indexes.resize(num_results);
        out_dists_sqr.resize(num_results);

        resultSet = nanoflann::KNNResultSet<real>(num_results);
        resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
        
        tree.index->findNeighbors(resultSet, &testp[0], nanoflann::SearchParams(10,0,true));

        t2 = high_resolution_clock::now();

        cout << "Time taken to get neighbours: " << 
        duration_cast<microseconds>(t2-t1).count() << " microseconds." << endl;

        ii = 0;
        for( auto const& cell:ret_indexes)
        {
            if(CheckCell(cell,cells,testp,pert))
            {
                cout << "Successfully found the cell: " << cell << endl;
                cout << "Index in the search: " << ii << endl;
                found = true;
                break;
            }
            ++ii;
        }    
        
        t3 = high_resolution_clock::now();
        cout << "Time taken to search: " << duration_cast<microseconds>(t3-t2).count() << " microseconds." << endl;

    }

    if(found == false)
    {
        num_results = 500;
        cout << "Last ditch effort..." << endl;
        cout << "Expanding search to " << num_results << endl;

        t1 = high_resolution_clock::now();

        ret_indexes.resize(num_results);
        out_dists_sqr.resize(num_results);

        resultSet = nanoflann::KNNResultSet<real>(num_results);
        resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
        
        tree.index->findNeighbors(resultSet, &testp[0], nanoflann::SearchParams(10,0,true));

        t2 = high_resolution_clock::now();

        cout << "Time taken to get neighbours: " << 
        duration_cast<microseconds>(t2-t1).count() << " microseconds." << endl;

        ii = 0;
        for( auto const& cell:ret_indexes)
        {
            if(CheckCell(cell,cells,testp,pert))
            {
                cout << "Successfully found the cell: " << cell << endl;
                cout << "Index in the search: " << ii << endl;
                found = true;
                break;
            }
            ++ii;
        }    
        
        t3 = high_resolution_clock::now();
        cout << "Time taken to search: " << duration_cast<microseconds>(t3-t2).count() << " microseconds." << endl;

    }

}

int main (int argc, char **argv)
{
    Eigen::initParallel();

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
	high_resolution_clock::time_point t2;

    SIM svar;
    FLUID fvar;
    AERO avar;
    MESH cells;

    StateVecD testp = StateVecD::Zero();
    
    read_params(argc,argv,svar,testp);
    
    t2 = high_resolution_clock::now();

    cout << "Time taken to read settings: " << duration_cast<microseconds>(t2-t1).count() << " microseconds." << endl; 

    #if SIMDIM == 3
        if(svar.CDForFOAM == 0)
        {
            TAU::Read_BMAP(svar);
            TAU::Read_TAUMESH_FACE(svar,cells,fvar,avar);
        }
        else
        {
            FOAM::Read_FOAM(svar,cells);
        }
    #else
        if (svar.CDForFOAM == 0)
        {
            vector<uint> used_verts;
            TAU::Read_BMAP(svar);
            TAU::Read_TAUMESH_EDGE(svar,cells,fvar,avar,used_verts);
        }
        else
        {
            cout << "OpenFOAM mesh input is currently unsupported in two-dimensions" << endl;
            exit(-1);
        }
    #endif

    high_resolution_clock::time_point t3 = high_resolution_clock::now();

    cout << "Time taken to read mesh: " << duration_cast<microseconds>(t3-t2).count() << " microseconds." << endl; 

    Vec_Tree tree(SIMDIM,cells.cCentre,20);
	
    high_resolution_clock::time_point t4 = high_resolution_clock::now();

	tree.index->buildIndex();

    high_resolution_clock::time_point t5 = high_resolution_clock::now();

    cout << "Time taken to build Vec_Tree: " << duration_cast<microseconds>(t4-t3).count() << 
    " KDTree: " << duration_cast<microseconds>(t5-t4).count() << " microseconds." << endl; 

    test_containtment(cells,tree,testp);

    return 0;
}