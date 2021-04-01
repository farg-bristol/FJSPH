
#include "Var.h"
#include "BinaryIO.h"
#include "IOFunctions.h"
#include "IO.h"

void ReadLatestTime(SIM& svar, FLUID& fvar, AERO& avar, State& pn)
{
#ifdef DEBUG
	dbout << "Reading frame info file for restart information." << endl;
#endif

	// Re-read the settings and fluid file
	// cout << svar.outfolder << endl;

	string outdir = svar.outfolder;
	if(outdir.back() != '/')
    	outdir.append("/");


  	Read_SIM_Var(outdir,svar,fvar,avar);

  	Read_FLUID_Var(outdir,svar,fvar,avar);

  	svar.outfolder = outdir;

	// Now get the data from the files. Start with the boundary
	if(svar.outtype == 0)
	{
		State boundary, fuel;
		string file;	

		// Read the fuel
		void* fuelHandle = NULL;
		string fuelf = outdir;
		int32_t fuelFrames, boundFrames;
		double fuelTime, boundTime;
		fuelf.append("Fuel.szplt");

		if(tecFileReaderOpen(fuelf.c_str(),&fuelHandle))
		{
			cout << "Error opening szplt file. Path:" << endl;
			cout << file << endl;
			exit(-1);
		}

		cout << "Checking Fuel file..." << endl;
		CheckContents(fuelHandle,svar,fuelFrames,fuelTime);

		if (svar.Bcase != 4 && svar.Bcase !=0)
		{
			void* boundHandle = NULL;
			string boundf = outdir;
			boundf.append("Boundary.szplt");

			if(tecFileReaderOpen(boundf.c_str(),&boundHandle))
			{
				cout << "Error opening szplt file. Path:" << endl;
				cout << file << endl;
				exit(-1);
			}

			cout << "Checking Boundary file..." << endl;
			CheckContents(boundHandle,svar,boundFrames,boundTime);

			if(fuelFrames!= boundFrames)
			{
				cout << "Caution! Number of frames is not consistent between fuel and boundary files." << endl;
			}


			if(fuelTime != boundTime)
			{
				cout << "Caution! Frame times are not consistent between fuel and boundary files." << endl;

				if(fuelTime > boundTime)
				{
					double time = 0.0;
					for(int32_t frame = fuelFrames-1; frame > 1; frame--)
					{
						if(tecZoneGetSolutionTime(fuelHandle, frame, &time))
						{
							cout << "Failed to get time data for frame : " << frame << " from fuel file." << endl;
							continue;
						}
						
						if(time == boundTime)
						{
							cout << "Found the correct frame" << endl;
							fuelFrames = frame;
							fuelTime = time;
							break;
						}
					}
				}
				else
				{
					double time = 0.0;
					for(int32_t frame = boundFrames-1; frame > 1; frame--)
					{
						if(tecZoneGetSolutionTime(boundHandle, frame, &time))
						{
							cout << "Failed to get time data for frame : " << frame << " from boundary file." << endl;
							continue;
						}
						

						if(time == fuelTime)
						{
							cout << "Found the correct frame" << endl;
							boundFrames = frame;
							boundTime = time;
							break;
						}
					}
				}

				if(fuelTime != boundTime)
				{
					cout << "Could not find a consistent time in each file. Stopping." << endl;
					exit(-1);
				}
			}

			cout << "Attempting to read the boundary..." << endl;
			Read_Binary_Timestep(boundHandle,svar,boundFrames,boundary);
		}

	    svar.frame = fuelFrames-1;

	    // Read the actual data.
		cout  << "Attempting to read the fuel..." << endl;
		Read_Binary_Timestep(fuelHandle,svar,fuelFrames,fuel);


		pn = boundary;
		svar.bndPts = boundary.size();
		svar.simPts = fuel.size();
		pn.insert(pn.end(),fuel.begin(),fuel.end());	
		svar.totPts = pn.size();
		
		if(svar.simPts + svar.bndPts != svar.totPts)
		{
			cout << "Mismatch of array sizes. Total array is not the sum of sim and boundary arrays" << endl;
			exit(-1);
		}

		// if(svar.totPts != totPts || svar.bndPts != bndPts || svar.simPts != simPts)
		// {
		// 	cout << "Mismatch of the particle numbers in the frame file and data file." << endl;
		// }
	}
	else
	{
		// TODO: ASCII Restart.
		// Particle numbers can be found from frame file.
		// Find EOF, then walk back from there how many particles.

		cout << "ASCII file collation not yet implemented." << endl;
		exit(-1);
	}
}

void WriteLatestTime(SIM& svar, State& pn)
{
	string fout = svar.outfolder;
	fout.append("FuelLatest.szplt");

	string bout = svar.outfolder;
	bout.append("BoundLatest.szplt");

	void* fuelFile;
	void* boundFile;


	Init_Binary_PLT(svar,fout,"Simulation Particles",fuelFile);

	Write_Binary_Timestep(svar,pn,svar.bndPts,svar.totPts,"Fuel",1,fuelFile);


	Init_Binary_PLT(svar,fout,"Boundary Particles",boundFile);

	Write_Binary_Timestep(svar,pn,0,svar.bndPts,"Boundary",2,boundFile);

}


int main(int argc, char *argv[])
{
	SIM svar;
	FLUID fvar;
	AERO avar;
	State pn;

	string file = argv[1];

	if(file.back() != '/')
    	file.append("/");

  	Read_SIM_Var(file,svar,fvar,avar);

  	/*Get fluid properties from fluid.dat*/
	Read_FLUID_Var(file,svar,fvar,avar);

	char cCurrentPath[FILENAME_MAX];
	if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
	{
		cerr << "Failed to get current working directory." << endl;
		exit(-1);
	}


	/*Get output absolute path*/
	string pathname = cCurrentPath;
  	pathname.append("/");
  	pathname.append(svar.infolder);
  	pathname.append(svar.outfolder);
  	pathname.append("/");
  
  	/*Check for output file name*/		
	check_folder(pathname);
	
	/*Check if there is a slash at the end.*/
  	if (pathname.back() != '/')
  	{
  		pathname.append("/");
  	}

  	svar.outfolder = pathname;
	cout << "Output folder: " << svar.outfolder << endl << endl;

	ReadLatestTime(svar,fvar,avar,pn);

	WriteLatestTime(svar,pn);

	cout << "Lastest frame written successfully. Exiting." << endl;

	return 0;

}