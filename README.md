# FJSPH ReadMe 

# What is this repository about?
WCSPH, or Weakly Compressible Smoothed Particle Hydrodynamics, is a great method for the modelling of fluid flows - particularly liquid motion. Typical applications are scenarios such as tank sloshing, wave breaking and many others. The aim for this specific repository is to apply this to situations with a high Weber number, like jets in a crossflow or coflow, typically found in fuel injection or jettison. 

The specific scenario of the PhD associated with this repository is aircraft fuel jettison. For the initial breakup of the fuel where it's acting very much like a liquid, the current plan is to use SPH for the model. After the initial breakup the SPH will transition to using a particle tracking method, that is far more efficient, until the particles are out of the simulation space. 

# What is WCSPH?
Weakly Compressible Smoothed Particle Hydrodynamics (WCSPH) is one formulation of the SPH equations, which were first introduced by Gingold and Monaghan (1977), and independently by Lucy (1977), for astrophyiscal simulations. It was realised some years later that this formulation is also applicable to fluid dynamics with reasonable accuracy. Most notable is the ease with which liquid motion can be modelled. Unlike many other methods, SPH doesn't require any special treatment of the boundary between liquid and gas phases in order to track the interface. Of course, for good physics surface tension must be considered, but the nature of the particle formulation For a good summary on SPH equations, right from the base equations to all the additions made through the years, Monaghan (2005) has the most notable contributions, excepting surface tension models. For a review of surface tension models, Huber et al. (2015) provides a review of three primary models. 

Recently work has been done to massively improve the performance of SPH in regions where negative pressures exist, and to smooth the density field. These are the delta-SPH formulations and particle shifting methods, and further still the ALE (Arbitrary Lagrangian Eulerian) formulation. The ALE formulation combines the shifting more completely into the SPH equations, minimising the loss of conservation, at the gain of particle distribution. 

# What does this code do?
At present the code is still in large development, and many features are not present, that I would like to be. IO is a large area where work has been avoided. This repository is hugely designed for simulating fuel exiting a pipe, into an aerodynamic solution provided using a TAU mesh format. Examples of these are provided to give an idea how to interact with the program. Equations wise, the current code has the following set of formulations (detailed in the references stated below):

Base SPH: Monaghan (1992) + delta-SPH continuity correction (Marrone,2018). Currently undergoing testing for ALE-SPH (Sun,2019).

Laminar Viscosity: Linking of artificial viscosity to actual viscosity (1997).

Surface Tension: Continuum Surface Force (Brackbill,1997) with interactions restricted to surface particles (Ordoubadi,2015).

Aerodynamic coupling: Upwind droplet continuum correction (Gissler et al.,2017).

# Dependencies 
The code has *four* dependent libraries at present:

[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) - Tested with 3.4-rc1, 3.3.9, 3.3.8.

[NanoFLANN](https://github.com/jlblancoc/nanoflann) - Jose Luis Blanco-Claraco (University of Almería).

[NetCDF](https://www.unidata.ucar.edu/downloads/netcdf/) - Tested with NetCDF 4.7.3. 

[TECIO](https://www.tecplot.com/products/tecio-library/) - Tested with 2016 R2, 2018 R2, 2019 R1, 2020 R2 versions of TECIO.

The code has been built such that the first two libraries need to be in local subfolder "Third_Party", then in  "Eigen", and "NanoFLANN". NetCDF is assumed to be installed on the system path, and the TECIO files are in the standard install location for linux tecplot (/usr/local/tecplot/360ex_2020r2/include) and included in the makefile. 
For your system, the include directories in the makefile will need to be modified for your tecplot version and install location to correctly locate the TECIO files.
It is possible to download only the TECIO library, via the [SU2 code github](https://github.com/su2code/SU2/tree/master/externals/tecio).
The bare minimum of Eigen has been provided in the repository, and additional download of the library shouldn't be necessary. 


# Building
A makefile is provided to ease compliation. 
There are several options in the build file, but because there are so few dependencies and files associated at present, it is still a single compile line.
In the makefile are four primary options: `2D`, `debug2`, `3D`, and `debug3`. 

These are build options for two and three dimensions, and with or without debug flags. 
Additionally are options to build associated programs Cell2Edge (`c2e`), for building 2D Edge based mesh files from TAU files, and Cell2Face (`c2f`) for 3D TAU meshes to face based data. 

Lastly are some options that can be passed to the compiler, using the `FLAGS` option, e.g. `make 2D FLAGS='-DALE'`. Currently the two accepted flags are `-DALE` and `-DCUBIC`. The ALE flag turns on the particle shifting and ALE momentum and continuity equations, and is recommended to improve distribution of the particles. The cubic flag changes the kernel to the cubic spline, instead of the Wendland C2. While provided, it is not recommended, as in the tests I've done it quickly becomes unstable.

## C++ Standards
The NanoFLANN code requires C++11 standard support, and so will not build on Visual Studio 13, due to an incomplete support of the C++11 standard. If you are on windows and don't have access to VS15, then MinGW is the next best bet (unless you are on windows 10, and have WSL enabled).

Syntax has been used that requires **g++ 9 or later**, so the code will not compile with g++ 8 or older. 

# Input
On the command line either one, or two options exists; a clean run or run from a current solution. If no inputs are provided, the code exits without performing anything. Comments are provided in the input file that should be self-explanatory as to their function.

Option 1: `./WCSPH Droplet2D`
Provide the input folder for the simulation. This starts the simulation from time 0, and removes any current output files to start writing anew.

Option 2: `./WCSPH -r Droplet2D`
This restarts the simulation, assuming that the output data exists and can be restarted from.

More than two arguments simply states a warning that they will be ignored. 

In the input folder should be all files required for the simulation.
The two files that are required for any and all simulations is the "settings" and "fluid" files. 
If the simulation involves a mesh, then the mesh and solution file are expected in the same folder, and these file names are expected in the fluid file, after the solution file. For an example of a mesh based simulation, look at the **RAE2822** folder. 

# Output
The code uses the TECIO library to output szplt files, because of it's highly efficient storage, and great ability to append to an existing file, and ability to view the results with the simulation still running. 
There are three files that are output in this format: "Boundary.szplt", "Fuel.szplt", and "Restart.szplt". 
The first two files provide the point data for all the SPH particles. The boundary and fuel are separate if one wished to save space and not output the boundary files during the simulation.
If this option is selected however, restart is not possible. 

Lastly, the restart file is written to provide a lossless storage of parameters to restart the simulation.
For the restart to work, certain output options need to be selected to have restart data. Three options exist, two of which are for mesh-based rruns.
Output option 2 is restart data for non mesh based runs, i.e. with uniform flow. 
Output option 5 is restart data for a mesh based run, with less information about the cell, whereas option 3 is all cell data is also written for the point.

# Future Plans
A great amount more work could be done in this field. 
Section needs to be added to later.


# References 
Colagrossi, A. & Landrini, M.
*Numerical simulation of interfacial flows by smoothed particle hydrodynamics*,
Journal of Computational Physics, Elsevier BV, 2003 , 191 , 448-475

Gingold, R. A. & Monaghan, J. J.
*Smoothed particle hydrodynamics: theory and application to non-spherical stars*, 
Monthly Notices of the Royal Astronomical Society, Oxford University Press (OUP), 1977 , 181 , 375-389

Huber, M.; Reinhardt, S.; Weiskopf, D. & Eberhardt, B.
*Evaluation of Surface Tension Models for SPH-Based Fluid Animations Using a Benchmark Test*, 
Workshop on Virtual Reality Interaction and Physical Simulation, The Eurographics Association, 2015

Lucy, L. B.
*A numerical approach to the testing of the fission hypothesis*, 
The Astronomical Journal, IOP Publishing, 1977 , 82 , 1013

Monaghan, J. J.
*Smoothed particle hydrodynamics*, 
Annual Review of Astronomy and Astrophysics, Annual Reviews, 1992 , 30 , 543-574

Monaghan, J. J.
*Smoothed particle hydrodynamics*, 
Reports on progress in physics, IOP Publishing, 2005 , 68 , 1703

Morris, J. P.; Fox, P. J. & Zhu, Y.
*Modeling Low Reynolds Number Incompressible Flows Using SPH*, 
Journal of Computational Physics, Elsevier BV, 1997 , 136 , 214-226

Nair, P. & Pöschel, T.
*Dynamic capillary phenomena using Incompressible SPH*,
Chemical Engineering Science, Elsevier BV, 2018 , 176 , 192-204

Schechter, H. & Bridson, R.
*Ghost SPH for animating water*, 
ACM Transactions on Graphics, Association for Computing Machinery (ACM), 2012 , 31 , 1-8

