# FJSPH ReadMe

# What is this repository about?

WCSPH, or Weakly Compressible Smoothed Particle Hydrodynamics, is a great method for the modelling of fluid flows - particularly liquid motion. Typical applications are scenarios such as tank sloshing, wave breaking and many others. The aim for this specific repository is to apply this to situations with a high Weber number, like jets in a crossflow or coflow, typically found in fuel injection or jettison.

The specific scenario of the PhD associated with this repository is aircraft fuel jettison. For the initial breakup of the fuel where it's acting very much like a liquid, the current plan is to use SPH for the model. After the initial breakup the SPH will transition to using a particle tracking method, that is far more efficient, until the particles are out of the simulation space.

# What is WCSPH?

Weakly Compressible Smoothed Particle Hydrodynamics (WCSPH) is one formulation of the SPH equations, which were first introduced by Gingold and Monaghan (1977), and independently by Lucy (1977), for astrophyiscal simulations. It was realised some years later that this formulation is also applicable to fluid dynamics with reasonable accuracy. Most notable is the ease with which liquid motion can be modelled. Unlike many other methods, SPH doesn't require any special treatment of the boundary between liquid and gas phases in order to track the interface. Of course, for good physics surface tension must be considered, but the nature of the particle formulation For a good summary on SPH equations, right from the base equations to all the additions made through the years, Monaghan (2005) has the most notable contributions, excepting surface tension models. For a review of surface tension models, Huber et al. (2015) provides a review of three primary models.

Recently work has been done to massively improve the performance of SPH in regions where negative pressures exist, and to smooth the density field. These are the delta-SPH formulations and particle shifting methods, and further still the ALE (Arbitrary Lagrangian Eulerian) formulation. The ALE formulation combines the shifting more completely into the SPH equations, minimising the loss of conservation, at the gain of particle distribution.

# What does this code do?

At present the code is still in large development, and many features are not present, that I would like to be. IO is a large area where work has been avoided. This repository is hugely designed for simulating fuel exiting a pipe, into an aerodynamic solution provided using a TAU mesh format. Examples of these are provided to give an idea how to interact with the program. Equations wise, the current code has the following set of formulations (detailed in the references stated below):

Base SPH: Monaghan (1992) + delta-SPH continuity correction (Marrone,2018) + ALE-SPH (Sun,2019).

Laminar Viscosity: Linking of artificial viscosity to actual viscosity (1997).

Surface Tension: Pairwise force (Nair and Poeschel, 2018).

Aerodynamic coupling: Upwind droplet continuum correction (Gissler et al.,2017).

# Dependencies

The code has _four_ dependent libraries at present:

[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) - Tested with 3.4-rc1, 3.3.9, 3.3.8.

[NanoFLANN](https://github.com/jlblancoc/nanoflann) - Jose Luis Blanco-Claraco (University of Almería).

[nlohmann JSON for Modern C++](https://json.nlohmann.me/) - Tested with v3.7.3.

[NetCDF](https://www.unidata.ucar.edu/downloads/netcdf/) - Tested with NetCDF 4.7.3.

[TECIO](https://www.tecplot.com/products/tecio-library/) - Tested with 2016 R2, 2018 R2, 2019 R1, 2020 R2 versions of TECIO.

The code has been built such that the first three libraries need to be in local subfolder "Third_Party". NetCDF is assumed to be installed on the system path, and the TECIO files are in the standard install location for linux tecplot (e.g. /usr/local/tecplot/360ex_2020r2/include) and included in the makefile.
For your system, the include directories in the makefile will need to be modified for your tecplot version and install location to correctly locate the TECIO files.
It is possible to download only the TECIO library, via the [SU2 code github](https://github.com/su2code/SU2/tree/master/externals/tecio).
The bare minimum of Eigen has been provided in the repository, and additional download of the library shouldn't be necessary.

# Building

A makefile is provided to ease compliation.
There are several options in the build file, but because there are so few dependencies and files associated at present, it is still a single compile line.
The target directory for compilation has been taken as `$(HOME)/bin` for the compilation, and this will need adjusting to wherever you wish the program to exist.
In the makefile are four primary options: `2D`, `debug2`, `3D`, and `debug3`.

These are build options for two and three dimensions, and with or without debug flags.
Additionally are options to build associated programs Cell2Edge (`c2e`), for building 2D Edge based mesh files from TAU files, and Cell2Face (`c2f`) for 3D TAU meshes to face based data. This does not need to be run for OpenFOAM meshes, as it is already a face based file format.

Lastly are some options that can be passed to the compiler, using the `FLAGS` option, e.g. `make 2D FLAGS='-DNOFROZEN'`.
Currently the accepted flags are:
`-DCUBIC` - Use the cubic spline kernel (not recommended).
`-DCSF` - Vernaud et al. surface tension model.
`-DHEST` - He et al. surface tension model.
`-DPAIRWISE` - Use the pairwise surface tension model. Current default.
`-DNOST` - No surface tension model at all.
`-DISOEOS` - use Isothermal equation of state.
`-DCOLEEOS` - Use cole equation of state. Current default.
`-DTIC` - Engage tensile instability control. Can help for certain things but does slow computation down.
`-DLINEAR` - Use the Marrone et al. linearised pressure momentum relationship.
`-DNOFROZEN` - Don't use frozen dissipation terms.
`-DNODSPH` - Turn delta-SPH off.
`-DFOD=0` - Turns the real data type from double to float.

## C++ Standards

A minimum of c++17 is required for filesystem header.
The NanoFLANN code requires C++11 standard support, and so will not build on Visual Studio 13, due to an foam_is_incomplete support of the C++11 standard. If you are on windows and don't have access to VS15, then MinGW is the next best bet (unless you are on windows 10, and have WSL enabled).

Syntax has been used that requires **g++ 9 or later**.

# Input

On the command line, only one argument is accepted; the input para file.
The executable is assumed to be on the system path, but can be put in the same folder to run.
Inside the input file, it is hoped that the title of variables give a good indication of their function, but work will need to be done to create documentation of the options available.
A number of 2D and 3D examples are provided in the examples folder, that should be possible to adapt to the required task.

Example run: `FJSPH.2d para2D`
Provide the input folder for the simulation. This starts the simulation from time 0, and removes any current output files to start writing anew.
Once the simulation writes the first timestep, it will append a line to the end of the para file to say there is data available to restart from.
If a restart is not desired this line will need to be deleted.

If the simulation involves a mesh, then the mesh and solution file are expected in the same folder; look at the **RAE2822** folder for an example.

# Output

The code uses the TECIO library to output szplt files, because of it's highly efficient storage, and great ability to append to an existing file, and ability to view the results with the simulation still running.
There are two files that are output, appending "\_boundary.szplt" and "\_fuel.szplt" to the output prefix defined in the para file.  
If there are no boundary points, the boundary file will not be created.

For successful restart certain outputs are required to give enough information. For any simulation, any option above 1 is possible to restart from. However for mesh based simuations, only 3 or 5 can be restarted from, as these contain the cell information.

# Future Plans

Implement Artificial Compressibility Smoothed Particle Hydrodynamics (ACSPH).

# References

Colagrossi, A. & Landrini, M.
_Numerical simulation of interfacial flows by smoothed particle hydrodynamics_,
Journal of Computational Physics, Elsevier BV, 2003 , 191 , 448-475

Gingold, R. A. & Monaghan, J. J.
_Smoothed particle hydrodynamics: theory and application to non-spherical stars_,
Monthly Notices of the Royal Astronomical Society, Oxford University Press (OUP), 1977 , 181 , 375-389

Huber, M.; Reinhardt, S.; Weiskopf, D. & Eberhardt, B.
_Evaluation of Surface Tension Models for SPH-Based Fluid Animations Using a Benchmark Test_,
Workshop on Virtual Reality Interaction and Physical Simulation, The Eurographics Association, 2015

Lucy, L. B.
_A numerical approach to the testing of the fission hypothesis_,
The Astronomical Journal, IOP Publishing, 1977 , 82 , 1013

Monaghan, J. J.
_Smoothed particle hydrodynamics_,
Annual Review of Astronomy and Astrophysics, Annual Reviews, 1992 , 30 , 543-574

Monaghan, J. J.
_Smoothed particle hydrodynamics_,
Reports on progress in physics, IOP Publishing, 2005 , 68 , 1703

Morris, J. P.; Fox, P. J. & Zhu, Y.
_Modeling Low Reynolds Number foam_is_incompressible Flows Using SPH_,
Journal of Computational Physics, Elsevier BV, 1997 , 136 , 214-226

Nair, P. & Pöschel, T.
_Dynamic capillary phenomena using foam_is_incompressible SPH_,
Chemical Engineering Science, Elsevier BV, 2018 , 176 , 192-204
