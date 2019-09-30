# WCSPH ReadMe 
Latest update: 30/09/2019


# What is this repository about?
WCSPH, or Weakly Compressible Smoothed Particle Hydrodynamics, is a great method for the modelling of fluid flows - particularly liquid motion. Typical applications are scenarios such as tank sloshing, wave breaking and many others. The aim for this specific repository is to apply this to situations with a high gas to liquid velocity ratio, like jets in a crossflow or coflow, typically found in fuel injection or jettison. 

The specific scenario of the PhD associated with this repository is aircraft fuel jettison. For the initial breakup of the fuel where it's acting very much like a liquid, the current plan is to use SPH for the model. After the initial breakup the SPH will transition to using a particle tracking method, that is far more efficient, until the particles are out of the simulation space. 

# What is WCSPH?
Weakly Compressible Smoothed Particle Hydrodynamics (WCSPH) is one formulation of the SPH equations, which were first introduced by Gingold and Monaghan (1977), and independently by Lucy (1977), for astrophyiscal simulations. It was realised a few years later that this formulation is also applicable to fluid dynamics with reasonable accuracy. Most notable is the ease with which liquid motion can be modelled. Unlike many other methods, SPH doesn't require any special treatment of the boundary between liquid and gas phases. The equations handle this implicitly. For a good summary on SPH equations, right from the base equations to all the additions made through the years, Monaghan (2005) has the most notable contributions, excepting surface tension models. For a review of surface tension models, Huber et al. (2015) provides a review of three primary models. 

# What does this code do?
At present the code is still in large development, so functionality isn't anywhere near where it should be to be a functional tool in your project, at least without modifications of your own. SPH has the brilliant feature of essentially being able to plug in equations for your particular need. At some point it may be worth analysing which models are most appropriate for the intended application of this repository, however accuracy is unlikely to be a priority compared to speed. As such most of the code here is intended for simplicity. The current code has the following set of formulations (detailed in the references stated above):

Base SPH: Monaghan (1992).

Laminar Viscosity: Morris (1997).

Surface Tension: Nair & Poeschel (2018). (Optional - Currently inactive)

Aerodynamic coupling: Gissler et al. (2017).

# Dependencies 
The code only has *four* dependent libraries at present:

[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) 

[NanoFLANN](https://github.com/jlblancoc/nanoflann) - Jose Luis Blanco-Claraco (University of Almería).

[TECIO](https://github.com/su2code/SU2/tree/master/externals/tecio) - TecPlot.

[GMPXX](https://gmplib.org/manual/index.html#Top)

The code has been built such that Eigen and NanoFLANN need to be in local subfolders "Eigen", "NanoFLANN" respectively. The bare minimum of each has been provided in the repository, so download from other sources shouldn't be necessary. For exactly what needs to be in these folders (particularly for NanoFLANN), see the include headers in the source. 
TECIO and and the GMP library are assumed to be on the system path. GMP is normally very easy to install on linux or MinGW, and the makefile includes the necessary library links. The hyperlink also provides info on how to install. Please also make sure you have installed the c++ library as well, GMPXX, otherwise it won't work. TECIO on the other hand, is a little more involved. The hyperlink brings you to the version distributed with SU2. This has a makefile and a Cmake file, allowing for installation either in the system path or locally. These will make a library .a file, which will need to be linked to either on the system or in the compile command (not included). It may also be advisable to have boost-1.55 or later available, however I don't believe it is needed. 

# Building
A makefile is provided to ease compliation. There are several options in the build file, but because there are so few dependencies and files associated at present, it is still a single compile line, so even if you dont have a make environment, it's easy to see what command is needed. To build, all the user needs to write is `make 2D` or `make 3D` and the code should then compile. 

## C++ 11
A word of warning that the code requires C++11 standard support, and so will not build on Visual Studio 13, due to an incomplete support of the C++11 standard. If you are on windows and don't have access to VS15, then MinGW is the next best bet (unless you are on windows 10, and have WSL enabled. This is also the development platform).

# Input
On the command line, either zero, one, or two arguments can be provided. If no inputs are provided, the program assumes a standard set of settings for a simulation. The first argument is the input folder for settings, that will then overwrite the defaults (more detail on the options later). In this folder should be a file named `settings.dat`, `fluid.dat`, and depending on the options chosen, a third file for the mesh or VLM settings. A sample input is provided in `/Inputs/`. Comments are provided in the input file that should be self-explanatory as to their function. The second argument is the output folder name. If not provided it defaults to the same name as the input folder. More than two arguments simply states a warning that they will be ignored. 

In time, I'll add a standalone file describing in better detail what the parameters in the .dat files do (for now they are annotated with their function in the file), as some are more obscure, but obviously because they are still being developed. 

# Output
The code outputs a set of .plt tecplot360 files. This has strands in it, so should automatically load ready to animate when brought into tecplot. Unfortunately I haven't figured out how to get tecplot to open straight into scatter, so some steps are needed to get the data to be visible. 

When first loaded into tecplot, it'll show as a `XY Line` plot in the top left. Change this to a `2D Cartesian`. Then tick `scatter`, and the data should then become visible. To make it more pretty and presentable, you can change the particle properties in `Zone Style`. The boundary and simulation particles are separate zones, so each can be modified appropriately. To play the simulation, simply press the play button on the left. 

An option may be added to allow writing without strands, as this severely imacts load times into tecplot, and prevents you from retaining the view style between files. 

# Future Plans
The future is somewhat uncertain regarding this project. At present higher accuracy simulations are needed for the code to be run on before an evaluation of accuracy can be done. I have reservations that the fidelity is sufficient, however further testing is needed to have confidence in that. 

# References 
Gingold, R. A. & Monaghan, J. J.
*Smoothed particle hydrodynamics: theory and application to non-spherical stars*, 
Monthly Notices of the Royal Astronomical Society, Oxford University Press (OUP), 1977 , 181 , 375-389

Gissler, C. & Band, S. & Peer, A. & Ihmsen, M. & Teschner, M.
*Approximate Air-Fluid Interactions for SPH*, 
Workshop on Virtual Reality Interaction and Physical Simulation, The Eurographics Association, 2017,

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
