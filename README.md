# WCSPH ReadMe 
Latest update: 19/12/18


# What is this repository about?
WCSPH, or Weakly Compressible Smoothed Particle Hydrodynamics, is a great method for the modelling of fluid flows - particularly liquid motion. Typical applications are scenarios such as tank sloshing, wave breaking and many others. The aim for this specific repository is to apply this to situations with a high Weber number, like jets in a crossflow or coflow, typically found in fuel injection or jettison. 

The specific scenario of the PhD associated with this repository is aircraft fuel jettison. For the initial breakup of the fuel where it's acting very much like a liquid, the current plan is to use SPH for the model. After the initial breakup the SPH will transition to using a particle tracking method, that is far more efficient, until the particles are out of the simulation space. 

# What is WCSPH?
Weakly Compressible Smoothed Particle Hydrodynamics (WCSPH) is one formulation of the SPH equations, which were first introduced by Gingold and Monaghan (1977), and independently by Lucy (1977), for astrophyiscal simulations. It was realised a few years later that this formulation is also applicable to fluid dynamics with reasonable accuracy. Most notable is the ease with which liquid motion can be modelled. Unlike many other methods, SPH doesn't require any special treatment of the boundary between liquid and gas phases. The equations handle this implicitly. For a good summary on SPH equations, right from the base equations to all the additions made through the years, Monaghan (2005) has the most notable contributions, excepting surface tension models. For a review of surface tension models, Huber et al. (2015) provides a review of three primary models. 

# What does this code do?
At present the code is still in large development, so functionality isn't anywhere near where it should be to be a functional tool in your project, at least without modifications of your own. SPH has the brilliant feature of essentially being able to plug in equations for your particular need. At some point it may be worth analysing which models are most appropriate for the intended application of this repository, however accuracy is unlikely to be a priority compared to speed. As such most of the code here is intended for simplicity. The current code has the following set of formulations (detailed in the references stated above):

Base SPH: Monaghan (1992) + XSPH correction.

Laminar Viscosity: Morris (1997).

Surface Tension: Nair & Poeschel (2018).

Density Reinitialisation: Colagrossi & Landrini (2003).

# Dependencies 
The code only has *three* dependent libraries at present:

[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) 

[NanoFLANN](https://github.com/jlblancoc/nanoflann) - Jose Luis Blanco-Claraco (University of Almería).

[Poisson Disk Sampling](https://github.com/thinks/poisson-disk-sampling) - Tommy Hinks.

The code has been built such that these need to be in local subfolders "Eigen", "NanoFLANN", and "PDS" respectively. The bare minimum of each has been provided in the repository, so download from other sources shouldn't be necessary. For exactly what needs to be in these folders (particularly for NanoFLANN), see the include headers in the source. 

# Building
A makefile is provided to ease compliation. There are several options in the build file, but because there are so few dependencies and files associated at present, it is still a single compile line, so even if you dont have a make environment, it's easy to see what command is needed. To build, all the user needs to write is `make build` and the code should then compile. To run, only one input file has been provided as a sample, and it is `make cross` in order to run the program with that input. 

## C++ 11
A word of warning that the code requires C++11 standard support, and so will not build on Visual Studio 13, due to an incomplete support of the C++11 standard. If you are on windows and don't have access to VS15, then MinGW is the next best bet (unless you are on windows 10, and have WSL enabled).

# Input
On the command line, either zero, one, or two arguments can be provided. If no inputs are provided, the program assumes a standard set of settings for a simulation. The first argument is the input file for settings, that will then overwrite the defaults (aore detail on the options later). Comments are provided in the input file that should be self-explanatory as to their function. The second argument is the output file name. If not provided it defaults to 'Test.plt'. More than two arguments simply states a warning that they will be ignored. 

In time, I'll add a standalone file describing in better detail what the parameters in the .dat files do (for now they are annotated with their function in the file), as some are more obscure, but obviously because they are still being developed. 

# Output
The code outputs a tecplot360 file, "Test.plt". This has strands in it, so should automatically load ready to animate when brought into tecplot. Unfortunately I haven't figured out how to get tecplot to open straight into scatter, so some steps are needed to get the data to be visible. 

When first loaded into tecplot, it'll show as a `XY Line` plot in the top left. Change this to a `2D Cartesian`. Then tick `scatter`, and the data should then become visible. To make it more pretty and presentable, you can change the particle properties in `Zone Style`. The boundary and simulation particles are separate zones, so each can be modified appropriately. To play the simulation, simply press the play button on the left. 

An option may be added to allow writing without strands, as this severely imacts load times into tecplot, and prevents you from retaining the view style between files. 

# Future Plans
The plan for the future is to continue developing code for the crossflow case and try to evaluate which way is best to apply the force to the fluid. At present the methods wishing to be tested are:

* Force on all upstream particles. (works pretty well)
* Force on all particles defined as a surface. (Terrible)
* Force on all particles defined as a surface, proportional to the surface tension. (Works really well)
* Enact force and surface tension on particles using generated 'ghost' particles. (complicated)

Natually a lot of work needs to be done in order to successfully test these, and evaluate each of their traits. At present, some testing of the first method has been done, with reassuring success. Work is afoot to incorporate the other options into a code to allow easy comparison. The list is generally arranged in what I anticipate being low to high accuracy. What I am unsure about is which is the highest computation, particularly for a given accuracy. Using ghost particles is likely to be the most computationally heavy, however could well achieve the highest accuracy for the least tuning work. Finding the relationship between surface tension and aerodynamic force is likely to require lots of tuning, which may never even finish. Therefore part of me leans to that to save the uncertainty. 

Another point of attraction for the use of ghost particles is that there appears to be a scarcity in research in this direction of SPH, at least a second phase other than the boundary (see Schechter & Brison, 2012). If implementation is effective, then there is possibility for it to be used in several other applications where surface tension is critical to the scenario, but only a single phase is of interest. 

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
