# WCSPH ReadMe 
Latest update: 14/12/18

# What is this repository about?
WCSPH, or Weakly Compressible Smoothed Particle Hydrodynamics, is a great method for the modelling of fluid flows - particularly liquid motion. Typical applications are scenarios such as tank sloshing, wave breaking and many others. The aim for this specific repository is to apply this to situations with a high Weber number, like jets in a crossflow or coflow, typically found in fuel injection or jettison. 
The specific scenario of the PhD associated with this repository is aircraft fuel jettison. For the initial breakup of the fuel where it's acting very much like a liquid, the current plan is to use SPH for the model.

# What is WCSPH?
Weakly Compressible Smoothed Particle Hydrodynamics (WCSPH) is one formulation of the SPH equations, which were first introduced by Gingold and Monaghan (1977), and independently by Lucy (1977), for astrophyiscal simulations. It was realised a few years later that this formulation is also applicable to fluid dynamics with reasonable accuracy. Most notable is the ease with which liquid motion can be modelled. Unlike many other methods, SPH doesn't require any special treatment of the boundary between liquid and gas phases. The equations handle this implicitly. For a good summary on SPH equations, right from the base equations to all the additions made through the years, Monaghan (2005) has the most notable contributions, excepting surface tension models. For a review of surface tension models, Huber et al. (2015) provides a review of three primary models. 

# What does this code do?
At present the code is still in large development, so functionality isn't anywhere near where it should be to be a functional tool in your project, at least without modifications of your own. SPH has the brilliant feature of essentially being able to plug in equations for your particular need. At some point it may be worth analysing which models are most appropriate for the intended application of this repository, however accuracy is unlikely to be a priority compared to speed. As such most of the code here is intended for simplicity. The current code has the following set of formulations (detailed in the references stated above):

Base SPH: Monaghan (1992) + XSPH correction.

Laminar Viscosity: Morris (1997).

Surface Tension: Nair & Poeschel (2018).

Density Reinitialisation: Colagrossi & Landrini (2003).

# Dependencies 
The code only has two dependent libraries at present:

[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). 

[NanoFLANN](https://github.com/jlblancoc/nanoflann) - Jose Luis Blanco-Claraco (University of Almería).

The code has been built such that these need to be in local subfolders "Eigen" and "NanoFLANN" respectively. The bare minimum of Eigen and NanoFLANN have been provided in the repository, so download from other sources shouldn't be necessary. 

# Building
A makefile is provided to ease compliation. There are several options in the build file, but because there are so few dependencies and files associated at present, it is still a single compile line. To build, all the user needs to write is `make build` and the code should then compile. To run, only one input file has been provided as a sample, and it is `make cross` in order to run the program (on windows) with that input. 

## VC++ 13
A word of warning that the code will not build on Visual Studio 13, due to an incomplete support of the C++11 standard that is required for the NanoFLANN library. If you are on windows and don't have access to VS15, then MinGW is the next best bet (unless you are on windows 10, and have WSL enabled).


# References 
Colagrossi, A. & Landrini, M.
*Numerical simulation of interfacial flows by smoothed particle hydrodynamics*
Journal of Computational Physics, Elsevier BV, 2003 , 191 , 448-475

Gingold, R. A. & Monaghan, J. J.
*Smoothed particle hydrodynamics: theory and application to non-spherical stars* 
Monthly Notices of the Royal Astronomical Society, Oxford University Press (OUP), 1977 , 181 , 375-389

Huber, M.; Reinhardt, S.; Weiskopf, D. & Eberhardt, B.
*Evaluation of Surface Tension Models for SPH-Based Fluid Animations Using a Benchmark Test* 
Workshop on Virtual Reality Interaction and Physical Simulation, The Eurographics Association, 2015

Lucy, L. B.
*A numerical approach to the testing of the fission hypothesis* 
The Astronomical Journal, IOP Publishing, 1977 , 82 , 1013

Monaghan, J. J.
*Smoothed particle hydrodynamics* 
Annual Review of Astronomy and Astrophysics, Annual Reviews, 1992 , 30 , 543-574

Monaghan, J. J.
*Smoothed particle hydrodynamics* 
Reports on progress in physics, IOP Publishing, 2005 , 68 , 1703

Morris, J. P.; Fox, P. J. & Zhu, Y.
*Modeling Low Reynolds Number Incompressible Flows Using SPH* 
Journal of Computational Physics, Elsevier BV, 1997 , 136 , 214-226

Nair, P. & Pöschel, T.
*Dynamic capillary phenomena using Incompressible SPH*
Chemical Engineering Science, Elsevier BV, 2018 , 176 , 192-204
