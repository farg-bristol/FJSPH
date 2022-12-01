# Compiler to use. In C++, so g++.
CC:=gcc
CXX:=g++-10

NPROCS = $(shell grep -c 'processor' /proc/cpuinfo)
MAKEFLAGS += -j$(NPROCS)

# Libraries to include
TECLIB:=-ltecio -lopenblas
NETLIB:=-lnetcdf #-lH5hut -lhdf5

#-lgmpxx -lgmp
# Compiler flags. 
CXXFLAGS:=-std=c++17 -Wall -ffast-math -funroll-loops -O3 -fopenmp -march=native -mtune=native -DHAS_BLAS

# TECINC:=-I/usr/local/tecplot/360ex_2020r2/include/
# TECLINK:=-L/usr/local/tecplot/360ex_2020r2/bin/

# TECINC:=-I/usr/local/tecplot/360ex_2021r2/include/
# TECLINK:=-L/usr/local/tecplot/360ex_2021r2/bin/

TECINC:=-I/usr/local/tecplot/360ex_2022r1/include/
TECLINK:=-L/usr/local/tecplot/360ex_2022r1/bin/

H5INC=-I/usr/include/hdf5/serial/
H5LIB=-L/usr/lib/x86_64-linux-gnu/hdf5/serial/

EIGENINC=-I$(HOME)/Eigen

TARGETDIR:=$(HOME)/bin

DIM=-DSIMDIM=2

SRC := src/Add.cpp src/AsciiIO.cpp src/BinaryIO.cpp src/CDFIO.cpp src/Containment.cpp\
	 	src/Droplet.cpp src/FJSPH.cpp src/FOAMIO.cpp src/Geometry.cpp src/Init.cpp src/Integration.cpp src/IO.cpp\
	    src/IPT.cpp src/Neighbours.cpp src/Newmark_Beta.cpp src/Resid.cpp src/Runge_Kutta.cpp\
	    src/Shifting.cpp src/Speedtest.cpp src/VLM.cpp\
		src/shapes/arc.cpp src/shapes/circle.cpp src/shapes/cylinder.cpp src/shapes/inlet.cpp src/shapes/line.cpp src/shapes/square.cpp

SRCDIR := src
# Make the object directories in the excutable, so that it's not uploaded to cloud all the time.
2DOBJROOT := $(TARGETDIR)/obj_2d
3DOBJROOT := $(TARGETDIR)/obj_3d

DSPHEXT := dsph
2DDSPHDIR := $(2DOBJROOT)/$(DSPHEXT)
3DDSPHDIR := $(3DOBJROOT)/$(DSPHEXT)

ALEEXT := ALE
2DALEDIR := $(2DOBJROOT)/$(ALEEXT)
3DALEDIR := $(3DOBJROOT)/$(ALEEXT)

# Create a list of object files to go alongside the source files
# SRC := $(shell find $(SRCDIR) -type f -name *.cpp)
2DOBJ := $(patsubst $(SRCDIR)/%,$(2DDSPHDIR)/%,$(SRC:.cpp=.o))
3DOBJ := $(patsubst $(SRCDIR)/%,$(3DDSPHDIR)/%,$(SRC:.cpp=.o))
2DDEP := $(2DOBJ:.o=.d)
3DDEP := $(3DOBJ:.o=.d)

# DSPHEXT := dsph
# 2DDSPHDIR:= $(2DOBJDIR)/$(DSPHEXT)
# 3DDSPHDIR:= $(3DOBJDIR)/$(DSPHEXT)
# 2DDSPHOBJ := $(patsubst $(2DOBJDIR)/%,$(2DDSPHDIR)/%,$(2DOBJ))
# 3DDSPHOBJ := $(patsubst $(3DOBJDIR)/%,$(3DDSPHDIR)/%,$(3DOBJ))
# 2DDEPD := $(2DDSPHOBJ:.o=.d)
# 3DDEPD := $(3DDSPHOBJ:.o=.d)

2DALEOBJ := $(patsubst $(2DOBJROOT)/$(DSPHEXT)/%,$(2DALEDIR)/%,$(2DOBJ))
3DALEOBJ := $(patsubst $(3DOBJROOT)/$(DSPHEXT)/%,$(3DALEDIR)/%,$(3DOBJ))
2DDEPA := $(2DALEOBJ:.o=.d)
3DDEPA := $(3DALEOBJ:.o=.d)

.PHONY: all clean

all: 2D 3D

2ddirs:
	@mkdir -p $(2DDSPHDIR)
	@mkdir -p $(2DDSPHDIR)/shapes
	@mkdir -p $(2DALEDIR)
	@mkdir -p $(2DALEDIR)/shapes

3ddirs:
	@mkdir -p $(3DDSPHDIR)
	@mkdir -p $(3DDSPHDIR)/shapes
	@mkdir -p $(3DALEDIR)
	@mkdir -p $(3DALEDIR)/shapes

# dsphdirs : 
# 	@mkdir -p $(2DDSPHDIR)
# 	@mkdir -p $(2DDSPHDIR)/shapes
# 	@mkdir -p $(3DDSPHDIR)
# 	@mkdir -p $(3DDSPHDIR)/shapes


# SET2D:
# 	$(eval DIM=-DSIMDIM=2)

# SET3D:
# 	$(eval DIM=-DSIMDIM=3)
	
print-%  : ; @echo $* = $($*)

debug:
	$(eval CXXFLAGS=-g $(CXXFLAGS) -DDEBUG)

clean: 
	$(RM) $(2DOBJ)
	$(RM) $(3DOBJ)
	$(RM) $(2DDEP)
	$(RM) $(3DDEP)
	$(RM) $(TARGETDIR)/2DSPH
	$(RM) $(TARGETDIR)/3DSPH
	$(RM) $(2DALEOBJ)
	$(RM) $(3DALEOBJ)
	$(RM) $(2DDEPA)
	$(RM) $(3DDEPA)
	$(RM) $(TARGETDIR)/2DSPH.$(ALEEXT)
	$(RM) $(TARGETDIR)/3DSPH.$(ALEEXT)

cleanDSPH:
	$(RM) $(2DOBJ)
	$(RM) $(3DOBJ)
	$(RM) $(2DDEP)
	$(RM) $(3DDEP)
	$(RM) $(TARGETDIR)/2DSPH
	$(RM) $(TARGETDIR)/3DSPH

cleanALE:
	$(RM) $(2DALEOBJ)
	$(RM) $(3DALEOBJ)
	$(RM) $(2DDEPA)
	$(RM) $(3DDEPA)
	$(RM) $(TARGETDIR)/2DSPH.$(ALEEXT)
	$(RM) $(TARGETDIR)/3DSPH.$(ALEEXT)
	
2D: 2ddirs 2DDSPH 2DALE

3D: 3ddirs 3DDSPH 3DALE

-include $(2DDEP)
2DDSPH : $(2DOBJ)
	$(CXX) -DSIMDIM=2 $(FLAGS) $(CXXFLAGS) $(TECINC) $(EIGENINC) $(TECLINK) -o $(TARGETDIR)/2DSPH $^ $(TECLIB) $(NETLIB)

-include $(2DDEPA)
2DALE : $(2DALEOBJ)
	$(CXX) -DSIMDIM=2 -DALE $(FLAGS) $(CXXFLAGS) $(TECINC) $(EIGENINC) $(TECLINK) -o $(TARGETDIR)/2DSPH.$(ALEEXT)  $^ $(TECLIB) $(NETLIB)

-include $(3DDEP)
3DDSPH : $(3DOBJ)
	$(CXX) -DSIMDIM=3 $(FLAGS) $(CXXFLAGS) $(TECINC) $(EIGENINC) $(TECLINK) -o $(TARGETDIR)/3DSPH $^ $(TECLIB) $(NETLIB)

-include $(3DDEPA)
3DALE : $(3DALEOBJ)
	$(CXX) -DSIMDIM=3 -DALE $(FLAGS) $(CXXFLAGS) $(TECINC) $(EIGENINC) $(TECLINK) -o $(TARGETDIR)/3DSPH.$(ALEEXT) $^ $(TECLIB) $(NETLIB)


$(2DDSPHDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) -DSIMDIM=2 $(FLAGS) $(CXXFLAGS) $(TECINC) $(EIGENINC) $(TECLINK) -MMD -MP -c $< -o $@ $(TECLIB) $(NETLIB)	

$(2DALEDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) -DSIMDIM=2 -DALE $(FLAGS) $(CXXFLAGS) $(TECINC) $(EIGENINC) $(TECLINK) -MMD -MP -c $< -o $@ $(TECLIB) $(NETLIB)

$(3DDSPHDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) -DSIMDIM=3 $(FLAGS) $(CXXFLAGS) $(TECINC) $(EIGENINC) $(TECLINK) -MMD -MP -c $< -o $@ $(TECLIB) $(NETLIB)	

$(3DALEDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) -DSIMDIM=3 -DALE $(FLAGS) $(CXXFLAGS) $(TECINC) $(EIGENINC) $(TECLINK) -MMD -MP -c $< -o $@ $(TECLIB) $(NETLIB)

c2f:
	$(CXX) $(FLAGS) $(CXXFLAGS) $(INC) $(LLINK) $(TECINC) $(TECLINK) -o $(TARGETDIR)/Cell2Face src/Cell2Face.cpp $(TECLIB) $(NETLIB)

c2fD:
	$(CXX) -g -DDEBUG $(FLAGS) $(CXXFLAGS) $(INC) $(LLINK) -o $(TARGETDIR)/Cell2Face src/Cell2Face.cpp $(NETLIB)

c2e:
	$(CXX) $(FLAGS) $(CXXFLAGS) $(INC) $(LLINK) $(TECINC) $(TECLINK) -o $(TARGETDIR)/Cell2Edge src/Cell2Edge.cpp $(TECLIB) $(NETLIB)

c2eD:
	$(CXX) -g -DDEBUG $(FLAGS) $(CXXFLAGS) $(INC) $(LLINK) $(TECINC) $(TECLINK) -o $(TARGETDIR)/Cell2Edge src/Cell2Edge.cpp $(TECLIB) $(NETLIB)

tau2zcfd:
	$(CXX) -g -DDEBUG $(FLAGS) $(CXXFLAGS) $(H5INC) -o $(TARGETDIR)/TAU2zCFD src/TAUtozCFD.cpp $(NETLIB) $(H5LIB) -lhdf5

# mesh:
# 	$(CXX) $(FLAGS) $(CXXFLAGS) $(INC) -o $(TARGETDIR)/MakeMesh src/MakeMesh.cpp -lnetcdf

# meshD:
# 	$(CXX) -g -DDEBUG $(FLAGS) $(CXXFLAGS) $(INC) -o MakeMesh src/MakeMesh.cpp -lnetcdf