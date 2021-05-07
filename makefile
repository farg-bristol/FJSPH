# Compiler to use. In C++, so g++.
CC=gcc
CXX=g++

CX8=g++-8

# Libraries to include
LIBS=-ltecio -lnetcdf #-lH5hut -lhdf5
#-lgmpxx -lgmp
# Compiler flags. 
CXXFLAGS=-std=c++11 -Wall -ffast-math -funroll-loops -O3 -fopenmp

INC=-I/usr/local/tecplot/360ex_2020r2/include
LLINK=-L/usr/local/tecplot/360ex_2020r2/bin
ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
# Find the OS to execute correctly
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
#TARGET=$(ROOT_DIR)/WCSPH# Target executable
TARGET=WCSPH
else
#TARGET=$(ROOT_DIR)/WCSPH.exe
TARGET=WCSPH.exe
endif

#SOURCE=$(ROOT_DIR)/src/WCSPH.cpp
SOURCE=src/WCSPH.cpp


# Compile and Build
2D:
	$(CXX) -DSIMDIM=2 ${FLAGS} $(CXXFLAGS) $(INC) $(LLINK) -o $(TARGET) $(SOURCE) $(LIBS)

3D:
	$(CXX) -DSIMDIM=3 ${FLAGS} $(CXXFLAGS) $(INC) $(LLINK) -o $(TARGET) $(SOURCE) $(LIBS)

NOALE2D:
	$(CXX) -DSIMDIM=2 -DNOALE ${FLAGS} $(CXXFLAGS) $(INC) $(LLINK) -o $(TARGET) $(SOURCE) $(LIBS)

NOALE3D:
	$(CXX) -DSIMDIM=3 -DNOALE ${FLAGS} $(CXXFLAGS) $(INC) $(LLINK) -o $(TARGET) $(SOURCE) $(LIBS)
#add debug flag
debug2:
	$(CXX) -g -DDEBUG -DSIMDIM=2 ${FLAGS} $(CXXFLAGS) $(INC) $(LLINK) -o $(TARGET) $(SOURCE) $(LIBS)

debug3:
	$(CXX) -g -DDEBUG -DSIMDIM=3 ${FLAGS} $(CXXFLAGS) $(INC) $(LLINK) -o $(TARGET) $(SOURCE) $(LIBS)

post:
	$(CXX) ${FLAGS} $(CXXFLAGS) $(INC) $(LLINK) src/Post.cpp -o SPHPost $(LIBS)

debugP:
	$(CXX) -g -DDEBUG ${FLAGS} $(CXXFLAGS) $(INC) $(LLINK) -o SPHPost src/Post.cpp $(LIBS)

c2f:
	$(CXX) $(FLAGS) $(CXXFLAGS) $(INC) $(LLINK) -o Cell2Face src/Cell2Face.cpp $(LIBS)

c2fD:
	$(CXX) -g -DDEBUG $(FLAGS) $(CXXFLAGS) $(INC) $(LLINK) -o Cell2Face src/Cell2Face.cpp $(LIBS)

c2e:
	$(CXX) $(FLAGS) $(CXXFLAGS) $(INC) $(LLINK) -o Cell2Edge src/Cell2Edge.cpp $(LIBS)

c2eD:
	$(CXX) -g -DDEBUG $(FLAGS) $(CXXFLAGS) $(INC) $(LLINK) -o Cell2Edge src/Cell2Edge.cpp $(LIBS)

mesh:
	$(CXX) $(FLAGS) $(CXXFLAGS) $(INC) -o MakeMesh src/MakeMesh.cpp -lnetcdf_c++4

meshD:
	$(CXX) -g -DDEBUG $(FLAGS) $(CXXFLAGS) $(INC) -o MakeMesh src/MakeMesh.cpp -lnetcdf_c++4

clean:
	$(RM) $(TARGET)

