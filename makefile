# Compiler to use. In C++, so g++.
CC=gcc
CXX=g++

# Libraries to include
LIBS=-ltecio -lnetcdf 

# Compiler flags. 
CXXFLAGS=-std=c++11 -Wall -ffast-math -funroll-loops -O3 -fopenmp

INC=-I/usr/local/tecplot/360ex_2020r2/include
LLINK=-L/usr/local/tecplot/360ex_2020r2/bin

TARGETDIR=~/bin

SOURCE=src/WCSPH.cpp


# Compile and Build
2D:
	$(CXX) -DSIMDIM=2 ${FLAGS} $(CXXFLAGS) $(INC) $(LLINK) -o $(TARGETDIR)/2DFJSPH $(SOURCE) $(LIBS)

3D:
	$(CXX) -DSIMDIM=3 ${FLAGS} $(CXXFLAGS) $(INC) $(LLINK) -o $(TARGETDIR)/3DFJSPH $(SOURCE) $(LIBS)

#add debug flag
debug2:
	$(CXX) -g -DDEBUG -DSIMDIM=2 ${FLAGS} $(CXXFLAGS) $(INC) $(LLINK) -o $(TARGETDIR)/2DFJSPH $(SOURCE) $(LIBS)

debug3:
	$(CXX) -g -DDEBUG -DSIMDIM=3 ${FLAGS} $(CXXFLAGS) $(INC) $(LLINK) -o $(TARGETDIR)/3DFJSPH $(SOURCE) $(LIBS)

post:
	$(CXX) ${FLAGS} $(CXXFLAGS) $(INC) $(LLINK) src/Post.cpp -o $(TARGETDIR)/SPHPost $(LIBS)

debugP:
	$(CXX) -g -DDEBUG ${FLAGS} $(CXXFLAGS) $(INC) $(LLINK) -o $(TARGETDIR)/SPHPost src/Post.cpp $(LIBS)

c2f:
	$(CXX) $(FLAGS) $(CXXFLAGS) $(INC) $(LLINK) -o $(TARGETDIR)/Cell2Face src/Cell2Face.cpp $(LIBS)

c2fD:
	$(CXX) -g -DDEBUG $(FLAGS) $(CXXFLAGS) $(INC) $(LLINK) -o $(TARGETDIR)/Cell2Face src/Cell2Face.cpp $(LIBS)

c2e:
	$(CXX) $(FLAGS) $(CXXFLAGS) $(INC) $(LLINK) -o $(TARGETDIR)/Cell2Edge src/Cell2Edge.cpp $(LIBS)

c2eD:
	$(CXX) -g -DDEBUG $(FLAGS) $(CXXFLAGS) $(INC) $(LLINK) -o $(TARGETDIR)/Cell2Edge src/Cell2Edge.cpp $(LIBS)

clean:
	$(RM) $(TARGETDIR)/2DFJSPH
	$(RM) $(TARGETDIR)/3DFJSPH
	$(RM) $(TARGETDIR)/SPHPost
	$(RM) $(TARGETDIR)/Cell2Face
	$(RM) $(TARGETDIR)/Cell2Edge

