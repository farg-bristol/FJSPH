# Compiler to use. In C++, so g++.
CXX=g++

CX8=g++-8

# Libraries to include
LIBS=-ltecio -lnetcdf_c++4 #-lH5hut -lhdf5
#-lgmpxx -lgmp
# Compiler flags. If desired add -g for debugging info.
CFLAGS=-std=c++11 -Wall -ffast-math -funroll-loops -O3 -fopenmp

INC=#-I/usr/include/hdf5/serial/
LLINK=#-L/usr/lib/x86_64-linux-gnu/hdf5/serial/

# Find the OS to execute correctly
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
TARGET=WCSPH# Target executable
else
TARGET=WCSPH.exe
endif

SOURCE=src/WCSPH.cpp


# Compile and Build
2D:
	$(CXX) -DSIMDIM=2 ${FLAGS} $(CFLAGS) $(INC) $(LLINK) $(SOURCE) $(LIBS) -o $(TARGET)

3D:
	$(CXX) -DSIMDIM=3 ${FLAGS} $(CFLAGS) $(INC) $(LLINK) $(SOURCE) $(LIBS) -o $(TARGET)

#add debug flag
debug2:
	$(CXX) -g -DDEBUG -DSIMDIM=2 ${FLAGS} $(CFLAGS) $(INC) $(LLINK) $(SOURCE) $(LIBS) -o $(TARGET)

debug3:
	$(CXX) -g -DDEBUG -DSIMDIM=3 ${FLAGS} $(CFLAGS) $(INC) $(LLINK) $(SOURCE) $(LIBS) -o $(TARGET)

post:
	$(CXX) ${FLAGS} $(CFLAGS) $(INC) $(LLINK) src/Post.cpp $(LIBS) -o SPHPost

debugP:
	$(CXX) -g -DDEBUG ${FLAGS} $(CFLAGS) $(INC) $(LLINK) src/Post.cpp $(LIBS) -o SPHPost

c2f:
	$(CXX) $(FLAGS) $(CFLAGS) $(INC) $(LLINK) src/Cell2Face.cpp $(LIBS) -o Cell2Face

c2fD:
	$(CXX) -g -DDEBUG $(FLAGS) $(CFLAGS) $(INC) $(LLINK) src/Cell2Face.cpp $(LIBS) -o Cell2Face

c2e:
	$(CXX) $(FLAGS) $(CFLAGS) $(INC) $(LLINK) src/Cell2Edge.cpp $(LIBS) -o Cell2Edge

c2eD:
	$(CXX) -g -DDEBUG $(FLAGS) $(CFLAGS) $(INC) $(LLINK) src/Cell2Edge.cpp $(LIBS) -o Cell2Edge

mesh:
	$(CXX) $(FLAGS) $(CFLAGS) $(INC) src/MakeMesh.cpp -lnetcdf_c++4 -o MakeMesh

meshD:
	$(CXX) -g -DDEBUG $(FLAGS) $(CFLAGS) $(INC) src/MakeMesh.cpp -lnetcdf_c++4 -o MakeMesh

clean:
	$(RM) $(TARGET)

new: 
	$(RM) $(TARGET)
	$(CXX) $(CFLAGS) -o $(TARGET) $(SOURCE) $(LIBS)
