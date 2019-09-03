# Compiler to use. In C++, so g++.
CXX=g++

CX8=g++-8

# Libraries to include
LIBS=-ltecio -lgmpxx -lgmp
#-lnetcdf_c++4

# Compiler flags. If desired add -g for debugging info.
CFLAGS=-std=c++11 -Wall -Wextra -ffast-math -funroll-loops -O3 -fopenmp

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
	$(CXX) -DSIMDIM=2 $(CFLAGS) $(SOURCE) $(LIBS) -o $(TARGET)

3D:
	$(CXX) -DSIMDIM=3 $(CFLAGS) $(SOURCE) $(LIBS) -o $(TARGET)

#add debug flag
debug2:
	$(CXX) -g -DSIMDIM=2 $(CFLAGS) $(SOURCE) $(LIBS) -o $(TARGET)

debug3:
	$(CXX) -g -DSIMDIM=3 $(CFLAGS) $(SOURCE) $(LIBS) -o $(TARGET)

clean:
	$(RM) $(TARGET)

new: 
	$(RM) $(TARGET)
	$(CXX) $(CFLAGS) -o $(TARGET) $(SOURCE) $(LIBS)
