# Compiler to use. In C++, so g++.
CXX=g++

CX8=g++-8

# Libraries to include
INC=

# Compiler flags. If desired add -g for debugging info.
CFLAGS=-std=c++11 -Wall -Wextra -ffast-math -funroll-loops -O3 -fopenmp

# Find the OS to execute correctly
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
TARGET=WCSPH# Target executable
3DTARGET=WCSPH
else
TARGET=WCSPH.exe
3DTARGET=WCSPH.exe
endif

SOURCE=src/2D/WCSPH.cpp
#3DTARGET=WCXSPH
3DSOURCE=src/3D/WCSPH.cpp


# Compile and Build
2D:
	$(CXX) $(INC) $(CFLAGS) -o $(TARGET) $(SOURCE)

3D:
	$(CXX) $(INC) $(CFLAGS) -o $(3DTARGET) $(3DSOURCE)

#add debug flag
debug2:
	$(CXX) $(INC) -g $(CFLAGS) -o $(TARGET) $(SOURCE)

debug3:
	$(CXX) $(INC) -g $(CFLAGS) -o $(3DTARGET) $(3DSOURCE)

clean:
	$(RM) $(TARGET)

new: 
	$(RM) $(TARGET)
	$(CXX) $(INC) $(CFLAGS) -o $(TARGET) $(SOURCE)
