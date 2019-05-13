# Compiler to use. In C++, so g++.
CXX=g++

CX8=g++-8

# Libraries to include
INC=

# Compiler flags. If desired add -g for debugging info.
CFLAGS=-std=c++11 -Wall -Wextra -ffast-math -funroll-loops -O3

# Find the OS to execute correctly
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
TARGET=WCXSPH# Target executable
3DTARGET=WCXSPH
else
TARGET=WCXSPH.exe
3DTARGET=WCXSPH.exe
endif

SOURCE=src/2D/WCXSPH.cpp
#3DTARGET=WCXSPH
3DSOURCE=src/3D/WCXSPH.cpp


# Compile and Build
2D:
	$(CXX) $(INC) $(CFLAGS) -o $(TARGET) $(SOURCE)

3D:
	$(CXX) $(INC) -g $(CFLAGS) -o $(3DTARGET) $(3DSOURCE)

#add debug flag
debug:
	$(CXX) $(INC) -g $(CFLAGS) -o $(TARGET) $(SOURCE)

test:
	$(CX8) $(INC) -g $(CFLAGS) -o $(TARGET) $(SOURCE)

clean:
	$(RM) $(TARGET)

new: 
	$(RM) $(TARGET)
	$(CXX) $(INC) $(CFLAGS) -o $(TARGET) $(SOURCE)
