# Compiler to use. In C++, so g++.
CXX=g++

# Libraries to include
INC=

# Compiler flags. If desired add -g for debugging info.
CFLAGS=-g -std=c++11 -Wall -Wextra -lm -ffast-math -funroll-loops -O3

# Find the OS to execute correctly
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
TARGET=WCXSPH# Target executable
3DTARGET=WCXSPH3D
else
TARGET=WCXSPH.exe
3DTARGET=WCXSPH3D.exe
endif

SOURCE=src/WCXSPH.cpp
3DTARGET=WCXSPH3D
3DSOURCE=src/WCXSPH3D.cpp



# Compile and Build
build: $(SOURCE)
	$(CXX) $(INC) $(CFLAGS) -o $(TARGET) $(SOURCE)

#Crossflow run with different input file
cross:
	./$(TARGET) Cross.dat SurfaceLeft.plt

# Linux run without .exe
linux:
	./$(TARGET) Input.dat

clean:
	$(RM) $(TARGET)

new: WCXSPH.cpp
	$(RM) $(TARGET)
	$(CXX) $(INC) $(CFLAGS) -o $(TARGET) $(SOURCE)

3D:
	$(CXX) $(INC) $(CFLAGS) -o $(3DTARGET) $(3DSOURCE)

3Drun:
	./$(3DTARGET).exe 3DInput.dat
