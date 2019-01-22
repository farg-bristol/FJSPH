# Compiler to use. In C++, so g++.
CXX=g++

# Libraries to include
INC= 

# Compiler flags. If desired add -g for debugging info.
CFLAGS=-g -std=c++11 -Wall -Wextra -fopenmp -lm -ffast-math -funroll-loops -O3

# Target executable
TARGET=WCXSPH
SOURCE=WCXSPH.cpp
3DTARGET=WCXSPH3D
3DSOURCE=WCXSPH3D.cpp

# Compile and Build
build: WCXSPH.cpp
	$(CXX) $(INC) $(CFLAGS) -o $(TARGET) $(SOURCE) 

# Windows run .exe
win:
	./$(TARGET).exe Input.dat 

#Crossflow run with different input file
cross:
	./$(TARGET).exe Cross.dat SurfaceLeft.plt

# Linux run without .exe
linux: 
	./$(TARGET) Input.dat

clean:
	$(RM) $(TARGET).exe

new: WCXSPH.cpp
	$(RM) $(TARGET)
	$(CXX) $(INC) $(CFLAGS) -o $(TARGET) $(SOURCE) 

3D: 
	$(CXX) $(INC) $(CFLAGS) -o $(3DTARGET) $(3DSOURCE)

3Drun:
	./$(3DTARGET).exe 3DInput.dat





