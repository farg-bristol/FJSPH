# Compiler to use. In C++, so g++.
CXX=g++

# Libraries to include
INC= 

# Compiler flags. If desired add -g for debugging info.
CFLAGS= -g -std=c++11 -Wall -Wextra -fopenmp -lm -ffast-math -funroll-loops -O3

# Target executable
TARGET = WCXSPH

# Compile and Build
build: WCXSPH.cpp
	$(CXX) $(INC) $(CFLAGS) -o $(TARGET) $< 

# Windows run .exe
win:
	./$(TARGET).exe Input.dat 

#Crossflow run with different input file
cross:
	./$(TARGET).exe Cross.dat

# Linux run without .exe
linux: 
	./$(TARGET) Input.dat

clean:
	$(RM) $(TARGET).exe

new: WCXSPH.cpp
	$(RM) $(TARGET)
	$(CXX) $(INC) $(CFLAGS) -o $(TARGET) $< 

3D	: WCXSPH3D.cpp
	$(CXX) $(INC) $(CFLAGS) -o $(TARGET) $< 

3Drun:
	./$(TARGET).exe 3DInput.dat





