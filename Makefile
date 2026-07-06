# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++17 -Os -fopenmp

# Linker flags
LIBS = -lboost_iostreams -fopenmp

# Target executable
TARGET = RecombMix

# Source file
SRCS = ./src/RecombMix.cpp

# Object file
OBJS = $(SRCS:.cpp=.o)

# Link
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LIBS)

# Compile
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Run program
run: $(TARGET)
	./$(TARGET)

# Clean build
clean:
	rm -f $(OBJS) $(TARGET)