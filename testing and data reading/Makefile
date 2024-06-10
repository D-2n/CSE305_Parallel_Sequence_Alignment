# Compiler
CXX = g++
CXXFLAGS = -Wall -std=c++11

# Executable names
TARGET = pull_data
TEST_TARGET = testing

# Source and header files
SRCS = pull_data.cpp
TEST_SRCS = testing.cpp
HEADERS = read_test_data.h

# Libraries
LIBS = -lcurl

# Default target
all: $(TEST_TARGET)

# Compile pull_data.cpp to create object file
pull_data.o: pull_data.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c pull_data.cpp $(LIBS)

# Compile testing.cpp and link with pull_data.o to create testing executable
$(TEST_TARGET): $(TEST_SRCS) pull_data.o $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(TEST_TARGET) $(TEST_SRCS) pull_data.o $(LIBS)

# Clean target
clean:
	rm -f $(TARGET) $(TEST_TARGET) pull_data.o
