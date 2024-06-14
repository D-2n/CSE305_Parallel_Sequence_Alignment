CXX = g++

CXXFLAGS = -Wall -Iclasses -Ifunctions

SRCDIRS = . test_functions alignment_algorithm sequence_alignment


SOURCES := $(wildcard $(addsuffix /*.cpp,$(SRCDIRS)))

OBJECTS := $(SOURCES:.cpp=.o)


TARGET = testing

all: $(TARGET)


$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^


%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<


clean:
	rm -f $(TARGET) $(OBJECTS)


.PHONY: all clean
