# Simple Makefile for TreeReader project
CXX = g++
CXXFLAGS = -g -O2 -Wall $(shell root-config --cflags)
LIBS = $(shell root-config --libs)

#Executable Names
TARGET1 = ParticleFluxCalculation
TARGET2 = DoseCalculation
TARGET3 = InteractionAnalysis
# Source Files
SRC1 = ParticleFluxCalculation.cxx
SRC2 = DoseCalculation.cxx
SRC3 = InteractionAnalysis.cxx

#Default target
all: $(TARGET1) $(TARGET2) $(TARGET3)
#all: DoseCalculation

$(TARGET1): $(SRC1) TreeReader.cxx TreeReader.h
	$(CXX) $(CXXFLAGS) -o $(TARGET1) $(SRC1) TreeReader.cxx $(LIBS)

$(TARGET2): $(SRC2) TreeReader.cxx TreeReader.h
	$(CXX) $(CXXFLAGS) -o $(TARGET2) $(SRC2) TreeReader.cxx $(LIBS)
	
$(TARGET3): $(SRC2) TreeReader.cxx TreeReader.h
	$(CXX) $(CXXFLAGS) -o $(TARGET3) $(SRC3) TreeReader.cxx $(LIBS)


#DoseCalculation: DoseCalculation.cxx TreeReader.cxx TreeReader.h
#	$(CXX) $(CXXFLAGS) -o DoseCalculation DoseCalculation.cxx TreeReader.cxx $(LIBS)

clean:
	#rm -f DoseCalculation *.o
	rm -f $(TARGET1) $(TARGET2) $(TARGET3) *.o
	

test: $(TARGET2)
	./$(TARGET2)

.PHONY: all clean test
