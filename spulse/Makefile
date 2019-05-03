#OBJ = DSSD2Tr_March TSCorrSiGe DSSDGeMerge separate_March BuildDecay DSSDGeCorr PLASTIC2Tr
#OBJ = online_DSSD WASABI2Tr
#OBJ = enthe 
#sourcefile = main.cpp enthe.cpp enthe.h 

EXE = mapping
sourcefile = main.cpp mapping.h mapping.cpp
OBJ = $(subst,.cpp,.o,$(filter %.cpp,$(sourcefile)))

ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs)

CFLAGS = -Wall -O2 -I$(ROOTSYS)/include -lXMLParser -lSpectrum
#-lanacore

GXX = g++ 
#-std=c++0x

all:$(EXE)

$(EXE):main.o mapping.o
	$(GXX) $(ROOTCFLAGS) $(ROOTLIBS) -o $(EXE) main.o mapping.o
	

main.o: main.cpp mapping.h
	$(GXX)  $(ROOTCFLAGS) $(ROOTLIBS) -c main.cpp

mapping.o: mapping.cpp mapping.h
	$(GXX) $(ROOTCFLAGS) $(ROOTLIBS) -c mapping.cpp
clean:
	rm -f *~ *.o $(EXE)
