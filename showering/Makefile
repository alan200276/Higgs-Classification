# --------------------------------------------- #
# Makefile for myexample code                   #
# Pascal Nef, updated by Ben Nachman (2017)     #
#                                               #
# Note: source setup.sh before make             #
# --------------------------------------------- #

CXXFLAGS = -g -std=c++0x  -O2 -Wall 
LDFLAGS = -g

.PHONY: clean debug all

all: myexample myexample_ghost Ghost_Associated_B_hadron

myexample:  myexample.so  myexampleAnalysis.so 
	$(CXX) myexample.so myexampleAnalysis.so -o $@.exe \
	$(CXXFLAGS) -Wno-shadow  \
	`root-config --glibs`  \
	-L$(FASTJETLOCATION)/lib `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins ` -lRecursiveTools -lQCDAwarePlugin \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -ldl \
		

myexample.so: myexample.C  myexampleAnalysis.so   
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags` 

myexampleAnalysis.so : myexampleAnalysis.cc myexampleAnalysis.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins`  \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags`   
    
    
    
myexample_ghost:  myexample_ghost.so  myexampleAnalysis_ghost.so 
	$(CXX) myexample_ghost.so myexampleAnalysis_ghost.so -o $@.exe \
	$(CXXFLAGS) -Wno-shadow  \
	`root-config --glibs`  \
	-L$(FASTJETLOCATION)/lib `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins ` -lRecursiveTools -lQCDAwarePlugin \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -ldl \
		

myexample_ghost.so: myexample_ghost.C  myexampleAnalysis_ghost.so   
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags` 

myexampleAnalysis_ghost.so : myexampleAnalysis_ghost.cc myexampleAnalysis_ghost.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins`  \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags` 
    

Ghost_Associated_B_hadron:  Ghost_Associated_B_hadron.so  myexampleAnalysis_ghost.so 
	$(CXX) Ghost_Associated_B_hadron.so Ghost_Associated_B_hadron_Analysis.so -o $@.exe \
	$(CXXFLAGS) -Wno-shadow  \
	`root-config --glibs`  \
	-L$(FASTJETLOCATION)/lib `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins ` -lRecursiveTools -lQCDAwarePlugin \
	-L$(PYTHIA8LOCATION)/lib -lpythia8 -ldl \
		

Ghost_Associated_B_hadron.so: Ghost_Associated_B_hadron.C  Ghost_Associated_B_hadron_Analysis.so   
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags` 

Ghost_Associated_B_hadron_Analysis.so : Ghost_Associated_B_hadron_Analysis.cc Ghost_Associated_B_hadron_Analysis.h
	$(CXX) -o $@ -c $<  \
	$(CXXFLAGS) -Wno-shadow -fPIC -shared \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins`  \
	-I$(PYTHIA8LOCATION)/include \
	`root-config --cflags` 

    
    

clean:
	rm -rf *.exe
	rm -rf *.o
	rm -rf *.so
	rm -f *~
