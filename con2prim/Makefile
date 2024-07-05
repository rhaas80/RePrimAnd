REPRIMAND_HOME := $(PWD)/../../reprimand-build

CXX=g++

CXXFLAGS = -std=c++11 -Wall -Wextra -fopenmp
CPPFLAGS = -I$(REPRIMAND_HOME)/include
LIBS = -lRePrimAnd

# multi-arch stuff
MACHINE := $(shell dpkg-architecture -qDEB_HOST_MULTIARCH 2>/dev/null || gcc -print-multiarch 2>/dev/null || gcc -dumpmachine 2>/dev/null)

LDFLAGS = -L$(REPRIMAND_HOME)/lib/$(MACHINE) -Wl,--rpath,$(REPRIMAND_HOME)/lib/$(MACHINE) -fopenmp

.PHONY: clean

minimal: minimal.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) minimal.cc -o minimal $(LIBS)

clean:
	rm minimal
