REPRIMAND_HOME := $(PWD)/../../reprimand-build

CXX=g++

CXXFLAGS = -std=c++11 -Wall -Wextra -fopenmp -g
CPPFLAGS = -I$(REPRIMAND_HOME)/include
LIBS = -lRePrimAnd -lgsl -lgslcblas

# multi-arch stuff
MACHINE := $(shell dpkg-architecture -qDEB_HOST_MULTIARCH 2>/dev/null || gcc -print-multiarch 2>/dev/null || gcc -dumpmachine 2>/dev/null)

LDFLAGS = -fopenmp -g
# add both lib and lib64 directories, one never knows what the OS / meson used
LDFLAGS += -L$(REPRIMAND_HOME)/lib/$(MACHINE) -Wl,--rpath,$(REPRIMAND_HOME)/lib/$(MACHINE)
LDFLAGS += -L$(REPRIMAND_HOME)/lib64/$(MACHINE) -Wl,--rpath,$(REPRIMAND_HOME)/lib64/$(MACHINE)
# GSL may not set all required paths
LDFLAGS += -L$(GSL_ROOT_DIR)/lib -Wl,--rpath,$(GSL_ROOT_DIR)/lib

.PHONY: clean

minimal: minimal.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) minimal.cc -o minimal $(LIBS)

clean:
	rm minimal
