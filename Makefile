INSTALLDIR=/usr/local/bin

src = $(wildcard *.cpp *.c *.cxx sweep/*.cc)
csrc = $(wildcard *.c)
obj_cpp = $(src:.cpp=.o)
	obj = $(obj_cpp:.c=.o)
	hdr = $(wildcard *.h *.hpp)

CC=g++
COPT=-g -lm

CXX=g++
CXXOPT = -g
#CXXOPT = -O3
CXXLIB = 

%.o: %.c %.h %.cc sweep/%.cc
	$(CC) $(COPT) -c $< -o $@ -lm

%.o: %.cpp %.hpp
	$(CXX) $(CXXOPT) -c $< -o $@

%.o: %.cxx %.hxx
	$(CC) $(COPT) -c $< -o $@ -lm


weakpwh: clipper.cpp sweep/advancing_front.cc sweep/cdt.cc sweep/sweep.cc sweep/sweep_context.cc common/shapes.cc weakpwh.cpp
	$(CXX) $(CXXOPT) $^  -o $@ $(CXXLIBS)

.PHONY: debug
debug: $(obj) $(hdr)
	$(CXX) -g -o weakpwh_debug $^ $(CXXLIBS)

.PHONY: clean
clean:
	rm -f *.o weakpwh 

.PHONY: install
install: weakpwh
	mkdir -p $(INSTALLDIR)
	cp $< $(INSTALLDIR)

