CXX := g++

CXXFLAGS := -Wall -g -Wno-sign-compare -O3 -I./Eigen/

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin) # Mac
	LDFLAGS := -framework GLUT -framework OpenGL
else # Linux?
	LDFLAGS := -lglut -lGLU -lGL
endif

OBJ := main.o fluid.o grid.o sparse.o

run: $(OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

# Nicked from http://www.gnu.org/software/make/manual/make.html#Automatic-Prerequisites
%.d: %.cpp
	@set -e; rm -f $@; \
	$(CXX) -MM $(CXXFLAGS) $< -o $@.tmp; \
	sed 's,\($*\)\.o[ :]*,\1.o: ,g' < $@.tmp > $@; \
	rm -f $@.tmp

-include $(OBJ:.o=.d)

clean:
	rm -f *.o *.d run
