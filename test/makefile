
# OS Name (Linux or Darwin)
OSUPPER = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])
OSLOWER = $(shell uname -s 2>/dev/null | tr [:upper:] [:lower:])

# Flags to detect 32-bit or 64-bit OS platform
OS_SIZE = $(shell uname -m | sed -e "s/i.86/32/" -e "s/x86_64/64/")
OS_ARCH = $(shell uname -m | sed -e "s/i386/i686/")

# Flags to detect either a Linux system (linux) or Mac OSX (darwin)
DARWIN = $(strip $(findstring DARWIN, $(OSUPPER)))

# OS-specific build flags
ifneq ($(DARWIN),) 
      CCFLAGS   = -arch $(OS_ARCH)
      ARMA_INC  = -I/opt/local/include
      ARMA_LIB  = -L/opt/local/lib -larmadillo -framework Accelerate -llapack
else
      ARMA_INC  = -I/usr/include
  ifeq ($(OS_SIZE),32)
      CCFLAGS   = -m32
      ARMA_LIB = -L/usr/lib
  else
      CCFLAGS   = -m64
      ARMA_LIB = -L/usr/lib64
  endif
      ARMA_LIB += -larmadillo -lblas -llapack
endif

# Debug build flags
ifeq ($(DEBUGVERSION),1)
      CCFLAGS   += -DDEBUG=1
else
      CCFLAGS   += -DDEBUG=0
endif 


all: build

CXX = g++
OPTIM_CXXFLAGS   = -O2 -funroll-loops
WARN_CXXFLAGS   = -Wall -pedantic
CXX_FLAGS =  -std=c++11 $(CCFLAGS) $(OPTIM_CXXFLAGS) $(WARN_CXXFLAGS)
USEGNUPLOTINTERP := 0
ifneq ($(USEGNUPLOTINTERP),0)
   CXX_FLAGS += -DUSEGNUPLOTINTERP=$(USEGNUPLOTINTERP)
endif

build: reconstructor


mockwavefrontgenerator.o: ../src/mockwavefrontgenerator.cpp ../src/mockwavefrontgenerator.h
	$(CXX) $(CXX_FLAGS) -o $@ -c $<

wavefrontreconstructor.o: ../src/wavefrontreconstructor.cpp \
	                       ../src/wavefrontreconstructor.h
	$(CXX) $(CXX_FLAGS) $(ARMA_INC) -o $@ -c $<

legsqwfreconstructor.o: ../src/legsqwfreconstructor.cpp \
	                           ../src/legsqwfreconstructor.h \
										../src/wavefrontreconstructor.h
	$(CXX) $(CXX_FLAGS) $(ARMA_INC) -o $@ -c $<

zernikewfreconstructor.o: ../src/zernikewfreconstructor.cpp \
	                           ../src/zernikewfreconstructor.h \
										../src/wavefrontreconstructor.h
	$(CXX) $(CXX_FLAGS) $(ARMA_INC) -o $@ -c $<

legendrepolynomial.o: ../src/legendrepolynomial.cpp ../src/legendrepolynomial.h
	$(CXX) $(CXX_FLAGS) -o $@ -c $<

zernikepolynomials.o: ../src/zernikepolynomials.cpp ../src/zernikepolynomials.h
	$(CXX) $(CXX_FLAGS) -o $@ -c $<

reconstructor.o: reconstructor.cpp
	$(CXX) $(CXX_FLAGS) $(ARMA_INC) -o $@ -c $<

reconstructor: reconstructor.o mockwavefrontgenerator.o zernikewfreconstructor.o \
	        zernikepolynomials.o legsqwfreconstructor.o legendrepolynomial.o \
			  wavefrontreconstructor.o
	$(CXX) $(ARMA_LIB) -o $@ $+

runtest: build
	./reconstructor

clean:
	rm -f *.o
	rm -f reconstructor



