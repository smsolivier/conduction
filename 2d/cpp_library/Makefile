SRCS = 
SRCS += helper.cc 
SRCS += VectorMath.cc 
SRCS += LinearSolver.cc 
SRCS += table.cc 
SRCS += Timer.cc 
SRCS += WriteCurve.cc
SRCS += GaussQuad.cc
SRCS += Polynomial.cc

# SRCS_GPU = MatrixVectorProduct.cu 

OBJ = $(subst .cc,.o,$(SRCS))
OBJ_GPU = $(subst .cu,.o,$(SRCS_GPU))

MYLIB = libMyLib.a 

all : $(OBJ) $(OBJ_GPU) $(MYLIB)

%.o : %.cc %.hh
	g++ -c -fopenmp -std=c++14 $<

%.o : %.cu
	nvcc -c -std=c++11 $<

$(MYLIB) : $(OBJ) $(OBJ_GPU)
	ar rcs $(MYLIB) $(OBJ) $(OBJ_GPU)

.PHONY : clean 
clean :
	rm -f $(OBJ) $(OBJ_GPU) $(MYLIB)
