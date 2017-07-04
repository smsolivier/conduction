CXX=g++
CFLAGS = -std=c++14 \
	-fopenmp 

LIBLOC = /home/sam/cpp
LIBS = MyLib
LIBS := $(addprefix -l, $(LIBS))
CLIBS = -I $(LIBLOC) -L $(LIBLOC) $(LIBS)


RM = rm -f 

SRC = 1d.cc \
	Basis.cc \
	Element.cc

OBJ = $(subst .cc,.o,$(SRC))

OBJDIR = obj/

OBJ := $(addprefix $(OBJDIR), $(OBJ))

exe = x

all : $(OBJDIR) $(exe) 

$(OBJDIR) :
	mkdir $(OBJDIR) 

$(exe) : $(OBJ) 
	$(CXX) $(CFLAGS) -o $(exe) $(OBJ) $(CLIBS)

$(OBJDIR)%.o : %.cc
	$(CXX) $(CFLAGS) -c $< -o $@ $(CLIBS)

.PHONY : clean 
clean : 
	$(RM) $(OBJDIR) -r 
	$(RM) $(exe) 
