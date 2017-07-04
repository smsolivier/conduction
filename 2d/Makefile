CXX=g++
CFLAGS = -std=c++14 \
	-fopenmp 

LIBLOC = /home/sam/cpp
LIBS = MyLib
LIBS := $(addprefix -l, $(LIBS))
CLIBS = -I $(LIBLOC) -L $(LIBLOC) $(LIBS)


RM = rm -f 

SRC = 2d.cc \
	Element.cc \
	Basis.cc

OBJ = $(subst .cc,.o,$(SRC))

OBJDIR = obj/

DATADIR = data/

OBJ := $(addprefix $(OBJDIR), $(OBJ))

exe = x

all : $(OBJDIR) $(DATADIR) $(exe) 

$(OBJDIR) :
	mkdir $(OBJDIR) 

$(DATADIR) :
	mkdir $(DATADIR)

$(exe) : $(OBJ) 
	$(CXX) $(CFLAGS) -o $(exe) $(OBJ) $(CLIBS)

$(OBJDIR)%.o : %.cc
	$(CXX) $(CFLAGS) -c $< -o $@ $(CLIBS)

.PHONY : clean 
clean : 
	$(RM) -r $(OBJDIR) 
	$(RM) -r $(DATADIR)
	$(RM) $(exe) 
