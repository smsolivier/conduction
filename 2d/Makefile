CXX=g++
CFLAGS = -std=c++14 \
	-fopenmp -g

LIBLOC = ./cpp_library
LIBS = MyLib
LIBS := $(addprefix -l, $(LIBS))
CLIBS = -I $(LIBLOC) -L $(LIBLOC) $(LIBS)


RM = rm -f 

SRC = $(wildcard src/*.cc)

OBJ = $(subst .cc,.o,$(SRC))

OBJDIR = obj/

DATADIR = data/

OBJ := $(addprefix $(OBJDIR), $(notdir $(OBJ)))

exe = x

all : $(OBJDIR) $(DATADIR) $(exe)

$(OBJDIR) :
	mkdir $(OBJDIR) 

$(DATADIR) :
	mkdir $(DATADIR)

$(exe) : $(OBJ)
	$(CXX) $(CFLAGS) -o $(exe) $(OBJ) $(CLIBS)

# $(OBJDIR)%.o : %.cc
# 	$(CXX) $(CFLAGS) -c $< -o $@ $(CLIBS)

$(OBJDIR)%.o:src/%.cc Makefile
	$(CXX) $(CFLAGS) -c $< -o $@ $(CLIBS)
	$(CXX) -MM $(CFLAGS) $< $(CLIBS) > src/$*.d

.PHONY : clean 
clean : 
	$(RM) -r $(OBJDIR) 
	$(RM) -r $(DATADIR)
	$(RM) $(exe) 

-include $(OBJ:.o=.d)