CXX=g++
CFLAGS = -std=c++14 -fopenmp 

# source file locations 
SRC = src

# include from source locations 
CFLAGS += -I$(SRC)

# cpp library 
LIBLOC = ./cpp_library
LIBS = MyLib
LIBS := $(addprefix -l, $(LIBS))
CLIBS = -I $(LIBLOC) -L $(LIBLOC) $(LIBS)

# look for .cc files in these folders 
VPATH = src

# all .cc files 
SRCFILES = $(notdir $(wildcard $(SRC)/*.cc $(SPMAT)/*.cc))

# where to put object and dependency files 
OBJDIR = obj
DEPDIR = dep

# convert .cc -> objdir/.o 
OBJ = $(patsubst %.cc,$(OBJDIR)/%.o,$(SRCFILES))

# output directory 
DATADIR = data

# convert .o -> depdir/.d 
DEPS = $(patsubst $(OBJDIR)/%.o, $(DEPDIR)/%.d, $(OBJ))

# name of output executable 
exe = x

all : $(exe)

# create executable 
# depends on cpp_library and object files 
$(exe) : $(LIBLOC)/libMyLib.a $(OBJ)  
	mkdir -p $(DATADIR)
	$(CXX) $(CFLAGS) $(OBJ) $(CLIBS) -o $(exe)

# build all object and dependency files (uses VPATH to find .cc files)
$(OBJDIR)/%.o : %.cc Makefile
	mkdir -p $(OBJDIR)
	$(CXX) $(CFLAGS) -c $< -o $@ $(CLIBS)
	mkdir -p $(DEPDIR)
	$(CXX) -MM $(CFLAGS) $< $(CLIBS) | sed '1s/^/$(OBJDIR)\//' > $*.d; mv $*.d $(DEPDIR)

# run the make file for the library 
$(LIBLOC)/libMyLib.a : $(wildcard $(LIBLOC)/*.cc $(LIBLOC)/*.hh) 
	cd $(LIBLOC); make clean; make CXX=$(CXX)

# clean up 
.PHONY : clean 
clean : 
	rm -f -r $(OBJDIR) 
	rm -f -r $(DATADIR)
	rm -f $(exe)
	rm -f -r $(OBJDIR) $(DEPDIR)
	cd $(LIBLOC); make clean 

listsrc :
	@echo $(SRCFILES)
	
# add in additional dependcy information 
-include $(DEPS)