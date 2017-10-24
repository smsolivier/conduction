CXX=g++
CFLAGS = -std=c++14 -fopenmp 

# source file locations 
SRC = src
SPMAT = spmat

# include from source locations 
CFLAGS += -I$(SRC) -I$(SPMAT)

# cpp library 
LIBLOC = ./cpp_library
LIBS = MyLib
LIBS := $(addprefix -l, $(LIBS))
CLIBS = -I $(LIBLOC) -L $(LIBLOC) $(LIBS)

# look for .cc files in these folders 
VPATH = src spmat

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

# create executable 
all : $(exe)

$(exe) : $(OBJ)
	mkdir -p $(DATADIR)
	$(CXX) $(CFLAGS) $(OBJ) $(CLIBS) -o $(exe)

# build all object and dependency files (uses VPATH to find .cc files)
$(OBJDIR)/%.o : %.cc Makefile
	mkdir -p $(OBJDIR)
	$(CXX) $(CFLAGS) -c $< -o $@ $(CLIBS)
	mkdir -p $(DEPDIR)
	$(CXX) -MM $(CFLAGS) $< $(CLIBS) | sed '1s/^/$(OBJDIR)\//' > $*.d; mv $*.d $(DEPDIR)

.PHONY : clean 
clean : 
	rm -f -r $(OBJDIR) 
	rm -f -r $(DATADIR)
	rm -f $(exe)
	rm -f -r $(OBJDIR) $(DEPDIR)

listsrc :
	@echo $(SRCFILES)
	
# add in additional dependcy information 
-include $(DEPS)