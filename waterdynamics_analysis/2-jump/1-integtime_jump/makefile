SOURCES=calc.f90 gser.f gcf.f gammq.f gammln.f fit.f main.f90 

PROJ=./integtime_jump

#COMPILER=ifort
COMPILER=gfortran

ifeq ($(COMPILER),ifort)
 WARNFLAGS= -traceback -check all -W1 -warn notruncated-source -fpe0 -ftrapuv -debug all 
 CODEOPTFLAGS= -O2
 LARGEARRAYFLAGS= -mcmodel=large -i-dynamic
# LARGEARRAYFLAGS= -i-dynamic
# LARGEARRAYFLAGS= -heap-arrays
# CODEOPTFLAGS= -O2 -p
# CODEOPTFLAGS= -g -pg
 PARAFLAGS= -openmp-stubs
else
 WARNFLAGS=  -fbounds-check -Wall -Wextra -Winline -Wno-unused-parameter -Wno-unused-variable
 #WARNFLAGS= -ffree-form -fbounds-check -Wall -Wextra -Winline -Wno-unused-parameter -Wno-unused-variable
 CODEOPTFLAGS= -fdec -O3 -funroll-loops
# CODEOPTFLAGS= -g
 PARAFLAGS= -fopenmp
endif

#FLAGS= $(CODEOPTFLAGS) $(PARAFLAGS)
FLAGS=$(WARNFLAGS) $(CODEOPTFLAGS) 
#FLAGS=$(WARNFLAGS) $(CODEOPTFLAGS) $(PARAFLAGS)
#FLAGS=$(WARNFLAGS) $(CODEOPTFLAGS) $(LARGEARRAYFLAGS)
#FLAGS=$(WARNFLAGS) $(CODEOPTFLAGS) $(LARGEARRAYFLAGS) $(PARAFLAGS)

$(PROJ): $(SOURCES)
	$(COMPILER) $(FLAGS) $(SOURCES) -o $@
	#$(COMPILER) -static $(FLAGS) $(SOURCES) -o $@

clean:
	@rm -rf *.o 

