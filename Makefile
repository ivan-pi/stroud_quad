#----------------------------------------------------------------------
# This section specifies the name of the linker and various
# compilation and linkage options.
#----------------------------------------------------------------------

FC = gfortran
FFLAGS = -Wall
FLFLAGS =

PROGRAMS = quad_driver
OBJECTS = quad_mod.o

#--------------------------------------------------------------
# This is the default target.
#--------------------------------------------------------------

.PHONY : all
all : $(PROGRAMS)

#-------------------------------------------------------------
# This section shows how to link each example program.
#-------------------------------------------------------------

quad_driver: quad_driver.o quad_mod.o
	$(FC) $(FFLAGS) -o $@ $(FLFLAGS) $^

#-------------------------------------------------------------
# This section contains default rule for creating Fortran
# object files.
#-------------------------------------------------------------

%.o: %.f
	$(FC) $(FFLAGS) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

#----------------------------------------------
# This section shows how to clean up afterward.
#----------------------------------------------

.PHONY : clean cleanobj cleanmod
clean : cleanobj cleanmod
	rm -f $(PROGRAMS)
cleanobj :
	rm -f *.o
cleanmod :
	rm -f *.mod *__genmod.f90