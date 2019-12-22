#----------------------------------------------------------------------
# This section specifies the name of the linker and various
# compilation and linkage options.
#----------------------------------------------------------------------

FC = gfortran
FFLAGS = -Wall -fbacktrace -std=legacy
FLFLAGS =

PROGRAMS = legacy_quad_driver
OBJECTS = 	flgama.o hermit.o hrecur.o hroot.o jacobi.o laguer.o \
			lgrecr.o lgroot.o recur.o root.o

#--------------------------------------------------------------
# This is the default target.
#--------------------------------------------------------------

.PHONY : all
all : $(PROGRAMS)

#-------------------------------------------------------------
# This section shows how to link each example program.
#-------------------------------------------------------------

legacy_quad_driver: legacy_quad_driver.o $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(FLFLAGS) $^

#-------------------------------------------------------------
# This section contains default rule for creating Fortran
# object files.
#-------------------------------------------------------------

%.o: %.f
	$(FC) $(FFLAGS) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

hermit.o: flgama.o hroot.o
hroot.o: hrecur.o
jacobi.o: flgama.o root.o
root.o: recur.o
laguer.o: lgroot.o
lgroot.o: lgrecr.o

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