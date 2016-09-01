#========================================================================
# Makefile to compile the routing model
#========================================================================

# Define compiler
FC90 = pgf90
FC77 = pgf77

#========================================================================
# PART 2: Assemble all of the sub-routines
#========================================================================
EXE = test1.exe

SRC90 = nrtype.f90 \
        public_var.f90 \
        data_type.f90 \
        ascii_util.f90 \
        var_lookup.f90 \
        stringmod.f90 \
        globalData.f90 \
        paramMaster.f90 \
        get_ixname.f90 \
        subset_meta.f90 \
        mo_nml.f90 \
        mo_xor4096.f90 \
        mo_opt_run.f90 \
        mo_dds.f90 \
        vic_routines.f90 \
        vic_subroutines_parallel.f90 \
        eval_obj.f90 \
        main_calibration.f90

OBJS = nrtype.o \
       public_var.o \
       data_type.o \
       ascii_util.o \
       var_lookup.o \
       stringmod.o \
       globalData.o \
       paramMaster.o \
       get_ixname.o \
       subset_meta.o \
       mo_nml.o \
       mo_xor4096.o \
       mo_opt_run.o \
       mo_dds.o \
       vic_routines.o \
       vic_subroutines_parallel.o \
       eval_obj.o \
       main_calibration.o

# Flag
FFLAGS90 = -O3 -Mbackslash -Mfreeform
FFLAGS77 = -O3 -c

#========================================================================
# PART 3: Rule 
#========================================================================
all:  $(EXE) clean

$(EXE): $(OBJS)
	$(FC90) -I./ -o $@ $(OBJS)

$(OBJS): $(SRC90) $(SRC77)

%.o: %.f90
	$(FC90) $(FFLAGS90) -c $@ $<

%.o: %.f
	$(FC77) $(FFLAGS77) -c $@ $<

# Remove object files
clean:
						rm -f *.o
						rm -f *.mod

# Copy the executable to the bin directory
#install:
#		mv $(EXE) $(EXE_PATH)

