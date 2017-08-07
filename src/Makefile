#========================================================================
# Makefile to compile MPR 
#========================================================================
#
# To troubleshoot your paths and setup, type 'make check'
# 
# You will need to check F_MASTER and FC. You will probably
# also have to specify the location of your netcdf and lapack libraries
# 
#========================================================================
# PART 0: User-configuration  
#========================================================================
# Define core directory.
# This is the parent directory of 'build' where all the source codes reside
# 'mpr-flex" under home directory is default.
F_MASTER = ~/mpr-flex/

# Define the Fortran Compiler. 
# ifort is default. The other options are gfortran and pgf90
FC  = ifort

# Define the NetCDF libraries and path to include files
# You may need to change the netCDF paths
ifeq "$(FC)" "gfortran"
  NCDF_PATH=/opt/netcdf4-gcc
  LIBNETCDF = -Wl,-rpath,$(NCDF_PATH)/lib \
    		      -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
  INCNETCDF = -I$(NCDF_PATH)/include
endif
ifeq "$(FC)" "ifort"
  NCDF_PATH = /opt/netcdf4-intel
  LIBNETCDF = -Wl,-rpath,$(NCDF_PATH)/lib \
    		      -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
  INCNETCDF = -I$(NCDF_PATH)/include
endif
ifeq "$(FC)" "pgf90"
  NCDF_PATH=/opt/netcdf4-pgi
  LIBNETCDF = -Wl,-rpath,$(NCDF_PATH)/lib \
              -L$(NCDF_PATH)/lib -lnetcdff -lnetcdf
  INCNETCDF = -I$(NCDF_PATH)/include
endif

#========================================================================
# PART 1: Define directory paths
#========================================================================
# Core directory that contains source code
F_KORE_DIR = $(F_MASTER)src/

# Location of the compiled modules
MOD_PATH = $(F_MASTER)src/

# Define the directory for the executables
EXE_PATH = $(F_MASTER)bin/

#========================================================================
# PART 1: Assemble all of the sub-routines
#========================================================================
# Define subroutines
DATATYPES = \
		nrtype.f90 \
    public_var.f90 \
		data_type.f90 \
		var_lookup.f90 \

UTILS = \
    ascii_util.f90 \
    stringmod.f90 \

IO = \
    globalData.f90 \
    popMeta.f90 \
    read_config.f90 \
    read_soildata.f90 \
    read_topodata.f90 \
    read_vegdata.f90 \
    read_mapdata.f90 \
    read_ncdata.f90 \
    write_param_nc.f90 \

CORE = \
    get_ixname.f90 \
    process_meta.f90 \
    mo_xor4096.f90 \
    mo_dds.f90 \
    sce.f90 \
    sac_routines.f90 \
    vic_routines.f90 \
    model_wrapper.f90 \
    tf.f90 \
    modelLayer.f90 \
    upscaling.f90 \
    mpr_routine.f90 \
    eval_model.f90 \
    mo_opt_run.f90 \

# concatanate model subroutines
TEMP_MODSUB = $(DATATYPES) $(UTILS) $(IO) $(CORE)

# insert appropriate directory name
MODSUB = $(patsubst %, $(F_KORE_DIR)%, $(TEMP_MODSUB))

#========================================================================
# PART 2: Define  driver programs, and executables
#========================================================================
# Define the driver routine
DRIVER = $(F_KORE_DIR)main_calibration.f90

# Define the executable
EXE = mpr-flex.exe

#========================================================================
# PART 3: Compile the puppy
#========================================================================
# Define flags
ifeq "$(FC)" "gfortran"
 FLAGS = -p -Wall -g -ffree-line-length-none -fmax-errors=0 -fbacktrace -fcheck=all
endif
ifeq "$(FC)" "ifort"
 #FLAGS1 = -warn all -noerror_limit -FR -O0 -auto -WB -traceback -g -fltconsistency
 FLAGS = -heap-arrays 10 -debug all -check all -warn all -FR -O2 -auto -WB -traceback -g -fltconsistency -fpe0
endif
ifeq "$(FC)" "pgf90"
 #FLAGS1 = -Bstatic -Mbackslash -g -Mchkptr -Mchkstk -Mpgicoff -Minform=inform -Ktrap=divz,inv -Mprof=lines,time
 FLAGS = -Mbackslash -Mnobounds -g -Mchkptr -Mchkstk -Mpgicoff -Minform=inform -Ktrap=divz,inv
endif

# tasks
all: compile clean install

# compile
compile:
	$(FC) $(FLAGS) $(MODSUB) $(DRIVER) \
	$(LIBNETCDF) $(INCNETCDF) -o $(EXE)

# Remove object files
clean:
						rm -f *.o
						rm -f *.mod
#
# Copy the executable to the bin directory
install:
		mv $(EXE) $(EXE_PATH)

check:
  $(info)
  $(info Displaying make variables:)
  $(info   F_MASTER   : $(F_MASTER))
  $(info   EXE_PATH   : $(EXE_PATH))
  $(info   FC         : $(FC))
  $(info   FLAGS      : $(FLAGS))
  $(info   NCDF_PATH  : $(NCDF_PATH))
  $(info   LIBNETCDF  : $(LIBNETCDF))
  $(info   INCNETCDF  : $(INCNETCDF))
  $(info Add $(NCDF_PATH)/lib to your LD_LIBRARY_PATH)
  $(info)
