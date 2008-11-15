CC = gcc

CFLAGS = -O3

LIB_DIRS = 

# This is the default on RedHat (and possibly other) Linux systems (METIS and ARPACK not available)
LIBS = -lgsl -llapack -lblas -lg2c -lm

# This is how it might look with ATLAS installed
# LIBS = -lgsl -llapack -lcblas -lf77blas -latlas -lm

DEFINES = 

MEX     = mex             # Matlab Mex executable
MEXEXT  = mexglx          # Matlab Mex extension (mexglx = Linux)
LIBSMEX = $(LIBS)         # Mex-specific libraries (if not LIBS from above)
