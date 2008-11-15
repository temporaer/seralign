bin_PROGRAMS = sdpa

sdpa_SOURCES = sdpa_chordal.cpp sdpa_dataset.cpp sdpa_dpotrf.cpp sdpa_io.cpp \
sdpa_jordan.cpp sdpa_linear.cpp sdpa_main.cpp sdpa_newton.cpp sdpa_parts.cpp \
sdpa_struct.cpp sdpa_tool.cpp sdpa_algebra.h sdpa_chordal.h sdpa_dataset.h \
sdpa_dpotrf.h sdpa_include.h sdpa_io.h sdpa_jordan.h sdpa_linear.h sdpa_newton.h \
sdpa_parts.h sdpa_right.h sdpa_struct.h sdpa_tool.h

if ENABLE_METIS	
sdpa_LDADD = $(pthread_lib) -L$(SPOOLES_LIBDIR) -lspooles -L$(METIS_LIBDIR) -lmetis $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) 
sdpa_CXXFLAGS = $(pthread_cflags) -I$(SPOOLES_INCLUDEDIR) -DUseMETIS=1 -I$(METIS_INCLUDEDIR)
else
sdpa_LDADD = $(pthread_lib) -L$(SPOOLES_LIBDIR) -lspooles $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) 
sdpa_CXXFLAGS = $(pthread_cflags) -I$(SPOOLES_INCLUDEDIR)  -DUseMETIS=0
endif

if BUILD_SPOOLES
SUBDIRS = spooles
endif
