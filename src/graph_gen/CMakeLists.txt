SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wall" )
ENDIF()

# configure triangle and the perl-wrapper
FIND_PROGRAM(TRIANGLE_BINARY triangle PATHS ${CMAKE_SOURCE_DIR}/third_party/triangle)
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/util/randpat/triangulate.pl ${SERALIGN_SCRIPT_DIR}/triangulate.pl @ONLY )
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/util/randpat/TriangleWrap.pm ${SERALIGN_SCRIPT_DIR}/TriangleWrap.pm @ONLY )
ADD_DEFINITIONS(-DTRIANGULATION_COMMAND="${SERALIGN_SCRIPT_DIR}/triangulate.pl")

include_directories( 
	. 
	${CMAKE_SOURCE_DIR}/configuration
	${CMAKE_SOURCE_DIR}/prob_dt
	${CMAKE_SOURCE_DIR}/tools
	)

SET(SRC 
		adjmat_gen.cpp        adjmat_gen.hpp
        random_adjmat_gen.cpp random_adjmat_gen.hpp
        mutagenesis.cpp       mutagenesis.hpp
		adjmat_gen_cfg.cpp 
		DegreeSort.cpp        DegreeSort.hpp
		random_pat_adjmat_gen.cpp random_pat_adjmat_gen.hpp
		jumbled_adjmat_gen.cpp jumbled_adjmat_gen.hpp
		)

add_library( graph_gen SHARED ${SRC} )
