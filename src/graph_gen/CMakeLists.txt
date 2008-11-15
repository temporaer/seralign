SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wall" )
ENDIF()

include_directories( 
	. 
	${CMAKE_SOURCE_DIR}/configuration
	)

SET(SRC 
		adjmat_gen.cpp        adjmat_gen.hpp
        random_adjmat_gen.cpp random_adjmat_gen.hpp
		adjmat_gen_cfg.cpp 
		)

add_library( graph_gen SHARED ${SRC} )
