SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wall" )
ENDIF()

include_directories( 
	. 
	${CMAKE_SOURCE_DIR}/configuration
	${CMAKE_SOURCE_DIR}/graph_gen
	${CMAKE_SOURCE_DIR}/seriation_gen
	${CMAKE_SOURCE_DIR}/prob_dt
	${CMAKE_SOURCE_DIR}/sdp_wrappers
	${CMAKE_SOURCE_DIR}/tools
	${CMAKE_SOURCE_DIR}/db
	${CMAKE_SOURCE_DIR}/postproc
	${CMAKE_SOURCE_DIR}/third_party/maart
	${CMAKE_SOURCE_DIR}/third_party/gnuplot-cpp
	)

add_library( actions SHARED action_cfg.cpp 
              action.hpp    action.cpp 
			  serialize.hpp serialize.cpp 
			  build_db.hpp  build_db.cpp
			  fastmap_distancemat.hpp fastmap_distancemat.cpp 
			  )

target_link_libraries( actions maartutils gnuplotcpp )
