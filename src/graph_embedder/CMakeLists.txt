SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wall" )
ENDIF()

include_directories ( 
	. 
	${CMAKE_SOURCE_DIR}/tools
	${CMAKE_SOURCE_DIR}/seriation_gen
	${CMAKE_SOURCE_DIR}/prob_dt
	${CMAKE_SOURCE_DIR}/graph_gen
	)

add_library( graphemb SHARED 
	heatkernel_embedder.hpp   heatkernel_embedder.cpp 
	spectral_embedder.hpp     spectral_embedder.cpp 
	fastmap_embedder.hpp      fastmap_embedder.cpp 
	feature_embedder.hpp      feature_embedder.cpp 
	graph_embedder.hpp        graph_embedder.cpp 
	graph_embedder_cfg.cpp
	)
target_link_libraries(graphemb  ${Boost_LIBRARIES} ${NANA_LIBRARIES} )
