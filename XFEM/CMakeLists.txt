cmake_minimum_required (VERSION 2.6)
project (XFEM_LIB)

######### GET HEADERS #########
set (XFEM_HEADERS
	include/CutModelTopology_CGAL.h  
	include/LinearTetrahedron.h  
	include/loadCutModelTopologyInVectorXFEMLinearTetrahedra.h  
	include/FEM_Elementtype.h  
	include/XFemLinearTetrahedron.h
	include/utils.h
	)

######### GET SOURCE FILES #########
set (XFEM_SOURCES
	src/CutModelTopology_CGAL.cpp 
	src/loadCutModelTopologyInVectorXFEMLinearTetrahedra.cpp
	src/LinearTetrahedron.cpp 
	src/XFemLinearTetrahedron.cpp
	)

######### ADD XFEM LIBRARY #########
add_library (${PROJECT_NAME} ${XFEM_SOURCES} ${XFEM_HEADERS}) #
