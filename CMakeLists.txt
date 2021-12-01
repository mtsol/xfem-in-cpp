cmake_minimum_required (VERSION 2.6)
project (XFEMExec)



################## GET LIBRARIES ##################

######### CMAKE #########
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
 
######### CGAL #########
find_package(CGAL REQUIRED)
include_directories(${CGAL_INCLUDE_DIRS})
set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${CGAL_LIBRARY})
# to avoid assertion fails
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -frounding-math)

######### VTK #########
find_package(VTK REQUIRED)
include_directories(${VTK_INCLUDE_DIRS})
set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${VTK_LIBRARIES})

######### EIGEN #########
find_package(EIGEN REQUIRED)
include_directories(${EIGEN_INCLUDE_DIRS})
#message("TEST " ${EIGEN_INCLUDE_DIRS})
set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${EIGEN_LIBRARIES})

######### BOOST #########
set(EXTERNAL_LIBS ${EXTERNAL_LIBS} -lboost_system )

######### TINYXML #########
set(EXTERNAL_LIBS ${EXTERNAL_LIBS} -ltinyxml )

######### XFEM #########
include_directories ("XFEM/include")
add_subdirectory (XFEM)
set(XFEM_LIBS XFEM_LIB)

##################



################## ADD EXECUTABLE AND LINK LIBRARIES ##################

add_executable (XFEMExec Main.cxx)

target_link_libraries (XFEMExec  ${XFEM_LIBS})
target_link_libraries (XFEMExec  ${EXTERNAL_LIBS})

##################