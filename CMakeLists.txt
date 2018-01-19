cmake_minimum_required (VERSION 2.6)
project (AleksaMasterThesis)

add_definitions("-std=c++14")
add_definitions("-lboost_system")
add_definitions("-lboost_filesystem")
add_definitions("-DBOOST_NO_CXX11_SCOPED_ENUMS")
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
set(CMAKE_BUILD_TYPE Release)
find_package(Boost 1.45.0 COMPONENTS filesystem system) 

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
endif()

#-------------------------------------------------------------------------------
IF (CMAKE_CXX_COMPILER_ID STREQUAL Clang)
	add_definitions("-Wno-tautological-compare")

ENDIF()

#-------------------------------------------------------------------------------
IF (APPLE)
	set (CMAKE_MODULE_LINKER_FLAGS "-Wl,-flat_namespace -Wl,-undefined -Wl,suppress")
ENDIF (APPLE)

#-------------------------------------------------------------------------------
IF (POLICY CMP0043)
	cmake_policy(SET CMP0043 OLD)
ENDIF()

#-------------------------------------------------------------------------------
# Control CMAKE_BUILD_TYPE, default: Debug
#-------------------------------------------------------------------------------
if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE "Debug" CACHE string "Build type" FORCE)
endif (NOT CMAKE_BUILD_TYPE)
message (STATUS "Build type: " ${CMAKE_BUILD_TYPE})

#-------------------------------------------------------------------------------
# Add library tools (tools)
#-------------------------------------------------------------------------------
aux_source_directory(tools TOOLS_FILES)
add_library (tools STATIC ${TOOLS_FILES})

list (INSERT LIBS 0 tools)

#-------------------------------------------------------------------------------
# Add library tools (models)
#-------------------------------------------------------------------------------

aux_source_directory(models MODELS_FILES)
add_library (models STATIC ${MODELS_FILES})

list (INSERT LIBS 0 models)

#-------------------------------------------------------------------------------
# Add library tools (methods)
#-------------------------------------------------------------------------------

aux_source_directory(methods METHODS_FILES)
add_library (methods STATIC ${METHODS_FILES})

list (INSERT LIBS 0 methods)

#-------------------------------------------------------------------------------
# Add library tools (methods)
#-------------------------------------------------------------------------------

aux_source_directory(mcmc MCMC_FILES)
add_library (mcmc STATIC ${MCMC_FILES})

list (INSERT LIBS 0 mcmc)

#-------------------------------------------------------------------------------
# Add library tools (methods)
#-------------------------------------------------------------------------------

aux_source_directory(smc SMC_FILES)
add_library (smc STATIC ${SMC_FILES})

list (INSERT LIBS 0 smc)

#-------------------------------------------------------------------------------
# Add the executable 'adaptive_timestepping'
#-------------------------------------------------------------------------------
add_executable (main main.cpp)
target_link_libraries (main ${Boost_LIBRARIES}  ${LIBS})

#-------------------------------------------------------------------------------
# Add the executable 'adaptive_timestepping'
#-------------------------------------------------------------------------------
add_executable (adaptive_timestepping adaptive_timestepping.cpp)
target_link_libraries (adaptive_timestepping ${Boost_LIBRARIES}  ${LIBS})

#-------------------------------------------------------------------------------
# Add the executable 'mcmc_main'
#-------------------------------------------------------------------------------
add_executable (mcmc_main mcmc_main.cpp)
target_link_libraries (mcmc_main  ${Boost_LIBRARIES}  ${LIBS})

#-------------------------------------------------------------------------------
# Add the executable 'inverse_finite'
#-------------------------------------------------------------------------------
add_executable (phmc phmc.cpp)
target_link_libraries (phmc ${Boost_LIBRARIES}  ${LIBS})


#-------------------------------------------------------------------------------
# Add the executable 'inverse_finite'
#-------------------------------------------------------------------------------
add_executable (pf_smc_main pf_smc_main.cpp)
target_link_libraries (pf_smc_main ${Boost_LIBRARIES}  ${LIBS})


#-------------------------------------------------------------------------------
# Add the executable 'smc_main'
#-------------------------------------------------------------------------------
add_executable (smc_main smc_main.cpp)
target_link_libraries (smc_main ${Boost_LIBRARIES}  ${LIBS})

#-------------------------------------------------------------------------------
# Add the executable 'pcn'
#-------------------------------------------------------------------------------
add_executable (pcn pcn.cpp)
target_link_libraries (pcn ${Boost_LIBRARIES}  ${LIBS})

