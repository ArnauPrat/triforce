#Triforce is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#Triforce is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.2)
PROJECT(triforce CXX)
IF( UNIX )
	SET(CMAKE_CXX_FLAGS_DEBUG "-O0 -std=c++11 -g -pg -fopenmp -DPROFILE -Wall")	
	SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -std=c++11 -fopenmp -DNDEBUG -Wall ")	
ENDIF(UNIX)

IF( WIN32 )
	SET(CMAKE_CXX_FLAGS_DEBUG "/Od /openmp /W3 /Zi")	
	SET(CMAKE_CXX_FLAGS_RELEASE "/O2 /openmp /Oi /Ot /W3 /Zi")	
    SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "/PROFILE")
    SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "/PROFILE")
    IF( MSVC )
        FOREACH( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
            STRING( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
            SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${PROJECT_BINARY_DIR} )
            SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${PROJECT_BINARY_DIR} )
            SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${PROJECT_BINARY_DIR} )
        ENDFOREACH( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )
    ENDIF( MSVC )
ENDIF( WIN32 )

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release)
endif()

INCLUDE_DIRECTORIES(./src)
FILE( GLOB_RECURSE SOURCE_FILES "src/*" )
ADD_EXECUTABLE(triforce ${SOURCE_FILES})	
