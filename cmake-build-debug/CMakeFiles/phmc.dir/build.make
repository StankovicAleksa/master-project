# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/aleksa/opt/clion/bin/cmake/bin/cmake

# The command to remove a file.
RM = /home/aleksa/opt/clion/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/aleksa/Dropbox/master/code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/aleksa/Dropbox/master/code/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/phmc.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/phmc.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/phmc.dir/flags.make

CMakeFiles/phmc.dir/phmc.cpp.o: CMakeFiles/phmc.dir/flags.make
CMakeFiles/phmc.dir/phmc.cpp.o: ../phmc.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/aleksa/Dropbox/master/code/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/phmc.dir/phmc.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/phmc.dir/phmc.cpp.o -c /home/aleksa/Dropbox/master/code/phmc.cpp

CMakeFiles/phmc.dir/phmc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/phmc.dir/phmc.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/aleksa/Dropbox/master/code/phmc.cpp > CMakeFiles/phmc.dir/phmc.cpp.i

CMakeFiles/phmc.dir/phmc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/phmc.dir/phmc.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/aleksa/Dropbox/master/code/phmc.cpp -o CMakeFiles/phmc.dir/phmc.cpp.s

CMakeFiles/phmc.dir/phmc.cpp.o.requires:

.PHONY : CMakeFiles/phmc.dir/phmc.cpp.o.requires

CMakeFiles/phmc.dir/phmc.cpp.o.provides: CMakeFiles/phmc.dir/phmc.cpp.o.requires
	$(MAKE) -f CMakeFiles/phmc.dir/build.make CMakeFiles/phmc.dir/phmc.cpp.o.provides.build
.PHONY : CMakeFiles/phmc.dir/phmc.cpp.o.provides

CMakeFiles/phmc.dir/phmc.cpp.o.provides.build: CMakeFiles/phmc.dir/phmc.cpp.o


# Object files for target phmc
phmc_OBJECTS = \
"CMakeFiles/phmc.dir/phmc.cpp.o"

# External object files for target phmc
phmc_EXTERNAL_OBJECTS =

phmc: CMakeFiles/phmc.dir/phmc.cpp.o
phmc: CMakeFiles/phmc.dir/build.make
phmc: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so
phmc: /usr/lib/x86_64-linux-gnu/libboost_system.so
phmc: libsmc.a
phmc: libmcmc.a
phmc: libmethods.a
phmc: libmodels.a
phmc: libtools.a
phmc: CMakeFiles/phmc.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/aleksa/Dropbox/master/code/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable phmc"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/phmc.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/phmc.dir/build: phmc

.PHONY : CMakeFiles/phmc.dir/build

CMakeFiles/phmc.dir/requires: CMakeFiles/phmc.dir/phmc.cpp.o.requires

.PHONY : CMakeFiles/phmc.dir/requires

CMakeFiles/phmc.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/phmc.dir/cmake_clean.cmake
.PHONY : CMakeFiles/phmc.dir/clean

CMakeFiles/phmc.dir/depend:
	cd /home/aleksa/Dropbox/master/code/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/aleksa/Dropbox/master/code /home/aleksa/Dropbox/master/code /home/aleksa/Dropbox/master/code/cmake-build-debug /home/aleksa/Dropbox/master/code/cmake-build-debug /home/aleksa/Dropbox/master/code/cmake-build-debug/CMakeFiles/phmc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/phmc.dir/depend
