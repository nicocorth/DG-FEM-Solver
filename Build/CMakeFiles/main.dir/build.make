# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/nicolas/Documents/Current_Work/TFE/DG_Solver

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/nicolas/Documents/Current_Work/TFE/DG_Solver/Build

# Include any dependencies generated for this target.
include CMakeFiles/main.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.dir/flags.make

CMakeFiles/main.dir/main.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nicolas/Documents/Current_Work/TFE/DG_Solver/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/main.cpp.o -c /home/nicolas/Documents/Current_Work/TFE/DG_Solver/main.cpp

CMakeFiles/main.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nicolas/Documents/Current_Work/TFE/DG_Solver/main.cpp > CMakeFiles/main.dir/main.cpp.i

CMakeFiles/main.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nicolas/Documents/Current_Work/TFE/DG_Solver/main.cpp -o CMakeFiles/main.dir/main.cpp.s

CMakeFiles/main.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/main.dir/main.cpp.o.requires

CMakeFiles/main.dir/main.cpp.o.provides: CMakeFiles/main.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/main.dir/main.cpp.o.provides

CMakeFiles/main.dir/main.cpp.o.provides.build: CMakeFiles/main.dir/main.cpp.o


CMakeFiles/main.dir/Mesh/Element.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/Mesh/Element.cpp.o: ../Mesh/Element.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nicolas/Documents/Current_Work/TFE/DG_Solver/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/main.dir/Mesh/Element.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/Mesh/Element.cpp.o -c /home/nicolas/Documents/Current_Work/TFE/DG_Solver/Mesh/Element.cpp

CMakeFiles/main.dir/Mesh/Element.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/Mesh/Element.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nicolas/Documents/Current_Work/TFE/DG_Solver/Mesh/Element.cpp > CMakeFiles/main.dir/Mesh/Element.cpp.i

CMakeFiles/main.dir/Mesh/Element.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/Mesh/Element.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nicolas/Documents/Current_Work/TFE/DG_Solver/Mesh/Element.cpp -o CMakeFiles/main.dir/Mesh/Element.cpp.s

CMakeFiles/main.dir/Mesh/Element.cpp.o.requires:

.PHONY : CMakeFiles/main.dir/Mesh/Element.cpp.o.requires

CMakeFiles/main.dir/Mesh/Element.cpp.o.provides: CMakeFiles/main.dir/Mesh/Element.cpp.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/Mesh/Element.cpp.o.provides.build
.PHONY : CMakeFiles/main.dir/Mesh/Element.cpp.o.provides

CMakeFiles/main.dir/Mesh/Element.cpp.o.provides.build: CMakeFiles/main.dir/Mesh/Element.cpp.o


CMakeFiles/main.dir/Solver/unknowns.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/Solver/unknowns.cpp.o: ../Solver/unknowns.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nicolas/Documents/Current_Work/TFE/DG_Solver/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/main.dir/Solver/unknowns.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.dir/Solver/unknowns.cpp.o -c /home/nicolas/Documents/Current_Work/TFE/DG_Solver/Solver/unknowns.cpp

CMakeFiles/main.dir/Solver/unknowns.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/Solver/unknowns.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nicolas/Documents/Current_Work/TFE/DG_Solver/Solver/unknowns.cpp > CMakeFiles/main.dir/Solver/unknowns.cpp.i

CMakeFiles/main.dir/Solver/unknowns.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/Solver/unknowns.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nicolas/Documents/Current_Work/TFE/DG_Solver/Solver/unknowns.cpp -o CMakeFiles/main.dir/Solver/unknowns.cpp.s

CMakeFiles/main.dir/Solver/unknowns.cpp.o.requires:

.PHONY : CMakeFiles/main.dir/Solver/unknowns.cpp.o.requires

CMakeFiles/main.dir/Solver/unknowns.cpp.o.provides: CMakeFiles/main.dir/Solver/unknowns.cpp.o.requires
	$(MAKE) -f CMakeFiles/main.dir/build.make CMakeFiles/main.dir/Solver/unknowns.cpp.o.provides.build
.PHONY : CMakeFiles/main.dir/Solver/unknowns.cpp.o.provides

CMakeFiles/main.dir/Solver/unknowns.cpp.o.provides.build: CMakeFiles/main.dir/Solver/unknowns.cpp.o


# Object files for target main
main_OBJECTS = \
"CMakeFiles/main.dir/main.cpp.o" \
"CMakeFiles/main.dir/Mesh/Element.cpp.o" \
"CMakeFiles/main.dir/Solver/unknowns.cpp.o"

# External object files for target main
main_EXTERNAL_OBJECTS =

bin/main: CMakeFiles/main.dir/main.cpp.o
bin/main: CMakeFiles/main.dir/Mesh/Element.cpp.o
bin/main: CMakeFiles/main.dir/Solver/unknowns.cpp.o
bin/main: CMakeFiles/main.dir/build.make
bin/main: ../gmsh-4.4.1-Linux64-sdk/lib/libgmsh.so.4.4.1
bin/main: CMakeFiles/main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/nicolas/Documents/Current_Work/TFE/DG_Solver/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable bin/main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.dir/build: bin/main

.PHONY : CMakeFiles/main.dir/build

CMakeFiles/main.dir/requires: CMakeFiles/main.dir/main.cpp.o.requires
CMakeFiles/main.dir/requires: CMakeFiles/main.dir/Mesh/Element.cpp.o.requires
CMakeFiles/main.dir/requires: CMakeFiles/main.dir/Solver/unknowns.cpp.o.requires

.PHONY : CMakeFiles/main.dir/requires

CMakeFiles/main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.dir/clean

CMakeFiles/main.dir/depend:
	cd /home/nicolas/Documents/Current_Work/TFE/DG_Solver/Build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/nicolas/Documents/Current_Work/TFE/DG_Solver /home/nicolas/Documents/Current_Work/TFE/DG_Solver /home/nicolas/Documents/Current_Work/TFE/DG_Solver/Build /home/nicolas/Documents/Current_Work/TFE/DG_Solver/Build /home/nicolas/Documents/Current_Work/TFE/DG_Solver/Build/CMakeFiles/main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main.dir/depend

