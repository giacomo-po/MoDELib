# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

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
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/giacomo/Documents/MODELIB2/tools/DDomp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/giacomo/Documents/MODELIB2/tools/DDomp/build

# Include any dependencies generated for this target.
include CMakeFiles/DDomp.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/DDomp.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/DDomp.dir/flags.make

CMakeFiles/DDomp.dir/main.cpp.o: CMakeFiles/DDomp.dir/flags.make
CMakeFiles/DDomp.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/giacomo/Documents/MODELIB2/tools/DDomp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/DDomp.dir/main.cpp.o"
	/opt/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/DDomp.dir/main.cpp.o -c /Users/giacomo/Documents/MODELIB2/tools/DDomp/main.cpp

CMakeFiles/DDomp.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DDomp.dir/main.cpp.i"
	/opt/local/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/giacomo/Documents/MODELIB2/tools/DDomp/main.cpp > CMakeFiles/DDomp.dir/main.cpp.i

CMakeFiles/DDomp.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DDomp.dir/main.cpp.s"
	/opt/local/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/giacomo/Documents/MODELIB2/tools/DDomp/main.cpp -o CMakeFiles/DDomp.dir/main.cpp.s

# Object files for target DDomp
DDomp_OBJECTS = \
"CMakeFiles/DDomp.dir/main.cpp.o"

# External object files for target DDomp
DDomp_EXTERNAL_OBJECTS =

DDomp: CMakeFiles/DDomp.dir/main.cpp.o
DDomp: CMakeFiles/DDomp.dir/build.make
DDomp: CMakeFiles/DDomp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/giacomo/Documents/MODELIB2/tools/DDomp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable DDomp"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DDomp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/DDomp.dir/build: DDomp

.PHONY : CMakeFiles/DDomp.dir/build

CMakeFiles/DDomp.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/DDomp.dir/cmake_clean.cmake
.PHONY : CMakeFiles/DDomp.dir/clean

CMakeFiles/DDomp.dir/depend:
	cd /Users/giacomo/Documents/MODELIB2/tools/DDomp/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/giacomo/Documents/MODELIB2/tools/DDomp /Users/giacomo/Documents/MODELIB2/tools/DDomp /Users/giacomo/Documents/MODELIB2/tools/DDomp/build /Users/giacomo/Documents/MODELIB2/tools/DDomp/build /Users/giacomo/Documents/MODELIB2/tools/DDomp/build/CMakeFiles/DDomp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/DDomp.dir/depend

