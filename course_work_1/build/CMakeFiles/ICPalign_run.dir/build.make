# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_SOURCE_DIR = /home/hu/Documents/compM080-compGV18-2019/course_work_1/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/hu/Documents/compM080-compGV18-2019/course_work_1/build

# Include any dependencies generated for this target.
include CMakeFiles/ICPalign_run.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ICPalign_run.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ICPalign_run.dir/flags.make

CMakeFiles/ICPalign_run.dir/main.cpp.o: CMakeFiles/ICPalign_run.dir/flags.make
CMakeFiles/ICPalign_run.dir/main.cpp.o: /home/hu/Documents/compM080-compGV18-2019/course_work_1/src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/hu/Documents/compM080-compGV18-2019/course_work_1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ICPalign_run.dir/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ICPalign_run.dir/main.cpp.o -c /home/hu/Documents/compM080-compGV18-2019/course_work_1/src/main.cpp

CMakeFiles/ICPalign_run.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ICPalign_run.dir/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/hu/Documents/compM080-compGV18-2019/course_work_1/src/main.cpp > CMakeFiles/ICPalign_run.dir/main.cpp.i

CMakeFiles/ICPalign_run.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ICPalign_run.dir/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/hu/Documents/compM080-compGV18-2019/course_work_1/src/main.cpp -o CMakeFiles/ICPalign_run.dir/main.cpp.s

CMakeFiles/ICPalign_run.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/ICPalign_run.dir/main.cpp.o.requires

CMakeFiles/ICPalign_run.dir/main.cpp.o.provides: CMakeFiles/ICPalign_run.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/ICPalign_run.dir/build.make CMakeFiles/ICPalign_run.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/ICPalign_run.dir/main.cpp.o.provides

CMakeFiles/ICPalign_run.dir/main.cpp.o.provides.build: CMakeFiles/ICPalign_run.dir/main.cpp.o


CMakeFiles/ICPalign_run.dir/mytools.cpp.o: CMakeFiles/ICPalign_run.dir/flags.make
CMakeFiles/ICPalign_run.dir/mytools.cpp.o: /home/hu/Documents/compM080-compGV18-2019/course_work_1/src/mytools.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/hu/Documents/compM080-compGV18-2019/course_work_1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/ICPalign_run.dir/mytools.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ICPalign_run.dir/mytools.cpp.o -c /home/hu/Documents/compM080-compGV18-2019/course_work_1/src/mytools.cpp

CMakeFiles/ICPalign_run.dir/mytools.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ICPalign_run.dir/mytools.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/hu/Documents/compM080-compGV18-2019/course_work_1/src/mytools.cpp > CMakeFiles/ICPalign_run.dir/mytools.cpp.i

CMakeFiles/ICPalign_run.dir/mytools.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ICPalign_run.dir/mytools.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/hu/Documents/compM080-compGV18-2019/course_work_1/src/mytools.cpp -o CMakeFiles/ICPalign_run.dir/mytools.cpp.s

CMakeFiles/ICPalign_run.dir/mytools.cpp.o.requires:

.PHONY : CMakeFiles/ICPalign_run.dir/mytools.cpp.o.requires

CMakeFiles/ICPalign_run.dir/mytools.cpp.o.provides: CMakeFiles/ICPalign_run.dir/mytools.cpp.o.requires
	$(MAKE) -f CMakeFiles/ICPalign_run.dir/build.make CMakeFiles/ICPalign_run.dir/mytools.cpp.o.provides.build
.PHONY : CMakeFiles/ICPalign_run.dir/mytools.cpp.o.provides

CMakeFiles/ICPalign_run.dir/mytools.cpp.o.provides.build: CMakeFiles/ICPalign_run.dir/mytools.cpp.o


# Object files for target ICPalign_run
ICPalign_run_OBJECTS = \
"CMakeFiles/ICPalign_run.dir/main.cpp.o" \
"CMakeFiles/ICPalign_run.dir/mytools.cpp.o"

# External object files for target ICPalign_run
ICPalign_run_EXTERNAL_OBJECTS =

ICPalign_run: CMakeFiles/ICPalign_run.dir/main.cpp.o
ICPalign_run: CMakeFiles/ICPalign_run.dir/mytools.cpp.o
ICPalign_run: CMakeFiles/ICPalign_run.dir/build.make
ICPalign_run: /usr/lib/x86_64-linux-gnu/libGL.so
ICPalign_run: imgui/libimgui.a
ICPalign_run: glfw/src/libglfw3.a
ICPalign_run: /usr/lib/x86_64-linux-gnu/librt.so
ICPalign_run: /usr/lib/x86_64-linux-gnu/libm.so
ICPalign_run: /usr/lib/x86_64-linux-gnu/libX11.so
ICPalign_run: glad/libglad.a
ICPalign_run: CMakeFiles/ICPalign_run.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/hu/Documents/compM080-compGV18-2019/course_work_1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable ICPalign_run"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ICPalign_run.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ICPalign_run.dir/build: ICPalign_run

.PHONY : CMakeFiles/ICPalign_run.dir/build

CMakeFiles/ICPalign_run.dir/requires: CMakeFiles/ICPalign_run.dir/main.cpp.o.requires
CMakeFiles/ICPalign_run.dir/requires: CMakeFiles/ICPalign_run.dir/mytools.cpp.o.requires

.PHONY : CMakeFiles/ICPalign_run.dir/requires

CMakeFiles/ICPalign_run.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ICPalign_run.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ICPalign_run.dir/clean

CMakeFiles/ICPalign_run.dir/depend:
	cd /home/hu/Documents/compM080-compGV18-2019/course_work_1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/hu/Documents/compM080-compGV18-2019/course_work_1/src /home/hu/Documents/compM080-compGV18-2019/course_work_1/src /home/hu/Documents/compM080-compGV18-2019/course_work_1/build /home/hu/Documents/compM080-compGV18-2019/course_work_1/build /home/hu/Documents/compM080-compGV18-2019/course_work_1/build/CMakeFiles/ICPalign_run.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ICPalign_run.dir/depend

