# COMPM080/GV18 Tutorial 1 - C++ 

In this tutorial, we introduce an example framework using libigl to access and visualize point cloud and mesh data.  
This tutorial contains two practices.  
1. compile an example viewer.  
2. use this viewer to show point cloud data and perform simple computations.  

---
# libigl 
libigl is a geometry processing library developed by https://libigl.github.io/  
It's a header-only library so you do not need to compile this library before you used it in your project.

Libigl heavily uses Eigen (https://eigen.tuxfamily.org/) to handle mesh data. Good news is it's also a header-only library. Here we list several useful references to help you get familiar with Libigl and Eigen's coding style.

* [Eigen quick references](http://eigen.tuxfamily.org/dox/group__QuickRefPage.html)
* [Eigen for matlab users](http://igl.ethz.ch/projects/libigl/matlab-to-eigen.html)
* [Libigl tutorials](https://libigl.github.io/tutorial/)
* [C++ template introduction](http://www.cplusplus.com/doc/oldtutorial/templates/)

---
# FAQ

* in pratice-2, I saw an error `[error] scan data not found`.
  please change input data's path at [main.cpp#L318](https://github.com/smartgeometry-ucl/compM080-compGV18-2019/blob/master/tutorial_001/cpp/pratices_002/src/main.cpp#L318) to absolute path such as `"/mydrive/compM080-compGV18-2019/tutorial_001/cpp/data/static_dock.sens"`.  
* cannot excute my built program in VS2017. You need to change start project to project `pratices_2_bin`. Please check *WINDOWS IDE setting:*. 


---

# Install:

## Windows 
download your preferred C++ IDE/Compiler. If you don't have a preference, we suggest you install VS2017.

[Visual Studio 2017](https://visualstudio.microsoft.com/thank-you-downloading-visual-studio/?sku=Community&rel=15)

And install git and cmake 
[git](https://git-scm.com/download/win)
[git-extension, gui for git](https://github.com/gitextensions/gitextensions/releases/download/v3.00.00/GitExtensions-3.00.00.4433.msi)
[cmake](https://github.com/Kitware/CMake/releases/download/v3.13.3/cmake-3.13.3-win64-x64.zip)

## OSX
use your package manager to install git, cmake, glew, and g++

Binary Install files:  
[cmake](https://github.com/Kitware/CMake/releases/download/v3.13.3/cmake-3.13.3-Darwin-x86_64.dmg)

Reference cmd:  
* [Homebrew](http://brew.sh/)  
`brew install git`  
`brew install cmake`  
`brew install glew`  
`brew install gcc`

* [MacPorts](https://www.macports.org/)  
`sudo port install cmake`  
`sudo port install git`  
`sudo port install glew`  
`sudo port install gcc`  

## Linux 
use your package manager to install git, cmake, opengl and glew  

Reference cmd:  
`sudo apt-get install git`  
`sudo apt-get install build-essential`  
`sudo apt-get install cmake`  
`sudo apt-get install cmake-qt-gui` or `sudo apt-get install cmake-gui`  
`sudo apt-get install libx11-dev`  
`sudo apt-get install mesa-common-dev libgl1-mesa-dev libglu1-mesa-dev xorg-dev libglew-dev`  

---
# Practice 1
In this practice, you will learn how to compile our example viewer via cmake.
After successfully compile the files, please check our code to see how to create viewer.

![p1](/tutorial_001/cpp/docimgs/p1.JPG "")

## i. get libigl via git submodule or git clone  
option 1) use submodule update:  
1. git clone this repository : `git clone https://github.com/smartgeometry-ucl/compM080-compGV18-2019.git`
2. `cd compM080-compGV18-2019`
3. `git submodule update --init --recursive`
4. check libigl directory in `/tutorial_001/cpp/libigl` is not empty  

option 2) download/clone libigl at `/tutorial_001/cpp/libigl`  

## ii. use cmake to build file.  

### option 1) use terminal (OSX/LINUX)
1. open cmd and switch to `/tutorial_001/cpp/pratices_001`  
2. build:  
`mkdir build; cd build`  
`cmake -DCMAKE_BUILD_TYPE=Release  ../src`  
`make`
3. run execution file via `./pratices_1_bin`

* For OSX users, if your building process is failed, please try to specify your compiler:  
`cmake -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ ../src`

example output:
````$cmake ..
-- The C compiler identification is AppleClang 7.3.0.7030031
-- The CXX compiler identification is AppleClang 7.3.0.7030031
-- Check for working C compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc
-- Check for working C compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc -- works
-- Detecting C compiler ABI info
.....
HEAD is now at a37e6e5... Update CMake.
-- Configuring done
-- Generating done
-- Build files have been written to: /Users/tirth/MainDataBox/1_Workspace/develop.project/2019.geometry_course/compM080-compGV18-2019/tutorial_001/cpp/pratices_001/src/build
$make
Scanning dependencies of target pratices_1_bin
...
[ 96%] Building CXX object CMakeFiles/pratices_1_bin.dir/main.cpp.o
[100%] Linking CXX executable pratices_1_bin
[100%] Built target pratices_1_bin
````

### option 2) use camke-gui (OSX/LINUX/WINDOWS)

For LINUX run `cmake-gui`  
For WINDOWS/OSX open cmake application  

#### References steps:
1. press configure
2. press generate
3. press Open

![configure](/tutorial_001/cpp/docimgs/cmake.JPG "")
![generate](/tutorial_001/cpp/docimgs/cmake_gen.JPG "")

#### OSX:
the default compiler is MAKE.  
After you did cmake you need to open a terminal to excute `make` in the build directory.  


#### WINDOWS IDE setting:  
**You need to set the default start project to execute program in VS2017**

![p2](/tutorial_001/cpp/docimgs/vs15.jpg "")

---
# Practice 2
In this practice, we showcase a scan data captured via using Kinect. We follow ScanNet's (https://github.com/ScanNet/ScanNet) approach using `.sens` file format to store scanned data. All IO functions are provided in this practice.  **Your job** is to implement a simple function `average_depth` in `tutorial_001\cpp\pratices_002\src\mytools.cpp`. This function will take a bunch of depth frames and compute an average depth image.  

**[What is depth image](https://www.quora.com/What-is-depth-image):**  
````Now depth image has values according how far is object. Pixel represent distance from camera.````


Note that depth images usually contains some invalid values (depth equal to zero) due to the artifact generated in the scanning process. You will need to bypass these invalid pixels when you are doing average. 

**comparsion**:  
![comparsion](/tutorial_001/cpp/docimgs/comp.jpg "")

**average depth example**:  
![average depth example](/tutorial_001/cpp/docimgs/p2.JPG "")




