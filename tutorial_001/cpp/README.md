# COMPM080/GV18 Tutorial 1 - C++ 

In this tutorial, we introduce a example framework using libigl to access and visualize point cloud and mesh data.

---
# libigl 
libigl is a geometry processing library developed by https://libigl.github.io/  
It's a header-only library so you do not need to compile this library before you used it in your project.

Libigl heavily using a useful linear algebra library Eigen, which is also a header-only library. (https://eigen.tuxfamily.org/)
Here we list several useful references for helping you familar with Libigl and Eigen's coding style.

*[Eigen quick references](http://eigen.tuxfamily.org/dox/group__QuickRefPage.html)
*[Eigen for matlab users](http://igl.ethz.ch/projects/libigl/matlab-to-eigen.html)
*[Libigl tutorials](https://libigl.github.io/tutorial/)

---

# Install:

## Windows 
download your preferred C++ IDE/Compiler. If you don't have preferrence, we suggest you install VS2017 (https://visualstudio.microsoft.com/zh-hant/vs/)

[Visual Studio 2017](https://visualstudio.microsoft.com/thank-you-downloading-visual-studio/?sku=Community&rel=15)

And install git and cmake 
[git](https://git-scm.com/download/win)
[cmake](https://github.com/Kitware/CMake/releases/download/v3.13.3/cmake-3.13.3-win64-x64.zip)

## OSX

use your package manager to install  
git, cmake and glew

You'll need to install GLEW using
* [Homebrew](http://brew.sh/) 
`brew install git`
`brew install cmake`
`brew install glew`

* [MacPorts](https://www.macports.org/) 
`sudo port install cmake`.
`sudo port install git`.
`sudo port install glew +universal`.

## Linux 
sudo apt-get install git
sudo apt-get install build-essential
sudo apt-get install cmake
sudo apt-get install libx11-dev
sudo apt-get install mesa-common-dev libgl1-mesa-dev libglu1-mesa-dev xorg-dev libglew-dev 

---
