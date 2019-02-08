# COMPM080/GV18 Tutorial 3 - C++ 

In this tutorial, we introduce how to recover a Signed Distance Field (SDF) by using RBF fitting method.  
[Reference paper](http://mesh.brown.edu/DGP/pdfs/Carr-sg2001.pdf)  
[Slide on moodle/topic4](https://moodle-1819.ucl.ac.uk/course/view.php)

---
# FAQ 

* There is an run-time error in windows env. : make sure you select "Visual Studio 201X - **Win64**" in CMake's Configuration (most of computers are 64bits now)
* CMake shows it cannot find libigl : library is missing. please read *STEP-2. get libigl* section at https://github.com/smartgeometry-ucl/compM080-compGV18-2019/tree/master/tutorial_002/cpp#step-2-get-libigl 

---
# C++ Dev. Cycle
1. do cmake to generate building file (configure->generate->Open IDE or `make`)
2. build your source code (via IDE's build or `make`)
3. run your program via `./practice_xxx.bin` or press Run in your IDE
4. modify source code in `src`
5. (optional) change to Debug build or RelWithDebug build to debug:  
[for WINDOWS] change your build type to RelWithDebug in VisualStudio  
[ref for OSX/LINUX](https://stackoverflow.com/questions/7724569/debug-vs-release-in-cmake)  
6. go back to **2 build your source code**

---
# Useful references 
* [Eigen quick references](http://eigen.tuxfamily.org/dox/group__QuickRefPage.html)
* [Eigen for matlab users](http://igl.ethz.ch/projects/libigl/matlab-to-eigen.html)
* [Libigl tutorials](https://libigl.github.io/tutorial/)
* [C++ template introduction](http://www.cplusplus.com/doc/oldtutorial/templates/)


---

# Install 
## STEP-1. get the code of tutorial 3:

### option (a) For those who are familiar with Git. 
1. cd to your previous `compM080-compGV18-2019` directory  
`cd ~/my_code/compM080-compGV18-2019`
2. commit all your change  
`git add .`  
`git commit -m 'added all my changed`  
3. get new code 
`git pull origin master`

### option (b) just download/clone this repository again.  
`git clone https://github.com/smartgeometry-ucl/compM080-compGV18-2019.git`

## STEP-2. get libigl (optinal)
If you can build tutorial-2 than you donot need to reset libigl.  
location: `compM080-compGV18-2019/libigl`   

---
# Practice - SDF
In this practice, you will learn how to recover a Signed Distance Field (SDF) by using RBF fitting method.  

## reference images
### 
![](/tutorial_003/cpp/doc_img/res1.JPG ) 

**Your job** is to implement the function `fit_SDF` and `calculate_SDF` at `tutorial_003\cpp\practices_sdf\src\mytools.cpp` 

**fit_SDF**  
````
void fit_SDF(Eigen::MatrixXd const & MODEL_PTS, Eigen::MatrixXd const & SIGN_DIST, Eigen::VectorXd & out_lambda_c)
{
	//
	// fit a RBF function on input data points
	//
	// MODEL_PTS: input data point, (N,3)
	// SIGN_DIST: SDF(MODEL_PTS): R^3 -> R^1 - the SDF value of each input vertex, (N,1)
	// out_lambda_c: the parameters of the recovered RBF function, (N+4,1)
	//
	// Hints: check few rows of the value you set in your matrix and make sure they are correct
	//
	const size_t N = MODEL_PTS.rows();

	Eigen::MatrixXf A(N, N);
	Eigen::MatrixXf P(N, 4);
	Eigen::MatrixXf b(N + 4, 1);

	// MX=b 
	Eigen::MatrixXf M(N + 4, N + 4);
  
  // ......
  
````

**calculate_SDF**  
````
void calculate_SDF(Eigen::MatrixXd const & MODEL_PTS,
	Eigen::VectorXd const & LAMBDA_C,
	Eigen::MatrixXd const & SAMPLE_PTS,
	Eigen::VectorXd & out_sign_dist)
{
	// given our RBF function and some input points
	// estimate the SDF value of each input points
	//
	// SAMPLE_PTS: input point, (N,3)
	// out_sign_dist: SDF value of each input point, (N,3)
	//

````
