# COMPM080/GV18 Tutorial 4 - C++ 

In this tutorial, we introduce how to retrieve boundary edges and vertices and compute Laplace-Beltrami operator using LibIgl.
 

---
# FAQ 

* There is an run-time error in windows env. : make sure you select "Visual Studio 201X - **Win64**" in CMake's Configuration (most of computers are 64bits now)
* CMake shows it cannot find libigl : library is missing. please read *STEP-2. get libigl* section at https://github.com/smartgeometry-ucl/compM080-compGV18-2019/tree/master/tutorial_002/cpp#step-2-get-libigl 

---
# Useful references 
* [Eigen quick references](http://eigen.tuxfamily.org/dox/group__QuickRefPage.html)
* [Eigen for matlab users](http://igl.ethz.ch/projects/libigl/matlab-to-eigen.html)
* [Libigl tutorials](https://libigl.github.io/tutorial/)
* [C++ template introduction](http://www.cplusplus.com/doc/oldtutorial/templates/)


---

# Install 
## STEP-1. get the code of tutorial 4:

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
# Boundary Vertices/Edges
In this practice, you will practice to use these two LibIgl function:  
1) igl::boundary_loop  
2) igl::boundary_facets

**Your job** is to implement the function `get_boundary_edges` and `get_boundary_vex` at `tutorial_004\cpp\practices_laplacian\src\mytools.cpp` 

Note that an edge is a boundary edge, if it has only a single face connected to it.

**get_boundary_edges**  
````

void get_boundary_edges( 
	Eigen::MatrixXi const & F,
	Eigen::MatrixXi & out_b_edge)
{
	// You job is to use igl::boundary_facets to find boundary edges
	//  
	// F : input faces
	//
	// out_bedge : output edges, K-by-2 (two vertices index)
	// 
	//  Hints:
	//   Eigen::MatrixXi b_edge, K-by-2
    //     
	//     igl::boundary_facets( F_in , b_edge )
    
  // ......
  
````


**get_boundary_vex**  
````
void get_boundary_vex(
	Eigen::MatrixXd const & V, 
	Eigen::MatrixXi const & F,
	Eigen::MatrixXd & out_bvex)
{
	// You job is to use igl::boundary_loop to find boundary vertices
	// 
	// V : input vertices, N-by-3
	// F : input faces
	//
	// out_bvex : output vertices, K-by-3
	// 
	//  Hints:
	//   Eigen::VectorXi b_bex_index
	//     igl::boundary_loop( F_in , b_bex_index )
	// 
````

## Results
### 
![](/tutorial_004/cpp/doc_img/res1.JPG ) 


---
# Laplace-Beltrami operator 
In this practice, you will practice to use these two LibIgl function:  
1) igl::cotmatrix(V, F, L)  
2) igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, Area )  
   
   L and Area are Eigen::SparseMatrix<double>  

![](/tutorial_004/cpp/doc_img/cota.jpg )  


**Your job** is to implement the function `compute_H` at `tutorial_004\cpp\practices_laplacian\src\mytools.cpp`  

````
void compute_H(
	Eigen::MatrixXd const & V,
	Eigen::MatrixXi const & F,
	Eigen::VectorXd & H
	)
{
	// You job is to use igl::cotmatrix, igl::massmatrix, igl::invert_diag to compute mean curvature H at each vertex
	// 
	// V : input vertices, N-by-3
	// F : input faces
	//
	// H : output vertices, K-by-1
	// 
````

## Results
### 
![](/tutorial_004/cpp/doc_img/res2.JPG ) 


Example of using SparseMatrix  
**Most block operation in Eigen::SparseMatrix are read-only.

[SparseQuickRefPage](https://eigen.tuxfamily.org/dox/group__SparseQuickRefPage.html)  

````
// copy L to A 

Eigen::SparseMatrix<double> L(N,N);
// ...

Eigen::SparseMatrix<double> A(N,N);
for (size_t i = 0; i < N; i++)
{
    for (size_t j = 0; j < N; j++)
    {
        // add new value 
        A.insert(i, j) = L.coeff(i, j);
        // update 
        A.coeff(i,j) += 1.0;
    }
}

````


````
//conversion
Eigen::SparseMatrix<double> A(N,N);
Eigen::MatrixXd dsA = A.toDense();

````
