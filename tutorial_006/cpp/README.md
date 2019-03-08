# COMPM080/GV18 Tutorial 6 - C++ 

In this tutorial, we introduce how to solve Poisson equation.

---
# FAQ 

* CMake shows it cannot find libigl : library is missing. please read *STEP-2. get libigl* section at https://github.com/smartgeometry-ucl/compM080-compGV18-2019/tree/master/tutorial_002/cpp#step-2-get-libigl 

---
# Useful references 
* [Eigen quick references](http://eigen.tuxfamily.org/dox/group__QuickRefPage.html)
* [Eigen for matlab users](http://igl.ethz.ch/projects/libigl/matlab-to-eigen.html)
* [Libigl tutorials](https://libigl.github.io/tutorial/)
* [Linear Solver](https://eigen.tuxfamily.org/dox/group__TopicLinearAlgebraDecompositions.html)
---

# Install 
## STEP-1. get the code of tutorial 6:

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
# Poisson Equation
![](/tutorial_006/cpp/doc_img/eq.JPG) 

## Results
### 
![](/tutorial_006/cpp/doc_img/input.JPG) 
![](/tutorial_006/cpp/doc_img/res.JPG) 


## Solving Linear System
* [Linear Solver](https://eigen.tuxfamily.org/dox/group__TopicLinearAlgebraDecompositions.html)

Input scale: 4000-by-4000 matrix  
DenseLU : time cost=5.003 sec  
SparseLU : time cost=0.021 sec  

### Dense Solver

solve Ax = b.
````
//#include <Eigen/LU>  
x = A.lu()  .solve(b));   
//#include <Eigen/Cholesky>  
x = A.ldlt().solve(b));  
//#include <Eigen/Cholesky>
x = A.llt() .solve(b));         
// #include <Eigen/QR>  
x = A.qr()  .solve(b));       
// #include <Eigen/SVD>       
x = A.svd() .solve(b));  
````

### Sparse Solver
* [Sparse Solver Concept](https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html)

solve Ax = b  
1. `Eigen::SparseMatrix<double> A;` and  `Eigen::VectorXd b, x;`
2. `SolverClassName<SparseMatrix<double> > solver;`
3. `solver.compute(A);`
4. `x = solver.solve(b);`

example:  
[SparseLU](https://eigen.tuxfamily.org/dox/classEigen_1_1SparseLU.html)  

````
Eigen::SparseMatrix<double> A = ...;
Eigen::VectorXd b = ...;

// solve Ax = b
Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

solver.analyzePattern(A_sp);
solver.factorize(A_sp);

if (solver.info() != Eigen::Success) { std::cout << "decomposition failed\n";exit(1); }

Eigen::VectorXd x = solver.solve(b);

if (solver.info() != Eigen::Success) {std::cout << "solving failed\n";exit(1);}

````

---
### **Your job** is to implement the function `get_boundary_vexIndex` and `solve_poission_eq` at [`tutorial_006\cpp\practices_poisson\src\mytools.cpp`](tutorial_006\cpp\practices_poisson\src\mytools.cpp)   

````
void get_boundary_vexIndex(Eigen::MatrixXi const & F, Eigen::VectorXi & boundary_i)
{
    // ***
    // You job is to call boundary_loop to get the indices of boundary vertices
    //
    // igl::boundary_loop( Eigen::FaceType/Matrix , Eigen::Index List Type/Vector/Maxtrix)


}
````

````
void solve_poission_eq(
        Eigen::MatrixXd const & V, 
        Eigen::MatrixXi const & F, 
        Eigen::MatrixXd & out_solved_vex)
{
    //
    // Solve for L * x = 0 s.t. boundary vertices remain fixed.
    //
    
    // get L
    Eigen::SparseMatrix<double> C, M, MInv;
    igl::cotmatrix(V, F, C);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    igl::invert_diag(M, MInv);

    // get the indices of boundary vertices  
    Eigen::VectorXi boundary_i;
    get_boundary_vexIndex(F, boundary_i);


    // List of all vertex indices
    Eigen::VectorXi all;
    Eigen::VectorXi interior_i;

    // I = low:step:hi
    // igl::colon(low,step,hi,I);
    igl::colon<int>(0, V.rows() - 1, all);

    // get the indices of interior vertices
    // interior_i = all_i - boundary_i
    Eigen::VectorXi IA;
    igl::setdiff(all, boundary_i, interior_i, IA);

    // get boundary vertices
    Eigen::MatrixXd bvex(boundary_i.rows(), 3);

    for (size_t i = 0; i < boundary_i.rows(); i++)
    {
        bvex.row(i) = V.row(boundary_i[i]);
    }

    //================================================= 
    // ***
    // You job is to slice Laplacian matrix to get C_ii, C_ib, and M_ii
    //
    // hints:
    // use igl::slice 
    // B = A(I,J)
    // igl::slice(A,I,J,B)
    //
    // ***
    
    Eigen::SparseMatrix<double>  C_ii, C_ib;
    Eigen::SparseMatrix<double>  M_ii;
    
    // ...
        
    //=================================================
    // ***
    // You job is to construct Ax=b 
    // where A = C_ii 
    // and   b = -1 * C_ib * x_b 
    // then solve Ax=b (e.g. use A.lu.solve(b) )
    // https://eigen.tuxfamily.org/dox/group__TopicLinearAlgebraDecompositions.html
    // 
    // Note that you will need to solve Ax=b for each dimension (xyz).
    // ***

    Eigen::MatrixXd new_V(V.rows(), V.cols());
    
    // set your A and b
    Eigen::MatrixXd b_xyz;
    Eigen::MatrixXd A;
    
	A = C_ii.toDense();
    
    // solve for each dimension (XYZ)
    for (int i = 0; i < 3; i++)
    {
        // solve Ax=b 
        // ...
        
        // restore vertices order
        // ...
        
    }
    
    out_solved_vex = new_V;
}

````

---

### get Laplace-Beltrami matrix
1) igl::cotmatrix(V, F, L)  
2) igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, Area )  
   
   L and Area are Eigen::SparseMatrix<double>  

### Example of using SparseMatrix  
**Eigen::SparseMatrix do NOT support many block operations such as ::head(3), ::topLeftCorner  
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
