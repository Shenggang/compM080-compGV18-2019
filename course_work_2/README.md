# Coursework 2 F&Q list


1. Gaussian curvature: Gaussian curvature only depends on angle deficit AND your A_i. the choices of A_i includes Barycentric area and Voronoi area. 
2. curvature visualization: just make sure you can see the difference between flat area and complex area. 

3. question 3, eigenvectors of the Laplace-Beltrami : please use cotangent version

4. Spectra and Eigen: Spectra requires Eigen 3.3.7. check your eigen's major version before use it:  
`	std::cout<<"eigen version.:"<<EIGEN_WORLD_VERSION<<","<<EIGEN_MAJOR_VERSION  << EIGEN_MINOR_VERSION<<"\n";`

5. question 3, reference paper: [Spectral Geometry Processing with Manifold Harmonics](http://www.cs.jhu.edu/~misha/ReadingSeminar/Papers/Vallet08.pdf).  
Section 2.3 and 3.1 provide the details related to question3.  
Equation (7) shows how to transform between canonical and MHB space (eigen space)  
`star_0` => our `M`  
`x` => canonical coordinates  
`\tilde{x}` => the coordinates in eigen space  

6. Sparse solver opotions in cpp - use **igl::eigs**  
[igl::eigs(A,B,EIGS_TYPE_SM ,x,s)](https://github.com/libigl/libigl/blob/508cb9940f4d1e8e54137d5afe2fd2eb9c4dc672/include/igl/eigs.h).  
It solves the generalized eigenvalue problem `A x = s B x`.  
STEPS:  
i. apply a hotfix at here [url](https://gist.github.com/yushiangw/c4579efca5a7ee7c5a0e3903dc8d165a) by replacing the file `libigl\include\igl\eigs.cpp`. The original implementation has numerical problems.  
ii. solve : `(-1*C) x = s M x`. `(-1*C)` is a trick to make this matrix become P.S.D.  
iii. `(-1*C) x = s M x` is equal to `M^(-1)(-1*C) x = s x` so now `x` is the eigenvectors of `M^(-1)*(-1*C)`. 

7. Sparse solver opotions in cpp - use **Spectra**  
STEPS:  
i. download the header only library Spectra at [https://spectralib.org/]   
ii. Installation (1/2): Update Eigen to [3.3.7](http://eigen.tuxfamily.org/index.php?title=Main_Page) by replacing the directory [libigl\external\eigen]().  
check your eigen's major version before use it:
`std::cout<<"eigen version.:"<<EIGEN_WORLD_VERSION<<","<<EIGEN_MAJOR_VERSION << EIGEN_MINOR_VERSION<<"\n";`  
iii. Installation (2/2): Update your include list in CmakeList.txt ([see L16](https://gist.github.com/yushiangw/3e81f5cdf0e4f6a6cfadd8740f4117e5#file-cmakelists-txt-L16)) and do cmake again  
iv. solve the eigenvectors of the matrix `((M^(-0.5))*-1*C*(M^(-0.5))` by using `Spectra::SymEigsSolver` ([example](https://spectralib.org/quick-start.html) ).  
v. multiply the extracted eigenvectors `x` with `M^(-0.5)` so that `y=M^(-0.5)*x`. `y` is the eigenvectors of the matrix `M^(-1)*(-1*C)`  
8. Sparse solver in python: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html
