# Coursework 2 F&Q list


1. Gaussian curvature: Gaussian curvature only depends on angle deficit AND your A_i. the choices of A_i includes Barycentric area and Voronoi area. 
2. curvature visualization: just make sure you can see the difference between flat area and complex area. 

3. question 3, eigenvectors of the Laplace-Beltrami : please use cotanget version

4. Spectra and Eigen: Spectra requires Eigen 3.3.7. check your eigen's major version before use it:  
`	std::cout<<"eigen version.:"<<EIGEN_WORLD_VERSION<<","<<EIGEN_MAJOR_VERSION  << EIGEN_MINOR_VERSION<<"\n";`

5. Sparse solver opotions in cpp:  
* Spectra https://spectralib.org/  (update Eigen to [3.3.7](http://eigen.tuxfamily.org/index.php?title=Main_Page) by replacing the directory libigl\external\eigen   
* [igl::eigs(A,B,EIGS_TYPE_SM ,u,s)](https://github.com/libigl/libigl/blob/508cb9940f4d1e8e54137d5afe2fd2eb9c4dc672/include/igl/eigs.h).   solves the generalized eigenvalue problem `A u = s B u`.  

6. question 3, reference paper: [Spectral Geometry Processing with Manifold Harmonics](http://www.cs.jhu.edu/~misha/ReadingSeminar/Papers/Vallet08.pdf). Section 2.3 and 3.1
