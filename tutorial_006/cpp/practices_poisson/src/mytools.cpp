#include "mytools.h"


void get_boundary_vexIndex(Eigen::MatrixXi const & F, Eigen::VectorXi & boundary_i)
{
    // ***
    // You job is to call boundary_loop to get the indices of boundary vertices
    //
    // igl::boundary_loop( Eigen::FaceType/Matrix , Eigen::Index List Type/Vector/Maxtrix)
	igl::boundary_loop(F, boundary_i);
}

void solve_poission_eq(Eigen::MatrixXd const & V, Eigen::MatrixXi const & F, Eigen::MatrixXd & out_solved_vex)
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

	if (boundary_i.size()==0)
	{
		std::cout << "no boundary vertex found\n";
		out_solved_vex = V;
		return;
	}

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
    
    igl::slice(C, interior_i, interior_i, C_ii);
    igl::slice(C, interior_i, boundary_i, C_ib);
    igl::slice(M, interior_i, interior_i, M_ii);
    
    
    //=================================================
    std::cout << "C_ii.row(0)=" << C_ii.row(0).head(3);
    std::cout << "expected C_ii.row(0)=" << "=-7.91335, 0, 0" << "\n\n";
    std::cout << "M_ii.row(0)=" << M_ii.row(0).head(3);
    std::cout << "expected M_ii.row(0)=" << "0.000121022, 0, 0" << "\n\n";
    
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
	
    b_xyz = -1 * C_ib * bvex;
	A = C_ii.toDense();

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(C_ii);
    solver.factorize(C_ii);
    if (solver.info() != Eigen::Success) { std::cout << "decomposition failed\n";exit(1); }

	// ... 
    for (int i = 0; i < boundary_i.rows(); ++i)
    {
        new_V.row(boundary_i[i]) = V.row(boundary_i[i]);
    }

    // now solve for each dimension (XYZ)
    for (int i = 0; i < 3; i++)
    {
        std::cout << "start solving .." << i << "\n";
        std::clock_t start = std::clock();
        
        // solve Ax=b 
        // ...
        Eigen::VectorXd solution = solver.solve(b_xyz.col(i));
        
        if (solver.info() != Eigen::Success) {std::cout << "solving failed\n";exit(1);}

        // restore vertices order   
        // ...
        for (int j = 0; j < interior_i.rows(); ++j)
        {
            new_V(interior_i[j], i) = solution[j];
        }
        
        double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        std::cout << "time cost="  << duration << '\n';
    }
    out_solved_vex = new_V;
}
