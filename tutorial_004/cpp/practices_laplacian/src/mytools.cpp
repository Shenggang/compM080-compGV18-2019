#include "mytools.h"

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

}

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
	//   Eigen::MatrixXi b_edge
	//     igl::boundary_facets( F_in , b_edge )
	//  

}

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
	//
	// Hints
	// Compute Laplace-Beltrami operator
	//	Eigen::SparseMatrix<double> L, Area, AreaInv;
	//	igl::cotmatrix(V, F, L);
	//	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, Area );
	//	

	
	//------------------------------------------
	// replace this 
	H.resize(V.rows());
	H.setZero();
	
    

}


void calculate_vertex_normal(Eigen::MatrixXd const & V, Eigen::MatrixXi const & F, Eigen::MatrixXd const & FN, Eigen::MatrixXd & out_VN)
{
	//
	// input:
	//   V: vertices
	//   F: face 
	//   FN: face normals
	// output:
	//   out_VN
	// 
	std::vector<std::vector<int> > VF;
	std::vector<std::vector<int> > VFi;
	igl::vertex_triangle_adjacency(V.rows(), F, VF, VFi);


	Eigen::MatrixXd VN(V.rows(), 3);

	for (int i = 0; i < V.rows(); i++)
	{
		Eigen::RowVector3d nv(0, 0, 0);

		for (int j = 0; j < VF[i].size(); j++)
		{
			nv = nv + FN.row(VF[i][j]);
		}

		nv = nv / VF[i].size();

		VN.row(i) = nv;
	}

	out_VN = VN;
}
