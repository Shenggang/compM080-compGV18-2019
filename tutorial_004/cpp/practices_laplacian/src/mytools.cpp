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
	Eigen::VectorXi b_bex_index;
	igl::boundary_loop(F, b_bex_index);
	out_bvex.resize(b_bex_index.rows(), 3);
	for (int i = 0; i < b_bex_index.rows(); ++i)
	{
		out_bvex.row(i) = V.row(b_bex_index(i));
	}

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
	igl::boundary_facets(F, out_b_edge);
}

double area_of_triangle(
	Eigen::MatrixXd const & vertex_1,
	Eigen::MatrixXd const & vertex_2,
	Eigen::MatrixXd const & vertex_3
	)
{
	// vertex_i : 1-by-3 matrix, a row of the vertices matrix
	//
	// return : double, area of triangle with vertices at v1,2,3
	double side[3];
	double s;
	side[0] = (vertex_2 - vertex_1).norm();
	side[1] = (vertex_3 - vertex_2).norm();
	side[2] = (vertex_1 - vertex_3).norm();
	s = 0.5*(side[0] + side[1] + side[2]);
	return sqrt(s*(s-side[0])*(s-side[1])*(s-side[2]));
}

void compute_area(
	Eigen::MatrixXd const & V,
	Eigen::MatrixXi const & F,
	Eigen::SparseMatrix<double> & Area
	)
{
	// V : input vertices, N-by-3
	// F : input faces
	//
	// Area : output diagonal matrix, N-by-N, each entry represents the normalising area around a vertex
	int faces = F.rows();
	int vertices = V.rows();
	double area_tri[faces];
	std::vector<std::vector<int>> neighbourFacesOf(vertices);
	for (int i = 0; i < faces; ++i)
	{
		neighbourFacesOf[F(i,0)].push_back(i);
		neighbourFacesOf[F(i,1)].push_back(i);
		neighbourFacesOf[F(i,2)].push_back(i);
		area_tri[i] = area_of_triangle(V.row(F(i,0)), V.row(F(i,1)), V.row(F(i,2)))/3;
	}

	Area.resize(vertices, vertices);
	Area.setZero();
	for (int i = 0; i < vertices; ++i)
	{
		double area = 0;
		for (int index:neighbourFacesOf[i])
		{
			area += area_tri[index];
		}
		Area.insert(i,i) = area;
	}
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

	Eigen::SparseMatrix<double> L, Area, AreaInv;
	igl::cotmatrix(V,F,L);
	//igl::massmatrix(V,F, igl::MASSMATRIX_TYPE_VORONOI, Area);
	compute_area(V, F, Area);

	AreaInv.resize(Area.rows(), Area.rows());
	AreaInv.setZero();
	for (int i = 0; i < Area.rows(); ++i)
	{
		AreaInv.insert(i,i) = 1/Area.coeff(i,i);
	}

	Eigen::MatrixXd LBV = AreaInv*L*V;
	//------------------------------------------
	// replace this 
	H.resize(LBV.rows());
	
	for (int i = 0; i < LBV.rows(); ++i)
	{
		H(i) = 0.5*LBV.row(i).norm();
	}
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
