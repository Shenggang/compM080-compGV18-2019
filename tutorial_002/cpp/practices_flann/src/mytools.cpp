#include "mytools.h"

void get_cube(Eigen::MatrixXd & out_V, Eigen::MatrixXi & out_F)
{
	out_V = (Eigen::MatrixXd(8, 3) <<
		0.0, 0.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, 1.0, 0.0,
		0.0, 1.0, 1.0,
		1.0, 0.0, 0.0,
		1.0, 0.0, 1.0,
		1.0, 1.0, 0.0,
		1.0, 1.0, 1.0).finished();

	out_F = (Eigen::MatrixXi(12, 3) <<
		1, 7, 5,
		1, 3, 7,
		1, 4, 3,
		1, 2, 4,
		3, 8, 7,
		3, 4, 8,
		5, 7, 8,
		5, 8, 6,
		1, 5, 6,
		1, 6, 2,
		2, 6, 8,
		2, 8, 4).finished().array() - 1;
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
	//   Your job is to implement vertex normal calculation
	//
    std::vector<std::vector<int> > VF;
	std::vector<std::vector<int> > VFi;
	igl::vertex_triangle_adjacency(V.rows(), F, VF, VFi);
	out_VN.resize(V.rows(),V.cols());

	int i = 0;
	for (i = 0; i < V.rows(); ++i)
	{
		std::vector<int> incident_faces = VF[i];
		Eigen::Matrix<double, 1, 3> vn(0,0,0);
		for (std::vector<int>::iterator it = incident_faces.begin(); it != incident_faces.end(); ++it)
		{
			vn += FN.row(*it);
		}
		vn /= incident_faces.size();
		out_VN.row(i) = vn.row(0);
	}
}

void calculate_vertex_normal_flann(Eigen::MatrixXd const & V, Eigen::MatrixXi const & F, Eigen::MatrixXd & out_VN)
{
	//
	// input:
	//   V: vertices
	//   F: face 
	// output:
	//   out_VN
	//
	// Your job is to implement vertex normal calculation vis using flann and igl:fitplane
	//  
	// igl::fit_plane(V, N, C);
	// Input:
	//   V #Vx3 matrix. The 3D point cloud, one row for each vertex.
	// Output: 
	//   N 1x3 Vector. The normal of the fitted plane.
	//   C 1x3 Vector. A point that lies in the fitted plane.
	//

	out_VN.resize(V.rows(), V.cols());
	out_VN.setZero();

	/*
	typedef nanoflann::KDTreeEigenMatrixAdaptor< Eigen::MatrixXd >  my_kd_tree_t;

	my_kd_tree_t   mat_index(V, 50);
	mat_index.index->buildIndex();

	Eigen::RowVector3d v_cen = V.colwise().sum() / V.rows();


	for (int idx = 0; idx < V.rows(); idx++) {

		Eigen::RowVector3d i_vex = V.row(idx);

		const size_t num_results = 8;
		std::vector<size_t>   ret_indexes(num_results);
		std::vector<double> out_dists_sqr(num_results);

		nanoflann::KNNResultSet<double> resultSet(num_results);

		resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
		mat_index.index->findNeighbors(resultSet, i_vex.data() , nanoflann::SearchParams(25));

		Eigen::MatrixXd Selectpoints(num_results, 3);
		Eigen::MatrixXd nn_vex(num_results, 3);

		for (size_t i = 0; i < num_results; i++) {
			//std::cout << "ret_index["<<i<<"]=" << ret_indexes[i] << " out_dist_sqr=" << out_dists_sqr[i] << endl;
			Selectpoints(i, 0) = V(ret_indexes[i], 0);
			Selectpoints(i, 1) = V(ret_indexes[i], 1);
			Selectpoints(i, 2) = V(ret_indexes[i], 2);
		}

		Eigen::RowVector3d pl_NV, Ct;

		igl::fit_plane(Selectpoints, pl_NV, Ct);

		if ((v_cen - i_vex).dot(pl_NV)  > 0) {
			pl_NV = pl_NV * -1;
		}

		out_VN.row(idx) = pl_NV;
	}
	*/
	
	// build kd-tree uses vertices V
	nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd> mat_index(V, 50);
	mat_index.index->buildIndex();

	Eigen::RowVector3d v_cen = V.colwise().sum() / V.rows();

	for (int i = 0; i < V.rows(); ++i)
	{
		Eigen::RowVector3d cur_vex = V.row(i);

		// set K nearest samples
		const size_t par_K = 8;
		
		// create results objects
		std::vector<size_t> indexes(par_K);
		std::vector<double> dists_sqr(par_K);

		// bind results
		nanoflann::KNNResultSet<double> res(par_K);
		res.init(indexes.data(), dists_sqr.data());

		// find KNN, 50 in SearchParams is ignored but kept for compatibility
		mat_index.index->findNeighbors(res, cur_vex.data(), nanoflann::SearchParams(50));

		// Find coordinates of KNN
		Eigen::MatrixXd nn_vex(indexes.size(), 3);
		for (size_t j = 0 ; j < indexes.size(); ++j)
		{
			nn_vex.row(j) = V.row(indexes[j]);
		} 

		// Fit to plane and get plane normal
		Eigen::RowVector3d normal, centre;
		igl::fit_plane(nn_vex, normal, centre);

		// Flip normal if it points towards inside
		if ((v_cen - cur_vex).dot(normal) > 0)
		{
			normal *= -1;
		}

		out_VN.row(i) = normal;
	}

}
