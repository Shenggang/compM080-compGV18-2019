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

void calculate_vertex_normal_flann(Eigen::MatrixXd const & V, Eigen::MatrixXd & out_VN)
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
	
	// build kd-tree uses vertices V
	nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd> mat_index(V, 50);
	mat_index.index->buildIndex();

	Eigen::RowVector3d v_cen = V.colwise().sum() / V.rows();

	for (int i = 0; i < V.rows(); ++i)
	{
		Eigen::RowVector3d cur_vex = V.row(i);

		// set K nearest samples
		const size_t par_K = 20;
		
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

void get_rotation_X(double const & degrees, Eigen::Matrix3d& mat)
{
	//
	// input:
	//   degrees: double. Degrees to be turned.
	// output:
	//   mat: 3x3 Matrix. Rotation matrix about X-axis
	//
	mat.setZero();
	double rad = degrees*PI/180;
	mat(1,1) = cos(rad);
	mat(2,2) = cos(rad);
	mat(1,2) = -sin(rad);
	mat(2,1) = sin(rad);
	mat(0,0) = 1;
}

void get_rotation_Y(double const & degrees, Eigen::Matrix3d& mat)
{
	//
	// input:
	//   degrees: double. Degrees to be turned.
	// output:
	//   mat: 3x3 Matrix. Rotation matrix about Y-axis
	//
	mat.setZero();
	double rad = degrees*PI/180;
	mat(0,0) = cos(rad);
	mat(2,2) = cos(rad);
	mat(0,2) = sin(rad);
	mat(2,0) = -sin(rad);
	mat(1,1) = 1;
}

void get_rotation_Z(double const & degrees, Eigen::Matrix3d& mat)
{
	//
	// input:
	//   degrees: double. Degrees to be turned.
	// output:
	//   mat: 3x3 Matrix. Rotation matrix about Z-axis
	//
	mat.setZero();
	double rad = degrees*PI/180;
	mat(0,0) = cos(rad);
	mat(1,1) = cos(rad);
	mat(0,1) = -sin(rad);
	mat(1,0) = sin(rad);
	mat(2,2) = 1;
}

void subsample(int const & total, float const & proportion, std::vector<int> & indices)
{
	//
	// input:
	//   total: 		int. Total number of items to sample from.
	//   proportion:    double. Probability of being selected.
	// output:
	//   indices:		int vector. List of selected indices.
	//
	std::default_random_engine generator;
  	std::uniform_real_distribution<double> distribution(0.0,1.0);

	indices.clear();
	for (int i = 0; i < total; ++i)
	{
		if (distribution(generator) <=  proportion)
		{
			indices.push_back(i);
		}
	}
}

void register_closest_point(Eigen::MatrixXd const & V_tar, Eigen::MatrixXd const & V_src, Eigen::MatrixXd const & N_tar,
									Eigen::MatrixXd & V_ref, Eigen::MatrixXd & V_NP, Eigen::MatrixXd & N_ref)
{
	//
	// input:
	//   V_tar: Nx3 Matrix. Target point cloud (being the reference)
	//   V_src: Mx3 Matrix. Source point cloud (to be aligned)
	//   N_tar: Nx3 Matrix. Normal at each target point in V_tar.
	// output:
	//   V_ref: subset of V_tar. The not rejected points in V_tar.
	//   V_NP: The nearest match in V_src for each point in V_ref.
	//   N_ref : subset of N_tar. similar to V_ref.
	//
	std::vector<int> selection_ref, selection_src;
	selection_ref.clear();
	selection_src.clear();

	nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd> mat_index(V_src, 50);
	mat_index.index->buildIndex();

	for (int i = 0; i < V_tar.rows(); ++i)
	{
		Eigen::RowVector3d cur_vex = V_tar.row(i);

		// set K nearest samples
		const size_t par_K = 1;
		
		// create results objects
		std::vector<size_t> indexes(par_K);
		std::vector<double> dists_sqr(par_K);

		// bind results
		nanoflann::KNNResultSet<double> res(par_K);
		res.init(indexes.data(), dists_sqr.data());

		// find NN, 50 in SearchParams is ignored but kept for compatibility
		mat_index.index->findNeighbors(res, cur_vex.data(), nanoflann::SearchParams(50));

		// store result
		selection_ref.push_back(i);
		selection_src.push_back(indexes[0]);
	}
	V_NP.resize(selection_src.size(), V_tar.cols());
	V_ref.resize(selection_ref.size(), V_tar.cols());
	N_ref.resize(selection_ref.size(), N_tar.cols());

	for (int i = 0; i < selection_ref.size(); ++i)
	{
		V_ref.row(i) = V_tar.row(selection_ref[i]);
		N_ref.row(i) = N_tar.row(selection_ref[i]);
		V_NP.row(i) = V_src.row(selection_src[i]);
	}
}

void register_closest_point(Eigen::MatrixXd const & V_tar, Eigen::MatrixXd const & V_src, Eigen::MatrixXd & V_ref, Eigen::MatrixXd & V_NP)
{
	//
	// input:
	//   V_tar: Nx3 Matrix. Target point cloud (being the reference)
	//   V_src: Mx3 Matrix. Source point cloud (to be aligned)
	// output:
	//   V_ref: subset of V_tar. The not rejected points in V_tar.
	//   V_NP: The nearest match in V_src for each point in V_ref.
	//
	Eigen::MatrixXd N_tar, N_ref;
	N_tar.resize(V_tar.rows(), V_tar.cols());
	register_closest_point(V_tar, V_src, N_tar, V_ref, V_NP, N_ref);
}

double estimate_rotation_translation(Eigen::MatrixXd const & V_tar, Eigen::MatrixXd const & V_src, Eigen::Matrix3d & R, Eigen::Vector3d & t)
{
	//
	// input:
	//   V_tar: Nx3 Matrix. Target point cloud (being the reference)
	//   V_src: Nx3 Matrix. Source point cloud, ordered such that each point is closest in V_src to the corresponding point in V_tar
	// output:
	//   R: 3x3 Matrix. The best rotation that aligns two clouds
	//	 t: 3x1 Vector. The best translation that aligns two cloud centroids.
	// return:
	//   double  quadratic loss
	//

	Eigen::Vector3d tar_mean, src_mean;

	// calculate centroids
	tar_mean = V_tar.colwise().sum()/V_tar.rows();
	src_mean = V_src.colwise().sum()/V_src.rows();

	// calculate deviation from centriods
	Eigen::MatrixXd dev_tar, dev_src;
	dev_tar = V_tar - tar_mean.transpose().replicate(V_tar.rows(), 1);
	dev_src = V_src - src_mean.transpose().replicate(V_src.rows(), 1);

	// calculate covariance matrix
	Eigen::Matrix3d cov;
	cov = dev_src.transpose()*dev_tar;

	// compute SVD
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);

	R = svd.matrixV() * svd.matrixU().transpose();
	if (R.determinant() < 0)
	{
		R.col(2) *= -1;
	}

	t = tar_mean - R*src_mean;

	//compute loss
	Eigen::MatrixXd residue;
	residue = (V_tar.transpose() - (R*V_src.transpose() + t.replicate(1, V_src.rows()))).transpose();

	double loss = 0;
	for (int i = 0; i < residue.rows(); ++i)
	{
		loss += residue.row(i).squaredNorm();
	}
	std::cout << "Loss = " << loss << std::endl;
	return loss;
}

double estimate_rotation_translation(Eigen::MatrixXd const & V_tar, Eigen::MatrixXd const & N_tar, Eigen::MatrixXd const & V_src, Eigen::Matrix3d & R, Eigen::Vector3d & t)
{
	//
	// input:
	//   V_tar: Nx3 Matrix. Target point cloud (being the reference)
	//   N_tar: Nx3 Matrix. Normal at each vertex in V_tar.
	//   V_src: Nx3 Matrix. Source point cloud, ordered such that each point is closest in V_src to the corresponding point in V_tar
	// output:
	//   R: 3x3 Matrix. The best rotation that aligns two clouds
	//	 t: 3x1 Vector. The best translation that aligns two cloud centroids.
	// return:
	//   double  quadratic loss
	//

	// Formulate Least square problem
	Eigen::MatrixXd A, b;
	int N = V_tar.rows();
	A.resize(N,6);
	b.resize(N,1);
	for (int i = 0; i < N; ++i)
	{
		Eigen::RowVector3d n,p,q;
		n = N_tar.row(i);
		p = V_tar.row(i);
		q = V_src.row(i);
		double a1,a2,a3;
		A(i,0) = n(2) * q(1) - n(1) * q(2);
		A(i,1) = n(0) * q(2) - n(2) * q(0);
		A(i,2) = n(1) * q(0) - n(0) * q(1);
		for (int j = 0; j < 3; ++ j)
		{
			A(i,j+3) = n(j);
			b(i) += n(j)*(p(j) - q(j));
		}
	}

	// // compute SVD
	// std::cout << "compute svd" << std::endl;
	// Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

	// // Construct pseudo-inverse
	// std::cout << "inverse" << std::endl;
	// Eigen::MatrixXd sigmaInv;
	// sigmaInv.resize(A.cols(), A.rows());
	// sigmaInv.setZero();
	// int size = (A.rows() < A.cols()) ? A.rows() : A.cols();
	// for (int i = 0; i < size; ++i)
	// {
	// 	auto sv = svd.singularValues()[i];
	// 	sigmaInv(i,i) = (sv == 0) ? 0 : 1/sv; 
	// }
	// Eigen::MatrixXd Ainv;
	// Ainv = svd.matrixV() * sigmaInv * (svd.matrixU().transpose());
	// std::cout << Ainv*A << std::endl;

	// construct pseudo-inverse
	Eigen::MatrixXd Ainv, ATA;
	ATA = A.transpose()*A;
	Ainv = ATA.inverse()*(A.transpose());

	// Solve for x
	Eigen::VectorXd x;
	x = Ainv * b;

	// Retrieve t from x
	t = x.bottomRows(3);

	// Reconstruct R from euler angles
	Eigen::Matrix3d rotx, roty, rotz;
	get_rotation_X(x(0)/PI*180, rotx); get_rotation_Y(x(1)/PI*180, roty); get_rotation_Z(x(2)/PI*180, rotz); 
	R = rotz * roty * rotx;

	// compute loss
	Eigen::MatrixXd residue;
	residue = (V_tar.transpose() - (R*V_src.transpose() + t.replicate(1, V_src.rows()))).transpose();

	double loss = 0;
	for (int i = 0; i < residue.rows(); ++i)
	{
		loss += pow(residue.row(i).dot(N_tar.row(i)), 2);
	}
	std::cout << "Loss = " << loss << std::endl;
	return loss;
}
