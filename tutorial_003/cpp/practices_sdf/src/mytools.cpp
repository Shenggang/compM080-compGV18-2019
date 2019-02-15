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

void create_input_pcd(Eigen::MatrixXd const & V, Eigen::MatrixXd const & VN, 
	Eigen::MatrixXd & out_pts, Eigen::VectorXd & out_sign_dist)
{
	// create input point cloud for SDF fitting problem
	// Your points should contains V itself and V+VN*randome_displacement and V-VN*randome_displacement 
	// 

	static std::default_random_engine generator;

	// select a random vertex
	std::uniform_real_distribution<double> distribution(0, 1);

	const size_t NUM_V = V.rows();

	out_pts.resize(NUM_V * 3, 3);
	out_sign_dist.resize(NUM_V * 3, 1);


	for (int i = 0; i < V.rows(); i++)
	{
		double rd_disp1 = distribution(generator);
		double rd_disp2 = distribution(generator);

		out_pts.row(i) = V.row(i);
		out_pts.row(i + NUM_V) = V.row(i) + VN.row(i) * rd_disp1;
		out_pts.row(i + NUM_V * 2) = V.row(i) - VN.row(i) * rd_disp2;

		out_sign_dist(i) = 0;
		out_sign_dist(i + NUM_V)=rd_disp1;
		out_sign_dist(i + NUM_V * 2)=-1 * rd_disp2;

	}


}


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
    // Hints: check "Solve Ax = b section" at https://eigen.tuxfamily.org/dox/AsciiQuickReference.txt
	//
	const size_t N = MODEL_PTS.rows();

	Eigen::MatrixXf A(N, N);
	Eigen::MatrixXf P(N, 4);
	Eigen::MatrixXf b(N + 4, 1);

	// MX=b 
	Eigen::MatrixXf M(N + 4, N + 4);
    
    
    // [remove this]
    out_lambda_c.resize(N+4,1);
    out_lambda_c.setZero();
    
}


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
     
	const size_t N = SAMPLE_PTS.rows();

	// [remove this]
	out_sign_dist.resize(N);
	out_sign_dist.setZero();
}
