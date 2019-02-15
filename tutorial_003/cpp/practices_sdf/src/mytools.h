#ifndef MYTOOLS
#define MYTOOLS

#define NOMINMAX

#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <igl/vertex_triangle_adjacency.h>
#include <random>
#include <Eigen/LU>
void get_cube(Eigen::MatrixXd & out_V , Eigen::MatrixXi & out_F);


void calculate_vertex_normal(
	Eigen::MatrixXd const & V, 
	Eigen::MatrixXi const & F, 
	Eigen::MatrixXd const & FN,
	Eigen::MatrixXd & out_VN);


void create_input_pcd(Eigen::MatrixXd const & V, Eigen::MatrixXd const & VN,
	Eigen::MatrixXd & out_pts, Eigen::VectorXd & out_sign_dist);

void fit_SDF(Eigen::MatrixXd const & MODEL_PTS, Eigen::MatrixXd const & SIGN_DIST, Eigen::VectorXd & out_lambda_c);

void calculate_SDF(Eigen::MatrixXd const & MODEL_PTS,
	Eigen::VectorXd const & LAMBDA_C,
	Eigen::MatrixXd const & SAMPLE_PTS,
	Eigen::VectorXd & out_sign_dist);


#endif


