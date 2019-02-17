#ifndef MYTOOLS
#define MYTOOLS

#define NOMINMAX

#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <igl/vertex_triangle_adjacency.h>
#include <nanoflann.hpp>
#include <igl/fit_plane.h>
#include <math.h>

#define PI 3.14159265
void get_cube(Eigen::MatrixXd & out_V , Eigen::MatrixXi & out_F);


void calculate_vertex_normal(
	Eigen::MatrixXd const & V,
	Eigen::MatrixXi const & F,
	Eigen::MatrixXd const & FN,
	Eigen::MatrixXd & out_VN);

void calculate_vertex_normal_flann(
	Eigen::MatrixXd const & V, 
	Eigen::MatrixXd & out_VN);

void get_rotation_X(
	double const & degrees,
	Eigen::Matrix3d& mat);

void get_rotation_Y(
	double const & degrees,
	Eigen::Matrix3d& mat);

void get_rotation_Z(
	double const & degrees,
	Eigen::Matrix3d& mat);

void subsample(
	int const & total,
	float const & proportion,
	std::vector<int> & indices);

void register_closest_point(
	Eigen::MatrixXd const & V_tar,
	Eigen::MatrixXd const & V_src,
	Eigen::MatrixXd const & N_tar,
	Eigen::MatrixXd & V_ref,
	Eigen::MatrixXd & V_NP,
	Eigen::MatrixXd & N_ref);

void register_closest_point(
	Eigen::MatrixXd const & V_tar,
	Eigen::MatrixXd const & V_src,
	Eigen::MatrixXd & V_ref,
	Eigen::MatrixXd & V_NP);

double estimate_rotation_translation(
	Eigen::MatrixXd const & V_tar,
	Eigen::MatrixXd const & V_src,
	Eigen::Matrix3d & R,
	Eigen::Vector3d & t);

double estimate_rotation_translation(
	Eigen::MatrixXd const & V_tar,
	Eigen::MatrixXd const & N_tar,
	Eigen::MatrixXd const & V_src,
	Eigen::Matrix3d & R,
	Eigen::Vector3d & t);



#endif


