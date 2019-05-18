#ifndef MYTOOLS
#define MYTOOLS

#define NOMINMAX

#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <random>
#include <igl/boundary_loop.h>
#include <igl/boundary_facets.h>
#include <igl/eigs.h>
#include <igl/setdiff.h> 
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <math.h>
#include <cmath>


void get_boundary_vex(Eigen::MatrixXd const & V_in, Eigen::MatrixXi const & F_in, Eigen::MatrixXd & out_bvex);
void get_boundary_edges(Eigen::MatrixXi const & F_in, Eigen::MatrixXi & out_bedge);

void compute_uniform_matrix(
	Eigen::MatrixXd const & V,
	Eigen::MatrixXi const & F,
	Eigen::SparseMatrix<double> & L);

double angle_between_vectors(
	Eigen::Vector3d const & v1,
	Eigen::Vector3d const & v2);

void compute_cotan_matrix(
	Eigen::MatrixXd const & V,
	Eigen::MatrixXi const & F,
	Eigen::SparseMatrix<double> & L);

void compute_mean_matrix(
	Eigen::MatrixXd const & V,
	Eigen::MatrixXi const & F,
	Eigen::SparseMatrix<double> & L);

double area_of_triangle(
	Eigen::MatrixXd const & vertex_1,
	Eigen::MatrixXd const & vertex_2,
	Eigen::MatrixXd const & vertex_3);

void compute_area(
	Eigen::MatrixXd const & V,
	Eigen::MatrixXi const & F,
	Eigen::SparseMatrix<double> & Area);

void compute_H(Eigen::MatrixXd const & V, Eigen::MatrixXi const & F, Eigen::VectorXd & H);
void compute_H(Eigen::MatrixXd const & V, Eigen::SparseMatrix<double> const & L, Eigen::VectorXd & H);

void compute_gaussian_curvature(Eigen::MatrixXd const & V, Eigen::MatrixXi const & F, Eigen::VectorXd & K);

void compute_eigendecomposition(
	Eigen::MatrixXd const & V, Eigen::MatrixXi const & F, 
	int const n, Eigen::MatrixXd & eigen_mesh, Eigen::SparseMatrix<double> &Area);

void tutte_flattening_uniform(
	Eigen::MatrixXd const & V, Eigen::MatrixXi const & F, Eigen::MatrixXd & out_V);

void tutte_flattening_LB(
	Eigen::MatrixXd const & V, Eigen::MatrixXi const & F, Eigen::MatrixXd & out_V);

void tutte_flattening_mean(
	Eigen::MatrixXd const & V, Eigen::MatrixXi const & F, Eigen::MatrixXd & out_V);

void compute_face_angle(
	Eigen::MatrixXd const & V_original, Eigen::MatrixXi const & F, Eigen::MatrixXd & FA);

double compute_distortion(Eigen::MatrixXd const & FA_original, Eigen::MatrixXd const & FA_cmpr);

void calculate_vertex_normal(
	Eigen::MatrixXd const & V, 
	Eigen::MatrixXi const & F, 
	Eigen::MatrixXd const & FN,
	Eigen::MatrixXd & out_VN);



#endif


