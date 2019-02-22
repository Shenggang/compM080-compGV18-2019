#ifndef MYTOOLS
#define MYTOOLS

#define NOMINMAX

#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <igl/vertex_triangle_adjacency.h>
#include <random>
#include <Eigen/LU>
#include <igl/boundary_loop.h>
#include <igl/boundary_facets.h>
#include <igl/setdiff.h> 
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>


void get_boundary_vex(Eigen::MatrixXd const & V_in, Eigen::MatrixXi const & F_in, Eigen::MatrixXd & out_bvex);
void get_boundary_edges(Eigen::MatrixXi const & F_in, Eigen::MatrixXi & out_bedge);

void compute_H(Eigen::MatrixXd const & V, Eigen::MatrixXi const & F, Eigen::VectorXd & H);

void calculate_vertex_normal(
	Eigen::MatrixXd const & V, 
	Eigen::MatrixXi const & F, 
	Eigen::MatrixXd const & FN,
	Eigen::MatrixXd & out_VN);



#endif


