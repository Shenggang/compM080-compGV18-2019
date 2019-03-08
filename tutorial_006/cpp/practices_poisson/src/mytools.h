#ifndef MYTOOLS
#define MYTOOLS

#define NOMINMAX

#include <vector>
#include <iostream>
#include <ctime>

#include <Eigen/Dense>
#include <Eigen/LU>

#include <igl/boundary_loop.h>
#include <igl/boundary_facets.h>
#include <igl/setdiff.h> 
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>

void get_boundary_vexIndex( Eigen::MatrixXi const & F , Eigen::VectorXi & boundary_i);

void solve_poission_eq(Eigen::MatrixXd const & V, Eigen::MatrixXi const & F, Eigen::MatrixXd & out_solved_vex);

#endif


