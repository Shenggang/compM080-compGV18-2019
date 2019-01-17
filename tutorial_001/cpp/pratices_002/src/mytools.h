#ifndef MYTOOLS
#define MYTOOLS

#define NOMINMAX

#include <vector>
#include <iostream>
#include <Eigen/Dense>

void average_depth(std::vector<Eigen::MatrixXf> const & in_depth,
	int const NUM_FRAMES,
	Eigen::MatrixXf & out_depth);


#endif


