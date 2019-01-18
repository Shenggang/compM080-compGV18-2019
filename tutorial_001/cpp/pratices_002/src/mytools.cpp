#include "mytools.h"

void average_depth(std::vector<Eigen::MatrixXf> const & in_depth,
	int const NUM_FRAMES,
	Eigen::MatrixXf & out_depth)
{
	// average depth:
	// This function averages depth value over K frames (K=NUM_FRAMES)
	// Note that you need to bypass the invalid depth value, whose depth pixel value equal to 0
	//
	// syntax hints:
	// Eigen 
	//  R = (Q.array() == 0).select(P, R) 
	// Numpy equivalent
	//  R[Q==0]=P[Q==0]
	// 
    
    
    out_depth = in_depth[0];
    
}

