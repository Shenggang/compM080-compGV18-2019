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
    int rows = in_depth[0].rows();
	int columns = in_depth[0].cols();
	out_depth.resize(rows, columns);

	int i,j,k;
	for (i = 0; i < rows; ++i)
	{
		for (j = 0; j < columns; ++j)
		{
			float sum = 0;
			int count = 0;
			for (k = 0; k < NUM_FRAMES; ++k)
			{
				if (in_depth[k](i,j) == 0)
				{
					continue;
				}
				sum += in_depth[k](i,j);
				count++;
			}
			out_depth(i,j) = sum/count;
		}
	}

}

