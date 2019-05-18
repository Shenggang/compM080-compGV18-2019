#ifndef SEAMSTER
#define SEAMSTER

#define NOMINMAX

#include <iostream>
#include <random>
#include <vector>
#include <Eigen/Dense>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/boundary_loop.h>
#include <igl/boundary_facets.h>
#include <igl/cut_mesh.h>
#include <igl/eigs.h>
#include <igl/edge_topology.h>
#include <igl/setdiff.h> 
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/sort.h>
#include <math.h>
#include "mytools.h"

using namespace std;
using namespace Eigen;

class EdgeTopology
{
public:
    MatrixXi EV;
    MatrixXi EF;
    MatrixXi FE;

    EdgeTopology(MatrixXd &,MatrixXi &);
    EdgeTopology(MatrixXd,MatrixXi,MatrixXi,MatrixXi,MatrixXi);

    void computeDirectedEdgeCorrespondence();

private:
    MatrixXd V;
    MatrixXi F;
    bool isDirected;
};

template <typename T>
bool isAinVec(T A, vector<T> vec);

int findEdge(MatrixXi const & EV, int source, int end);

void compute_seam(MatrixXd const & V, MatrixXi const & F, EdgeTopology const & et, vector<int> const & Terminals, vector<int> & out_SEi);

void cut_mesh(MatrixXd const & V, MatrixXi const & F, EdgeTopology const & et, vector<int> const & SEi, MatrixXd & out_V, MatrixXi & out_F); 

void find_seam_terminals(MatrixXd const & V, MatrixXi const & F, vector<int> & Terminals);

#endif