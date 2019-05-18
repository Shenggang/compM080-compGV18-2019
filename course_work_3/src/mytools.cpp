#include "mytools.h"

void get_boundary_vex(
	Eigen::MatrixXd const & V, 
	Eigen::MatrixXi const & F,
	Eigen::MatrixXd & out_bvex)
{
	// You job is to use igl::boundary_loop to find boundary vertices
	// 
	// V : input vertices, N-by-3
	// F : input faces
	//
	// out_bvex : output vertices, K-by-3
	// 
	//  Hints:
	//   Eigen::VectorXi b_bex_index
	//     igl::boundary_loop( F_in , b_bex_index )
	// 
	Eigen::VectorXi b_bex_index;
	igl::boundary_loop(F, b_bex_index);
	out_bvex.resize(b_bex_index.rows(), 3);
	for (int i = 0; i < b_bex_index.rows(); ++i)
	{
		out_bvex.row(i) = V.row(b_bex_index(i));
	}

}

void get_boundary_edges( 
	Eigen::MatrixXi const & F,
	Eigen::MatrixXi & out_b_edge)
{
	// You job is to use igl::boundary_facets to find boundary edges
	//  
	// F : input faces
	//
	// out_bedge : output edges, K-by-2 (two vertices index)
	// 
	//  Hints:
	//   Eigen::MatrixXi b_edge
	//     igl::boundary_facets( F_in , b_edge )
	//  
	igl::boundary_facets(F, out_b_edge);
}

void compute_uniform_matrix(
	Eigen::MatrixXd const & V,
	Eigen::MatrixXi const & F,
	Eigen::SparseMatrix<double> & L
)
{
	// V : input vertices, N-by-3
	// F : input faces
	//
	// L : output uniform laplace operator, N-by-N.
	int nv = V.rows();
	std::vector<std::vector<int> > A;
	igl::adjacency_list(F,A);
	L.resize(nv,nv);
	L.reserve(10*nv);

	std::vector<Eigen::Triplet<double>> Ltemp;
	for (int i = 0; i < nv; ++i)
	{
		int n = A[i].size();
		Ltemp.push_back(Eigen::Triplet<double>(i,i,-1));
		for (int j:A[i])
		{
			Ltemp.push_back(Eigen::Triplet<double>(i,j,double(1)/n));
		}
	}
	L.setFromTriplets(Ltemp.begin(), Ltemp.end());
}

double angle_between_vectors(
	Eigen::Vector3d const & v1,
	Eigen::Vector3d const & v2
	)
{
	// returns angle between two vectors
	double adotb = v1.dot(-v2);
	if (adotb == 0)
		return 0;
	adotb = adotb/v1.norm()/v2.norm();
	if (adotb >=1)
		return 0;
	if (adotb <= -1)
		return M_PI;
	return acos(adotb);
}

void compute_cotan_matrix(
	Eigen::MatrixXd const & V,
	Eigen::MatrixXi const & F,
	Eigen::SparseMatrix<double> & L
	)
{
	// V : input vertices, N-by-3
	// F : input faces
	//
	// L : output cotangent matrix, N-by-N.
	int nv = V.rows();
	int nf = F.rows();
	std::vector<std::vector<int> > A;
	igl::adjacency_list(F,A);
	
	std::vector<Eigen::Triplet<double>> Ltemp;

	// Code from libigl, reduces computation time by preallocation
	L.resize(nv,nv);
	L.reserve(10*nv);
	Eigen::Matrix<int, -1, 2> edges;
	edges.resize(3,2);
	edges <<
		1,2,
		2,0,
		0,1;

	// compute and store (cot(a_ij)+cot(b_ij))/2 for each edge ij
	//std::cout << "Computing cot ij" << std::endl;
	for (int i = 0; i < nf; ++i)
	{
		int v[3];
		for (int j = 0; j < 3; ++j){
			v[j] = F(i,j);
		}
		std::vector<Eigen::MatrixXd> vertices, sides;
		for (int j = 0; j < 3; ++j)
		{
			Eigen::Vector3d vertex(V(v[j],0), V(v[j],1), V(v[j],2));
			vertices.push_back(vertex);
		}
		sides.push_back(vertices[1] - vertices[2]);
		sides.push_back(vertices[2] - vertices[0]);
		sides.push_back(vertices[0] - vertices[1]);
		for (int j = 0; j < 3; ++j)
		{
			int source = edges(j,0);
			int dest = edges(j,1);
			double cot = 0.5/tan(angle_between_vectors(sides[source], sides[dest]));
			Ltemp.push_back(Eigen::Triplet<double>(v[source], v[dest], cot));
			Ltemp.push_back(Eigen::Triplet<double>(v[dest], v[source], cot));
			Ltemp.push_back(Eigen::Triplet<double>(v[source], v[source], -cot));
			Ltemp.push_back(Eigen::Triplet<double>(v[dest], v[dest], -cot));
		}
	}
	L.setFromTriplets(Ltemp.begin(), Ltemp.end());
}

void compute_mean_matrix(
	Eigen::MatrixXd const & V,
	Eigen::MatrixXi const & F,
	Eigen::SparseMatrix<double> & L
)
{
	// V : input vertices, N-by-3
	// F : input faces
	//
	// L : output mean matrix, N-by-N.
	int nv = V.rows();
	int nf = F.rows();
	std::vector<std::vector<int> > A;
	igl::adjacency_list(F,A);
	
	std::vector<Eigen::Triplet<double>> Ltemp;

	// Code from libigl, reduces computation time by preallocation
	L.resize(nv,nv);
	L.reserve(10*nv);
	Eigen::Matrix<int, -1, 2> edges;
	edges.resize(3,2);
	edges <<
		1,2,
		2,0,
		0,1;
	
	// compute and store (tan(a_ij/2)+tan(b_ij/2))/||v_i-v_j|| for each edge ij
	//std::cout << "Computing cot ij" << std::endl;
	for (int i = 0; i < nf; ++i)
	{
		int v[3];
		for (int j = 0; j < 3; ++j){
			v[j] = F(i,j);
		}
		std::vector<Eigen::Vector3d> vertices, sides;
		for (int j = 0; j < 3; ++j)
		{
			Eigen::Vector3d vertex(V(v[j],0), V(v[j],1), V(v[j],2));
			vertices.push_back(vertex);
		}
		sides.push_back(vertices[1] - vertices[2]);
		sides.push_back(vertices[2] - vertices[0]);
		sides.push_back(vertices[0] - vertices[1]);
		for (int j = 0; j < 3; ++j)
		{
			int v1 = edges(j,0);
			int v2 = edges(j,1);
			double mean = tan(angle_between_vectors(sides[v1], sides[v2])/2);
			double d1,d2;
			d1 = sqrt(sides[v1].dot(sides[v1]));
			d2 = sqrt(sides[v2].dot(sides[v2]));
			Ltemp.push_back(Eigen::Triplet<double>(v[j], v[v2], mean/d1));
			Ltemp.push_back(Eigen::Triplet<double>(v[j], v[v1], mean/d2));
			Ltemp.push_back(Eigen::Triplet<double>(v[j], v[j], -mean/d1-mean/d2));
		}
	}
	L.setFromTriplets(Ltemp.begin(), Ltemp.end());
}

double area_of_triangle(
	Eigen::MatrixXd const & vertex_1,
	Eigen::MatrixXd const & vertex_2,
	Eigen::MatrixXd const & vertex_3
	)
{
	// vertex_i : 1-by-3 matrix, a row of the vertices matrix
	//
	// return : double, area of triangle with vertices at v1,2,3
	double side[3];
	double s;
	side[0] = (vertex_2 - vertex_1).norm();
	side[1] = (vertex_3 - vertex_2).norm();
	side[2] = (vertex_1 - vertex_3).norm();
	s = 0.5*(side[0] + side[1] + side[2]);
	return sqrt(s*(s-side[0])*(s-side[1])*(s-side[2]));
}

void compute_area(
	Eigen::MatrixXd const & V,
	Eigen::MatrixXi const & F,
	Eigen::SparseMatrix<double> & Area
	)
{
	// V : input vertices, N-by-3
	// F : input faces
	//
	// Area : output diagonal matrix, N-by-N, each entry represents the normalising area around a vertex
	int faces = F.rows();
	int vertices = V.rows();
	double area_tri[faces];
	std::vector<std::vector<int>> neighbourFacesOf(vertices);
	for (int i = 0; i < faces; ++i)
	{
		neighbourFacesOf[F(i,0)].push_back(i);
		neighbourFacesOf[F(i,1)].push_back(i);
		neighbourFacesOf[F(i,2)].push_back(i);
		area_tri[i] = area_of_triangle(V.row(F(i,0)), V.row(F(i,1)), V.row(F(i,2)))/3;
	}

	Area.resize(vertices, vertices);
	Area.reserve(vertices);
	std::vector<Eigen::Triplet<double>> temp;
	for (int i = 0; i < vertices; ++i)
	{
		double area = 0;
		for (int index:neighbourFacesOf[i])
		{
			area += area_tri[index];
		}
		temp.push_back(Eigen::Triplet<double>(i,i, area));
	}
	Area.setFromTriplets(temp.begin(), temp.end());
}

void compute_H(
	Eigen::MatrixXd const & V,
	Eigen::MatrixXi const & F,
	Eigen::VectorXd & H
	)
{
	// You job is to use igl::cotmatrix, igl::massmatrix, igl::invert_diag to compute mean curvature H at each vertex
	// 
	// V : input vertices, N-by-3
	// F : input faces
	//
	// H : output vertices, K-by-1
	// 
	//
	// Hints
	// Compute Laplace-Beltrami operator
	//	Eigen::SparseMatrix<double> L, Area, AreaInv;
	//	igl::cotmatrix(V, F, L);
	//	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, Area );
	//	

	Eigen::SparseMatrix<double> L, Area, AreaInv;
	compute_area(V, F, Area);

	AreaInv.resize(Area.rows(), Area.cols());
	AreaInv.reserve(Area.rows());
	std::vector<Eigen::Triplet<double>> temp;
	for (int i = 0; i < Area.rows(); ++i)
	{
		temp.push_back(Eigen::Triplet<double>(i,i, double(1)/Area.coeff(i,i)));
	}
	AreaInv.setFromTriplets(temp.begin(), temp.end());

	//Eigen::MatrixXd LBV = AreaInv*L*V;
	Eigen::SparseMatrix<double> UL;
	compute_uniform_matrix(V,F,UL);

	compute_H(V, UL, H);
}


void compute_H(Eigen::MatrixXd const & V, Eigen::SparseMatrix<double> const & L,
				Eigen::VectorXd & H)
{
	// Input
	// V 	vertices
	// L 	Laplace operator
	//
	// Ouput
	// H 	Mean curvature
	Eigen::MatrixXd DV = L*V;

	H.resize(V.rows());

	for (int i = 0; i < V.rows(); ++i)
	{
		H(i) = 0.5*DV.row(i).norm();
	}
}

void compute_gaussian_curvature(Eigen::MatrixXd const & V, Eigen::MatrixXi const & F,
							Eigen::VectorXd & K)
{
	// Input
	// V 	vertices
	// F 	faces
	//
	// Ouput
	// K 	Gaussian curvature
	int nv = V.rows();
	int nf = F.rows();
	Eigen::Matrix<int, 3, 2> edges;
	edges.resize(3,2);
	edges <<
		1,2,
		2,0,
		0,1;

	K.resize(nv);
	K.setConstant(2*M_PI);

	for (int i = 0; i < nf; ++i)
	{
		int v[3];
		for (int j = 0; j < 3; ++j){
			v[j] = F(i,j);
		}
		std::vector<Eigen::MatrixXd> vertices, sides;
		for (int j = 0; j < 3; ++j)
		{
			Eigen::Vector3d vertex(V(v[j],0), V(v[j],1), V(v[j],2));
			vertices.push_back(vertex);
		}
		sides.push_back(vertices[1] - vertices[2]);
		sides.push_back(vertices[2] - vertices[0]);
		sides.push_back(vertices[0] - vertices[1]);

		for (int vertex = 0; vertex < 3; ++vertex)
		{
			double angle =angle_between_vectors(sides[edges(vertex,0)], sides[edges(vertex,1)]);
			K(v[vertex]) -= angle;
		}
	}
	
}


void compute_eigendecomposition(Eigen::MatrixXd const & V, Eigen::MatrixXi const & F, int const n, 
									Eigen::MatrixXd & eigen_mesh, Eigen::SparseMatrix<double> &Area)
{
	//
	// Input
	// V 	Nx3 matrix, verticese
	// F    faces
	// n	Number of eigenvectors needed
	//
	// Output
	// eigen_mesh  Nxn matrix, stores n eigenvectors of LB with largest eigenvalue
	Eigen::SparseMatrix<double> C;
	compute_cotan_matrix(V,F,C);
	compute_area(V,F,Area);

	C = -C;
	Eigen::VectorXd eigvalues;
	igl::eigs(C, Area, n, igl::EIGS_TYPE_SM, eigen_mesh, eigvalues);
}

void tutte_flattening_uniform(Eigen::MatrixXd const & V, Eigen::MatrixXi const & F, Eigen::MatrixXd & out_V)
{
	//
	// input:
	//   V: vertices
	//   F: face 
	// output:
	//   out_V
	// 

	// Check if it contains boundary
	Eigen::VectorXi b_bex_index;
	igl::boundary_loop(F, b_bex_index);
	if (b_bex_index.rows() == 0)
	{
		// If not, use the first face as boundary
		b_bex_index.resize(3);
		b_bex_index(0) = F(0,0);
		b_bex_index(1) = F(0,1);
		b_bex_index(2) = F(0,2);
	}

	// Set polygon boundary
	int n = b_bex_index.rows();

	Eigen::MatrixXd b_bex;
	b_bex.resize(n, 3);
	for (int i = 0; i < n; ++i)
	{
		double theta = 2*M_PI*i/n;
		b_bex(i,0) = 3*cos(theta);
		b_bex(i,1) = 3*sin(theta);
		b_bex(i,2) = 0;
	}

	// Compute Laplacian 
	Eigen::SparseMatrix<double> L;
	compute_uniform_matrix(V, F, L);

	// Compute mapped vertices
	// List of all vertex indices
    Eigen::VectorXi all;
    Eigen::VectorXi interior_i;

    // I = low:step:hi
    // igl::colon(low,step,hi,I);
    igl::colon<int>(0, V.rows() - 1, all);

    // get the indices of interior vertices
    // interior_i = all_i - boundary_i
    Eigen::VectorXi IA;
    igl::setdiff(all, b_bex_index, interior_i, IA);

    Eigen::SparseMatrix<double>  L_ii, L_ib;
    
    igl::slice(L, interior_i, interior_i, L_ii);
    igl::slice(L, interior_i, b_bex_index, L_ib);

    Eigen::MatrixXd new_V(V.rows(), V.cols());
    
	// set your A and b
    Eigen::MatrixXd b_xy;
	Eigen::MatrixXd A;
	
    b_xy = -1 * L_ib * b_bex;

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(L_ii);
    solver.factorize(L_ii);
    if (solver.info() != Eigen::Success) { std::cout << "decomposition failed\n";exit(1); }

	// ... 
    for (int i = 0; i < n; ++i)
    {
        new_V.row(b_bex_index[i]) = b_bex.row(i);
    }
	new_V.col(2).setZero();

    // now solve for each dimension (XYZ)
    for (int i = 0; i < 2; i++)
    {
        std::cout << "start solving .." << i << "\n";
        std::clock_t start = std::clock();
        
        // solve Ax=b 
        // ...
        Eigen::VectorXd solution = solver.solve(b_xy.col(i));
        
        if (solver.info() != Eigen::Success) {std::cout << "solving failed\n";exit(1);}

        // restore vertices order   
        // ...
        for (int j = 0; j < interior_i.rows(); ++j)
        {
            new_V(interior_i[j], i) = solution[j];
        }
        
        double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        std::cout << "time cost="  << duration << '\n';
    }
    out_V = new_V;
}

void tutte_flattening_LB(Eigen::MatrixXd const & V, Eigen::MatrixXi const & F, Eigen::MatrixXd & out_V)
{
	//
	// input:
	//   V: vertices
	//   F: face 
	// output:
	//   out_V
	// 

	// Check if it contains boundary
	Eigen::VectorXi b_bex_index;
	igl::boundary_loop(F, b_bex_index);
	if (b_bex_index.rows() == 0)
	{
		// If not, use the first face as boundary
		b_bex_index.resize(3);
		b_bex_index(0) = F(0,0);
		b_bex_index(1) = F(0,1);
		b_bex_index(2) = F(0,2);
	}

	// Set polygon boundary
	int n = b_bex_index.rows();

	Eigen::MatrixXd b_bex;
	b_bex.resize(n, 3);
	for (int i = 0; i < n; ++i)
	{
		double theta = 2*M_PI*i/n;
		b_bex(i,0) = 3*cos(theta);
		b_bex(i,1) = 3*sin(theta);
		b_bex(i,2) = 0;
	}

	// Compute Laplacian 
	Eigen::SparseMatrix<double> C;
	compute_cotan_matrix(V,F,C);

	int negative= 0;
	for (int k = 0; k < C.outerSize(); ++k)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(C, k); it; ++it)
		{
			if (it.value() < 0)
				negative++;
		}
	}
	std::cout << "Num of negative entries is " << negative << std::endl; 

	// Compute mapped vertices
	// List of all vertex indices
    Eigen::VectorXi all;
    Eigen::VectorXi interior_i;

    igl::colon<int>(0, V.rows() - 1, all);

    // get the indices of interior vertices
    Eigen::VectorXi IA;
    igl::setdiff(all, b_bex_index, interior_i, IA);
    
    Eigen::SparseMatrix<double>  C_ii, C_ib;
    
    igl::slice(C, interior_i, interior_i, C_ii);
    igl::slice(C, interior_i, b_bex_index, C_ib);
    
    Eigen::MatrixXd new_V(V.rows(), V.cols());
    
	// set your A and b
    Eigen::MatrixXd b_xy;
	Eigen::MatrixXd A;
	
    b_xy = -1 * C_ib * b_bex;

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(C_ii);
    solver.factorize(C_ii);
    if (solver.info() != Eigen::Success) { std::cout << "decomposition failed\n";exit(1); }

	// ... 
    for (int i = 0; i < n; ++i)
    {
        new_V.row(b_bex_index[i]) = b_bex.row(i);
    }
	new_V.col(2).setZero();

    // now solve for each dimension (XYZ)
    for (int i = 0; i < 2; i++)
    {
        std::cout << "start solving .." << i << "\n";
        std::clock_t start = std::clock();
        
        // solve Ax=b 
        // ...
        Eigen::VectorXd solution = solver.solve(b_xy.col(i));
        
        if (solver.info() != Eigen::Success) {std::cout << "solving failed\n";exit(1);}

        // restore vertices order   
        // ...
        for (int j = 0; j < interior_i.rows(); ++j)
        {
            new_V(interior_i[j], i) = solution[j];
        }
        
        double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        std::cout << "time cost="  << duration << '\n';
    }
    out_V = new_V;
}

void tutte_flattening_mean(Eigen::MatrixXd const & V, Eigen::MatrixXi const & F, Eigen::MatrixXd & out_V)
{
	//
	// input:
	//   V: vertices
	//   F: face 
	// output:
	//   out_V
	// 

	// Check if it contains boundary
	Eigen::VectorXi b_bex_index;
	igl::boundary_loop(F, b_bex_index);
	if (b_bex_index.rows() == 0)
	{
		// If not, use the first face as boundary
		b_bex_index.resize(3);
		b_bex_index(0) = F(0,0);
		b_bex_index(1) = F(0,1);
		b_bex_index(2) = F(0,2);
	}

	// Set polygon boundary
	int n = b_bex_index.rows();

	Eigen::MatrixXd b_bex;
	b_bex.resize(n, 3);
	for (int i = 0; i < n; ++i)
	{
		double theta = 2*M_PI*i/n;
		b_bex(i,0) = 3*cos(theta);
		b_bex(i,1) = 3*sin(theta);
		b_bex(i,2) = 0;
	}

	// Compute Laplacian 
	Eigen::SparseMatrix<double> C;
	compute_mean_matrix(V,F,C);

	// Compute mapped vertices
	// List of all vertex indices
    Eigen::VectorXi all;
    Eigen::VectorXi interior_i;

    // I = low:step:hi
    // igl::colon(low,step,hi,I);
    igl::colon<int>(0, V.rows() - 1, all);

    // get the indices of interior vertices
    // interior_i = all_i - boundary_i
    Eigen::VectorXi IA;
    igl::setdiff(all, b_bex_index, interior_i, IA);
    
    Eigen::SparseMatrix<double>  C_ii, C_ib;
    
    igl::slice(C, interior_i, interior_i, C_ii);
    igl::slice(C, interior_i, b_bex_index, C_ib);

    Eigen::MatrixXd new_V(V.rows(), V.cols());
	new_V.col(2).setZero();
    
	// set your A and b
    Eigen::MatrixXd b_xy;
	Eigen::MatrixXd A;
	
    b_xy = -1 * C_ib * b_bex;

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(C_ii);
    solver.factorize(C_ii);
    if (solver.info() != Eigen::Success) { std::cout << "decomposition failed\n";exit(1); }

	// ... 
    for (int i = 0; i < n; ++i)
    {
        new_V.row(b_bex_index[i]) = b_bex.row(i);
    }
	new_V.col(2).setZero();

    // now solve for each dimension (XYZ)
    for (int i = 0; i < 2; i++)
    {
        std::cout << "start solving .." << i << "\n";
        std::clock_t start = std::clock();
        
        // solve Ax=b 
        // ...
        Eigen::VectorXd solution = solver.solve(b_xy.col(i));
        
        if (solver.info() != Eigen::Success) {std::cout << "solving failed\n";exit(1);}

        // restore vertices order   
        // ...
        for (int j = 0; j < interior_i.rows(); ++j)
        {
            new_V(interior_i[j], i) = solution[j];
        }
        
        double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        std::cout << "time cost="  << duration << '\n';
    }
    out_V = new_V;
}

void compute_face_angle(Eigen::MatrixXd const & V_original, Eigen::MatrixXi const & F, Eigen::MatrixXd & FA)
{
	int nf = F.rows();
	Eigen::Matrix<int, 3, 2> edges;
	edges.resize(3,2);
	edges <<
		1,2,
		2,0,
		0,1;
	FA.resize(nf, 3);
	FA.setZero();
	for (int i = 0; i < nf; ++i)
	{
		int v[3];
		for (int j = 0; j < 3; ++j){
			v[j] = F(i,j);
		}
		std::vector<Eigen::Vector3d> vertices, sides;
		for (int j = 0; j < 3; ++j)
		{
			Eigen::Vector3d vertex(V_original(v[j],0), V_original(v[j],1), V_original(v[j],2));
			vertices.push_back(vertex);
		}
		sides.push_back(vertices[1] - vertices[2]);
		sides.push_back(vertices[2] - vertices[0]);
		sides.push_back(vertices[0] - vertices[1]);

		for (int vertex = 0; vertex < 3; ++vertex)
		{
			FA(i,vertex) = angle_between_vectors(sides[edges(vertex,0)], sides[edges(vertex,1)]);
			if (std::isnan(FA(i,vertex)))
			{
				std::cout << "(" << i << "," << vertex << ") -- " <<  sides[edges(vertex,0)].transpose() << "  |||   " << sides[edges(vertex,1)].transpose() << std::endl;
				std::cout << sides[edges(vertex,0)].dot(-sides[edges(vertex,1)]) << "   " << sides[edges(vertex,0)].norm() << "   " << sides[edges(vertex,1)].norm() << std::endl;
				std::cout << sides[edges(vertex,0)].dot(-sides[edges(vertex,1)])/sides[edges(vertex,0)].norm()/sides[edges(vertex,1)].norm() << std::endl;
				std::cout << acos(-1) << "   "<< acos(1) << std::endl;
			}
		}
	}
}

double compute_distortion(Eigen::MatrixXd const & FA_original, Eigen::MatrixXd const & FA_cmpr)
{
	double sum = 0;
	double PI2 = 2*M_PI;
	Eigen::MatrixXd diff = FA_original - FA_cmpr;
	for (int i = 0; i < FA_original.rows(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			while (std::abs(diff(i,j)) > M_PI)
				if (diff(i,j) > M_PI)
					diff(i,j) -= PI2;
				else
					diff(i,j) += PI2;
			sum += std::abs(diff(i,j));
		}
	}
	return sum;
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
