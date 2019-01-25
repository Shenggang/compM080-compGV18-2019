# COMPM080/GV18 Tutorial 2 - C++ 

In this tutorial, we introduce how to retrieve neighboring information via using `igl::vertex_triangle_adjacency` and demonstrate two approaches to calculate vertex normal. 1) using adjacent face normal. 2) using FLANN to retrieve K-nearest neighbors for each vertex and fit a plane on it.

---
# FAQ 

* ..

---
# Access neighboring vertices 

## Adjacent list (today's example)
````
    // a list of list of int. storage the index of adjacent face and the corresponding local index
    adjacent_list, vertex_local_index = build_adjacent_list
    
    v_i = V[i] // a vertex  
    neighboring_faces_indices = adjacent_list[v_i]
    local_vertex_indices = vertex_local_index[v_i]
    
    //loop over all neighbors
    for i = 1..N :
    	face_index = neighboring_faces_indices[i]
    	v_local_index = local_vertex_indices[i]
	
	for k in 1...3 :
	    if v_local_index != k : 
	        V[k] // a neighboring vertex
    
````
## Halfedge style neighboring access  
````
    v_0 = V[i] // a vertex 
    halfedge_0 = v_0.outgoing_hfedge()
    halfedge_i = halfedge_0
    
    while(1):
        k = halfedge_i.vertex_index
	V[k] // neighboring vertex
        halfedge_i = halfedge_i.opposite_hfedge()
	halfedge_i = halfedge_i.next_hfedge()
	
	if halfedge_i == halfedge_0:
	    break
````
---
# Useful references 
* [libigl's vertex_triangle_adjacency](https://github.com/libigl/libigl/blob/master/include/igl/vertex_triangle_adjacency.cpp)
* [OpenMesh's Halfedge Introduction](https://www.openmesh.org/media/Documentations/OpenMesh-6.3-Documentation/a00010.html)
* [CGAL 4.13 - Halfedge Data Structures](https://doc.cgal.org/latest/HalfedgeDS/index.html)
* [libigl's Fake Halfedge](https://github.com/libigl/libigl/blob/master/include/igl/HalfEdgeIterator.h)
* [Eigen quick references](http://eigen.tuxfamily.org/dox/group__QuickRefPage.html)
* [Eigen for matlab users](http://igl.ethz.ch/projects/libigl/matlab-to-eigen.html)
* [Libigl tutorials](https://libigl.github.io/tutorial/)
* [C++ template introduction](http://www.cplusplus.com/doc/oldtutorial/templates/)


---

# Install 
## STEP-1. get the code of tutorial 2:

### option (a) For those who are familiar with Git. 
1. cd to your previous `compM080-compGV18-2019` directory  
`cd ~/my_code/compM080-compGV18-2019`
2. commit all your changed you added in the previous tutorial   
`git add .`  
`git commit -m 'added all my changed'   
`git pull origin master'

### option (b) just download/clone this repository again.  
`git clone https://github.com/smartgeometry-ucl/compM080-compGV18-2019.git`

## STEP-2. get libigl
In this tutorial we add another libigl directory at the top level `compM080-compGV18-2019-private\libigl` in order to make library sharing become easier. Due to this you will need to clone libigl again. 

### option (a) use submodule update:  
1. at this repository's directory `compM080-compGV18-2019`
3. do `git submodule update --init --recursive`
4. check libigl directory in `compM080-compGV18-2019/libigl` is not empty  

### option (b) download/clone libigl at `compM080-compGV18-2019/libigl`  


---
# Practice 1
In this practice, you will learn how to access neighboring vertics and use it to calculate vertex normals. 
**Your job** is to implement the function `calculate_vertex_normal` at `tutorial_002\cpp\practices_access_neighbors\src\mytools.h`
## access neighboring vertices by using igl::vertex_triangle_adjacency

![p1](/tutorial_002/cpp/docimgs/adj_vex.JPG )

````
// build adjacent matrix 
// Inputs:
//  number of vertices
//  face def.
//
// Outputs:
//   VF   #V list of lists of incident faces (adjacency list) 
//   VFi  #V list of lists    of  index of incidence     within incident faces listed
//        vector<vector<int>>        VI                    triangle_face 

   std::vector<std::vector<int> > VF;
   std::vector<std::vector<int> > VFi;
   igl::vertex_triangle_adjacency(m_V.rows(), m_F, VF, VFi);

    //-----------------------------------------------------------
    // visualize adjacent vertices of 'm_VF[sel_vidx]' at here

    Eigen::MatrixXd adj_vex;
    adj_vex.resize(m_VF[sel_vidx].size() * 2, 3);

    int count = 0;
    // get adjacent faces & vertices
    for (size_t i = 0; i < m_VF[sel_vidx].size(); i++)
    {
      int face_idx = m_VF[sel_vidx][i];
      int v_local_idx = m_VFi[sel_vidx][i];

      for (int iv = 0; iv < 3; iv++)
      {
        if (iv != v_local_idx)
        {
          Eigen::RowVector3d const vex = m_V.row(m_F(face_idx, iv));
          adj_vex.row(count) = vex;
          count++;
        }
      }
    } 

````

## vertex normal reference image:

Hints:  
1. first calculate face normal (`igl::per_face_normals`)
2. access each vertex's neighboring faces
3. average face normal to get a vertex normal 

![p1](/tutorial_002/cpp/docimgs/vn2.JPG )


---
# Practice 2
In practice 2, we show an example to retrieve K-nearest neighbor of an input vertex.  
[nanoflann: https://github.com/jlblancoc/nanoflann](https://github.com/jlblancoc/nanoflann)   

**Your job** is to implement the function `calculate_vertex_normal_flann` at `tutorial_002\cpp\calculate_vertex_normal_flann\src\mytools.h`

## example

![](/tutorial_002/cpp/docimgs/knn.JPG )

````
// FLANN demo
			
typedef nanoflann::KDTreeEigenMatrixAdaptor< Eigen::MatrixXd >  my_kd_tree_t;

// build kd-tree uses vertices m_V
my_kd_tree_t mat_index( m_V, 50 /* max leaf */);
mat_index.index->buildIndex();

// select a random vertex
std::uniform_int_distribution<int> distribution(0, m_V.rows());
int rd_vidx = distribution(generator);
Eigen::RowVector3d rd_pts = m_V.row(rd_vidx);

// set K nearest samples
const size_t par_K = 20;

// create results objects
std::vector<size_t> indexes(par_K);
std::vector<double> dists_sqr(par_K);

// bind results
nanoflann::KNNResultSet<double> res(par_K);
res.init(indexes.data(), dists_sqr.data());

// [find KNN]
//   SearchParams Note: The first argument (checks_IGNORED_) is ignored, but kept for compatibility
//   checks_IGNORED_ was used to specify maximum # of leaf to be visited 
mat_index.index->findNeighbors(res, rd_pts.data(), nanoflann::SearchParams(50));

// check result
Eigen::MatrixXd nn_vex(indexes.size(),3);
for (size_t i=0 ; i < indexes.size() ;i++)
{ 
	nn_vex.row(i) = m_V.row(indexes[i]);
}

````
Hints:  
1. follow the above example, construct kdTree structure
2. access each vertex to retrieve its k-nearest neighbors
3. use `igl:fit_plane` to fit a plane
4. flip normal vector if needed (to get outer-pointing normal)

## reference images
### average face normal  
![](/tutorial_002/cpp/docimgs/vex_nv.JPG )

### fit a plane on KNN
![](/tutorial_002/cpp/docimgs/flann_nv.JPG )

