#include "seamster.h"

EdgeTopology::EdgeTopology(MatrixXd &V,MatrixXi &F)
{
    this->V = V;
    this->F = F;
    igl::edge_topology(this->V, this->F, this->EV, this->FE, this->EF);  
}

EdgeTopology::EdgeTopology(MatrixXd V,MatrixXi F,MatrixXi EV,MatrixXi FE,MatrixXi EF)
{
    this->V = V;
    this->F = F;
    this->EV = EV;
    this->FE = FE;
    this->EF = EF;
}

void EdgeTopology::computeDirectedEdgeCorrespondence()
{
    if (isDirected)
        return;
    for (int face = 0; face < FE.rows(); ++face)
    {
        for (int edge = 0; edge < 3; ++edge)
        {
            int v1 = EV(FE(face, edge), 0);
            FE(face, edge)++;
            if (v1 != F(face, edge))
                FE(face, edge) *= -1;
        }
    }
    isDirected = true;
}

template <typename T>
bool isAinVec(T A, vector<T> vec)
{
    for (int i = 0; i < vec.size(); ++i)
    {
        if (A == vec[i])
        {
            return true;
        }
    }
    return false;
}

int findEdge(MatrixXi const & EV, int const source, int const end)
{
    for (int i = 0; i < EV.rows(); ++i)
    {
        if (((EV(i,0) == source) && (EV(i,1) == end)) ||
        ((EV(i,0) == end) && (EV(i,1) == source)))
        {
            return i;
        }
    }
    return -1;
}

double Dijkstra(vector<vector<int>> const & adj, vector<vector<double>> const & adj_dist, int const source, vector<int> const & target, vector<int> & path)
{
    vector<int> unvisited, prev;
    vector<double> dist;
    for (int i = 0; i < adj.size(); ++i)
    {
        dist.push_back(999);
        prev.push_back(-1);
        unvisited.push_back(i);
    }
    dist[source] = 0;
    double min_dist = 99999;
    int c_tar = 0;
    while (unvisited.size() > 0)
    {
        int min = 9999;
        int u = -1;
        int index = 0;
        for (int i = 0; i < unvisited.size(); ++i)
        {
            if (dist[unvisited[i]] < min)
            {
                u = unvisited[i];
                min = dist[unvisited[i]];
                index = i;
            }
        }
        unvisited[index] = unvisited[unvisited.size()-1];
        unvisited.pop_back();
        for (int i = 0; i < adj[u].size(); ++i)
        {
            int v = adj[u][i];
            double alt = dist[u] + adj_dist[u][i];
            if (alt < dist[v])
            {
                dist[v] = alt;
                prev[v] = u;
            }
        }
    }
    for (int i = 0; i < target.size(); ++i)
    {
        int t = target[i];
        if (dist[t] < min_dist)
        {
            min_dist = dist[t];
            c_tar = t;
        } 
    }
    path.clear(); 
    //cout << "path = ";
    while (prev[c_tar] >=0 && c_tar != source)
    {
    //    cout << c_tar <<" - ";
        path.push_back(c_tar);
        c_tar = prev[c_tar];
    }
    //cout << c_tar << endl;
    path.push_back(c_tar);
    return min_dist;
}

void save_path(MatrixXi const & EV, vector<int> const & path, vector<int> & out_SEi)
{
    for (int i = 1; i < path.size(); ++i)
    {
        out_SEi.push_back(findEdge(EV, path[i-1], path[i]));
    }
}

void compute_seam(MatrixXd const & V, MatrixXi const & F, EdgeTopology const & et, vector<int> const & Terminals, vector<int> & out_SEi)
{
    //
	// V : input vertices, N-by-3
	// F : input faces
    // Terminals : input, terminal vertices
	//
	// out_SEi : output egdes in seam
    // 

    // Initalise adjacency list and edge lists
    vector<vector<int>> adj;
    vector<vector<double>> adj_dist;
    igl::adjacency_list(F, adj);

    // Find boundary vertices
    Eigen::VectorXi b_bex_index, b_bex;
	igl::boundary_loop(F, b_bex_index);
    b_bex.resize(V.rows());
    b_bex.setZero();
    for (int i = 1; i < b_bex_index.rows(); ++i)
    {
        b_bex(b_bex_index(i)) = -100;
    }

    // Initialise used vertex store
    vector<vector<int>> current_T(0), front(0);
    vector<VectorXi> used_vertices;
    for (int i = 0; i < Terminals.size(); ++i)
    {
        vector<int> t,f;
        t.push_back(Terminals[i]);
        f.push_back(Terminals[i]);
        current_T.push_back(t);
        front.push_back(f);
        VectorXi vec;
        vec.resize(V.rows());
        vec.setZero();
        vec(Terminals[i]) = 1;
        used_vertices.push_back(vec+b_bex);
    }

    // Compute edge length
    for (int i = 0; i < adj.size(); ++i)
    {
        vector<int> cur_adj = adj[i];
        vector<double> dists;
        for (int j = 0; j < cur_adj.size(); ++j)
        {
            if (b_bex[i] + b_bex[cur_adj[j]] == -200)
            {
                dists.push_back(99999);
            } else
            {
                dists.push_back(sqrt((V.row(i)-V.row(cur_adj[j])).squaredNorm()));
            }
        }
        adj_dist.push_back(dists);
    }

    int n = current_T.size();
    while (n > 1)
    {
        // Compute cost for each terminal, avoid reusing vertex within one terminal set
        int min_S, min_T;
        min_S = min_T = 0;
        double min_dist = 999;
        for (int i = 0; i < n; ++i)
        {
            for (int t = 0; t < front[i].size(); ++t)
            {
                int cur_terminal = front[i][t];
                for (int j = 0; j < adj[cur_terminal].size(); ++j)
                {
                    int index = adj[cur_terminal][j];
                    if (used_vertices[i](index) == 0)
                    {
                        double dist = adj_dist[cur_terminal][j];
                        if (dist < min_dist)
                        {
                            min_dist = dist;
                            min_S = i;
                            min_T = index;
                        }
                    }
                }   
            }
        }
        // Update terminal, collapse terminal if meet
        front[min_S].push_back(min_T);
        used_vertices[min_S](min_T) = 1;
        for (int i = 0; i < n; ++i)
        {
            if (i == min_S) continue;
            if (used_vertices[i](min_T) >= 1)
            {
                // meeting at min_T, find shortest path from min_T to both terminal sets and store the pathes
                // fuse terminal i and min_S by poping terminal i and fuse used_vertices
                vector<int> path1,path2;
                Dijkstra(adj, adj_dist, min_T, current_T[min_S], path1);
                Dijkstra(adj, adj_dist, min_T, current_T[i], path2);
                save_path(et.EV, path1, out_SEi);
                save_path(et.EV, path2, out_SEi);
                // fuse seed sets
                for (int j = 0; j < current_T[i].size(); ++j)
                {
                    current_T[min_S].push_back(current_T[i][j]);    
                }
                current_T[i] = current_T[n-1];
                current_T.pop_back();
                // fuse front
                front[min_S].pop_back();
                for (int j = 0; j < front[i].size(); ++j)
                {
                    front[min_S].push_back(front[i][j]);
                }
                front[i] = front[n-1];
                front.pop_back();
                // fuse unusable vertices
                used_vertices[min_S] += used_vertices[i];
                used_vertices[i] = used_vertices[n-1];
                used_vertices.pop_back();
                n--;
            }
        }
    }
}

void cut_mesh(MatrixXd const & V, MatrixXi const & F, EdgeTopology const & et, vector<int> const & SEi, MatrixXd & out_V, MatrixXi & out_F) 
{
    MatrixXi cuts;
    cuts.resize(F.rows(), 3);
    cuts.setZero();
    for (int i = 0; i < SEi.size(); ++i)
    {
        int f1,f2;
        f1 = et.EF(SEi[i],0);
        f2 = et.EF(SEi[i],1);
        for (int j = 0; j < 3; ++j)
        {
            if (et.FE(f1,j) == SEi[i])
            {
                cuts(f1,j) = 1;
            }
            if (et.FE(f2,j) == SEi[i])
            {
                cuts(f2,j) = 1;
            }
        }
    }
    igl::cut_mesh(V,F,cuts,out_V,out_F);
}

void find_seam_terminals(MatrixXd const & V, MatrixXi const & F, vector<int> & Terminals)
{
    VectorXd K;
    compute_gaussian_curvature(V, F, K);
    K = K.cwiseAbs();
    struct pair
    {
        int index;
        double k;

        static int cmp(const void *a, const void *b)
        {
            struct pair *ia = (struct pair *)a;
            struct pair *ib = (struct pair *)b;
            return ia->k < ib->k ? 1 : -1;
        }
    };
    vector<pair> KI(0);
    for (int i =0; i < K.rows(); ++i)
    {
        pair p;
        p.index = i;
        p.k = K(i);
        KI.push_back(p);
    }
    qsort(KI.data(), KI.size(), sizeof KI[0], pair::cmp);
    double sum_D = K.sum();
    double D = 0;
    int index = 0;
    Eigen::VectorXi b_bex_index, b_bex;
	igl::boundary_loop(F, b_bex_index);
    b_bex.resize(V.rows());
    b_bex.setZero();
    for (int i = 0; i < b_bex_index.rows(); ++i)
    {
        b_bex(b_bex_index(i)) = 1;
    }
    while (D < 0.05*sum_D)
    {
        if (b_bex(index) != 1)
        {
            D += KI[index].k;
            Terminals.push_back(KI[index].index);
        }
        index++;
    }
}