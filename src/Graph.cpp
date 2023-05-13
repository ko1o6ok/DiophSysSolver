#include "Graph.h"
#include <utility>
#include <cmath>
Graph::~Graph() {
    vertices.clear();
    adj_matrix.clear();
}
Graph::Graph() {
    num_vertices = 0;
}
void Graph::connect_vertices(unsigned long vert1, unsigned long vert2) {
    // Assume vert1 <= vert2
    if(vert1 > vert2)
        std::swap(vert1,vert2);
    adj_matrix[vert1][vert2] = 1; // Используется только верхняя половина матрицы
}
Graph::Graph(const vector<vector<double>>& pnt_cld) {
    num_vertices = pnt_cld.size();
    point_cloud = pnt_cld;
    // Заполним матрицу смежности
    for (int i = 0; i < num_vertices; ++i) {
        vector<short> app;
        app.reserve(num_vertices);
        for (int j = 0; j < num_vertices; ++j)
            app.push_back(0);
        adj_matrix.push_back(app);
    }
}
double euclidean_distance(vector<double> point_1, vector<double> point_2){
    double S = 0.0;
    for (int i = 0; i < point_1.size(); ++i) {
        auto t = point_1[i]-point_2[i];
        S += t*t;
    }
    return sqrt(S);
}
Graph Graph::connect_eps_neighbours(double eps) {
    Graph diff_gr(point_cloud);
    // Проходим по облаку точек
    for (int i = 0; i < num_vertices; ++i){
        vertices.push_back(i);
        diff_gr.vertices.push_back(i);
        for (int j = i+1; j < num_vertices; ++j) {
            auto dist = euclidean_distance(point_cloud[i],point_cloud[j]);
            if((adj_matrix[i][j] == 0)&&(dist <= eps)){
                connect_vertices(i,j);
                diff_gr.connect_vertices(i,j);
            }

        }
    }
    return diff_gr;
}
vector<unsigned long> Graph::lower_neighbours(unsigned long vertex,Graph difference) {
    vector<unsigned long> res;
    for(auto& v:vertices)
        if((v < vertex)&&(adj_matrix[v][vertex] != 0)&(difference.adj_matrix[v][vertex]!=0))
            res.push_back(v);
    return res;
}
Graph::Graph(const Graph &another_g) {
    vertices = another_g.vertices;
    num_vertices = another_g.num_vertices;
    adj_matrix = another_g.adj_matrix;
    point_cloud = another_g.point_cloud;
}
bool Graph::is_empty() {
    for (int i = 0; i < num_vertices; ++i)
        for (int j = i+1; j < num_vertices; ++j)
            if(adj_matrix[i][j]!=0)
                return false;
    return true;
}
