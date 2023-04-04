
#include "Graph.h"

#include <utility>
#include <cmath>
#include <iostream>
Graph::~Graph() {
    vertices.clear();
    adj_matrix.clear();
}

Graph::Graph(list<unsigned long> verts) {
    num_vertices = verts.size();
    vertices = std::move(verts);
    for (int i = 0; i < num_vertices; ++i) {
        vector<short> app;
        app.reserve(num_vertices);
        for (int j = 0; j < num_vertices; ++j)
            app.push_back(0);
        adj_matrix.push_back(app);
    }
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
void Graph::connect_eps_neighbours(double eps) {

    // Проходим по облаку точек
    for (int i = 0; i < num_vertices; ++i)
        for (int j = i+1; j < num_vertices; ++j) {
            auto dist = euclidean_distance(point_cloud[i],point_cloud[j]);
            //cout << "I've computed the distance "<< dist << endl;
            if((adj_matrix[i][j] == 0)&&(dist <= eps)){
                connect_vertices(i,j);
                //cout << "Connecting vertices "<< i<<" and "<< j << endl;
            }

        }
}

void Graph::disconnect_vertices(unsigned long vert1, unsigned long vert2) {
    // Assume vert1 <= vert2
    if(vert1 > vert2)
        std::swap(vert1,vert2);
    adj_matrix[vert1][vert2] = 0; // Используется только верхняя половина матрицы
}

void Graph::print_adj_matrix() {
    for (int i = 0; i < num_vertices; ++i){
        for (int j = 0; j < num_vertices; ++j)
            cout<< adj_matrix[i][j] << ", ";
        cout << endl;
    }
}
