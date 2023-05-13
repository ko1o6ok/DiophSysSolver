#ifndef FASTHNF_GRAPH_H
#define FASTHNF_GRAPH_H
#include <vector>
#include <list>
using namespace std;
class Graph {
public:
    list<unsigned long> vertices; // Упорядоченные по возрастанию номера вершин
    unsigned long num_vertices; // Число вершин
    vector<vector<short>> adj_matrix; // Матрица смежности. ВАЖНО: используется только верхняя половина. m <= n
    vector<vector<double>> point_cloud; // Его множество точек

    Graph();// Конструктор по умолчанию
    Graph(Graph const &another_g); // Конструктор копирования
    void connect_vertices(unsigned long vert1,unsigned long vert2); // Соединить две вершины ребром
    explicit Graph(const vector<vector<double>>& pnt_cld); // Создание ПУСТОГО графа из облака точек
    Graph connect_eps_neighbours(double eps); // Соединить только вершины с ЕВКЛИДОВЫМ расстоянием < eps
    vector<unsigned long> lower_neighbours(unsigned long vertex,Graph difference); // Все вершины, предшествующие данной
    ~Graph(); // Деструктор
    bool is_empty();
};
double euclidean_distance(vector<double> point_1, vector<double> point_2); // Евклидово расстояние между двумя точками
#endif //FASTHNF_GRAPH_H
