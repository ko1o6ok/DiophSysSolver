//
// Created by ф on 01.04.2023.
//

#ifndef FASTHNF_SIMPLEX_H
#define FASTHNF_SIMPLEX_H

#include <vector>
#include <list>
#include "Matrix.h"
using namespace std;

class Vertex{
public:
    const int dimension = 0;
    unsigned long index;
    explicit Vertex(unsigned long ind);
};

class Simplex {
public:
    unsigned long dimension;
    //vector<Simplex> faces;
    list<Vertex> vertices; // Вершины должны составлять упорядоченный по индексам список
    vector<unsigned long> make_word(); // Вернуть упорядоченную по возрастанию последовательность номеров вершин
    Simplex();
    ~Simplex();
    explicit Simplex(vector<unsigned long> v); // По вершинам
    Simplex(Simplex const &s); // Конструктор копирования
    void add_vertex(unsigned long ind); // Добавить новую вершину
    friend ostream& operator<<(ostream& stream,const Simplex& s); // Вывод

    vector<vector<unsigned long>> get_border_words(); // Слова, соответствующие оператору границы ; - + - + ...
};

struct MyNode{
    vector<MyNode*> children;
    unsigned long last_index;
};
MyNode* add_child(MyNode* node, unsigned long last_index);
MyNode* make_a_node(unsigned long last_index);
// Симплексное дерево
class SimplexTree{
public:
    MyNode* root; // Корень
    vector<vector<double>> point_cloud; // Порождающее множество точек
    unsigned long max_dimension; // Максимальная размерность входящего симплекса
    unsigned long num_vertices; // Число вершин
    void insert_simplex(Simplex s);
    void insert_simplex(vector<unsigned long> v);
    vector<Simplex> all_simplexes_of_dim(int k) const; // Строит по дереву все симплексы размерности k
    void construct_from_point_cloud(unsigned long max_dimension,double eps); // Строит симплекс-дерево на основе графа ближайших соседей
    Matrix<long int> border_operator_matrix(int dimension) const; // Матрица оператора границы
    vector<int> betti_numbers() const; // Выписать числа Бэтти данного комплекса

    // Базовый конструктор
    SimplexTree();
    explicit SimplexTree(vector<vector<double>>& pnt_cld); // Просто множество симплексов размерности 0 на основе множества точек
    void print() const;
    // Базовый деструктор
    ~SimplexTree();
};
void print_tree(const string& prefix,MyNode* rt);// Выписать дерево
void simplexes_of_dim(MyNode* current_node,int k,vector<Simplex>& simplexes,const Simplex& current_simplex);
#endif //FASTHNF_SIMPLEX_H
