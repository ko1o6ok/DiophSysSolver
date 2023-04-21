#ifndef FASTHNF_SIMPLEX_H
#define FASTHNF_SIMPLEX_H

#include <vector>
#include <list>
#include <chrono>

#include "Graph.h"
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
    Graph g; // Его граф
    unsigned long max_dimension; // Максимальная размерность входящего симплекса
    unsigned long num_vertices; // Число вершин
    void insert_simplex(Simplex s);
    void insert_simplex(vector<unsigned long> v);
    vector<Simplex> all_simplexes_of_dim(int k) const; // Строит по дереву все симплексы размерности k
    //Graph construct_from_point_cloud(unsigned long max_dimension,double eps); // Строит симплекс-дерево на основе графа ближайших соседей
    Matrix<long int> border_operator_matrix(int dimension,int& adds,bool& no_simplexes) const; // Матрица оператора границы
    Matrix<long int> border_operator_matrix(int dimension) const; // Матрица оператора границы
    vector<int> betti_numbers() const; // Выписать числа Бэтти данного комплекса

    // Базовый конструктор
    SimplexTree();
    explicit SimplexTree(vector<vector<double>>& pnt_cld); // Просто множество симплексов размерности 0 на основе множества точек
    void print() const;
    // Базовый деструктор
    ~SimplexTree();

    Graph eps_upgrade(double eps);
};
vector<vector<double>> read_to_pnt_cld(const string& filename);// Чтение облака точек из файла
void print_tree(const string& prefix,MyNode* rt);// Выписать дерево
void simplexes_of_dim(MyNode* current_node,int k,vector<Simplex>& simplexes,const Simplex& current_simplex);
// Запись чисел Бэтти для данного комплекса в файл
void write_betti_num_to_file(double max_eps,double step,const string& filename,vector<vector<double>> pnt_cld,int max_dim);

vector<vector<unsigned int>> distribute_load(vector<unsigned int> info,unsigned int N);

#endif //FASTHNF_SIMPLEX_H
