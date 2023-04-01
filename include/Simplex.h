//
// Created by ф on 01.04.2023.
//

#ifndef FASTHNF_SIMPLEX_H
#define FASTHNF_SIMPLEX_H

#include <vector>
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
    int dimension;
    vector<Simplex> faces;
    vector<Vertex> vertices; // Вершины должны составлять упорядоченный по индексам список
    vector<unsigned long> make_word(); // Вернуть упорядоченную по возрастанию последовательность номеров вершин
};

struct MyNode{
    vector<MyNode*> children;
    unsigned long last_index;
    static struct MyNode* make_a_node(unsigned long last_index);
    void add_child(MyNode* node, unsigned long last_index);
};

// Симплексное дерево
class SimplexTree{
public:
    MyNode* root;
    void insert_simplex(Simplex s);
    void make_a_root();


public:
    // Базовый конструктор
    SimplexTree();
};


#endif //FASTHNF_SIMPLEX_H
