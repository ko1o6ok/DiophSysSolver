#ifndef FASTHNF_SIMPLICIALCOMPLEX_H
#define FASTHNF_SIMPLICIALCOMPLEX_H

#include "Simplex.h"

class SimplicialComplex {
    int max_dimension;
    vector<Vertex> vertices; // Вершины
    vector<Simplex> simplexes; // Все симплексы
    Matrix<int> border_operator_matrix(int k); // Матрица линейного оператора границы для размерности k
    int k_th_betti_number(int k); // Число Бэтти размерности k
    SimplexTree reduce_to_simplex_tree();// Сжимаем комплекс в дерево
    vector<Simplex> simplexes_of_dim(int k);
};


#endif //FASTHNF_SIMPLICIALCOMPLEX_H
