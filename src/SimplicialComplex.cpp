

#include "SimplicialComplex.h"

#include <utility>

SimplexTree SimplicialComplex::reduce_to_simplex_tree() {
    auto tree = SimplexTree();
    for(auto& simplex:simplexes)
        tree.insert_simplex(simplex);
    return tree;
}

vector<Simplex> SimplicialComplex::simplexes_of_dim(int k) {
    vector<Simplex> res;
    res.reserve(vertices.size()); // Здесь можно много подумать и даже использовать какой-нибудь тервер
    for(auto& simplex:simplexes)
        if(simplex.dimension == k)
            res.push_back(simplex);
    return res;
}

SimplicialComplex::SimplicialComplex() {
    max_dimension = 0;
}
