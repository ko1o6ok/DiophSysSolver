
#include <algorithm>
#include <utility>
#include "Simplex.h"

Vertex::Vertex(unsigned long ind) {
    index = ind;
}

vector<unsigned long> Simplex::make_word() {
    vector<unsigned long> word;
    word.reserve(vertices.size());
    for (auto& vert:vertices) {
        word.push_back(vert.index);
    }
    return word;
}

Simplex::Simplex(vector<int> v) {
    dimension = v.size() -1;
    sort(v.begin(), v.end());
    for(auto t:v)
        vertices.emplace_back(t);
}

ostream &operator<<(ostream &stream, const Simplex& s) {
    cout << "This simplex is "<< endl;
    cout << "- ";
    for (auto v:s.vertices)
        stream << v.index << " - ";
    return stream;
}

Simplex::~Simplex() {
    vertices.clear();
}

void Simplex::add_vertex(unsigned long ind) {
    vertices.emplace_back(ind);
}

Simplex::Simplex(Simplex const &s) {
    dimension = s.dimension;
    for(auto v:s.vertices)
        vertices.push_back(v);
}

Simplex::Simplex() {
    dimension = 0;
}

void SimplexTree::insert_simplex(Simplex s) {
    MyNode* current_node = root;

    auto word = s.make_word();
    for (int i = 0;i < word.size();i++) {
        auto ind = word[i];
        bool exists = false;
        for(auto son:current_node->children)
            if(son->last_index == ind){
                current_node = son;
                exists = true;
                break;
            }
        if(!exists){
            for (int j = i; j < word.size(); j++) {
                auto my_ind = word[j];
                current_node = add_child(current_node,my_ind);
                num_vertices++;
            }
            break;
        }
    }
}


SimplexTree::SimplexTree() {
    root = make_a_node(0);
    num_vertices = 0;
}

//vector<Simplex> SimplexTree::simplexes_of_dim(int k) {
//    vector<Simplex> res;
//    // Симплексы размерности k соотносятся с путями из корня длины k+1
//
//    return res;
//}

void SimplexTree::insert_simplex(vector<int> v) {
    auto s = Simplex(std::move(v));
    insert_simplex(s);
}

vector<Simplex> SimplexTree::all_simplexes_of_dim(int k) const {
    vector<Simplex> simplexes;
    simplexes_of_dim(root,k+1,simplexes,Simplex());
    return simplexes;
}

void print_tree(const string& prefix,MyNode* rt) {
    if(rt != nullptr){
        cout << prefix;
        cout << "|-->";
        cout << rt->last_index << endl;
        for(auto& child:rt->children)
            print_tree(prefix + "|   ",child);
    }
}

void simplexes_of_dim(MyNode* current_node,int k,vector<Simplex>& simplexes,const Simplex& current_simplex){
    if(k<0)
        return;
    if(current_node!= nullptr){
        if(k == 0){
            //cout << "Pushing "<< current_simplex << endl;
            simplexes.push_back(current_simplex);
        }

        for(auto& child:current_node->children){
            Simplex new_simplex (current_simplex);
            new_simplex.add_vertex(child->last_index);
            simplexes_of_dim(child,k-1,simplexes,new_simplex);
        }
    }
}

SimplexTree::~SimplexTree() = default;

struct MyNode* make_a_node(unsigned long last_ind) {
    auto* node = new MyNode();
    node->last_index = last_ind;
    return node;
}

MyNode* add_child(MyNode *node, unsigned long last_ind) {
    auto* my_node = make_a_node(last_ind);
    node->children.push_back(my_node);
    return my_node;
}
