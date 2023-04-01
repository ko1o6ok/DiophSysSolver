
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

void SimplexTree::insert_simplex(Simplex s) {
    make_a_root();
    MyNode* current_node = root;
    auto word = s.make_word();
    for (int i = 0;i < word.size();i++) {
        auto ind = word[i];
        bool exists = false;
        for(auto& son:current_node->children)
            if(son->last_index == ind){
                current_node = son;
                exists = true;
                break;
            }
        if(!exists){
            for (int j = i; j < word.size(); j++) {
                auto my_ind = word[j];
                auto add_node = MyNode::make_a_node(my_ind);
                current_node->add_child(add_node,my_ind);
                current_node = add_node;
            }
            break;
        }
    }
}

void SimplexTree::make_a_root() {
    root = MyNode::make_a_node(0);
}

SimplexTree::SimplexTree() {
    root = MyNode::make_a_node(0);
}

struct MyNode *MyNode::make_a_node(unsigned long last_ind) {
    auto* node = new MyNode();
    node->last_index = last_ind;
    return node;
}

void MyNode::add_child(MyNode *node, unsigned long last_ind) {
    auto* my_node = make_a_node(last_ind);
    node->children.push_back(my_node);
}
