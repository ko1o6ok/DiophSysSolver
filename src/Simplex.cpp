
#include <algorithm>
#include <utility>
#include "Simplex.h"
#include "Graph.h"

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

Simplex::Simplex(vector<unsigned long> v) {
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

vector<vector<unsigned long>> Simplex::get_border_words() {
    auto v = make_word();
    vector<vector<unsigned long>> res;
    for (int i = 0; i < v.size(); ++i) {
        vector<unsigned long> add;
        for (int j = 0; j < v.size(); j++) {
            if(j != i)
                add.push_back(j);
        }
        res.push_back(add);
    }
    return res;
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
    root = make_a_node(-1);
    num_vertices = 0;
}

//vector<Simplex> SimplexTree::simplexes_of_dim(int k) {
//    vector<Simplex> res;
//    // Симплексы размерности k соотносятся с путями из корня длины k+1
//
//    return res;
//}

void SimplexTree::insert_simplex(vector<unsigned long> v){
    auto s = Simplex(std::move(v));
    insert_simplex(s);
}

vector<Simplex> SimplexTree::all_simplexes_of_dim(int k) const {
    vector<Simplex> simplexes;
    simplexes_of_dim(root,k+1,simplexes,Simplex());
    return simplexes;
}

void SimplexTree::construct_from_point_cloud(unsigned long max_dimension,double eps) {
    // Создали граф ближайших соседей
    Graph g(point_cloud);
    g.connect_eps_neighbours(eps);

//    for(auto& row:g.vertices){
//            cout << row << ", ";
//
//    }cout << endl;
    // На его основе создаём комплекс Вьеториса-Рипса:

    // Сначала заносим просто все точки
    for (int i = 0; i < num_vertices; ++i) {
        add_child(root,i);
    }

    // Индуктивно добавляем симплексы размерности <= max_dimension
    for (int k = 1; k < max_dimension+1; ++k) {

        auto simplexes = all_simplexes_of_dim(k-1); // Все симплексы размерности k-1

        if(!simplexes.empty()){

            vector<vector<unsigned long>> words; // Все "слова", соответствующие этим симплексам
            words.reserve(simplexes.size());
            for (auto& simplex:simplexes) {
                words.push_back(simplex.make_word());
            }
            for (auto& word:words) {
//                cout << "Current word is "<< endl;
//                for(auto& t:word)
//                    cout << t << ", ";
//                cout << endl;
//                if(word == vector<unsigned long>({1,2}))
//                    int klmn = 5;
                // Здесь важно учесть случай, когда буква в слове лишь одна!!
                unsigned long last, prev_last;
                if(k == 1)
                {
                    last = word[0];
                    prev_last = last;
                }else{
                    /*if(k-1>= word.size())
                        continue;*/
                    last = word[k-1];
                    prev_last = word[k-2];
                }



                auto s = last + 1;
                auto row = g.adj_matrix[prev_last];
                bool pushed = false;
                while (s < g.num_vertices){
                    if(row[s] != 0){
                        word.push_back(s);
                        pushed = true;
                        break;
                    }
                    s++;
                }
                if(pushed){
                    insert_simplex(word);
                    //cout << "Inserting "<<endl;
//                    for(auto& t:word)
//                        cout << t << ", ";
//                    cout << endl;
                }

        }
            words.clear();
            simplexes.clear();
        }
        else{
            break;
        }
    }


}

SimplexTree::SimplexTree(vector<vector<double>> &pnt_cld) {
    point_cloud = pnt_cld;
    root = make_a_node(0);
    num_vertices = pnt_cld.size();
}

void SimplexTree::print() const {
    print_tree("",root);
}

Matrix<long int> SimplexTree::border_operator_matrix(int dimension) const {
    auto simplexes = all_simplexes_of_dim(dimension); // Вытащили все симплексы нужной размерности
    if(simplexes.empty())
        return Matrix<long int>(1);

    unsigned int num_simplexes = simplexes.size();
    if(dimension == 0)
       return Matrix<long int>(num_simplexes);


    TDynamicVector<TDynamicVector<long int>> res;

//    cout << num_simplexes << endl;
//    m.print();
    vector<vector<unsigned long>> all_border_words;
    int  j = 0;

    bool is_first = true;

    for(auto& simplex:simplexes){
        // Применяем к нему оператор границы
//        cout << "His vertices "<< endl;
//        for(auto& v:simplex.vertices)
//            cout << v.index << ", ";
//        cout << endl;
        auto his_words = simplex.get_border_words();
        cout << "His words "<< endl;
        for(auto& w:his_words){
            for(auto& l:w)
                cout << l << ", ";
            cout << endl;
        }

        cout << endl;
        // Для каждого отдельного слова проверяем - было ли оно уже в матрице
        int coeff = -1;

        for (auto& word:his_words) {

            bool found = false;
            int i = 0;
            for (auto& w:all_border_words) {
                if(w == word){
                    found = true;
                    auto& tmp = res[i];
                    for (unsigned int k = 0; k < j-tmp.size()+1; ++k) {
                        tmp.push(0);
                    }
                    cout << "THUHUFDH "<<tmp << endl;
                    res[i][j] = coeff;

//                    cout << "Pushing word "<<endl;
//                    for(auto& t:word)
//                        cout << t << ", ";
//                    cout << endl;
//                    cout << " With coeff: "<<coeff << endl;
                    break;
                }
                i++;
            }
            if(!found){
                all_border_words.push_back(word);

                auto p = TDynamicVector<long>(num_simplexes);
                if(is_first){
                    res[0][j] = coeff;
                    is_first = false;
                }else{
                    p[j] = coeff;
                    //cout << p << endl;
                    //m.print();
                    //cout << endl;

                    res.push(p);
                    //m.print();
//                    cout << endl;
//                    cout << "Pushing word "<<endl;
//                    for(auto& t:word)
//                        cout << t << ", ";
//                    cout << endl;
//                    cout << " With coeff: "<<coeff << endl;
                }

            }

            coeff = - coeff;
        }
        j++;
    }

    // Далее матрицу нужно дополнить до квадратной
    unsigned int a = res.size();
    unsigned int b = res[0].size();
    cout << "HHHHHHHHHH" <<endl;
    for (int i = 0; i < res.size(); ++i) {
        cout << res[i] << endl;
    }
//    for (int i = 1; i < res.size(); ++i) {
//        if(res[i].size() > b)
//            b = res[i].size();
//    }
    if(a > b){
        // Высота больше ширины
        for (int i = 0; i < a; ++i) {
            auto& v = res[i];// Строка
            for (unsigned int k = v.size(); k < a; ++k) {
                v.push(0);
            }
        }

    }
    if(a < b){
        // Ширина больше высоты
        auto p = TDynamicVector<long>(b);
        for (int i = 0; i < b-a; ++i) {
            res.push(p);
        }
    }
    Matrix<long> m (res);
    //m.val.pMem+=m.val[0].size();
    return m;
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
