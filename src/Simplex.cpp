#include <algorithm>
#include <utility>
#include <fstream>
#include <sstream>
#include "Simplex.h"
#include "Graph.h"
vector<unsigned long> Simplex::make_word() {
    vector<unsigned long> word;
    word.reserve(vertices.size());
    for (auto& vert:vertices) {
        word.push_back(vert);
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
        stream << v << " - ";
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
                add.push_back(v[j]);
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
void SimplexTree::insert_simplex(vector<unsigned long> v){
    auto s = Simplex(std::move(v));
    insert_simplex(s);
}
vector<Simplex> SimplexTree::all_simplexes_of_dim(int k) const {
    vector<Simplex> simplexes;
    simplexes_of_dim(root,k+1,simplexes,Simplex());
    return simplexes;
}
SimplexTree::SimplexTree(vector<vector<double>> &pnt_cld) {
    point_cloud = pnt_cld;
    root = make_a_node(0);
    num_vertices = pnt_cld.size();
    g = Graph(pnt_cld);
    max_dimension = num_vertices;
    for (int i = 0; i < num_vertices; ++i) {
        add_child(root,i);
    }
}
Matrix<long int> SimplexTree::border_operator_matrix(int dimension,int& adds,bool& no_simplexes) const {
    adds = 0;
    auto simplexes = all_simplexes_of_dim(dimension); // Вытащили все симплексы нужной размерности
    if(simplexes.empty()){
        no_simplexes = true;
        return Matrix<long int>(1);
    }
    unsigned int num_simplexes = simplexes.size();
    if(dimension == 0){
        auto M = Matrix<long int>(num_simplexes);
        for (int i = 0; i < num_simplexes; ++i) {
            M[0][i] = 0;
        }
        return M;
    }
    TDynamicVector<TDynamicVector<long int>> res;
    vector<vector<unsigned long>> all_border_words;
    int  j = 0;
    bool is_first = true;
    for(auto& simplex:simplexes){
        // Применяем к нему оператор границы
        auto his_words = simplex.get_border_words();
        // Для каждого отдельного слова проверяем - было ли оно уже в матрице
        int coeff = 1;
        for (auto& word:his_words) {
            bool found = false;
            int i = 0;
            for (auto& w:all_border_words) {
                if(w == word){
                    found = true;
                    res[i][j] = coeff;
                    break;
                }
                i++;
            }
            if(!found){
                all_border_words.push_back(word);
                auto p = TDynamicVector<long>(num_simplexes);
                if(is_first){
                    res[0] = p;
                    res[0][j] = coeff;
                    is_first = false;
                }else{
                    p[j] = coeff;
                    res.push(p);
                }
            }
            coeff = - coeff;
        }
        j++;
    }
    // Далее матрицу нужно дополнить до квадратной
    unsigned int a = res.size();
    unsigned int b = res[0].size();
    if(a > b){
        adds+=a-b;
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
    return m;
}
vector<int> SimplexTree::betti_numbers() const {
    vector<int> res;
    auto prev = Matrix<long int>(num_vertices);
    for (int i = 0; i < num_vertices; ++i) {
        prev[0][i] = 0;
    }
    int adds;
    int prev_adds = 0;
    bool no_simplexes = false;
    bool prev_no_simplexes = false;
    unsigned long pnt_cld_sz = point_cloud.size();
    for (int i = 1; i < min(max_dimension,pnt_cld_sz)+2; ++i) {
        // Вычисляем матрицу оператора i-мерной границы
        auto M = border_operator_matrix(i,adds,no_simplexes);
        to_SNF(M); // Приводим её к смитовой нормальной форме
        int inserted_betti_number;
        if(no_simplexes){
            if(prev_no_simplexes)
                inserted_betti_number = 0;
            else{
                inserted_betti_number = prev.nullity()-prev_adds;
            }
        }else{
            inserted_betti_number = prev.nullity()-prev_adds - M.rank();
        }
        res.push_back(inserted_betti_number);
        prev = M;
        prev_adds = adds;
        prev_no_simplexes = no_simplexes;
    }
    return res;
}
Graph SimplexTree::eps_upgrade(double eps) {
    // Здесь важно учесть предыдущее состояние, чтобы не считать одно и то же дважды
    auto g_empty = g.is_empty();
    Graph difference = g.connect_eps_neighbours(eps);
    if((difference.is_empty())&&(!g_empty))
        return g;
    // На его основе создаём комплекс Вьеториса-Рипса:
    if(max_dimension == 0)
        return g;
    // Индуктивно добавляем симплексы размерности <= max_dimension
    for (int k = 1; k < max_dimension+1; ++k) {
        auto simplexes = all_simplexes_of_dim(k-1); // Все симплексы размерности k-1
        if(!simplexes.empty()){
            for(auto& simplex:simplexes){
                // Рассмотрим последовательность вершин данного симплекса (она упорядочена по возрастанию!)
                auto his_word = simplex.make_word();
                vector<unsigned long> N; // Вершины, которые будут добавлены
                for (auto& u:his_word) {
                    auto l_n = g.lower_neighbours(u,difference); // Предшествующие вершины
                    // Предшествующие вершины должны входить в граф разности иначе они уже были добавлены
                    // Вставляем, сортируем и удаляем дубликаты
                    N.insert(N.end(),l_n.begin(),l_n.end());
                    sort(N.begin(),N.end());
                    N.erase(unique(N.begin(),N.end()),N.end());
                }
                for(auto& vertex:N){
                    // Наращиваем новый симплекс из предыдущего
                    vector<unsigned long> ins({vertex});
                    // Вставляем, сортируем и удаляем дубликаты
                    bool not_connected = false;
                    for (auto& his_vert:his_word) {
                        if(difference.adj_matrix[vertex][his_vert] == 0){
                            not_connected = true;
                            break;
                        }
                    }
                    if(not_connected)
                        continue;
                    ins.insert(ins.end(),his_word.begin(),his_word.end());
                    sort(ins.begin(),ins.end());
                    ins.erase(unique(ins.begin(),ins.end()),ins.end());
                    if(ins.size()==k+1){
                        insert_simplex(ins);
                    }
                }
            }
        }
        else{
            break;
        }
    }
    // !!!!!!!
    num_vertices = g.num_vertices;
    // !!!!!
    return g;
}
void simplexes_of_dim(MyNode* current_node,int k,vector<Simplex>& simplexes,const Simplex& current_simplex){
    if(k<0)
        return;
    if(current_node!= nullptr){
        if(k == 0){
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
void write_betti_num_to_file(double max_eps,double step,const string& filename,vector<vector<double>> pnt_cld,const int max_dim){
    ofstream out;
    out.open(filename);
    double eps = 0.0;

    SimplexTree tree(pnt_cld);
    tree.max_dimension = max_dim;

    double time_taken = 0.0;
    double progress ;
    cout << "Progress: "<<endl;
    while(eps < max_eps){
        progress = eps/max_eps;
        int barWidth = 70;
        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();
        auto start = chrono::high_resolution_clock::now();
        tree.eps_upgrade(eps);
        auto b_n = tree.betti_numbers();
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        time_taken += duration.count();
        out << eps << " ";
        for(auto& t:b_n)
            out << t << " ";
        out << endl;
        eps += step;
    }
    cout << endl << "Time taken "<<time_taken<<" ms"<<endl;
    system("pause");
}
vector<vector<double>> read_to_pnt_cld(const string& filename){
    fstream in;
    in.open(filename,ios::in);
    vector<vector<double>> res;
    if(in.is_open()){
        string tp;
        while(getline(in,tp)){
            vector<double> add;
            stringstream ss(tp);
            string word;
            while(ss >> word)
                add.push_back(atof(word.c_str()));
            res.push_back(add);
        }
        in.close();
    }
    return res;
}


