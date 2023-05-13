
#include <algorithm>
#include <utility>
#include "Simplex.h"
#include "Graph.h"
#include <fstream>
#include <sstream>
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

//Graph SimplexTree::construct_from_point_cloud(unsigned long max_dim,double eps) {
//    //root->children.clear();
//    // Создали граф ближайших соседей
//    Graph gr(point_cloud);
//    max_dimension = max_dim;
//    gr.connect_eps_neighbours(eps);
////    cout << "HERE !!!"<<endl;
////    g.print_adj_matrix();
////    cout << "HERE !!!"<<endl;
////    for(auto& row:g.vertices){
////            cout << row << ", ";
////
////    }cout << endl;
//    // На его основе создаём комплекс Вьеториса-Рипса:
//
//    // Сначала заносим просто все точки
//    for (int i = 0; i < num_vertices; ++i) {
//        add_child(root,i);
//    }
//    if(max_dimension == 0)
//        return g;
//    // Индуктивно добавляем симплексы размерности <= max_dimension
//    for (int k = 1; k < max_dimension+1; ++k) {
//
//        auto simplexes = all_simplexes_of_dim(k-1); // Все симплексы размерности k-1
//
//        if(!simplexes.empty()){
//            for(auto& simplex:simplexes){
//                // Рассмотрим последовательность вершин данного симплекса (она упорядочена по возрастанию!)
//                auto his_word = simplex.make_word();
//                vector<unsigned long> N; // Вершины, которые будут добавлены
//                for (auto& u:his_word) {
//                    auto l_n = g.lower_neighbours(u); // Предшествующие вершины
////                    if(k == 2){
////                        cout << "SEARCHING LOWER NEIGHBOURS of "<<u << endl;
////
////                        for(auto& t:l_n)
////                            cout << t << ", ";
////                        cout << endl << "------------------"<<endl;
////                    }
//
//                    // Вставляем, сортируем и удаляем дубликаты
//                    N.insert(N.end(),l_n.begin(),l_n.end());
//                    sort(N.begin(),N.end());
//                    N.erase(unique(N.begin(),N.end()),N.end());
////                    cout << "N is"<<endl;
////                    for(auto& t:N)
////                        cout << t << ", ";
////                    cout << endl;
//                }
//                for(auto& vertex:N){
//                    // Наращиваем новый симплекс из предыдущего
//                    vector<unsigned long> ins({vertex});
//                    // Вставляем, сортируем и удаляем дубликаты
//                    bool not_connected = false;
//                    for (auto& his_vert:his_word) {
//                        if(gr.adj_matrix[vertex][his_vert] == 0){
//                            not_connected = true;
//                            break;
//                        }
//                    }
//                    if(not_connected)
//                        continue;
//                    ins.insert(ins.end(),his_word.begin(),his_word.end());
//                    sort(ins.begin(),ins.end());
//                    ins.erase(unique(ins.begin(),ins.end()),ins.end());
//                    if(ins.size()==k+1){
////                        if( k == 2){
////                            cout << "Inserting with k = "<< k <<endl;
////                            cout << ins.size() << " < - "<<endl;
////                            for(auto& t:ins)
////                                cout << t << ", ";
////                            cout << endl;
////                        }
//
//                        insert_simplex(ins);
//                    }
//
//                }
//            }
//        }
//        else{
//            break;
//        }
//    }
//
//    return gr;
//}

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

void SimplexTree::print() const {
    print_tree("",root);
}

Matrix<long int> SimplexTree::border_operator_matrix(int dimension,int& adds,bool& no_simplexes) const {
    adds = 0;
    auto simplexes = all_simplexes_of_dim(dimension); // Вытащили все симплексы нужной размерности

//    if(dimension == 2){
//        cout << "SIMPLEXES" << endl;
//        cout << "-----------------------"<< endl;
//        for (auto& simplex:simplexes) {
//            cout << simplex <<endl;
//        }
//        cout << "-----------------------"<< endl;
//    }

    //cout << simplexes.size() << endl;
    if(simplexes.empty()){
        no_simplexes = true;
        return Matrix<long int>(1);
    }


    unsigned int num_simplexes = simplexes.size();
    if(dimension == 0){
       // cout << "FOUND "<<num_simplexes <<" OF DIM 0"<<endl;
        auto M = Matrix<long int>(num_simplexes);
        for (int i = 0; i < num_simplexes; ++i) {
            M[0][i] = 0;
        }
        return M;
    }



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

        //TDynamicVector<long> ins(num_simplexes); // Каждый вставляемый вектор имеет одинаковую длину
//        cout << "Initial simplex "<<endl;
//        cout << simplex << endl;
//        cout << "His words "<< endl;
//        for(auto& w:his_words){
//            for(auto& l:w)
//                cout << l << ", ";
//            cout << endl;
//        }
//
//        cout << endl;
        // Для каждого отдельного слова проверяем - было ли оно уже в матрице
        int coeff = 1;

        for (auto& word:his_words) {
            //cout << "COEFF "<< coeff<<endl;
            bool found = false;
            int i = 0;
            for (auto& w:all_border_words) {
//                cout << "Compared :"<<endl;
//                for(auto& t:w)
//                    cout << t << ", ";
//                cout << endl<< "And: ";
//                for(auto& z:word)
//                    cout << z << ", ";
//                cout << endl << "result is "<< (w == word) << endl;
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
                    //cout << p << endl;
                    //m.print();
                    //cout << endl;

                    res.push(p);
//                    cout << "MMMMM" <<endl;
//                    for (int m = 0; m < res.size(); ++m) {
//                        cout << res[m] << endl;
//                    }
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
        //res.push(ins);
        j++;
    }

    // Далее матрицу нужно дополнить до квадратной
    unsigned int a = res.size();
    unsigned int b = res[0].size();
//    cout << "HHHHHHHHHH" <<endl;
//    for (int i = 0; i < res.size(); ++i) {
//        cout << res[i] << endl;
//    }
//    cout << "HHHHHHHHHH" <<endl;
//    for (int i = 1; i < res.size(); ++i) {
//        if(res[i].size() > b)
//            b = res[i].size();
//    }
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
//    auto m_t = m.transpose();
//    for (int i = 0; i < m.GetSize(); ++i) {
//        int counter = 1;
//        auto& row = m_t[i];
//        for (int k = 0; k < m.GetSize(); ++k) {
//            if(row[k] == 1){
//                row[k] = counter;
//                counter = -counter;
//            }
//        }
//    }
    //m.val.pMem+=m.val[0].size();
    return m;
}
Matrix<long int> SimplexTree::border_operator_matrix(int dimension) const {

    auto simplexes = all_simplexes_of_dim(dimension); // Вытащили все симплексы нужной размерности
    //cout << simplexes.size() << endl;
    if(simplexes.empty())
        return Matrix<long int>(1);

    unsigned int num_simplexes = simplexes.size();
    if(dimension == 0){
        auto M = Matrix<long int>(num_simplexes);
        for (int i = 0; i < num_simplexes; ++i) {
            M[0][i] = 0;
        }
        return M;
    }



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

        //TDynamicVector<long> ins(num_simplexes); // Каждый вставляемый вектор имеет одинаковую длину
//        cout << "Initial simplex "<<endl;
//        cout << simplex << endl;
//        cout << "His words "<< endl;
//        for(auto& w:his_words){
//            for(auto& l:w)
//                cout << l << ", ";
//            cout << endl;
//        }
//
//        cout << endl;
        // Для каждого отдельного слова проверяем - было ли оно уже в матрице
        int coeff = 1;

        for (auto& word:his_words) {

            bool found = false;
            int i = 0;
            for (auto& w:all_border_words) {
//                cout << "Compared :"<<endl;
//                for(auto& t:w)
//                    cout << t << ", ";
//                cout << endl<< "And: ";
//                for(auto& z:word)
//                    cout << z << ", ";
//                cout << endl << "result is "<< (w == word) << endl;
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
                    //cout << p << endl;
                    //m.print();
                    //cout << endl;

                    res.push(p);
//                    cout << "MMMMM" <<endl;
//                    for (int m = 0; m < res.size(); ++m) {
//                        cout << res[m] << endl;
//                    }
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
        //res.push(ins);
        j++;
    }

    // Далее матрицу нужно дополнить до квадратной
    unsigned int a = res.size();
    unsigned int b = res[0].size();
//    cout << "HHHHHHHHHH" <<endl;
//    for (int i = 0; i < res.size(); ++i) {
//        cout << res[i] << endl;
//    }
//    cout << "HHHHHHHHHH" <<endl;
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
vector<int> SimplexTree::betti_numbers() const {
    vector<int> res;
    auto prev = border_operator_matrix(0);
    to_SNF(prev);
    int adds;
    int prev_adds = 0;
    bool no_simplexes = false;
    bool prev_no_simplexes = false;
    unsigned long pnt_cld_sz = point_cloud.size();
    //cout << "MAX DIM IS "<<max_dimension <<endl;
    for (int i = 1; i < min(max_dimension,pnt_cld_sz)+2; ++i) {
        // Вычисляем матрицу оператора i-мерной границы
        //cout << i << " - th border matrix"<<endl;
        auto M = border_operator_matrix(i,adds,no_simplexes);


        //cout << "Was this matrix "<<endl;
        //M.print();
        to_SNF(M); // Приводим её к смитовой нормальной форме
//        cout << "Current matrix "<<endl;
//        M.print();
//        cout << endl << "Previous matrix "<<endl;
//        prev.print();
//        cout << endl;
//        cout << "Her number of adds "<<prev_adds<<endl;
        int inserted_betti_number;
        if(no_simplexes){
            //cout << "FOUND NO SIMPLEXES OF DIM "<<i<<endl;
            if(prev_no_simplexes)
                inserted_betti_number = 0;
            else{
                inserted_betti_number = prev.nullity()-prev_adds;
            }
        }else{
            inserted_betti_number = prev.nullity()-prev_adds - M.rank();
        }
//        if(prev.nullity() - prev_adds == 0){
//            inserted_betti_number = 0;
//        }
//        if(inserted_betti_number == 2){
//
//        }

//        cout << "Prev nullity is "<<prev.nullity()<<endl;
//        cout << "Current rank is "<<M.rank()<<endl;
//        cout << "Inserting Betti number "<<inserted_betti_number <<endl;
        res.push_back(inserted_betti_number);
        prev = M;
        prev_adds = adds;
        prev_no_simplexes = no_simplexes;

//        if(i == 2)
//            break;
    }
    return res;
}

Graph SimplexTree::eps_upgrade(double eps) {


    //max_dimension = max_dim;

    // Здесь важно учесть предыдущее состояние, чтобы не считать одно и то же дважды
    //Graph prev_g(g);
    auto g_empty = g.is_empty();
    Graph difference = g.connect_eps_neighbours(eps);

    if((difference.is_empty())&&(!g_empty))
        return g;
//    cout << "GRAPH !!!"<<endl;
//    g.print_adj_matrix();
//    cout << "DIFFERENCE !!!"<<endl;
//    difference.print_adj_matrix();
//    cout << "HERE !!!"<<endl;
//    for(auto& row:g.vertices){
//            cout << row << ", ";
//
//    }cout << endl;
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
//                    if(k == 2){
//                        cout << "SEARCHING LOWER NEIGHBOURS of "<<u << endl;
//
//                        for(auto& t:l_n)
//                            cout << t << ", ";
//                        cout << endl << "------------------"<<endl;
//                    }

                    // Вставляем, сортируем и удаляем дубликаты
                    N.insert(N.end(),l_n.begin(),l_n.end());
                    sort(N.begin(),N.end());
                    N.erase(unique(N.begin(),N.end()),N.end());
//                    cout << "N is"<<endl;
//                    for(auto& t:N)
//                        cout << t << ", ";
//                    cout << endl;
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
//                        if( (eps > 0.32)&&(eps<0.35)){
//                            cout << "Inserting with k = "<< k <<endl;
//                            cout << ins.size() << " < - "<<endl;
//                            for(auto& t:ins)
//                                cout << t << ", ";
//                            cout << endl;
//                        }
//                        if(ins == vector<unsigned long>({1,2,3,9})){
//                            cout << "IS IT REALLY CONNECTED??" << endl;
//                            cout << g.adj_matrix[1][2] << endl;
//                            cout << g.adj_matrix[1][3] << endl;
//                            cout << g.adj_matrix[1][9] << endl;
//                            cout << g.adj_matrix[2][3] << endl;
//                            cout << g.adj_matrix[2][9] << endl;
//                            cout << g.adj_matrix[3][9] << endl;
//                            cout << "IS IT REALLY CONNECTED??" << endl;
//                        }
                        //if(ins != vector<unsigned long>({1,2,3,9}))
                        insert_simplex(ins);
                    }

                }
            }
        }
        else{
            break;
        }
    }

    return g;
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

void write_betti_num_to_file(double max_eps,double step,const string& filename,vector<vector<double>> pnt_cld,const int max_dim){
    ofstream out;
    out.open(filename);
    double eps = 0.0;

    SimplexTree tree(pnt_cld);
    tree.max_dimension = max_dim;

    //double time_taken = 0.0;
    while(eps < max_eps){
        // Естественно, лучше не каждый раз его заново создавать, а наращивать существующее
        // Здесь можно ускорить!!


        // Вот это мы замеряем

        //auto start = chrono::high_resolution_clock::now();
        tree.eps_upgrade(eps);
        auto b_n = tree.betti_numbers();
        //auto stop = chrono::high_resolution_clock::now();
        //auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        //time_taken += duration.count();
//        cout <<eps<<","<<duration.count()<< endl;

//        if((eps > 0.32)&&(eps<0.35)){
//            cout << "Current eps is "<<eps << endl;
//            tree.print();cout <<endl;
//            gr.print_adj_matrix();
//            cout << endl;
//
//        }

            //cout <<
//        if(b_n[0] == 2){

//            break;
//        }
        out << eps << " ";
        for(auto& t:b_n)
            out << t << " ";
        out << endl;
        //break;
        eps += step;
    }
    //cout << "Time taken "<<time_taken << " ms"<<endl;
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

vector<vector<unsigned int>> distribute_load(vector<unsigned int> info,unsigned int N){
    // В info на i - м месте количество симплексов i+1 - го ранга
    unsigned int max_rank = info.size();

    vector<vector<unsigned int>> work_arr(N);
    vector<unsigned int> sums_array(N);

    for (int i = 0; i < N; ++i) {
        vector<unsigned int> v(max_rank);
        for(auto& t:v)
            t = 0;
        work_arr[i] = v;

        sums_array[i] = 0;
    }



    // Заполняем
    for(int r = 0;r<max_rank;r++){
        unsigned int num = info[r];
        long j = 0;

        while(num!=0){
            if(j>= N)
                j -= N;
            work_arr[j][r] += 1;
            num --;
            j ++;
        }
    }

    for (int i = 0; i < N; ++i) {
        unsigned int S = 0;
        for (int r = 0; r < work_arr[i].size(); ++r) {
            S += work_arr[i][r] * (r+1);
        }
        sums_array[i] = S;
    }
    // Далее берём максимум и минимум: уменьшаем максимум за счёт увеличения минимума до значения меньше текущего максимума
    // повторяем такую процедуру до тех пор, когда она уже ничего не меняет
    vector<vector<unsigned int>> prev_arr(N);

    unsigned int min;
    int min_pos;
    unsigned int max;
    int max_pos;

    while(prev_arr != work_arr){
//        cout << "-------------" << endl;
//        for (int i = 0; i < N; ++i) {
//            cout << "DEVICE - "<<i<<endl;
//            for (int j = 0; j < work_arr[i].size(); ++j) {
//                cout << "Rank "<<j + 1 <<" num is "<<work_arr[i][j]<< endl;
//            }
//        }
//        cout << "-------------" << endl;

        prev_arr = work_arr;

        min = sums_array[0];
        max = sums_array[0];
        min_pos = 0;
        max_pos = 0;

        for (int i = 1; i < N; ++i) {
            auto t = sums_array[i];
            if(t < min){
                min = t;
                min_pos = i;
            }
            if(t>max){
                max = t;
                max_pos = i;
            }
        }
        if(min_pos == max_pos)
            return work_arr;

        auto current_max = max;
        //unsigned int add = 0;

        auto& a = work_arr[min_pos];
        auto& b = work_arr[max_pos];


        unsigned int ind = b.size()-1;

        while(true){
            if((b[ind] > 0)&&(min + ind +1 < max)){
                b[ind] --;
                a[ind] ++;
                 max -= ind +1;
                 min += ind +1;
                 sums_array[max_pos] -= ind+1;
                 sums_array[min_pos] += ind +1;
            }
            else{
                ind --;
            }

            if(ind == -1)
                break;
        }
        cout << "SUMS ARRAY: ";
        for(auto& s:sums_array)
            cout << s <<", ";
        cout << endl;

    }
    return work_arr;
}
