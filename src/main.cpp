#include "Matrix.h"
#include "Simplex.h"
#include "Graph.h"
#include <chrono>
#include <fstream>
using namespace std::chrono;

int main() {
    auto pnt_cld = read_to_pnt_cld(R"(C:\C++_proj\FastHNF\data.txt)");
//    for (auto& v:pnt_cld) {
//        for(auto& t:v)
//            cout << t<<", ";
//        cout << endl;
//    }
//    double eps = 6.7;
//    SimplexTree tree(pnt_cld);
//    Graph gr = tree.construct_from_point_cloud(5,eps);
//    tree.print();
//    cout << endl;
//    gr.print_adj_matrix();
//    auto b_n = tree.betti_numbers();
    write_betti_num_to_file(1.5,0.01,R"(C:\C++_proj\FastHNF\res.txt)",pnt_cld,5);
//    vector<double> a ={0.0,0.0,0.2};
//    vector<double> b ={1.0,0.0,0.0};
//    vector<double> c ={1.5,0.0,0.0};
//    //vector<double> d = {0.0,2.0,0.0};
//    vector<vector<double>> cloud = {a,b,c};
//
//    SimplexTree my_tree(cloud);
////    my_tree.insert_simplex({0,1,2});
////    my_tree.insert_simplex({0,2});
////    my_tree.insert_simplex({1,2});
////    my_tree.insert_simplex({2});
//    my_tree.construct_from_point_cloud(1,1.7);
//    my_tree.print();
//    for(auto& s:my_tree.all_simplexes_of_dim(1)){
//        cout << s << endl;
//    }
    //auto M = my_tree.border_operator_matrix(2);
    //M.print();
    //cout << "-----------------------"<<endl;
    //to_SNF(M).print();
//    cout << "Betti numbers :"<<endl;
//    for(auto& b_n:my_tree.betti_numbers())
//        cout << b_n << endl;
    /*TDynamicVector<int> a({1,2,3});
    cout << a << endl;
    a.push(4);
    cout << a;*/
//    Graph g(cloud);
//    g.connect_eps_neighbours(2.2);
//    //g.connect_vertices(1,2);
//    g.print_adj_matrix();
//    auto tree = SimplexTree();
//
//    tree.insert_simplex({3,2,1});
//    tree.insert_simplex({4,2,1});
//    //tree.insert_simplex({2,3});
////    cout << s << endl;
////    add_child(tree.root,1);
//    print_tree("",tree.root);
//    vector<Simplex> simplexes = tree.all_simplexes_of_dim(2);
//    //cout << simplexes.size();
//    for (const auto& s:simplexes) {
//        cout <<s<< endl;
//    }
//    cout << tree.num_vertices ;
    // Алиса
//    auto rand_matrix = gen_public_key(5,20);
//    // Сохранили
//    auto public_key = rand_matrix;
//    auto pr = decompose(rand_matrix);
//    rand_matrix.print();
//    //rand_matrix.print();
//    //((pr.second * A) * pr.first).print();
//    //A.print();
//
//
//    vector<int> msg = {-1,1,0,1,1};
//
//    TDynamicVector<int> message(msg);
//
//    cout << "Initial message is "<< message<< endl;
//
//    auto encrypted_message = encrypt(message,public_key,1);
//
//    cout << "Encrypted message is "<< encrypted_message<< endl;
//    //rand_matrix.print();
//    //public_key.print();
//
//    auto decrypted_message = decrypt(encrypted_message,pr.first,rand_matrix,pr.second,1);
//
//    cout << "Decrypted message is "<< decrypted_message<< endl;
//    Matrix<int> A(10);
//    TDynamicVector<int> a1({1, 1,1});
//    TDynamicVector<int> a2 ({-1,-1,-1});
//    TDynamicVector<int> a3({1,-1,1});
//    A[0] = a1;
//    A[1] = a2;
//    A[2] = a3;
//    A.randomize(1);
//    A[0] = A[1];
//    A[2] = A[1];
//    cout<< "----------------------"<<endl;
//    cout << "Testing on start matrix:"<<endl;
//    A.print() ;
    //TDynamicVector<int> b;

//    auto start = high_resolution_clock::now();
//    compute_SNF(A);
//    auto stop = high_resolution_clock::now();
//    auto duration = duration_cast<milliseconds>(stop - start);
//    cout <<"v3: "<<duration.count()<< " ms"<< endl;
//    cout << "The result is:"<< endl;
//    A.print();
//    auto B = A.transpose();
//    compute_SNF(B); // Если сделать этот же алгоритм от трансп. матрицы, то он вычисляет SNF
//    cout << endl;
//    B.print();
    //cout << endl << "Is it in SNF? -> "<<check_SNF(A)<< endl;
    //indexes_wrong_part(A);
    //print_wrong_part(A);
//    auto start1 = high_resolution_clock::now();
//    solve_SLDE_v4(C,b);
//    auto stop1 = high_resolution_clock::now();
//    auto duration1 = duration_cast<milliseconds>(stop1 - start1);
//    cout <<"v4: "<<duration1.count()<< " ms"<< endl;
//    ofstream out;
//    out.open(R"(C:\C++_proj\FastHNF\out2.txt)");
//    if(out.is_open()){
//        for (int n = 2; n < 300; ++n) {
//            Matrix a(n);
//            a.randomize(35);
//            TDynamicVector<int> v(n);
//            v[0] = 1; v[1] = 3;
//
//            auto start = high_resolution_clock::now();
//            solve_SLDE_v3(a,v);
//            auto stop = high_resolution_clock::now();
//            auto duration = duration_cast<milliseconds>(stop - start);
//
//            // To get the value of duration use the count()
//// member function on the duration object
//            out <<  n <<","<<duration.count()<< endl;
//        }
//    }

//    int size = 100;
//
//    Matrix a(size);
//    a.randomize(13);
//    // cout << a << endl;
//    auto v = TDynamicVector<int>(20);
//
//    auto start1 = high_resolution_clock::now();
//    solve_SLDE_v1(a,v);
//    auto stop1 = high_resolution_clock::now();
//    auto duration1 = duration_cast<milliseconds>(stop1 - start1);
//    cout <<duration1.count()<< " ms"<<endl;
//
//    Matrix a2(size);
//    a2.randomize(35);
//
//    auto start2 = high_resolution_clock::now();
//    solve_SLDE_v2(a2,v);
//    auto stop2 = high_resolution_clock::now();
//    auto duration2 = duration_cast<milliseconds>(stop2 - start2);
//    cout <<duration2.count()<< " ms" <<endl;
//
//    Matrix a3(size);
//    a3.randomize(35);
//    auto start3 = high_resolution_clock::now();
//    solve_SLDE_v3(a3,v);
//    auto stop3 = high_resolution_clock::now();
//    auto duration3 = duration_cast<milliseconds>(stop3 - start3);
//    cout <<duration3.count()<< " ms" <<endl;
    //TDynamicVector<int> v(20);
    //solve_SLDE_v2(a,v);
//    Matrix a(70);
//
//    a.randomize(40);
//
//    //cout << "Start matrix" << endl << a;
//    TDynamicVector<int> v(5);
//    v[0] = 1; v[1] = 3;
//    cout  << "Type before start"<<endl;
//    string s;
//    cin >> s;
//    //auto start = high_resolution_clock::now();
//    solve_SLDE(a,v);
    return 0;
}
