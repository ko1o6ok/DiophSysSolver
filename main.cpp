#include "Matrix.h"
#include <chrono>
#include <fstream>
using namespace std::chrono;

int main() {
    Matrix<int> A(200);

    A.randomize(20);
    Matrix<int> C(A);
    //cout << A << "----------------------"<<endl;
    //cout << C;
//    a[0][0] = 7;
//    a[1][0] = 7;
//    a[2][0] = 7;
    //cout << "Testing on start matrix:"<<endl<<a;
    //cout << a ;
    TDynamicVector<int> b;

    auto start = high_resolution_clock::now();
    decompose(A);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout <<"v3: "<<duration.count()<< " ms"<< endl;

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
