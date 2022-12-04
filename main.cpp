#include "Matrix.h"
#include <chrono>
#include <fstream>
using namespace std::chrono;

int main() {
//    ofstream out;
//    out.open(R"(C:\C++_proj\FastHNF\out.txt)");
//    if(out.is_open()){
//        for (int n = 2; n < 300; ++n) {
//            Matrix a(n);
//
//            a.randomize(40);
////    a[0][0] = -7;
////    a[0][1] = -10;
////    a[0][2] = -2;
////
////    a[1][0] = -5;
////    a[1][1] = 5;
////    a[1][2] = 0;
////
////    a[2][0] = -10;
////    a[2][1] = -1;
////    a[2][2] = -6;
//            //cout << "Start matrix" << endl << a;
//            TDynamicVector<int> v(n);
//            v[0] = 1; v[1] = 3;
//
//            auto start = high_resolution_clock::now();
//            solve_SLDE(a,v);
//            auto stop = high_resolution_clock::now();
//            auto duration = duration_cast<milliseconds>(stop - start);
//
//            // To get the value of duration use the count()
//// member function on the duration object
//            out <<  n <<","<<duration.count()<< endl;
//        }
//    }

    Matrix a(70);

    a.randomize(40);

    //cout << "Start matrix" << endl << a;
    TDynamicVector<int> v(5);
    v[0] = 1; v[1] = 3;
    cout  << "Type before start"<<endl;
    string s;
    cin >> s;
    //auto start = high_resolution_clock::now();
    solve_SLDE(a,v);
    return 0;
}
