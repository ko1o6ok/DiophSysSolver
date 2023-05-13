#include "Matrix.h"
#include "Simplex.h"
#include <fstream>
int main() {
    cout << "<<PERSISTENT HOMOLOGY>>" <<endl;
    double max_eps,step;
    int max_dim;
    cout << "Input max epsilon (for the dot, please, use '.' ):"<<endl;
    cin >> max_eps;
    cout << "Input step size (for the dot, please, use '.' ):"<<endl;
    cin >> step;
    cout << "Input max dimension:"<<endl;
    cin >> max_dim;
    auto pnt_cld = read_to_pnt_cld(R"(C:\C++_proj\FastHNF\data.txt)");
    write_betti_num_to_file(max_eps,step,R"(C:\C++_proj\FastHNF\res.txt)",pnt_cld,max_dim);
    return 0;
}
