

#include "Matrix.h"
#include<cstdlib>
#include <cmath>
#include <cstring>
#include <algorithm>
using namespace std;
template <typename T>
Matrix<T>::Matrix(unsigned int n) {
    N = n;
    val = TDynamicVector<TDynamicVector<T>>(N);
    for (int i = 0; i < N; ++i)
            val[i] = TDynamicVector<T>(N);
}
template<> Matrix<long>::Matrix(unsigned int n) {
    N = n;
    val = TDynamicVector<TDynamicVector<long>>(N);
    for (int i = 0; i < N; ++i)
        val[i] = TDynamicVector<long>(N);
}
template <typename T>
Matrix<T>::Matrix(const TDynamicVector<TDynamicVector<T>>& v) {
    N = v.size();
    val = v;
}
template<> Matrix<long>::Matrix(const TDynamicVector<TDynamicVector<long>>& v) {
    N = v.size();
    val = v;
}
template <typename T>
TDynamicVector<T>& Matrix<T>::operator[](unsigned int i) {
    return val[i];
}
template<> TDynamicVector<long>& Matrix<long>::operator[](unsigned int i) {
    return val[i];
}
template <typename T>
TDynamicVector<T> Matrix<T>::operator*(const TDynamicVector<T>& vec) {
    // Предполагается, что размер vec равен N
    TDynamicVector<T> res(N);
    for (int i = 0; i < N; ++i) {
        res[i] = val[i] * vec;
    }
    return res;
}

template<typename T>
ostream &operator<<(ostream &ostr, Matrix<T>& m) {
    for (int i = 0; i <m.N; ++i) {
        ostr << m[i] ;
        ostr << endl;
    }
    return ostr;
}
template<typename T>
unsigned int Matrix<T>::GetSize() const {
    return N;
}
template<> vector<vector<double> >
Matrix<double>::add_matrix(vector<vector<double> > matrix_A,
           vector<vector<double> > matrix_B, int split_index,
           int multiplier)
{
    for (auto i = 0; i < split_index; i++)
        for (auto j = 0; j < split_index; j++)
            matrix_A[i][j]
                    = matrix_A[i][j]
                      + (multiplier * matrix_B[i][j]);
    return matrix_A;
}
// Strassen's algorithm
template<> vector<vector<double> >
Matrix<double>::multiply_matrix(vector<vector<double> > matrix_A,
                vector<vector<double> > matrix_B)
{
    int col_1 = matrix_A[0].size();
    int row_1 = matrix_A.size();
    int col_2 = matrix_B[0].size();
    int row_2 = matrix_B.size();

    if (col_1 != row_2) {
        cout << "\nError: The number of columns in Matrix "
                "A  must be equal to the number of rows in "
                "Matrix B\n";
        return {};
    }

    vector<double> result_matrix_row(col_2, 0);
    vector<vector<double> > result_matrix(row_1,
                                       result_matrix_row);

    if (col_1 == 1)
        result_matrix[0][0]
                = matrix_A[0][0] * matrix_B[0][0];
    else {
        int split_index = col_1 / 2;

        vector<double> row_vector(split_index, 0);

        vector<vector<double> > a00(split_index, row_vector);
        vector<vector<double> > a01(split_index, row_vector);
        vector<vector<double> > a10(split_index, row_vector);
        vector<vector<double> > a11(split_index, row_vector);
        vector<vector<double> > b00(split_index, row_vector);
        vector<vector<double> > b01(split_index, row_vector);
        vector<vector<double> > b10(split_index, row_vector);
        vector<vector<double> > b11(split_index, row_vector);

        for (auto i = 0; i < split_index; i++)
            for (auto j = 0; j < split_index; j++) {
                a00[i][j] = matrix_A[i][j];
                a01[i][j] = matrix_A[i][j + split_index];
                a10[i][j] = matrix_A[split_index + i][j];
                a11[i][j] = matrix_A[i + split_index]
                [j + split_index];
                b00[i][j] = matrix_B[i][j];
                b01[i][j] = matrix_B[i][j + split_index];
                b10[i][j] = matrix_B[split_index + i][j];
                b11[i][j] = matrix_B[i + split_index]
                [j + split_index];
            }

        vector<vector<double> > p(multiply_matrix(
                a00, add_matrix(b01, b11, split_index, -1)));
        vector<vector<double> > q(multiply_matrix(
                add_matrix(a00, a01, split_index,1), b11));
        vector<vector<double> > r(multiply_matrix(
                add_matrix(a10, a11, split_index,1), b00));
        vector<vector<double> > s(multiply_matrix(
                a11, add_matrix(b10, b00, split_index, -1)));
        vector<vector<double> > t(multiply_matrix(
                add_matrix(a00, a11, split_index,1),
                add_matrix(b00, b11, split_index,1)));
        vector<vector<double> > u(multiply_matrix(
                add_matrix(a01, a11, split_index, -1),
                add_matrix(b10, b11, split_index,1)));
        vector<vector<double> > v(multiply_matrix(
                add_matrix(a00, a10, split_index, -1),
                add_matrix(b00, b01, split_index,1)));

        vector<vector<double> > result_matrix_00(add_matrix(
                add_matrix(add_matrix(t, s, split_index,1), u,
                           split_index,1),
                q, split_index, -1));
        vector<vector<double> > result_matrix_01(
                add_matrix(p, q, split_index,1));
        vector<vector<double> > result_matrix_10(
                add_matrix(r, s, split_index,1));
        vector<vector<double> > result_matrix_11(add_matrix(
                add_matrix(add_matrix(t, p, split_index,1), r,
                           split_index, -1),
                v, split_index, -1));

        for (auto i = 0; i < split_index; i++)
            for (auto j = 0; j < split_index; j++) {
                result_matrix[i][j]
                        = result_matrix_00[i][j];
                result_matrix[i][j + split_index]
                        = result_matrix_01[i][j];
                result_matrix[split_index + i][j]
                        = result_matrix_10[i][j];
                result_matrix[i + split_index]
                [j + split_index]
                        = result_matrix_11[i][j];
            }

        a00.clear();
        a01.clear();
        a10.clear();
        a11.clear();
        b00.clear();
        b01.clear();
        b10.clear();
        b11.clear();
        p.clear();
        q.clear();
        r.clear();
        s.clear();
        t.clear();
        u.clear();
        v.clear();
        result_matrix_00.clear();
        result_matrix_01.clear();
        result_matrix_10.clear();
        result_matrix_11.clear();
    }
    return result_matrix;
}
// Адаптивный алгоритм Штрассена
// pivot point = 1_000
template<> Matrix<double> Matrix<double>::operator*(Matrix<double>& m) {
    if(N > 1000){
        size_t s = val.size();
        vector<vector<double>> my_matrix;
        for(int i = 0;i<s;i++){
            vector<double> add;
            auto temp = val[i];
            for(int j = 0;j<s;j++)
                add.push_back(temp[j]);
            my_matrix.push_back(add);
        }
        vector<vector<double>> his_matrix;
        for(int i = 0;i<s;i++){
            vector<double> add;
            auto temp = m.val[i];
            for(int j = 0;j<s;j++)
                add.push_back(temp[j]);
            his_matrix.push_back(add);
        }
        auto a(multiply_matrix(my_matrix,his_matrix));
        Matrix<double> res(s);
        for(int i = 0;i<s;i++)
            for(int j = 0;j<s;j++)
                res[i][j] = a[i][j];
        return res;
    }
    Matrix<double> help = m.transpose();
    Matrix<double> r(N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            r[i][j] = val[i] * help[j];
    return r;
}
template<> Matrix<int> Matrix<int>::operator*(Matrix<int>& m) {

    Matrix<int> help = m.transpose();
    Matrix<int> r(N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            r[i][j] = val[i] * help[j];
    return r;
}
template <typename T>
Matrix<T> Matrix<T>::transpose() {
    Matrix r(N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            r.val[i][j] = val[j][i];
        }
    }
    return r;
}

template<> void Matrix<int>::randomize(unsigned int range) {
    random_device rd;
    mt19937_64 gen(rd());
    uniform_int_distribution<> distr(-range,range);
    //srand((unsigned) time(nullptr));
    //auto t = range * 2;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            val[i][j] = distr(gen);
        }
    }
}


template<> void Matrix<int>::print() {
    for (int i = 0; i <N; ++i) {
        cout << val[i] ;
        cout << endl;
    }
}
//template<> void Matrix<short>::print() {
//    for (int i = 0; i <N; ++i) {
//        cout << val[i] ;
//        cout << endl;
//    }
//}
template<typename T>
void Matrix<T>::SetSize(unsigned int s) {
    N = s;
}

template<>
long Matrix<long>::nullity() const {
    int S = 0;
    for (int i = 0; i < N; ++i)
        if(val[i][i] == 0)
            S++;
    return S;
}

template<>
long Matrix<long>::rank() const {
    return N-nullity();
}

template<> void Matrix<long>::SetSize(unsigned int s) {
    N = s;
}
template<> void Matrix<long>::print() {
    for (int i = 0; i <N; ++i) {
        cout << val[i] ;
        cout << endl;
    }
}
template<> void Matrix<double>::print() {
    for (int i = 0; i <N; ++i) {
        cout << val[i] ;
        cout << endl;
    }
}
//vector<int> indexes_wrong_part(Matrix<int>& A){
//    unsigned int n = A.GetSize();
//    // Проходимся по каждой строке
//    // По столбцам не нужно из-за особенностей алгоритма
//   // int num = 0 ;
//   vector<int> indexes;
//    for (int i = 0; i < n; ++i) {
//        auto& row = A[i];
//        for (int j = i+1; j < n; ++j) {
//            if(row[j]!=0){
//                indexes.push_back(i);
//                //cout << row << endl;
//                //num++;
//                break;
//            }
//
//        }
//    }
//    return indexes;
//    //cout << "There's "<<num<<" of them"<< endl;
//}


void compute_SNF(Matrix<long>& A) {
    unsigned int n = A.GetSize();

    // (A | U)
    // U A R = D - почти диагональная
    //Разрешённые действия (строки/столбцы a и b):
    // a = a +k *b
    // a = -a
    // a < - > b (поменять местами)

    // Проход по столбцам - многомерный аналог алгоритма Евклида: сложность O(1)*n = O(n)
    // Умножение вектора на скаляр - O(n)
    // Итого n * O(n^2) = O(n^3)

    // Пока нет квантовых компьютеров, можно лишь понижать константу

    int num_nonzero = 10; // Число ненулевых элементов в текущем столбце в данный момент
    int pos_min_nonzero; // Позиция минимального ненулевого элемента
    int min_nz; // Минимальный по модулю ненулевой элемент в k-м столбце
    int second_min_nz; // Второй минимум
    int pos_sec_nz; // Позиция второго минимума

    int m; // Здесь будет целая часть отношения в алгоритме Евклида (div в общем)
    int multiplier;

    int saved_zero;

    bool is_all_zero; // Состоит ли столбец только из нулей
    for (int k = 0; k < n; k++) {
        is_all_zero = false;
        auto A_T = A.transpose(); // Чтобы компу было удобнее ходить в память

//        TDynamicVector<int> col = A_T[k];// Берём k-й столбец
        TDynamicVector<long> col(n);
        for (int p = 0; p < n; p++) {
            col[p] = A[p][k];
        }
        // Расширенный алгоритм Евклида над столбцом

        while (true) {
            num_nonzero = 1;
            // Нули нам явно не нужны. They are the good guys.
            int j = k;
            // Здесь ошибка: не учтено, что все могут быть нулевыми!!
            while (col[j] == 0){
                j++;
                if(j==n){
                    is_all_zero = true;
                    break;
                }

            }
            if(is_all_zero)
                break;
            // Нули в начале пройдены
            min_nz = col[j];
            pos_min_nonzero = j;

            second_min_nz = col[j];
            pos_sec_nz = j;
            bool flag = false;
            for (int i = j + 1; i < n; i++) {
                int t = col[i];
                if (t != 0) {
                    if (abs(t) < abs(min_nz)) {
                        second_min_nz = min_nz;
                        pos_sec_nz = pos_min_nonzero;

                        min_nz = t;
                        pos_min_nonzero = i;
                    } else {
                        if (!flag) {
                            second_min_nz = t;
                            pos_sec_nz = i;
                            flag = true;
                        }
                        if ((t != min_nz) && (abs(t) < abs(second_min_nz))) {
                            second_min_nz = t;
                            pos_sec_nz = i;
                        }
                    }
                    num_nonzero++;
                }
            }
            if (num_nonzero <= 1) break; // Алгоритм Евклида гарантирует, что мы точно придём сюда
            // Теперь совершаем алгоритм Евклида над двумя мин. элементами
            while ((min_nz != 0) && (second_min_nz != 0)) {
                m = floor(second_min_nz / min_nz);
                A[pos_sec_nz] = A[pos_sec_nz] - A[pos_min_nonzero] * m;
                second_min_nz = second_min_nz % min_nz;
                col[pos_sec_nz] = col[pos_sec_nz] % min_nz;
                if (second_min_nz == 0)
                    break;
                m = floor(min_nz / second_min_nz);
                A[pos_min_nonzero] = A[pos_min_nonzero] - A[pos_sec_nz] * m;
                min_nz = min_nz % second_min_nz;
                col[pos_min_nonzero] = col[pos_min_nonzero] % second_min_nz;

            }
//            cout << "Current stage: "<<endl;
//            A.print();
            if (min_nz == 0) {
                min_nz = second_min_nz;
                pos_min_nonzero = pos_sec_nz;
            }
            // Вертикальный ход
//            cout << "Current stage: "<<endl;
//            A.print();

            for (int l = 0; l < n; l++)
                if ((l != pos_min_nonzero) && (col[l] != 0)) {

                    multiplier = floor(col[l] / min_nz);
                    col[l] = col[l] - min_nz * multiplier;
                    A[l] = A[l] - A[pos_min_nonzero] * multiplier;
//                    cout << "Current stage: "<<endl;
//                    A.print();
                }
//            cout << "Current stage: "<<endl;
//            A.print();
        }
        if(is_all_zero)
            continue;
        // Апосля того, как был зачищен стобец, меняем местами k-ю и min_pos строки
        if (k != n - 1)
            swap(A[k], A[pos_min_nonzero]);
        // А теперь, пользуясь оставшимся ненулевым элементом, уменьшаем элементы строки
        int el = A[k][k];
        // Проходим по k-й строке
        auto T = A.transpose();
        for (int i = k + 1; i < n; i++) {
            long &t = A[k][i];
            if (abs(el) <= abs(t)) {
                T[i] = T[i] - T[k] * floor(t/el);
                //t = t % el;
            }
        }
        A = T.transpose();
    }
    // Финальная полировка
    for (int k = 0; k < n; k++) {
        int el = A[k][k];
        if(el != 0){
            // Проходим по k-й строке
            for (int i = k + 1; i < n; i++) {
                long &t = A[k][i];
                if (t != 0)
                    if (abs(el) <= abs(t))
                        t = t % el;
            }
        }

    }
}
Matrix<long> to_SNF(Matrix<long>& A){
    compute_SNF(A);
    //compute_SNF(A);
//    cout << "Computed SNF of "<<endl;
//    A.print();
//    cout << endl;
    A =  A.transpose();
    //A.print();
    //cout << endl;
    compute_SNF(A);
    //A.print();
    //cout << endl;
    //compute_SNF(T);
    return A;
}
//pair<Matrix<int>,Matrix<int>> decompose(Matrix<int>& A) {
//    auto temp = A;
//    unsigned int n = A.GetSize();
////    auto B = A.transpose();
////    A = B*A;
//    Matrix<int> U(n); // Матрица, запоминающая действия со строками
//    Matrix<int> R(n); // Матрица, запоминающая действия со столбцами
//    // Сначала U - единичная матрица
//    // Она будет присоединена к исходной матрице
//    for (int i = 0; i < n; ++i){
//        U.val[i][i] = 1;
//        R.val[i][i] = 1;
//    }
//
//    // (A | U)
//    // U A R = D - почти диагональная
//    //Разрешённые действия (строки/столбцы a и b):
//    // a = a +k *b
//    // a = -a
//    // a < - > b (поменять местами)
//    //cout << "Matrix U on step "<< -1 <<": "<< endl << U ;
//
//    // Проход по столбцам - многомерный аналог алгоритма Евклида: сложность O(1)*n = O(n)
//    // Умножение вектора на скаляр - O(n)
//    // Итого n * O(n^2) = O(n^3)
//
//    // Пока нет квантовых компьютеров, можно лишь понижать константу
//
//    int num_nonzero = 10; // Число ненулевых элементов в текущем столбце в данный момент
//    int pos_min_nonzero; // Позиция минимального ненулевого элемента
//    int min_nz; // Минимальный по модулю ненулевой элемент в k-м столбце
//    int second_min_nz; // Второй минимум
//    int pos_sec_nz; // Позиция второго минимума
//
//    int m; // Здесь будет целая часть отношения в алгоритме Евклида (div в общем)
//    int multiplier;
//
//    int saved_zero;
//    for(int k = 0; k<n; k++){
//        auto A_T = A.transpose(); // Чтобы компу было удобнее ходить в память
//
////        TDynamicVector<int> col = A_T[k];// Берём k-й столбец
//        TDynamicVector<int> col(n);
//        for(int p=0;p<n;p++){
//            col[p] = A[p][k];
//        }
//        // Расширенный алгоритм Евклида над столбцом
//
//        while (true){
//            num_nonzero = 1;
//
//            //int j = 0;
//            // Нули нам явно не нужны. They are the good guys.
//            int j = k;
//            while (col[j] == 0)
//                j++;
//            // Нули в начале пройдены
//            min_nz = col[j];
//            pos_min_nonzero = j;
//
//            second_min_nz = col[j];
//            pos_sec_nz = j;
//            bool flag = false;
//            for (int i = j+1;i<n;i++) {
//                int t = col[i];
//                if(t!=0){
//                    if(abs(t) < abs(min_nz)){
//                        second_min_nz = min_nz;
//                        pos_sec_nz = pos_min_nonzero;
//
//                        min_nz = t;
//                        pos_min_nonzero = i;
//                    }
//                    else{
//                        if(!flag){
//                            second_min_nz = t;
//                            pos_sec_nz = i;
//                            flag = true;
//                        }
//                        if((t!=min_nz)&&(abs(t) < abs(second_min_nz))){
//                            second_min_nz = t;
//                            pos_sec_nz = i;
//                        }
//                    }
//                    num_nonzero ++;
//                }
//            }
//            //cout << col <<"|"<< min_nz<<","<< second_min_nz<< endl;
//            if(num_nonzero == 1) break; // Алгоритм Евклида гарантирует, что мы точно придём сюда
//            // Теперь совершаем алгоритм Евклида над двумя мин. элементами
////            cout << "current matrix: "<< endl;
////            A.print();
////            //cout << A[1] << endl;
////            cout << "also matrix: "<< endl;
////            (U * temp * R).print();
//            while ((min_nz!=0)&&(second_min_nz!=0)){
//                m  = floor(second_min_nz/min_nz);
//                A[pos_sec_nz] = A[pos_sec_nz] - A[pos_min_nonzero] * m;
//                U[pos_sec_nz] = U[pos_sec_nz] - U[pos_min_nonzero] * m;
//                second_min_nz = second_min_nz % min_nz ;
//                col[pos_sec_nz] = col[pos_sec_nz] % min_nz;
//                if(second_min_nz == 0)
//                    break;
//                m  = floor(min_nz/second_min_nz);
//                A[pos_min_nonzero] = A[pos_min_nonzero] - A[pos_sec_nz] * m;
//                U[pos_min_nonzero] = U[pos_min_nonzero] - U[pos_sec_nz] * m;
//                min_nz = min_nz % second_min_nz;
//                col[pos_min_nonzero] = col[pos_min_nonzero] % second_min_nz;
//            }
//            if(min_nz == 0){
//                min_nz = second_min_nz;
//                pos_min_nonzero = pos_sec_nz;
//            }
//            // Вертикальный ход
//            for (int l = 0;l<n;l++)
//                if((l != pos_min_nonzero)&&(col[l]!=0)){
//                    multiplier = floor(col[l]/min_nz);
//                    col[l] = col[l] - min_nz * multiplier;
//
//                    A[l] = A[l] - A[pos_min_nonzero] * multiplier;
//
//                    U[l] = U[l] - U[pos_min_nonzero] * multiplier;
//                }
//
//        }
////        cout << "current matrix: "<< endl;
////        A.print();
////        //cout << A[1] << endl;
////        cout << "also matrix: "<< endl;
////        (U * temp * R).print();
//        // Апосля того, как был зачищен стобец, меняем местами k-ю и min_pos строки
//        if(k!=n-1){
//            swap(A[k],A[pos_min_nonzero]);
//            swap(U[k],U[pos_min_nonzero]);
//        }
////        cout << "current matrix: "<< endl;
////        A.print();
////        //cout << A[1] << endl;
////        cout << "also matrix: "<< endl;
////        (U * temp * R).print();
//
//        //cout << A << endl;
//        // А теперь, пользуясь оставшимся ненулевым элементом, уменьшаем элементы строки
//        int el = A[k][k];
//        // Проходим по k-й строке
//        //auto T = A.transpose();
//        for (int i = k+1; i < n; i++) {
//            int& t = A[k][i];
//            if(abs(el) <= abs(t)){
////                A.print() ;
////                cout << "________" << R[i] << endl;
////                cout << floor(t/el) << endl;
////                cout << R[k] * floor(t/el) << endl;
//                R[i] = R[i] -  R[k] * floor(t/el);
//                t = t % el;
//                //T[i] = T[i] -  T[k] * floor(abs(t/el));
////                cout << R[i] << endl;
////                cout << R[i][k] << " | "<< R[k][k] << "|" << floor(t/el) << endl;
////                cout << "Получившийся элемент "<<R[i][k] << endl;
////                t = t % el;
////                cout << "Получившийся элемент "<<t << endl;
////                A.print() ;
//            }
//        }
//        //A = T.transpose();
////        cout << "current matrix: "<< endl;
////        A.print();
////        //cout << A[1] << endl;
////        cout << "also matrix: "<< endl;
////        (U * temp * R).print();
//        //cout << A << endl;
//    }
//    // Финальная полировка
//    for (int k = 0; k < n; k++) {
//        int el = A[k][k];
//        // Проходим по k-й строке
//        for (int i = k+1; i < n; i++) {
//            int& t = A[k][i];
//            if(t!=0)
//                if(abs(el) <= abs(t)){
//
//                    R[i] =R[i] -  R[k] * floor(t/el);
//                    t = t % el;
//                }
//        }
//    }
//    return {U,R.transpose()};
//}

// ax = 1 mod modulus
//int modulo_inverse(int a, int modulus)
//{
//    int m0 = modulus, t, q;
//    int x0 = 0, x1 = 1;
//
//    if (modulus == 1)
//    {
//        return 0;
//    }
//    while (a > 1) {
//        q = a / modulus;
//        t = modulus;
//        modulus = a % modulus, a = t;
//
//        t = x0;
//
//        x0 = x1 - q * x0;
//
//        x1 = t;
//    }
//
//    if (x1 < 0)
//    {
//        x1 += m0;
//    }
//
//    return x1;
//}

// Минимальное положительное x, удовлетворяющее системе из n сравнений:
// x = remainder_i mod number_i
//int chinese_remainder_theorem(const int* numbers, const int* remainders, int n)
//{
//    /*
//        Computing product of all numbers
//        from the num[] array
//    */
//    int prod = 1;
//    for (int i = 0; i < n; i++)
//    {
//        prod *= numbers[i];
//    }
//
//    // Initialize result
//    int result = 0;
//
//    // Apply the formula above
//    for (int i = 0; i < n; i++) {
//        int temp = prod / numbers[i];
//        result += remainders[i] * modulo_inverse(temp, numbers[i]) * temp;
//    }
//
//    return result % prod;
//}
//pair<Matrix<int>,Matrix<int>> solve_SLDE_mod_p(Matrix<int>& A, TDynamicVector<int>& b, int p){
//    unsigned int n = A.GetSize();
//    Matrix<int> U(n); // Матрица, запоминающая действия со строками
//    Matrix<int> R(n); // Матрица, запоминающая действия со столбцами
//    // Сначала U - единичная матрица
//    // Она будет присоединена к исходной матрице слева
//    for (int i = 0; i < n; ++i){
//        U.val[i][i] = 1;
//        R.val[i][i] = 1;
//    }
//    // O(n^2)
//    for (int i = 0; i < n; ++i) {
//        A[i] = A[i].mult_modulo(1,p);
//    }
//    // Идея такая: ищем минимальный по модулю остаток r (чтобы понизить константу)
//    // Затем находим к нему обратный по модулю m: r m = 1 mod p
//    // Пользуясь этим, зануляем остальные элементы: a r m  = a mod p . v - > v - (a * r * m) *v
//    // Так в каждом столбце
//    // Пройдя по столбцу, также пользуемся этим элементом, чтобы занулить элементы справа от него
//    int min, min_pos; // Минимальный ненулевой остаток в текущем столбце и его позиция
//    int inverse; // Обратный по модулю
//    for (int i = 0; i < n; ++i) {
//
//        // Текущий столбец
//        TDynamicVector<int> col(n);
//        for(int t=0;t<n;t++){
//            col[t] = A[t][i];
//        }
//        // Находим минимальный ненулевой остаток и его позицию
//        // Пока что проигнорировал случай, когда весь столбец нулевой
//        int j = 0;
//        while(col[j]==0){
//
//            j++;
//            if(j==n)
//                break;
//        }
//        if(j==n) continue;
//        min = col[j];
//        min_pos = j;
//        for (int k = j+1; k < n; k++) {
//            int temp = col[k];
//            if(temp!=0){
//                if(abs(temp)<abs(min)){
//                    min = temp;
//                    min_pos = k;
//                }
//            }
//        }
//        // Находим к нему обратный по модулю
//        inverse = modulo_inverse(min,p);
//
//        // Зануляем весь столбец
//        auto v = A[min_pos]; // Работаем этой строкой
//        auto v2 = U[min_pos];
//        for (int k = 0; k < n; ++k)
//            if(k!=min_pos){
//
//                int& t = col[k];
//                if(t!=0){
//                    //cout << A << endl;
//                    auto B = v.mult_modulo(-inverse *A[k][i],p);
//                    //cout << "It's adding "<< B<<endl<<"To "<<A[k] << endl;
//                    A[k] = A[k].add_modulo(B,p) ;
//                    U[k] = U[k].add_modulo(v2.mult_modulo(-inverse * U[k][i],p),p) ;
//                    //cout << A << endl;
//                }
//
//            }
//
//        // Зачищен столбец
//        // Проходим по строке A[min_pos] и зануляем элементы справа
//        for (int k = i+1; k < n; k++) {
//            int& t  = A[min_pos][k];
//            if(t!=0){
//                t = 0;
//                R[k] = R[k].add_modulo(R[i].mult_modulo(-inverse * t,p),p);
//            }
//        }
//        if(i!=n-1)
//            swap(A[i],A[min_pos]);
//
//    }
//    cout << A;
//    return {U,R.transpose()};
//}
//vector<int> SieveOfEratosthenes(int n)
//{
//    vector<int> result;
//    // Create a boolean array "prime[0..n]" and initialize
//    // all entries it as true. A value in prime[i] will
//    // finally be false if i is Not a prime, else true.
//    bool prime[n + 1];
//    std::memset(prime, true, sizeof(prime));
//
//    for (int p = 2; p * p <= n; p++) {
//        // If prime[p] is not changed, then it is a prime
//        if (prime[p]) {
//            // Update all multiples of p greater than or
//            // equal to the square of it numbers which are
//            // multiple of p and are less than p^2 are
//            // already been marked.
//            for (int i = p * p; i <= n; i += p)
//                prime[i] = false;
//        }
//    }
//
//    // Print all prime numbers
//    for (int p = 2; p <= n; p++)
//        if (prime[p])
//            result.push_back(p);
//    return result;
//}
//pair<Matrix<int>,Matrix<int>> solve_SLDE_modular_method(Matrix<int>& A, TDynamicVector<int>& b){
//    unsigned int n = A.GetSize();
//    // Во-первых, оцениваем, насколько большими могут стать элементы (a_i <= det A <= |A1| * |A2| * ... - Hadamard bound)
//    // [[ Этот шаг можно улучшить, воспользовавшись Леммой о двух векторах и определителе ]]
//    // В реальности в большинстве случайных матриц определитель мал
//    long HadamardBoundSq = 1;
//    for (int i = 0; i < n; ++i) {
//        HadamardBoundSq *= (A[i] * A[i]);
//    }
//    long double HadamardBound = sqrt(HadamardBoundSq);
//    auto primes = SieveOfEratosthenes(floor(HadamardBound));
//    unsigned int k = primes.size();
//    // По факту необходимое количество простых НАМНОГО МЕНЬШЕ!! Этот шаг можно оптимизировать
//    //Для каждого элемента матриц U и R создаётся массив остатков по модулю простых
//    // Во-вторых, решаем k * n * n систем по модулю p_i
//    // В-третьих, собираем всё вместе через CRT
//
//}
//Matrix<int> gen_public_key(int size,int range){
//    Matrix<int> m(size);
//    m.randomize(range);
//    return m;
//}
//vector<Matrix<int>> gen_private_key(Matrix<int>& public_key){
//    auto L = decompose(public_key);
//    vector<Matrix<int>> RES(3);
//    RES[0] = L.first;
//    RES[1] = public_key;
//    RES[2] = L.second;
//    return RES;
//}
//TDynamicVector<int> encrypt(const TDynamicVector<int>& message, Matrix<int>& public_key,short sigma){
//    auto length = message.size();
//    random_device rd;
//    mt19937_64 gen(rd());
//    uniform_int_distribution<> distr(-sigma,sigma);
//
//    TDynamicVector<int> r(length);
//    for (int i = 0;i<length;i++) {
//        r[i] = distr(gen);
//    }
//
//    return   public_key * message + r;
//}

//TDynamicVector<int> decrypt(const TDynamicVector<int>& encrypted_message, Matrix<int> U,Matrix<int> private_key,Matrix<int> R,short sigma){
//    // Быстро инвертируем приватный ключ. Это легко сделать, потому что векторы почти ортогональны
//    unsigned int s = private_key.GetSize();
//    // Сделаем из нашей матрицы матрицу из даблов
//    Matrix<double> tmp(s);
//    for (int i = 0; i < s; ++i)
//        for (int j = 0; j < s; ++j)
//            tmp[i][j] = private_key[i][j];
//    // Сделаем матрицу из даблов
//    Matrix<double> U_(s);
//    for (int i = 0; i < s; ++i)
//        for (int j = 0; j < s; ++j)
//            U_[i][j] = U[i][j];
//    // Сделаем матрицу из даблов
//    Matrix<double> R_(s);
//    for (int i = 0; i < s; ++i)
//        for (int j = 0; j < s; ++j)
//            R_[i][j] = R[i][j];
//    // Создадим единичную матрицу
//    Matrix<double> inverse(s);
//    for (int i = 0; i < s; ++i)
//        for (int j = 0; j < s; ++j)
//            inverse[i][j] = (i==j);
//    // Теперь посчитаем обратную
//    for (int i = 0; i < s; ++i) {
//        auto el = tmp[i][i];
//        tmp[i] = tmp[i] * (1/el);
//        inverse[i] = inverse[i] * (1/el);
//        for (int j = 0; j < i; ++j) {
//            auto p = tmp[j][i];
//            if(p!=0){
//                tmp[j] = tmp[j] - tmp[i] * p;
//                inverse[j] = inverse[j] - inverse[i] * p;
//            }
//        }
//    }
//    //inverse.print();
//    auto res_matrix =  R_ * inverse * U_;
//    //res_matrix.print();
//    size_t ms = encrypted_message.size();
//    TDynamicVector<double> encrypted_message_(ms);
//    for (int i = 0; i < ms;i++) {
//        encrypted_message_[i] = encrypted_message[i];
//    }
//    auto v_ = res_matrix * encrypted_message_;
//    cout << v_ << endl;
//    //cout <<"The resulting vector is "<< v_ << endl;
//    // Теперь используется алгоритм Бабая
//    // D x = v
//    auto x =  v_;
//    TDynamicVector<int> x_(x.size());
//    //round_to_nearest ;
//    for (int i = 0; i < x.size(); ++i) {
//        x_[i] = lround(x[i]);
//    }
//    auto decrypted_message = private_key * x_;
//    return decrypted_message;
//    }