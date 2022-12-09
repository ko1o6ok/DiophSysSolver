

#include "Matrix.h"
#include <numeric>
#include<cstdlib>
Matrix::Matrix(unsigned int n) {
    N = n;
    val = TDynamicVector<TDynamicVector<int>>(N);
    for (int i = 0; i < N; ++i)
            val[i] = TDynamicVector<int>(N);
}

Matrix::Matrix(const TDynamicVector<TDynamicVector<int>>& v) {
    N = v.size();
    val = v;
}

TDynamicVector<int>& Matrix::operator[](unsigned int i) {
    return val[i];
}

TDynamicVector<int> Matrix::operator*(const TDynamicVector<int>& vec) {
    // Предполагается, что размер vec равен N
    TDynamicVector<int> res(N);
    for (int i = 0; i < N; ++i) {
        res[i] = val[i] * vec;
    }
    return res;
}

ostream &operator<<(ostream &ostr, Matrix& m) {
    for (int i = 0; i <m.N; ++i) {
        ostr << m[i] ;
        ostr << endl;
    }
    return ostr;
}

unsigned int Matrix::GetSize() const {
    return N;
}
// Адаптивный алгоритм Штрассена
// pivot point = 1_000
Matrix Matrix::operator*(Matrix& m) {
//    if(N > 1000){
//        Matrix a(multiply_matrix(val,m.val));
//        return a;
//    }
    Matrix help = m.transpose();
    Matrix r(N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            r[i][j] = val[i] * help[j];
    return r;
}



vector<vector<int> >
Matrix::add_matrix(vector<vector<int> > matrix_A,
           vector<vector<int> > matrix_B, int split_index,
           int multiplier = 1)
{
    for (auto i = 0; i < split_index; i++)
        for (auto j = 0; j < split_index; j++)
            matrix_A[i][j]
                    = matrix_A[i][j]
                      + (multiplier * matrix_B[i][j]);
    return matrix_A;
}
// Strassen's algorithm
vector<vector<int> >
Matrix::multiply_matrix(vector<vector<int> > matrix_A,
                vector<vector<int> > matrix_B)
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

    vector<int> result_matrix_row(col_2, 0);
    vector<vector<int> > result_matrix(row_1,
                                       result_matrix_row);

    if (col_1 == 1)
        result_matrix[0][0]
                = matrix_A[0][0] * matrix_B[0][0];
    else {
        int split_index = col_1 / 2;

        vector<int> row_vector(split_index, 0);

        vector<vector<int> > a00(split_index, row_vector);
        vector<vector<int> > a01(split_index, row_vector);
        vector<vector<int> > a10(split_index, row_vector);
        vector<vector<int> > a11(split_index, row_vector);
        vector<vector<int> > b00(split_index, row_vector);
        vector<vector<int> > b01(split_index, row_vector);
        vector<vector<int> > b10(split_index, row_vector);
        vector<vector<int> > b11(split_index, row_vector);

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

        vector<vector<int> > p(multiply_matrix(
                a00, add_matrix(b01, b11, split_index, -1)));
        vector<vector<int> > q(multiply_matrix(
                add_matrix(a00, a01, split_index), b11));
        vector<vector<int> > r(multiply_matrix(
                add_matrix(a10, a11, split_index), b00));
        vector<vector<int> > s(multiply_matrix(
                a11, add_matrix(b10, b00, split_index, -1)));
        vector<vector<int> > t(multiply_matrix(
                add_matrix(a00, a11, split_index),
                add_matrix(b00, b11, split_index)));
        vector<vector<int> > u(multiply_matrix(
                add_matrix(a01, a11, split_index, -1),
                add_matrix(b10, b11, split_index)));
        vector<vector<int> > v(multiply_matrix(
                add_matrix(a00, a10, split_index, -1),
                add_matrix(b00, b01, split_index)));

        vector<vector<int> > result_matrix_00(add_matrix(
                add_matrix(add_matrix(t, s, split_index), u,
                           split_index),
                q, split_index, -1));
        vector<vector<int> > result_matrix_01(
                add_matrix(p, q, split_index));
        vector<vector<int> > result_matrix_10(
                add_matrix(r, s, split_index));
        vector<vector<int> > result_matrix_11(add_matrix(
                add_matrix(add_matrix(t, p, split_index), r,
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

Matrix Matrix::transpose() {
    Matrix r(N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            r.val[i][j] = val[j][i];
        }
    }
    return r;
}

void Matrix::randomize(unsigned int range) {
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
void solve_SLDE_v1(Matrix& A, TDynamicVector<int> &b) {
    unsigned int n = A.GetSize();
    Matrix U(n); // Матрица, запоминающая действия со строками
    // Сначала U - единичная матрица
    // Она будет присоединена к исходной матрице
    for (int i = 0; i < n; ++i){
        U.val[i][i] = 1;
    }

    // (A | U)
    // U A = H - верхне-треуольная
    //Разрешённые действия (строки a и b):
    // a = a +k *b
    // a = -a
    // a < - > b (поменять местами)
    //cout << "Matrix U on step "<< -1 <<": "<< endl << U ;

    // Проход по столбцам - многомерный аналог алгоритма Евклида: сложность O(1)*n = O(n)
    // Умножение вектора на скаляр - O(n)
    // Итого n * O(n^2) = O(n^3)

    unsigned int num_nonzero = 10; // Число ненулевых элементов в текущем столбце в данный момент
    unsigned int pos_min_nonzero; // Позиция минимального ненулевого элемента
    int min_nz; // Минимальный по модулю ненулевой элемент в k-м столбце

    for(int k = 0; k<n; k++){
        auto A_T = A.transpose(); // Чтобы компу было удобнее ходить в память

//        TDynamicVector<int> col = A_T[k];// Берём k-й столбец
        TDynamicVector<int> col(n);
        for(int p=0;p<n;p++){
            col[p] = A[p][k];
        }
        // Расширенный алгоритм Евклида над столбцом

        while (true){
            num_nonzero = 1;

            //int j = 0;
            // Нули нам явно не нужны. They are the good guys.
            int j = k;
            while (col[j] == 0)
                j++;
            // Нули в начале пройдены
            min_nz = col[j];
            pos_min_nonzero = j;
            for (int i = j+1;i<n;i++) {
                int t = col[i];
                if(t!=0){
                    if(abs(t) < abs(min_nz)){
                        min_nz = t;
                        pos_min_nonzero = i;
                    }
                    num_nonzero ++;
                }
            }

            if(num_nonzero == 1) break; // Алгоритм Евклида гарантирует, что мы точно придём сюда

            // Вертикальный ход
            for (int l = 0;l<n;l++)
                if((l != pos_min_nonzero)&&(col[l]!=0)){

                    int multiplier = floor(col[l]/min_nz);

                    col[l] = col[l] % min_nz;

                    A[l] = A[l] - A[pos_min_nonzero] * multiplier;

                    U[l] = U[l] - U[pos_min_nonzero] * multiplier;
                }

        }

        // Апосля того, как был зачищен стобец, меняем местами k-ю и min_pos строки
        if(k!=n-1)
            swap(A[k],A[pos_min_nonzero]);
    }
//    cout << "RES: "<<endl << A;
//    cout << "U is:" << endl << U;
}

void solve_SLDE_v2(Matrix& A, TDynamicVector<int> &b) {
    unsigned int n = A.GetSize();
//    auto B = A.transpose();
//    A = B*A;
    Matrix U(n); // Матрица, запоминающая действия со строками
    Matrix R(n); // Матрица, запоминающая действия со столбцами
    // Сначала U - единичная матрица
    // Она будет присоединена к исходной матрице
    for (int i = 0; i < n; ++i){
        U.val[i][i] = 1;
        R.val[i][i] = 1;
    }

    // (A | U)
    // U A R = D - почти диагональная
    //Разрешённые действия (строки/столбцы a и b):
    // a = a +k *b
    // a = -a
    // a < - > b (поменять местами)
    //cout << "Matrix U on step "<< -1 <<": "<< endl << U ;

    // Проход по столбцам - многомерный аналог алгоритма Евклида: сложность O(1)*n = O(n)
    // Умножение вектора на скаляр - O(n)
    // Итого n * O(n^2) = O(n^3)

    unsigned int num_nonzero = 10; // Число ненулевых элементов в текущем столбце в данный момент
    unsigned int pos_min_nonzero; // Позиция минимального ненулевого элемента
    int min_nz; // Минимальный по модулю ненулевой элемент в k-м столбце

    for(int k = 0; k<n; k++){
        auto A_T = A.transpose(); // Чтобы компу было удобнее ходить в память

//        TDynamicVector<int> col = A_T[k];// Берём k-й столбец
        TDynamicVector<int> col(n);
        for(int p=0;p<n;p++){
            col[p] = A[p][k];
        }
        // Расширенный алгоритм Евклида над столбцом

        while (true){
            num_nonzero = 1;

            //int j = 0;
            // Нули нам явно не нужны. They are the good guys.
            int j = k;
            while (col[j] == 0)
                j++;
            // Нули в начале пройдены
            min_nz = col[j];
            pos_min_nonzero = j;
            for (int i = j+1;i<n;i++) {
                int t = col[i];
                if(t!=0){
                    if(abs(t) < abs(min_nz)){
                        min_nz = t;
                        pos_min_nonzero = i;
                    }
                    num_nonzero ++;
                }
            }

            if(num_nonzero == 1) break; // Алгоритм Евклида гарантирует, что мы точно придём сюда

            // Вертикальный ход
            for (int l = 0;l<n;l++)
                if((l != pos_min_nonzero)&&(col[l]!=0)){

                    int multiplier = floor(col[l]/min_nz);

                    col[l] = col[l] % min_nz;

                    A[l] = A[l] - A[pos_min_nonzero] * multiplier;

                    U[l] = U[l] - U[pos_min_nonzero] * multiplier;
                }

        }

        // Апосля того, как был зачищен стобец, меняем местами k-ю и min_pos строки
        if(k!=n-1)
            swap(A[k],A[pos_min_nonzero]);
        // А теперь, пользуясь оставшимся ненулевым элементом, уменьшаем элементы строки
        int el = A[k][k];
        // Проходим по k-й строке
        for (int i = k+1; i < n; i++) {
            int& t = A[k][i];
            if(abs(el) < abs(t)){

                R[i] =R[i] -  R[k] * floor(t/el);
                t = t % el;
            }
        }
    }
//    cout << "RES: "<<endl << A;
//    cout << "U is:" << endl << U;
//    auto realR = R.transpose();
//    cout << "R is" << endl << realR;
}
void solve_SLDE_v3(Matrix& A, TDynamicVector<int> &b) {
    unsigned int n = A.GetSize();
//    auto B = A.transpose();
//    A = B*A;
    Matrix U(n); // Матрица, запоминающая действия со строками
    Matrix R(n); // Матрица, запоминающая действия со столбцами
    // Сначала U - единичная матрица
    // Она будет присоединена к исходной матрице
    for (int i = 0; i < n; ++i){
        U.val[i][i] = 1;
        R.val[i][i] = 1;
    }

    // (A | U)
    // U A R = D - почти диагональная
    //Разрешённые действия (строки/столбцы a и b):
    // a = a +k *b
    // a = -a
    // a < - > b (поменять местами)
    //cout << "Matrix U on step "<< -1 <<": "<< endl << U ;

    // Проход по столбцам - многомерный аналог алгоритма Евклида: сложность O(1)*n = O(n)
    // Умножение вектора на скаляр - O(n)
    // Итого n * O(n^2) = O(n^3)

    int num_nonzero = 10; // Число ненулевых элементов в текущем столбце в данный момент
    int pos_min_nonzero; // Позиция минимального ненулевого элемента
    int min_nz; // Минимальный по модулю ненулевой элемент в k-м столбце
    int second_min_nz; // Второй минимум
    int pos_sec_nz; // Позиция второго минимума

    int m; // Здесь будет целая часть отношения в алгоритме Евклида (div в общем)
    int multiplier;

    int saved_zero;
    for(int k = 0; k<n; k++){
        auto A_T = A.transpose(); // Чтобы компу было удобнее ходить в память

//        TDynamicVector<int> col = A_T[k];// Берём k-й столбец
        TDynamicVector<int> col(n);
        for(int p=0;p<n;p++){
            col[p] = A[p][k];
        }
        // Расширенный алгоритм Евклида над столбцом

        while (true){
            num_nonzero = 1;

            //int j = 0;
            // Нули нам явно не нужны. They are the good guys.
            int j = k;
            while (col[j] == 0)
                j++;
            // Нули в начале пройдены
            min_nz = col[j];
            pos_min_nonzero = j;

            second_min_nz = col[j];
            pos_sec_nz = j;
            bool flag = false;
            for (int i = j+1;i<n;i++) {
                int t = col[i];
                if(t!=0){
                    if(abs(t) < abs(min_nz)){
                        second_min_nz = min_nz;
                        pos_sec_nz = pos_min_nonzero;

                        min_nz = t;
                        pos_min_nonzero = i;
                    }
                    else{
                        if(!flag){
                            second_min_nz = t;
                            pos_sec_nz = i;
                            flag = true;
                        }
                        if((t!=min_nz)&&(abs(t) < abs(second_min_nz))){
                            second_min_nz = t;
                            pos_sec_nz = i;
                        }
                    }
                    num_nonzero ++;
                }
            }
            //cout << col <<"|"<< min_nz<<","<< second_min_nz<< endl;
            if(num_nonzero == 1) break; // Алгоритм Евклида гарантирует, что мы точно придём сюда
            // Теперь совершаем алгоритм Евклида над двумя мин. элементами
            while ((min_nz!=0)&&(second_min_nz!=0)){
                m  = floor(second_min_nz/min_nz);
                A[pos_sec_nz] = A[pos_sec_nz] - A[pos_min_nonzero] * m;
                U[pos_sec_nz] = U[pos_sec_nz] - U[pos_min_nonzero] * m;
                second_min_nz = second_min_nz % min_nz ;
                col[pos_sec_nz] = col[pos_sec_nz] % min_nz;
                if(second_min_nz == 0)
                    break;
                m  = floor(min_nz/second_min_nz);
                A[pos_min_nonzero] = A[pos_min_nonzero] - A[pos_sec_nz] * m;
                U[pos_min_nonzero] = U[pos_min_nonzero] - U[pos_sec_nz] * m;
                min_nz = min_nz % second_min_nz;
                col[pos_min_nonzero] = col[pos_min_nonzero] % second_min_nz;
            }
            if(min_nz == 0){
                min_nz = second_min_nz;
                pos_min_nonzero = pos_sec_nz;
            }
            // Вертикальный ход
            for (int l = 0;l<n;l++)
                if((l != pos_min_nonzero)&&(col[l]!=0)){
//                    double s = col[l]/min_nz;
//                    if(s>0){
//                        multiplier = floor(s);
//                    }
//                    else{
//                        multiplier = -floor(-s);
//                    }
                    multiplier = floor(col[l]/min_nz);
                    col[l] = col[l] - min_nz * multiplier;

                    A[l] = A[l] - A[pos_min_nonzero] * multiplier;

                    U[l] = U[l] - U[pos_min_nonzero] * multiplier;
                }
        }

        // Апосля того, как был зачищен стобец, меняем местами k-ю и min_pos строки
        if(k!=n-1)
            swap(A[k],A[pos_min_nonzero]);
        // А теперь, пользуясь оставшимся ненулевым элементом, уменьшаем элементы строки
        int el = A[k][k];
        // Проходим по k-й строке
        for (int i = k+1; i < n; i++) {
            int& t = A[k][i];
            if(abs(el) < abs(t)){

                R[i] =R[i] -  R[k] * floor(t/el);
                t = t % el;
            }
        }
        //cout << A << endl;
    }

}
//void solve_SLDE_v3(Matrix& A, TDynamicVector<int> &b) {
//    // Находим два мин. элемента, затем производим над ними алг. Евклида и берём оставшееся ненулевое число как минимальное
//
//    unsigned int n = A.GetSize();
////    auto B = A.transpose();
////    A = B*A;
//    Matrix U(n); // Матрица, запоминающая действия со строками
//    Matrix R(n); // Матрица, запоминающая действия со столбцами
//    // Сначала U - единичная матрица
//    // Она будет присоединена к исходной матрице
//    for (int i = 0; i < n; ++i) {
//        U.val[i][i] = 1;
//        R.val[i][i] = 1;
//    }
//
//
//
//    unsigned int num_nonzero = 10; // Число ненулевых элементов в текущем столбце в данный момент
//    int pos_min_nonzero1; // Позиция первого минимального по модулю ненулевого элемента
//    int pos_min_nonzero2; // Позиция второго минимального по модулю ненулевого элемента
//    int pos_min_nonzero; // Позиция итогового элемента, полученного из них алгоритмом Евклида
//    // Теперь находим два минимальных элемента, из них алгоритмом Евклида делаем ещё более маленький, а затем уже им работаем
//
//    int min_nz1; // Первый минимальный по модулю ненулевой элемент в k-м столбце
//    int min_nz2; // Второй минимальный по модулю ненулевой элемент в k-м столбце
//
//    int min_nz ; // Полученный из них алгоритмом Евклида маленький элемент
//    // Проход по столбцам
//    for (int k = 0; k < n; k++) {
//        auto A_T = A.transpose(); // Чтобы компу было удобнее ходить в память
//
////        TDynamicVector<int> col = A_T[k];// Берём k-й столбец
//        TDynamicVector<int> col(n);
//        for (int p = 0; p < n; p++) {
//            col[p] = A[p][k];
//        }
//        // Расширенный алгоритм Евклида над столбцом
//        while (true) {
//            //cout << num_nonzero << endl;
//            num_nonzero = 0;
//
//            //int j = 0;
//            // Нули нам явно не нужны. They are the good guys.
////            int j = k;
////            while (col[j] == 0)
////                j++;
//            // Нули в начале пройдены
//            // Можно оптимизировать!
//            for (int i = k; i < n; i++) {
//                int t = col[i];
//                if (t != 0) {
//                    num_nonzero++;
//                }
//            }
//            //cout << col << "| "<< num_nonzero << endl;
//            if (num_nonzero == 1) break;
//
//            auto pr = col.pos_two_min_nonzero_el_s();
//            pos_min_nonzero1 = pr.first;
//            pos_min_nonzero2 = pr.second;
//            min_nz1 = col[pos_min_nonzero1];
//            min_nz2 = col[pos_min_nonzero2];
//
//            //col[pos_min_nonzero1] = min_nz1;
//            //col[pos_min_nonzero2] = min_nz2;
//
//            min_nz = min_nz1;
//            if(min_nz == 0)
//                cout << "WTF"<< endl;
//            pos_min_nonzero = pos_min_nonzero1;
//
//            // Вертикальный ход
//            for (int l = 0; l < n; l++)
//                if ((l != pos_min_nonzero) && (col[l] != 0)) {
//
//                    int multiplier = floor(col[l] / min_nz);
//
//                    col[l] = col[l] % min_nz;
//
//                    A[l] = A[l] - A[pos_min_nonzero] * multiplier;
//
//                    U[l] = U[l] - U[pos_min_nonzero] * multiplier;
//                }
//
//        }
//
//        // Апосля того, как был зачищен стобец, меняем местами k-ю и min_pos строки
//        if (k != n - 1)
//            swap(A[k], A[pos_min_nonzero]);
//        // А теперь, пользуясь оставшимся ненулевым элементом, уменьшаем элементы строки
//        int el = A[k][k];
//
//        // Проходим по k-й строке
//        for (int i = k + 1; i < n; i++) {
//            int &t = A[k][i];
//            if (abs(el) < abs(t)) {
//
//                R[i] = R[i] - R[k] * floor(t / el);
//                t = t % el;
//            }
//        }
//    cout << A <<endl<<"-------------------------"<< endl;
//    }
//}






//void solve_SLDE(Matrix& A, TDynamicVector<int> &b) {
//    unsigned int n = A.GetSize();
//    Matrix U(n);
//    // Сначала U - единичная матрица
//    // Она будет присоединена к исходной матрице
//    for (int i = 0; i < n; ++i)
//        U.val[i][i] = 1;
//    // (A | U)
//    // U A = H - верхне-треуольная
//    //Разрешённые действия (строки a и b):
//    // a = a +k *b
//    // a = -a
//    // a < - > b (поменять местами)
//    //auto A_T = A.transpose();
//    // Проход по столбцам - многомерный алгоритм Евклида
//    for(int i = 0; i<n; i++){
//        TDynamicVector<int> col ;
//        unsigned int min_pos ;
//        int min_val;
//
//        // a % b = a - floor(a/b) * b;
//        unsigned int number_of_nonzero = 10; // Рандомное значение, не равное единице. Чтобы цикл запустился
//        while (number_of_nonzero != 1){
//            col = A.transpose()[i];
//            min_pos = col.find_min_nonzero();
//            min_val = col[min_pos];
//            for (int j = 0; j < n; ++j){
//                number_of_nonzero = 0;
//                if((col[j]!= 0 )&&(j != min_pos)){
//                    number_of_nonzero ++;
//                    int multiplier = floor(col[j]/min_val);
//                    A[j] = A[j] - A[min_pos] * multiplier;
//                    U[j] = U[j] - U[min_pos] * multiplier;
//                }
//            }
//        }
//        // Апосля того, как был зачищен стобец, меняем местами i-ю и min_pos строки
//        swap(A[i],A[min_pos]);
//    }
//}
int NdimGCD(vector<int>& vals){
    int n = vals.size();
    unsigned int num_nonzero = 10;
    unsigned int pos_min_nonzero;
    while (num_nonzero != 1){
        num_nonzero = 1;

        int j = 0;
        while (vals[j] == 0)
            j++;
        int min_nz = vals[j];
        pos_min_nonzero = j;
        for (int i = j+1;i<n;i++) {
            int t = vals[i];
            if(t!=0){
                if(t < min_nz){
                    min_nz = t;
                    pos_min_nonzero = i;
                }

                num_nonzero ++;
            }
        }
        for (int l = 0;l<n;l++)
            if(l != pos_min_nonzero){
                vals[l] = vals[l] % min_nz;
            }
//        cout << "Current state: ";
//        for(int& t:vals)
//            cout << t << ", ";
//        cout << endl;
    }
    int S = 0;
    for(int& t: vals)
        S+=t;
    return S;
}

//pair<Matrix, Matrix> Matrix::FHNF() {
//    Matrix U(N);
//    // Сначала U - единичная матрица
//    // Она будет присоединена к исходной матрице
//    for (int i = 0; i < N; ++i)
//        U.val[i][i] = 1;
//    Matrix H(N);
//    return {U,H};
//}
