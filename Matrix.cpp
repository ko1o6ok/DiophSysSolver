

#include "Matrix.h"
#include <numeric>
#include<cstdlib>
#include <cstring>
#include <algorithm>
using namespace std;
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

double Matrix::entropy() {
    unsigned int r = N*N;
    // Расплющим матрицу в вектор
    vector<int> v(r);
    vector<int> v_old(r);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            v[N*i+j] = val[i][j];
        }
    }
    vector<int>::iterator ip;
    v_old = v;
    // Sorting the array
    sort(v.begin(), v.end());
    // Now v becomes 1 1 2 2 3 3 3 3 7 7 8 10

    // Using std::unique
    ip = unique(v.begin(), v.begin() + v.size());
    // Now v becomes {1 2 3 7 8 10 * * * * * *}
    // * means undefined

    // Resizing the vector so as to remove the undefined terms
    v.resize(std::distance(v.begin(), ip));
//    TDynamicVector<int> n_v(v);
//    cout << n_v;
    double H = 0;
    for(auto& el:v){
        double c = count(v_old.begin(),v_old.end(),el);

        double p = c/r;
        H += -p*log(p);
    }
    return H;
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
        cout << A.entropy() << endl;
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

        // Апосля того, как был зачищен столбец, меняем местами k-ю и min_pos строки
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
        std::cout << A.entropy() << endl;
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

    // Пока нет квантовых компьютеров, можно лишь понижать константу

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
        //std::cout << A.entropy() << endl;
    }

}
int choose_row(Matrix& A, int k){
    int metric; // Максимальное значение второго минимального
    unsigned int s = A.GetSize();
    int res; // рекомендуемая строка
    for (int i = k; i < s; ++i) {
        auto row = A[i];
        // Нули нам явно не нужны. They are the good guys.
        int j = k;
        while (row[j] == 0)
            j++;
        // Нули в начале пройдены
        int min_nz = row[j];
        int second_min_nz = row[j];
        bool flag = false;
        for (int t = j+1;t<s;t++) {
            int temp = row[i];
            if(t!=0){
                if(abs(temp) < abs(min_nz)){
                    second_min_nz = min_nz;
                    min_nz = temp;
                }
                else{
                    if(!flag){
                        second_min_nz = temp;
                        flag = true;
                    }
                    if((temp!=min_nz)&&(abs(temp) < abs(second_min_nz))){
                        second_min_nz = temp;
                    }
                }
            }
        }
        if((i == k)||( abs(metric)<abs(second_min_nz) )){
            metric = second_min_nz;
            res = i;
        }
    }
    return res;
}
void solve_SLDE_v4(Matrix& A, TDynamicVector<int> &b) {
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
    // По сравнению с предыд. алгоритмом мы здесь столбец ВЫБИРАЕМ. Так, чтобы два минимальных элемента были минимальными!
    for(int k = 0; k<n; k++){
        auto A_T = A.transpose(); // Чтобы компу было удобнее ходить в память
        // Вот здесь мы сначала проходимся по матрице и выбираем наиболее удобный столбец (как улучшение можно искать столбец ЛИБО строку
        // и производить соответствующие действия)
        // Естественно, столбцы смотрим, лишь начиная с k-го
        int chosen = choose_row(A_T,k); // Выбираем строку трансп. матрицы = столбец исходной!
        // SWAP-аем понравившийся с k-м
        swap(A_T[chosen],A_T[k]);
        A = A_T.transpose();
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

        // Апосля того, как был зачищен столбец, меняем местами k-ю и min_pos строки
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
        //std::cout << A.entropy() << endl;
    }

}
//void Matrix::blind_entropy_shift(){
//    double entr = entropy(); // Текущая энтропия
//    // На каждом шаге выполняем случайное действие, чтобы уменьшилась энтропия
//    // Пока на оптимизацию времени я забью, сделаем неоптимальный вариант
//    for (int i = 0; i < 7; ++i) {
//        Matrix B(*this); // Копируем
//        for(int k=-10;k<10;k++){
//            B[i] = B[i] + B[(i+1)%7] * k;
//            auto t = B.entropy();
//            if( t< entr){
//                entr = t;
//                *this = B;
//                break;
//            }
//        }
//        cout << *this << endl<< "Entropy ="<<entr<< endl;
//    }
//
//}
// ax = 1 mod modulus
int modulo_inverse(int a, int modulus)
{
    int m0 = modulus, t, q;
    int x0 = 0, x1 = 1;

    if (modulus == 1)
    {
        return 0;
    }
    while (a > 1) {
        q = a / modulus;
        t = modulus;
        modulus = a % modulus, a = t;

        t = x0;

        x0 = x1 - q * x0;

        x1 = t;
    }

    if (x1 < 0)
    {
        x1 += m0;
    }

    return x1;
}

// Минимальное положительное x, удовлетворяющее системе из n сравнений:
// x = remainder_i mod number_i
int chinese_remainder_theorem(const int* numbers, const int* remainders, int n)
{
    /*
        Computing product of all numbers
        from the num[] array
    */
    int prod = 1;
    for (int i = 0; i < n; i++)
    {
        prod *= numbers[i];
    }

    // Initialize result
    int result = 0;

    // Apply above formula
    for (int i = 0; i < n; i++) {
        int temp = prod / numbers[i];
        result += remainders[i] * modulo_inverse(temp, numbers[i]) * temp;
    }

    return result % prod;
}
pair<Matrix,Matrix> solve_SLDE_mod_p(Matrix& A, TDynamicVector<int>& b, int p){
    unsigned int n = A.GetSize();
    Matrix U(n); // Матрица, запоминающая действия со строками
    Matrix R(n); // Матрица, запоминающая действия со столбцами
    // Сначала U - единичная матрица
    // Она будет присоединена к исходной матрице слева
    for (int i = 0; i < n; ++i){
        U.val[i][i] = 1;
        R.val[i][i] = 1;
    }
    // O(n^2)
    for (int i = 0; i < n; ++i) {
        A[i] = A[i].mult_modulo(1,p);
    }
    // Идея такая: ищем минимальный по модулю остаток r (чтобы понизить константу)
    // Затем находим к нему обратный по модулю m: r m = 1 mod p
    // Пользуясь этим, зануляем остальные элементы: a r m  = a mod p . v - > v - (a * r * m) *v
    // Так в каждом столбце
    // Пройдя по столбцу, также пользуемся этим элементом, чтобы занулить элементы справа от него
    int min, min_pos; // Минимальный ненулевой остаток в текущем столбце и его позиция
    int inverse; // Обратный по модулю
    for (int i = 0; i < n; ++i) {

        // Текущий столбец
        TDynamicVector<int> col(n);
        for(int t=0;t<n;t++){
            col[t] = A[t][i];
        }
        // Находим минимальный ненулевой остаток и его позицию
        // Пока что проигнорировал случай, когда весь столбец нулевой
        int j = 0;
        while(col[j]==0){

            j++;
            if(j==n)
                break;
        }
        if(j==n) continue;
        min = col[j];
        min_pos = j;
        for (int k = j+1; k < n; k++) {
            int temp = col[k];
            if(temp!=0){
                if(abs(temp)<abs(min)){
                    min = temp;
                    min_pos = k;
                }
            }
        }
        // Находим к нему обратный по модулю
        inverse = modulo_inverse(min,p);

        // Зануляем весь столбец
        auto v = A[min_pos]; // Работаем этой строкой
        auto v2 = U[min_pos];
        for (int k = 0; k < n; ++k)
            if(k!=min_pos){

                int& t = col[k];
                if(t!=0){
                    //cout << A << endl;
                    auto B = v.mult_modulo(-inverse *A[k][i],p);
                    //cout << "It's adding "<< B<<endl<<"To "<<A[k] << endl;
                    A[k] = A[k].add_modulo(B,p) ;
                    U[k] = U[k].add_modulo(v2.mult_modulo(-inverse * U[k][i],p),p) ;
                    //cout << A << endl;
                }

            }

        // Зачищен столбец
        // Проходим по строке A[min_pos] и зануляем элементы справа
        for (int k = i+1; k < n; k++) {
            int& t  = A[min_pos][k];
            if(t!=0){
                t = 0;
                R[k] = R[k].add_modulo(R[i].mult_modulo(-inverse * t,p),p);
            }
        }
        if(i!=n-1)
            swap(A[i],A[min_pos]);

    }
    cout << A;
    return {U,R.transpose()};
}
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

vector<int> SieveOfEratosthenes(int n)
{
    vector<int> result;
    // Create a boolean array "prime[0..n]" and initialize
    // all entries it as true. A value in prime[i] will
    // finally be false if i is Not a prime, else true.
    bool prime[n + 1];
    std::memset(prime, true, sizeof(prime));

    for (int p = 2; p * p <= n; p++) {
        // If prime[p] is not changed, then it is a prime
        if (prime[p]) {
            // Update all multiples of p greater than or
            // equal to the square of it numbers which are
            // multiple of p and are less than p^2 are
            // already been marked.
            for (int i = p * p; i <= n; i += p)
                prime[i] = false;
        }
    }

    // Print all prime numbers
    for (int p = 2; p <= n; p++)
        if (prime[p])
            result.push_back(p);
    return result;
}
pair<Matrix,Matrix> solve_SLDE_modular_method(Matrix& A, TDynamicVector<int>& b){
    unsigned int n = A.GetSize();
    // Во-первых, оцениваем, насколько большими могут стать элементы (a_i <= det A <= |A1| * |A2| * ... - Hadamard bound)
    // [[ Этот шаг можно улучшить, воспользовавшись Леммой о двух векторах и определителе ]]
    // В реальности в большинстве случайных матриц определитель мал
    long HadamardBoundSq = 1;
    for (int i = 0; i < n; ++i) {
        HadamardBoundSq *= (A[i] * A[i]);
    }
    long double HadamardBound = sqrt(HadamardBoundSq);
    auto primes = SieveOfEratosthenes(floor(HadamardBound));
    unsigned int k = primes.size();
    // По факту необходимое количество простых НАМНОГО МЕНЬШЕ!! Этот шаг можно оптимизировать
    //Для каждого элемента матриц U и R создаётся массив остатков по модулю простых
    // Во-вторых, решаем k * n * n систем по модулю p_i
    // В-третьих, собираем всё вместе через CRT

}