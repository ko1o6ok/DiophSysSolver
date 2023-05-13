#include<cstdlib>
#include <cmath>
#include "Matrix.h"
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
            if (min_nz == 0) {
                min_nz = second_min_nz;
                pos_min_nonzero = pos_sec_nz;
            }
            for (int l = 0; l < n; l++)
                if ((l != pos_min_nonzero) && (col[l] != 0)) {
                    multiplier = floor(col[l] / min_nz);
                    col[l] = col[l] - min_nz * multiplier;
                    A[l] = A[l] - A[pos_min_nonzero] * multiplier;
                }
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
    // Я так делаю потому, что мой алгоритм приводит почти к SNF
    compute_SNF(A);
    A =  A.transpose();
    compute_SNF(A);
    return A;
}