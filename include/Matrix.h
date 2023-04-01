#ifndef FASTHNF_MATRIX_H
#define FASTHNF_MATRIX_H
#include <vector>
#include <iostream>
#include <random>
using namespace std;

// Динамический вектор -
// шаблонный вектор на динамической памяти
template<typename T>
class TDynamicVector
{

protected:
    size_t sz{};
    T* pMem;
public:
    TDynamicVector(T* arr, size_t s) : sz(s)
    {
        assert(arr != nullptr && "TDynamicVector ctor requires non-nullptr arg");
        pMem = new T[sz];
        std::copy(arr, arr + sz, pMem);
    }
    explicit TDynamicVector(size_t size = 1) : sz(size)
    {
        if (sz <= 0)
            throw length_error("Vector size should be greater than zero");
        pMem = new T[sz]();// {}; // У типа T д.б. констуктор по умолчанию
    }
    explicit TDynamicVector(vector<T> v){
        sz = v.size();
        pMem = new T[sz]();
        for (int i = 0; i < sz; ++i)
            pMem[i] = v[i];
    }
    TDynamicVector(const TDynamicVector& v) // конструктор копирования
    {
        sz = v.sz;
        pMem = new T[sz](); // Выделили свежую
        for (int i = 0; i < sz; ++i)
            pMem[i] = v.pMem[i];
    }
    TDynamicVector(TDynamicVector&& v) noexcept // Знаменитый вектор векторов
    {
        sz = v.sz;
        if (sz == 0)
            throw length_error("Vector size should be greater than zero");
        pMem = new T[sz](); // Это перепроверить
        for (int i = 0; i < sz; ++i)
            pMem[i] = v.pMem[i];
    }
    ~TDynamicVector() // Очистка памяти. ВАЖНО!
    {
        delete[] pMem;
    }
    TDynamicVector& operator=(const TDynamicVector& v) // Оператор присваивания. Добавить обработку самоприсваивания!!!
    {
        if (this != &v){
            delete[] pMem;
            sz = v.sz;
            pMem = new T[sz]();
            for (int i = 0; i < sz; ++i){
                pMem[i] = v.pMem[i];
            }
        }
        return *this;
    }
    TDynamicVector& operator=(TDynamicVector&& v) noexcept // Присваивание вектора векторов
    {
        if (this != &v){
            delete[] pMem;
            sz = v.sz;
            pMem = new T[sz](); // Это перепроверить
            for (int i = 0; i < sz; ++i){
                pMem[i] = v.pMem[i];
            }

        }
        return *this;
    }

    size_t size() const noexcept { return sz; }

    // индексация
    T& operator[](size_t ind)
    {
        if((ind<0)||(ind>=sz))
            throw length_error("AAAA");
        return pMem[ind];
    }
    const T& operator[](size_t ind) const
    {
        return pMem[ind];
    }
    // индексация с контролем
    T& at(size_t ind)
    {
        return &pMem[ind];
    }
    const T& at(size_t ind) const
    {
        return &pMem[ind];
    }

    // сравнение
    bool operator==(const TDynamicVector& v) const noexcept
    {
        if(sz != v.sz)
            return false;
        for (int i = 0; i < sz; ++i)
            if(pMem[i] != v[i])
                return false;
        return true;
    }
    bool operator!=(const TDynamicVector& v) const noexcept
    {
        if(sz != v.sz)
            return true;
        for (int i = 0; i < sz; ++i)
            if(pMem[i] != v[i])
                return true;
        return false;
    }

    // скалярные операции
    TDynamicVector operator+(T val)
    {
        for (int i = 0; i < sz; ++i)
            pMem[i] += val;
        return *this;
    }
    TDynamicVector operator-(T val)
    {
        for (int i = 0; i < sz; ++i)
            pMem[i] -= val;
        return *this;
    }
    TDynamicVector operator*(T val)
    {
        TDynamicVector<T> temp(sz);
        for (int i = 0; i < sz; ++i)
             temp[i]= pMem[i] * val;
        return temp;
    }

    // векторные операции
    TDynamicVector operator+(const TDynamicVector& v)
    {
        if(sz != v.sz)
            throw length_error("Sizes don't match!");
        auto t = TDynamicVector<T>(*this);
        for (int i = 0; i < sz; ++i)
            t.pMem[i] += v[i];
        return t;
    }
    TDynamicVector<int> mult_modulo(int by, int modulus){
        TDynamicVector<int> a(sz);
        for (int i = 0; i < sz; ++i) {
            int t = (pMem[i] * by) % modulus;
            if(t<0)
                a[i] = modulus+t;
            else{
                a[i] = t;
            }
        }
        return a;
    }
    TDynamicVector<int> add_modulo(const TDynamicVector<int>& to, int modulus){
        TDynamicVector<int> a(sz);
        for (int i = 0; i < sz; ++i) {
            int t = (pMem[i] + to.pMem[i]) % modulus;
            a[i] = t;

        }
        return a;
    }
    TDynamicVector operator-(const TDynamicVector& v)
    {
        if(sz != v.sz)
            throw length_error("Sizes don't match!");
        auto t = TDynamicVector<T>(*this);
        for (int i = 0; i < sz; ++i)
            t.pMem[i] -= v[i];
        return t;
    }
    T operator*(const TDynamicVector& v) //noexcept(noexcept(T())) - вынужден это закомментить, чтобы бросать exception-ы для тестов
    {
        if(sz != v.sz)
            throw length_error("Sizes don't match!");
        T S = pMem[0]*v[0];
        for (int i = 1; i < sz; ++i)
            S += pMem[i]*v[i];
        return S;
    }

    friend void swap(TDynamicVector& lhs, TDynamicVector& rhs) noexcept
    {
        std::swap(lhs.sz, rhs.sz);
        std::swap(lhs.pMem, rhs.pMem);
    }

    // ввод/вывод
    friend istream& operator>>(istream& istr, TDynamicVector& v)
    {
        for (size_t i = 0; i < v.sz; i++)
            istr >> v.pMem[i]; // требуется оператор>> для типа T
        return istr;
    }
    friend ostream& operator<<(ostream& ostr, const TDynamicVector& v)
    {
        for (size_t i = 0; i < v.sz; i++)
            ostr << v.pMem[i] << ", "; // требуется оператор<< для типа T
        return ostr;
    }

};



template<typename T>
class Matrix {
private:
    unsigned int N; // Размер
public:
    TDynamicVector<TDynamicVector<T>> val;
    unsigned int GetSize() const;
    // Операции выбраны исходя из необходимости
    explicit Matrix(unsigned int n); // n * n zero matrix
    explicit Matrix(const TDynamicVector<TDynamicVector<T>>& v); // Преобразователь типа
    TDynamicVector<T>& operator[](unsigned int i); // i-я строка
    Matrix operator*(Matrix& m); // the most important function, it needs to be fast. Using adaptive Strassen's method
    TDynamicVector<T> operator*(const TDynamicVector<T>& vec);
    Matrix transpose();
    friend ostream& operator<<(ostream& ostr, Matrix<T>& m);
    void print(); // ВТорой вывод
    static vector<vector<double> >
    add_matrix(vector<vector<double> > matrix_A,
               vector<vector<double> > matrix_B, int split_index,
               int multiplier );
    static vector<vector<double> >
    multiply_matrix(vector<vector<double> > matrix_A,
                    vector<vector<double> > matrix_B);
    void randomize(unsigned int range); // random int values

};
int modulo_inverse(int a, int rem); // Обратный по модулю
int chinese_remainder_theorem(const int* numbers, const int* remainders, int n); // Мин. положит. решение китайской теоремы об остатках

pair<Matrix<int>,Matrix<int>> decompose(Matrix<int>& A); //  [ THE BEST SO FAR ]
void compute_SNF(Matrix<int>& A);
pair<Matrix<int>,Matrix<int>> solve_SLDE_mod_p(Matrix<int>& A, TDynamicVector<int>& b,int p); // Ax = b mod p
// В перспективе будут возвращаться две матрицы: U и R. U A R = S, S - смитова нормальная форма
vector<int> SieveOfEratosthenes(int n); // Возвращает вектор простых чисел, меньших n
pair<Matrix<int>,Matrix<int>> solve_SLDE_modular_method(Matrix<int>& A, TDynamicVector<int>& b); // Решаем достаточно систем по модулю, затем собираем всё вместе
//--------------------------------------------------------------
// GGH crypto-algorithm
//--------------------------------------------------------------
Matrix<int> gen_public_key(int size,int range);
vector<Matrix<int>> gen_private_key(Matrix<int>& public_key);

// e = x A + r
TDynamicVector<int> encrypt(const TDynamicVector<int>& message, Matrix<int>& public_key,short sigma);

TDynamicVector<int> decrypt(const TDynamicVector<int>& message, Matrix<int> U,Matrix<int> private_key,Matrix<int> A,short sigma);

#endif //FASTHNF_MATRIX_H
