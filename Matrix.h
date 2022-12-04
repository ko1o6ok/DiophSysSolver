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
    explicit TDynamicVector(size_t size = 1) : sz(size)
    {
        if (sz <= 0)
            throw length_error("Vector size should be greater than zero");
//        if (sz > MAX_VECTOR_SIZE)
//            throw length_error("Vector size shouldn't be bigger than 10^8");
        pMem = new T[sz]();// {}; // У типа T д.б. констуктор по умолчанию
    }
    TDynamicVector(T* arr, size_t s) : sz(s)
    {
        assert(arr != nullptr && "TDynamicVector ctor requires non-nullptr arg");
        pMem = new T[sz];
        std::copy(arr, arr + sz, pMem);
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
    TDynamicVector operator-(double val)
    {
        for (int i = 0; i < sz; ++i)
            pMem[i] -= val;
        return *this;
    }
    TDynamicVector operator*(double val)
    {
        TDynamicVector<int> temp(sz);
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
    TDynamicVector operator-(const TDynamicVector& v)
    {
        if(sz != v.sz)
            throw length_error("Sizes don't match!");
        auto t = TDynamicVector<T>(*this);
        for (int i = 0; i < sz; ++i)
            t.pMem[i] -= v[i];
        return t;
    }
    unsigned int find_min_nonzero(){
        if(sz == 1)
            return 0;
        int m = pMem[0];
        int i = 0;
        while (pMem[i]==0){
            i++;
        }
        int m_pos = i;
        while (i < sz){
            auto t = pMem[i];
            if((t!=0)&&(t < m)){
                m = t;
                m_pos = i;
            }
            i++;
        }
        return m_pos;
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




class Matrix {
private:
    unsigned int N; // Размер

public:
    TDynamicVector<TDynamicVector<int>> val;
    unsigned int GetSize() const;
    // Операции выбраны исходя из необходимости
    explicit Matrix(unsigned int n); // n * n zero matrix
    explicit Matrix(const TDynamicVector<TDynamicVector<int>>& v); // Преобразователь типа
    TDynamicVector<int>& operator[](unsigned int i); // i-я строка
    Matrix operator*(Matrix& m); // the most important function, it needs to be fast. Using adaptive Strassen's method
    TDynamicVector<int> operator*(const TDynamicVector<int>& vec);
    Matrix transpose();
    friend ostream& operator<<(ostream& ostr, Matrix& m);
    static vector<vector<int> >
    add_matrix(vector<vector<int> > matrix_A,
               vector<vector<int> > matrix_B, int split_index,
               int multiplier );
    static vector<vector<int> >
    multiply_matrix(vector<vector<int> > matrix_A,
                    vector<vector<int> > matrix_B);
    void randomize(unsigned int range); // random int values

    // The function I'm trying to optimize
    // Fast Hermite Normal form of the val matrix
    //pair<Matrix,Matrix> FHNF(); // Returns pair of U and H, such that det(U) = +-1 and U H = A (H is upper-triangular)

};
void solve_SLDE(Matrix& A, TDynamicVector<int>& b); // Soving A x = b

#endif //FASTHNF_MATRIX_H
