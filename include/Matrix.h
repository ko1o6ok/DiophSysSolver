#ifndef FASTHNF_MATRIX_H
#define FASTHNF_MATRIX_H
#include <vector>
#include <iostream>
using namespace std;
// Динамический вектор -
// шаблонный вектор на динамической памяти
template<typename T>
class TDynamicVector
{
protected:
    size_t sz{};
public:
    T* pMem;
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
    void push(T el){
        TDynamicVector<T> temp(*this);
        delete[] pMem;
        pMem = new T[sz+1](); // Это перепроверить
        for (int i = 0; i < sz; ++i){
            pMem[i] = temp.pMem[i];
        }
        pMem[sz] = el;
        sz = sz+1;
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
};
template<typename T>
class Matrix {
private:
    unsigned int N; // Размер
public:
    TDynamicVector<TDynamicVector<T>> val;
    unsigned int GetSize() const;
    long nullity() const; // Для приведённой матрицы!!!
    long rank() const; // Для приведённой матрицы!!!
    // Операции выбраны исходя из необходимости
    explicit Matrix(unsigned int n); // n * n zero matrix
    explicit Matrix(const TDynamicVector<TDynamicVector<T>>& v); // Преобразователь типа
    TDynamicVector<T>& operator[](unsigned int i); // i-я строка
    TDynamicVector<T> operator*(const TDynamicVector<T>& vec);
    Matrix transpose();
};
void compute_SNF(Matrix<long>& A);
Matrix<long> to_SNF(Matrix<long>& A);
#endif //FASTHNF_MATRIX_H
