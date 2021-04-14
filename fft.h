#define _USE_MATH_DEFINES

#include <complex>
#include <vector>
#include <array>
#include <iostream>
#include <cmath>
#include <omp.h>

using complex = std::complex<float>;

//蝶形操作
inline void butter_fly(complex *X, complex *X1, complex *X2, size_t n, size_t N, complex *W)
{
    complex w,w0;
    w = 1;
    for (int i = 0; i < n / 2; i++)
    {
        w = W[i * N / n];
        X[i] = X1[i] + w * X2[i];
        X[i + n / 2] = X1[i] - w * X2[i];
        w *= w0;
    }
}

//使用递归实现的dit
complex *dit_recursion(complex *x, size_t n);

//使用非递归实现的dit
complex *dit(complex *x, size_t n);

complex *dit_recursion(complex *x, size_t n, size_t N, complex *W)
{
    if (n == 1)
    {
        complex *X = new complex[1];
        X[0] = x[0];
        return X;
    }

    complex *x1 = new complex[n / 2];
    complex *x2 = new complex[n / 2];
    for (int i = 0; i < n / 2; i++)
    {
        x1[i] = x[i * 2];
        x2[i] = x[i * 2 + 1];
    }

    complex *X1;
    complex *X2;
    complex *X = new complex[n];
    X1 = dit_recursion(x1, n / 2, N, W);
    X2 = dit_recursion(x2, n / 2, N, W);
    delete[] x1;
    delete[] x2;

    butter_fly(X, X1, X2, n, N, W);

    delete[] X1;
    delete[] X2;
    return X;
}

//使用递归实现的dit
complex *dit_recursion(complex *x, size_t N)
{
    complex *W = new complex[N];
    for (int i = 0; i < N; i++)
        W[i] = complex(cos(2 * M_PI * i / N), sin(2 * M_PI * i / N));
    complex *X = dit_recursion(x, N, N, W);
    delete W;
    return X;
}

//使用非递归实现的dit
complex *dit(complex *x, size_t N)
{
    size_t r = 0;
    while ((1 << r) < N) //r=logN
        r++;
    size_t *rev = new size_t[N]; //
    rev[0] = 0;
    rev[1] = N / 2;
    for (int i = 2; i < N; i++)
        rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (r - 1)); //magic
#pragma omp parallel for num_threads(8)
    for (int i = 0; i < N; i += 2)
        swap(x[i], x[rev[i]]);

    complex *W = new complex[N];
#pragma omp parallel for num_threads(8)
    for (int i = 0; i < N; i++)
        W[i] = complex(cos(2 * M_PI * i / N), sin(2 * M_PI * i / N));

    complex *X = new complex[N];
    complex *X_ = new complex[N];
    for (int i = 0; i < N; i++)
        X[i] = x[i];
    for (size_t len = 2; len <= N; len *= 2)
    {
#pragma omp parallel for num_threads(8)
        for (int i = 0; i < N; i += len)
            butter_fly(X_ + i, X + i, X + i + len / 2, len, N, W);
        std::swap(X, X_);
    }

    delete[] X_;
    return X;
}
