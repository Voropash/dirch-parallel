#include <Std.h>
#include <Definitions.h>
#include <MathFunctions.h>

///////////////////////////////////////////////////////////////////////////////

NumericType LaplasOperator(const CMatrix &matrix, const CUniformGrid &grid, size_t x, size_t y) {
    const NumericType ldx = (matrix(x, y) - matrix(x - 1, y)) / grid.X.Step(x - 1); // производные численные
    const NumericType rdx = (matrix(x + 1, y) - matrix(x, y)) / grid.X.Step(x);
    const NumericType tdy = (matrix(x, y) - matrix(x, y - 1)) / grid.Y.Step(y - 1);
    const NumericType bdy = (matrix(x, y + 1) - matrix(x, y)) / grid.Y.Step(y);
    const NumericType dx = (ldx - rdx) / grid.X.AverageStep(x); // численные вторые производные
    const NumericType dy = (tdy - bdy) / grid.Y.AverageStep(y);
    return (dx + dy);
}

// Вычисление невязки rij во внутренних точках.
void CalcR(const CMatrix &p, const CUniformGrid &grid, CMatrix &r) {
#ifndef DIRCH_NO_OPENMP
    // omp_set_num_threads(3); // устанавливаем 3 потока
#pragma omp parallel for // https://mindhalls.ru/pragma-omp-directives-samples/
    for (long x = 1; x < r.SizeX() - 1; x++) {
        // omp_get_thread_num() для того чтобы узнать номер потока
#else
        for( size_t x = 1; x < r.SizeX() - 1; x++ ) {
#endif
        for (size_t y = 1; y < r.SizeY() - 1; y++) {
            r(x, y) = LaplasOperator(p, grid, x, y) - F(grid.X[x], grid.Y[y]);
        }
    }
}

// Вычисление значений gij во внутренних точках.
void CalcG(const CMatrix &r, const NumericType alpha, CMatrix &g) {
#ifndef DIRCH_NO_OPENMP
    // omp_set_num_threads(3); // устанавливаем 3 потока
#pragma omp parallel for // https://mindhalls.ru/pragma-omp-directives-samples/
    for (long x = 1; x < g.SizeX() - 1; x++) {
#else
        for( size_t x = 1; x < g.SizeX() - 1; x++ ) {
#endif
        for (size_t y = 1; y < g.SizeY() - 1; y++) {
            g(x, y) = r(x, y) - alpha * g(x, y);
        }
    }
}

// Вычисление значений pij во внутренних точках, возвращается евклидова норма.
NumericType CalcP(const CMatrix &g, const NumericType tau, CMatrix &p) {
    NumericType squares = 0;
    for (size_t x = 1; x < p.SizeX() - 1; x++) {
        for (size_t y = 1; y < p.SizeY() - 1; y++) {
            const NumericType newValue = p(x, y) - tau * g(x, y); // считаем новое значение
            squares += (newValue - p(x, y)) * (newValue - p(x, y)); // добавляем квадрат разницы компоненты
            p(x, y) = newValue; // присваиваем новое значение
        }
    }
    return static_cast<NumericType> (pow(squares, 0.5));
}

// Вычисление значений pij во внутренних точках, возвращается сумма квадратов.
NumericType CalcP_2(const CMatrix &g, const NumericType tau, CMatrix &p) {
    NumericType squares = 0;
    for (size_t x = 1; x < p.SizeX() - 1; x++) {
        for (size_t y = 1; y < p.SizeY() - 1; y++) {
            const NumericType newValue = p(x, y) - tau * g(x, y); // считаем новое значение
            squares += (newValue - p(x, y)) * (newValue - p(x, y)); // добавляем квадрат разницы компоненты
            p(x, y) = newValue; // присваиваем новое значение
        }
    }
    return static_cast<NumericType> (squares);
}

// Вычисление alpha.
CFraction CalcAlpha(const CMatrix &r, const CMatrix &g, const CUniformGrid &grid) {
    NumericType numerator = 0;
    NumericType denominator = 0;
#ifndef DIRCH_NO_OPENMP
    // omp_set_num_threads(3); // устанавливаем 3 потока
#pragma omp parallel for reduction( +:numerator, denominator ) // https://habrahabr.ru/company/intel/blog/88574/
    // суммируем числитель и знаменатель
    for (long x = 1; x < r.SizeX() - 1; x++) {
#else
        for( size_t x = 1; x < r.SizeX() - 1; x++ ) {
#endif
        for (size_t y = 1; y < r.SizeY() - 1; y++) {
            const NumericType common = g(x, y) * grid.X.AverageStep(x) * grid.Y.AverageStep(y);
            numerator += LaplasOperator(r, grid, x, y) * common;
            denominator += LaplasOperator(g, grid, x, y) * common;
        }
    }
    return CFraction(numerator, denominator);
}

// Вычисление tau.
CFraction CalcTau(const CMatrix &r, const CMatrix &g, const CUniformGrid &grid) {
    NumericType numerator = 0;
    NumericType denominator = 0;
#ifndef DIRCH_NO_OPENMP
    // omp_set_num_threads(3); // устанавливаем 3 потока
#pragma omp parallel for reduction( +:numerator, denominator ) //https://habrahabr.ru/company/intel/blog/88574/
    // суммируем числитель и знаменатель
    for (long x = 1; x < r.SizeX() - 1; x++) {
#else
        for( size_t x = 1; x < r.SizeX() - 1; x++ ) {
#endif
        for (size_t y = 1; y < r.SizeY() - 1; y++) {
            const NumericType common = g(x, y) * grid.X.AverageStep(x) * grid.Y.AverageStep(y);
            numerator += r(x, y) * common;
            denominator += LaplasOperator(g, grid, x, y) * common;
        }
    }
    return CFraction(numerator, denominator);
}

///////////////////////////////////////////////////////////////////////////////
