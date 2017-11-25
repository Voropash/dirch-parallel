#pragma once

// Файл содержит задание используемых типов и функций
// для решения задачи u''xx + u''yy + F = 0

typedef double NumericType;
#ifdef MPI_VERSION
const MPI_Datatype MpiNumericType = MPI_DOUBLE;
#endif
const NumericType DefaultEps = static_cast<NumericType>( 0.0001 );

#include <MathObjects.h> // чтобы использовать CArea
#include <math.h>

// Область решения задачи.
const CArea Area(0, 3, 0, 3);

// Правая часть.
inline NumericType F(NumericType x, NumericType y) {
    const NumericType f = (x * x + y * y) / ((1 + x * y) * (1 + x * y));
    return f;
}

// Граничная функция.
inline NumericType Phi(NumericType x, NumericType y) {
    const NumericType phi = log(1 + x * y);
    return phi;
}
