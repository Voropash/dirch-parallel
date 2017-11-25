#include <Std.h>
#include <math.h>
#include <Definitions.h>
#include <MathObjects.h>
#include <Errors.h>

///////////////////////////////////////////////////////////////////////////////

CMatrix &CMatrix::operator=(const CMatrix &other) {
    sizeX = other.sizeX;
    sizeY = other.sizeY;
    values = other.values;
    return *this;
}

void CMatrix::Init(const size_t _sizeX, const size_t _sizeY) {
    sizeX = _sizeX;
    sizeY = _sizeY;
    values.resize(sizeX * sizeY);
    fill(values.begin(), values.end(), static_cast<NumericType>( 0 ));
}

///////////////////////////////////////////////////////////////////////////////

ostream &operator<<(ostream &out, const CMatrixPart &matrixPart) {
    out << "[" << matrixPart.BeginX << ", " << matrixPart.EndX << ") x "
        << "[" << matrixPart.BeginY << ", " << matrixPart.EndY << ")";
    return out;
}

///////////////////////////////////////////////////////////////////////////////

NumericType CUniformPartition::BorderFunc(NumericType t) { // функция, f(t) из методички - как параметр сетки
    return static_cast<NumericType> (pow((1.0 + t), 1.5) - 1) / (pow(2.0, 1.5) - 1);
}

void CUniformPartition::PartInit(NumericType p0, NumericType pN, size_t size, size_t begin, size_t end) {
    if (!(p0 < pN)) {
        throw CException("CUniformPartition: bad interval");
    }
    if (!(size > 1)) {
        throw CException("CUniformPartition: inavalid size");
    }
    if (!(begin < end && end <= size)) {
        throw CException("CUniformPartition: invalid [begin, end)");
    }

    ps.clear();
    ps.reserve(end - begin);

    for (size_t i = begin; i < end; i++) { // инициализация сетки с использованием функции f(t), задающей неравномерность
        NumericType part = static_cast<NumericType>( i ) / (size - 1);
        part = CUniformPartition::BorderFunc(part);
        const NumericType p = part * pN + (1 - part) * p0;
        ps.push_back(p);
    }
}

///////////////////////////////////////////////////////////////////////////////
