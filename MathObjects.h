#pragma once

///////////////////////////////////////////////////////////////////////////////

struct CArea { // область, где считаем
	NumericType X0;
	NumericType Xn;
	NumericType Y0;
	NumericType Yn;

	CArea( NumericType x0, NumericType xn, NumericType y0, NumericType yn ) :
		X0( x0 ), Xn( xn ), Y0( y0 ), Yn( yn )
	{
	}
};

///////////////////////////////////////////////////////////////////////////////

struct CFraction { // дробь для взаимодействий через MPI
	NumericType Numerator; // числитель
	NumericType Denominator; // знаменатель

	explicit CFraction( NumericType numerator = 0, NumericType denominator = 1 ) :
		Numerator( numerator ),
		Denominator( denominator )
	{
	}

	NumericType Value() const
	{
		return ( Numerator / Denominator );
	}
};

///////////////////////////////////////////////////////////////////////////////

class CMatrix { // матрица
public:
	CMatrix() :
		sizeX( 0 ),
		sizeY( 0 )
	{
	}

	CMatrix( size_t sizeX, size_t sizeY )
	{
		Init( sizeX, sizeY );
	}

	CMatrix( const CMatrix& other ) { *this = other; };
	CMatrix& operator=( const CMatrix& other );


	void Init( const size_t _sizeX, const size_t _sizeY );

	NumericType& operator()( size_t x, size_t y )
	{
		return values[y * sizeX + x];
	}
	NumericType operator()( size_t x, size_t y ) const
	{
		return values[y * sizeX + x];
	}

	size_t SizeX() const { return sizeX; }
	size_t SizeY() const { return sizeY; }

private:
	size_t sizeX;
	size_t sizeY;
	vector<NumericType> values;
};

///////////////////////////////////////////////////////////////////////////////

struct CMatrixPart { // описание границы подматрицы
	size_t BeginX;
	size_t EndX;
	size_t BeginY;
	size_t EndY;

	CMatrixPart() : // изначально пустая
		BeginX( 0 ), EndX( 0 ),
		BeginY( 0 ), EndY( 0 )
	{
	}

	CMatrixPart( size_t beginX, size_t endX, size_t beginY, size_t endY ) : // конструктор
		BeginX( beginX ), EndX( endX ),
		BeginY( beginY ), EndY( endY )
	{
		assert( BeginX < EndX );
		assert( BeginY < EndY );
	}

	void SetRow( size_t beginX, size_t endX, size_t y ) // строка от-до
	{
		assert( beginX < endX );
		BeginX = beginX;
		EndX = endX;
		BeginY = y;
		EndY = y + 1;
	}

	void SetColumn( size_t x, size_t beginY, size_t endY ) // столбец от-до
	{
		assert( beginY < endY );
		BeginX = x;
		EndX = x + 1;
		BeginY = beginY;
		EndY = endY;
	}

	size_t SizeX() const // getter длины
	{
		return ( EndX - BeginX );
	}

	size_t SizeY() const // getter ширины
	{
		return ( EndY - BeginY );
	}

	size_t Size() const // кол-во узлов сетки
	{
		return SizeX() * SizeY();
	}
};

ostream& operator<<( ostream& out, const CMatrixPart& matrixPart );

///////////////////////////////////////////////////////////////////////////////

class CUniformPartition {
private:
	CUniformPartition( const CUniformPartition& );
	CUniformPartition& operator=( const CUniformPartition& );

public:
	CUniformPartition() {}

	NumericType BorderFunc( NumericType t );
	void PartInit( NumericType p0, NumericType pN, size_t size, size_t begin, size_t end );
	void Init( NumericType p0, NumericType pN, size_t N )
	{
		PartInit( p0, pN, N, 0, N );
	}

	size_t Size() const
	{
		return ps.size();
	}
	NumericType operator[]( size_t i ) const
	{
		return Point( i );
	}
	NumericType Point( size_t i ) const
	{
		return ps[i];
	}
	NumericType Step( size_t i ) const
	{
		return ( Point( i + 1 ) - Point( i ) );
	}
	NumericType AverageStep( size_t i ) const
	{
		return ( Point( i + 1 ) - Point( i - 1 ) ) / static_cast<NumericType>( 2 );
	}

private:
	vector<NumericType> ps;
};

///////////////////////////////////////////////////////////////////////////////

struct CUniformGrid { // описание грида - по сути 2 вектора с вспомогательными методами
	CUniformPartition X;
	CUniformPartition Y;

	CMatrixPart Column( size_t x, size_t decreaseTop = 0, size_t decreaseBottom = 0 ) const
	{
		CMatrixPart part;
		part.SetColumn( x, decreaseTop, Y.Size() - decreaseBottom );
		return part;
	}
	CMatrixPart Row( size_t y, size_t decreaseLeft = 0, size_t decreaseRight = 0 ) const
	{
		CMatrixPart part;
		part.SetRow( decreaseLeft, X.Size() - decreaseRight, y );
		return part;
	}
};

///////////////////////////////////////////////////////////////////////////////
