#include "..\include\Vector.h"
#include "..\include\Matrix.h"
#include <assert.h>
#include <math.h>
#include <cstring>
#include <cfloat>

/*
================================
SparseMat
================================
*/
template <typename T> SparseMat<T>::SparseMat()
{
	m_row = new unsigned int[16];
	m_col = new unsigned int[16];
	m_data = new T[16];
}


template <typename T> SparseMat<T>::~SparseMat()
{
	delete[] m_row;
	delete[] m_col;
	delete[] m_data;
}


template <typename T> void SparseMat<T>::AddValue( T value, unsigned int nRow, unsigned int nColumn )
{
	if ( m_nCurrentSize + 1 >= m_nAllocatedSize )
	{
		unsigned int *new_row = new unsigned int[m_nAllocatedSize * 2];
		memccpy( new_row, m_row, m_nAllocatedSize );
		delete[] m_row;
		m_row = new_row;

		unsigned int *new_col = new unsigned int[m_nAllocatedSize * 2];
		memccpy( new_col, m_col, m_nAllocatedSize );
		delete[] m_col;
		m_row = new_col;

		unsigned int *new_data = new unsigned int[m_nAllocatedSize * 2];
		memccpy( new_data, m_data, m_nAllocatedSize );
		delete[] m_data;
		m_row = new_data;
	}

	m_row[m_nCurrentSize] = nRow;
	m_col[m_nCurrentSize] = nColumn;
	m_data[m_nCurrentSize] = value;
	m_nCurrentSize += 1;
}


/*
================================
MatN
================================
*/
MatN::MatN( const unsigned int n )
{
	m_row = n;
	m_col = n;
	m_data = new float[n * n];
	for ( unsigned int i = 0; i < n * n; i++ )
	{
		m_data[i] = 0.0f;
	}
}

MatN::MatN( const unsigned int n, const float *data )
{
	m_row = n;
	m_col = n;
	m_data = new float[n * n];
	for ( unsigned int i = 0; i < n * n; i++ )
	{
		m_data[i] = data[i];
	}
}

MatN::MatN( const unsigned int row, const unsigned int col )
{
	m_row = row;
	m_col = col;
	m_data = new float[row * col];
	for ( unsigned int i = 0; i < row * col; i++ )
	{
		m_data[i] = 0.0f;
	}
}

MatN::MatN( const unsigned int row, const unsigned int col, const float *data )
{
	m_row = row;
	m_col = col;
	m_data = new float[row * col];
	for ( unsigned int i = 0; i < row * col; i++ )
	{
		m_data[i] = data[i];
	}
}

MatN::MatN( const MatN &mat )
{
	m_row = mat.m_row;
	m_col = mat.m_col;
	m_data = new float[m_row * m_col];
	for ( unsigned int i = 0; i < m_row * m_col; i++ )
	{
		m_data[i] = mat.m_data[i];
	}
}

void MatN::operator=( MatN other )
{
	if ( m_col != other.m_col || m_row != other.m_row )
	{
		delete[] m_data;
		m_data = new float[other.m_col * other.m_row];
	}

	for ( unsigned int i = 0; i < m_col * m_row; i++ )
	{
		m_data[i] = other.m_data[i];
	}
}

bool MatN::operator==( MatN other )
{
	if ( m_col != other.m_col || m_row != other.m_row )
	{
		return false;
	}

	for ( unsigned int i = 0; i < m_col * m_row; i++ )
	{
		if ( m_data[i] != other.m_data[i] )
		{
			return false;
		}
	}
	return true;
}

VecN MatN::GetRowVec( const unsigned int row ) const
{
	assert( row < m_row );
	VecN returnVec( m_col );
	for ( unsigned int i = 0; i < m_col; i++ )
	{
		returnVec[i] = m_data[row * m_col + i];
	}
	return returnVec;
}

void MatN::SetRowVec( const unsigned int row, const VecN *vec )
{
	assert( row < m_row );
	for ( unsigned int i = 0; i < m_col; i++ )
	{
		m_data[row * m_col + i] = ( *vec )[i];
	}
}

VecN MatN::GetColVec( const unsigned int col ) const
{
	assert( col < m_col );
	VecN returnVec( m_row );
	for ( unsigned int i = 0; i < m_row; i++ )
	{
		returnVec[i] = m_data[col + i * m_col];
	}
	return returnVec;
}

void MatN::SetColVec( const unsigned int col, const VecN *vec )
{
	assert( col < m_col );
	for ( unsigned int i = 0; i < m_row; i++ )
	{
		m_data[col + i * m_col] = ( *vec )[i];
	}
}

float MatN::GetComponent( unsigned int m, unsigned int n ) const
{
	assert( m < m_row &&n < m_col );
	return m_data[m * m_row + n];
}

void MatN::SetComponent( unsigned int m, unsigned int n, float val )
{
	assert( m < m_row &&n < m_col );
	m_data[m * m_row + n] = val;
}

MatN MatN::operator+( MatN other ) const
{
	assert( m_col == other.m_col && m_row == other.m_row );
	MatN returnMat( m_row, m_col );
	for ( unsigned int i = 0; i < m_col * m_row; i++ )
	{
		returnMat.m_data[i] = m_data[i] + other.m_data[i];
	}
	return returnMat;
}

void MatN::operator+=( MatN other )
{
	assert( m_col == other.m_col && m_row == other.m_row );
	for ( unsigned int i = 0; i < m_col * m_row; i++ )
	{
		m_data[i] += other.m_data[i];
	}
}

MatN MatN::operator-( MatN other ) const
{
	assert( m_col == other.m_col && m_row == other.m_row );
	MatN returnMat( m_row, m_col );
	for ( unsigned int i = 0; i < m_col * m_row; i++ )
	{
		returnMat.m_data[i] = m_data[i] - other.m_data[i];
	}
	return returnMat;
}

void MatN::operator-=( MatN other )
{
	assert( m_col == other.m_col && m_row == other.m_row );
	for ( unsigned int i = 0; i < m_col * m_row; i++ )
	{
		m_data[i] -= other.m_data[i];
	}
}

MatN MatN::operator*( float scalar ) const
{
	MatN returnMat( m_row, m_col );
	for ( unsigned int i = 0; i < m_col * m_row; i++ )
	{
		returnMat.m_data[i] = m_data[i] * scalar;
	}
	return returnMat;
}

void MatN::operator*=( float scalar )
{
	for ( unsigned int i = 0; i < m_col * m_row; i++ )
	{
		m_data[i] *= scalar;
	}
}

VecN MatN::operator*( VecN vec ) const
{
	assert( m_col == vec.size() );
	VecN returnVec( m_row );
	for ( unsigned int i = 0; i < m_row; i++ )
	{
		for ( unsigned int j = 0; j < m_col; j++ )
		{
			returnVec[i] += m_data[i * m_col + j] * vec[j];
		}
	}
	return returnVec;
}

void MatN::operator*=( VecN vec )
{
	assert( m_col == vec.size() );
	float *newData = new float[m_row];
	for ( unsigned int i = 0; i < m_row; i++ )
	{
		newData[i] = 0.0f;
		for ( unsigned int j = 0; j < m_col; j++ )
		{
			newData[i] += m_data[i * m_col + j] * vec[j];
		}
	}
	delete[] m_data;
	m_data = newData;
}

MatN MatN::operator*( MatN other ) const
{
	assert( m_col == other.m_row );
	MatN returnMat( m_row, other.m_col );
	for ( unsigned int i = 0; i < m_row; i++ )
	{
		for ( unsigned int j = 0; j < other.m_col; j++ )
		{
			for ( unsigned int k = 0; k < other.m_row; k++ )
			{
				returnMat.m_data[i * other.m_col + j] += m_data[i * m_col + k] * other[k * other.m_col + j];
			}
		}
	}
	return returnMat;
}

void MatN::operator*=( MatN other )
{
	assert( m_col == other.m_row );
	float *newData = new float[m_col * m_row];
	for ( unsigned int i = 0; i < m_row; i++ )
	{
		for ( unsigned int j = 0; j < other.m_col; j++ )
		{
			newData[i * other.m_col + j] = 0.0f;
			for ( unsigned int k = 0; k < other.m_row; k++ )
			{
				newData[i * other.m_col + j] += m_data[i * m_col + k] * other[k * other.m_col + j];
			}
		}
	}
	delete[] m_data;
	m_data = newData;
}

MatN MatN::operator/( float scalar ) const
{
	MatN returnMat( m_row, m_col );
	for ( unsigned int i = 0; i < m_col * m_row; i++ )
	{
		returnMat.m_data[i] /= scalar;
	}
	return returnMat;
}

void MatN::operator/=( float scalar )
{
	for ( unsigned int i = 0; i < m_col * m_row; i++ )
	{
		m_data[i] /= scalar;
	}
}

void GetCofactor( const MatN &mat, MatN &temp, int p, int q, int n )
{
	int i = 0;
	int j = 0;

	//Looping for each element of the matrix
	for ( int r = 0; r < n; r++ )
	{
		for ( int c = 0; c < n; c++ )
		{
//Copying into temporary matrix only those element which are not in given row and column
			if ( r != p && c != q )
			{
				temp[i * temp.RowCount() + j] = mat[r * mat.RowCount() + c];
				j += 1;

				//Row is filled, so increase row index and reset col index
				if ( j == n - 1 )
				{
					j = 0;
					i++;
				}
			}
		}
	}
}

float DeterminantRecursive( const MatN &mat, int n )
{
	float D = 0; // Initialize result

	//Base case : if matrix contains single element
	if ( n == 1 )
	{
		return mat[0];
	}

	MatN temp( mat.RowCount() ); //To store cofactors

	float sign = 1.0f; //To store sign multiplier

	 //Iterate for each element of first row
	for ( int f = 0; f < n; f++ )
	{
// Getting Cofactor of A[0][f]
		GetCofactor( mat, temp, 0, f, n );
		D += sign * mat[f] * DeterminantRecursive( temp, n - 1 );

		// terms are to be added with alternate sign
		sign *= -1.0f;
	}

	return D;
}

float MatN::Determinant()
{
	return DeterminantRecursive( *this, m_row );
}

void Adjoint( const MatN &mat, MatN &adj )
{
	const int size = mat.RowCount();
	if ( size == 1 )
	{
		adj[0] = 1.0;
		return;
	}

	//temp is used to store cofactors of this
	float sign = 1.0f;
	MatN temp( size );

	for ( int i = 0; i < size; i++ )
	{
		for ( int j = 0; j < size; j++ )
		{
			GetCofactor( mat, temp, i, j, size );

			//sign of adj[j][i] positive if sum of row and column indexes is even.
			sign = ( ( i + j ) % 2 == 0 ) ? 1 : -1;

			// Interchanging rows and columns to get the transpose of the cofactor matrix
			adj[j * size + i] = sign * ( DeterminantRecursive( temp, size - 1 ) );
		}
	}
}

bool MatN::Inverse( MatN *inv )
{
	float det = Determinant();
	if ( fabsf( det ) < FLT_EPSILON )
	{
		return false;
	}

	MatN adj( m_row );
	Adjoint( *this, adj ); //Find adjoint

	//inverse( A ) = adj( A ) / det( A )
	for ( int i = 0; i < m_row; i++ )
	{
		for ( int j = 0; j < m_row; j++ )
		{
			( *inv )[i * m_row + j] = adj[i * m_row + j] / float( det );
		}
	}

	return true;
}

MatN MatN::Transpose() const
{
	MatN returnMat( m_col, m_row );
	for ( unsigned int i = 0; i < m_row; i++ )
	{
		for ( unsigned int j = 0; j < m_col; j++ )
		{
			returnMat[j * m_row + i] = m_data[i * m_col + j];
		}
	}
	return returnMat;
}

void MatN::Transposed()
{
	unsigned int temp = m_row;
	m_row = m_col;
	m_col = temp;
	for ( unsigned int i = 0; i < m_row; i++ )
	{
		for ( unsigned int j = 0; j < m_col; j++ )
		{
			float temp = m_data[j * m_row + i];
			m_data[j * m_row + i] = m_data[i * m_col + j];
			m_data[i * m_col + j] = temp;
		}
	}
}

Mat2 MatN::as_Mat2() const
{
	Mat2 returnMat;
	for ( unsigned int i = 0; i < m_row; i++ )
	{
		for ( unsigned int j = 0; j < m_col; j++ )
		{
			if ( i >= 2 || j >= 2 )
			{
				returnMat[i * m_row + j] = 0.0f;
			}
			else
			{
				returnMat[i * m_row + j] = m_data[i * m_row + j];
			}
		}
	}
	return returnMat;
}

Mat3 MatN::as_Mat3() const
{
	Mat3 returnMat;
	for ( unsigned int i = 0; i < m_row; i++ )
	{
		for ( unsigned int j = 0; j < m_col; j++ )
		{
			if ( i >= 3 || j >= 3 )
			{
				returnMat[i * m_row + j] = 0.0f;
			}
			else
			{
				returnMat[i * m_row + j] = m_data[i * m_row + j];
			}
		}
	}
	return returnMat;
}

Mat4 MatN::as_Mat4() const
{
	Mat4 returnMat;
	for ( unsigned int i = 0; i < m_row; i++ )
	{
		for ( unsigned int j = 0; j < m_col; j++ )
		{
			if ( i >= 4 || j >= 4 )
			{
				returnMat[i * m_row + j] = 0.0f;
			}
			else
			{
				returnMat[i * m_row + j] = m_data[i * m_row + j];
			}
		}
	}
	return returnMat;
}


/*
================================
Mat2
================================
*/
Mat2::Mat2()
{
	m_data[0] = 1.0f;
	m_data[1] = 0.0f;
	m_data[2] = 1.0f;
	m_data[3] = 0.0f;
}

Mat2::Mat2( float v )
{
	m_data[0] = v;
	m_data[1] = v;
	m_data[2] = v;
	m_data[3] = v;
}

Mat2::Mat2( const float *data )
{
	m_data[0] = data[0];
	m_data[1] = data[1];
	m_data[2] = data[2];
	m_data[3] = data[3];
}

Mat2::Mat2( const Mat2 &mat )
{
	m_data[0] = mat.m_data[0];
	m_data[1] = mat.m_data[1];
	m_data[2] = mat.m_data[2];
	m_data[3] = mat.m_data[3];
}

void Mat2::operator=( Mat2 other )
{
	m_data[0] = other.m_data[0];
	m_data[1] = other.m_data[1];
	m_data[2] = other.m_data[2];
	m_data[3] = other.m_data[3];
}

bool Mat2::operator==( Mat2 other )
{
	for ( unsigned int i = 0; i < 4; i++ )
	{
		if ( m_data[i] != other.m_data[i] )
		{
			return false;
		}
	}
	return true;
}

Vec2 Mat2::GetRowVec( const unsigned int row ) const
{
	assert( row < 2 );
	Vec2 returnVec;
	returnVec[0] = m_data[row * 2];
	returnVec[1] = m_data[row * 2 + 1];
	return returnVec;
}

void Mat2::SetRowVec( const unsigned int row, const Vec2 *vec )
{
	assert( row < 2 );
	m_data[row * 2] = ( *vec )[0];
	m_data[row * 2 + 1] = ( *vec )[1];
}

Vec2 Mat2::GetColVec( const unsigned int col ) const
{
	assert( col < 2 );
	Vec2 returnVec;
	returnVec[0] = m_data[col * 2];
	returnVec[1] = m_data[col + 1 * 2];
	return returnVec;
}

void Mat2::SetColVec( const unsigned int col, const Vec2 *vec )
{
	assert( col < 2 );
	m_data[col * 2] = ( *vec )[0];
	m_data[col + 1 * 2] = ( *vec )[1];
}

float Mat2::GetComponent( unsigned int m, unsigned int n ) const
{
	assert( m < 2 && n < 2 );
	return m_data[m * 2 + n];
}

void Mat2::SetComponent( unsigned int m, unsigned int n, float val )
{
	assert( m < 2 && n < 2 );
	m_data[m * 2 + n] = val;
}

Mat2 Mat2::operator+( Mat2 other ) const
{
	Mat2 returnMat;
	returnMat.m_data[0] = m_data[0] + other.m_data[0];
	returnMat.m_data[1] = m_data[1] + other.m_data[1];
	returnMat.m_data[2] = m_data[2] + other.m_data[2];
	returnMat.m_data[3] = m_data[3] + other.m_data[3];
	return returnMat;
}

void Mat2::operator+=( Mat2 other )
{
	m_data[0] += other.m_data[0];
	m_data[1] += other.m_data[1];
	m_data[2] += other.m_data[2];
	m_data[3] += other.m_data[3];
}

Mat2 Mat2::operator-( Mat2 other ) const
{
	Mat2 returnMat;
	returnMat.m_data[0] = m_data[0] - other.m_data[0];
	returnMat.m_data[1] = m_data[1] - other.m_data[1];
	returnMat.m_data[2] = m_data[2] - other.m_data[2];
	returnMat.m_data[3] = m_data[3] - other.m_data[3];
	return returnMat;
}

void Mat2::operator-=( Mat2 other )
{
	m_data[0] -= other.m_data[0];
	m_data[1] -= other.m_data[1];
	m_data[2] -= other.m_data[2];
	m_data[3] -= other.m_data[3];
}

Mat2 Mat2::operator*( float scalar ) const
{
	Mat2 returnMat;
	returnMat.m_data[0] = m_data[0] * scalar;
	returnMat.m_data[1] = m_data[1] * scalar;
	returnMat.m_data[2] = m_data[2] * scalar;
	returnMat.m_data[3] = m_data[3] * scalar;
	return returnMat;
}

void Mat2::operator*=( float scalar )
{
	m_data[0] *= scalar;
	m_data[1] *= scalar;
	m_data[2] *= scalar;
	m_data[3] *= scalar;
}

Vec2 Mat2::operator*( Vec2 vec ) const
{
	Vec2 returnVec;
	returnVec[0] = m_data[0] * vec[0] + m_data[1] * vec[1];
	returnVec[1] = m_data[2] * vec[0] + m_data[3] * vec[1];
	return returnVec;
}

Mat2 Mat2::operator*( Mat2 other ) const
{
	Mat2 returnMat;
	returnMat[0] = m_data[0] * other.m_data[0] + m_data[1] * other.m_data[2];
	returnMat[1] = m_data[0] * other.m_data[1] + m_data[1] * other.m_data[3];
	returnMat[2] = m_data[1] * other.m_data[0] + m_data[2] * other.m_data[2];
	returnMat[3] = m_data[1] * other.m_data[1] + m_data[2] * other.m_data[3];
	return returnMat;
}

void Mat2::operator*=( Mat2 other )
{
	float a, b, c, d;
	a = m_data[0] * other.m_data[0] + m_data[1] * other.m_data[2];
	b = m_data[0] * other.m_data[1] + m_data[1] * other.m_data[3];
	c = m_data[1] * other.m_data[0] + m_data[2] * other.m_data[2];
	d = m_data[1] * other.m_data[1] + m_data[2] * other.m_data[3];
	m_data[0] = a;
	m_data[1] = b;
	m_data[2] = c;
	m_data[3] = d;
}

Mat2 Mat2::operator/( float scalar ) const
{
	Mat2 returnMat;
	returnMat.m_data[0] = m_data[0] / scalar;
	returnMat.m_data[1] = m_data[1] / scalar;
	returnMat.m_data[2] = m_data[2] / scalar;
	returnMat.m_data[3] = m_data[3] / scalar;
	return returnMat;
}

void Mat2::operator/=( float scalar )
{
	m_data[0] /= scalar;
	m_data[1] /= scalar;
	m_data[2] /= scalar;
	m_data[3] /= scalar;
}

float Mat2::Determinant()
{
	return m_data[0] * m_data[3] - m_data[1] * m_data[2];
}

bool Mat2::Inverse( Mat2 *inv )
{
	const float det = Determinant();

	if ( fabs( det ) < FLT_EPSILON )
	{
		return false;
	}

	inv->m_data[0] = m_data[3] / det;
	inv->m_data[1] = m_data[1] / det * -1.0f;
	inv->m_data[2] = m_data[2] / det * -1.0f;
	inv->m_data[3] = m_data[0] / det;

	return true;
}

Mat2 Mat2::Transpose() const
{
	Mat2 returnMat;
	returnMat.m_data[0] = m_data[0];
	returnMat.m_data[1] = m_data[2];
	returnMat.m_data[2] = m_data[1];
	returnMat.m_data[3] = m_data[3];
	return returnMat;
}

void Mat2::Transposed()
{
	for ( unsigned int i = 0; i < 2; i++ )
	{
		for ( unsigned int j = 0; j < 2; j++ )
		{
			float temp = m_data[j * 2 + i];
			m_data[j * 2 + i] = m_data[i * 2 + j];
			m_data[i * 2 + j] = temp;
		}
	}
}

MatN Mat2::as_MatN() const
{
	MatN returnMat( 2 );
	returnMat[0] = m_data[0];
	returnMat[1] = m_data[1];
	returnMat[2] = m_data[2];
	returnMat[3] = m_data[3];
	return returnMat;
}

Mat3 Mat2::as_Mat3() const
{
	Mat3 returnMat;
	returnMat[0] = m_data[0];
	returnMat[1] = m_data[1];
	returnMat[3] = m_data[2];
	returnMat[4] = m_data[3];
	return returnMat;
}

Mat4 Mat2::as_Mat4() const
{
	Mat4 returnMat;
	returnMat[0] = m_data[0];
	returnMat[1] = m_data[1];
	returnMat[4] = m_data[2];
	returnMat[5] = m_data[3];
	return returnMat;
}


/*
================================
Mat3
================================
*/
Mat3::Mat3()
{
	m_data[0] = 1.0f;
	m_data[1] = 0.0f;
	m_data[2] = 0.0f;
	m_data[3] = 0.0f;
	m_data[4] = 1.0f;
	m_data[5] = 0.0f;
	m_data[6] = 0.0f;
	m_data[7] = 0.0f;
	m_data[8] = 1.0f;
}

Mat3::Mat3( float v )
{
	m_data[0] = v;
	m_data[1] = v;
	m_data[2] = v;
	m_data[3] = v;
	m_data[4] = v;
	m_data[5] = v;
	m_data[6] = v;
	m_data[7] = v;
	m_data[8] = v;
}

Mat3::Mat3( const Vec3 &vA, const Vec3 &vB, const Vec3 &vC )
{
	m_data[0] = vA[0];
	m_data[1] = vA[1];
	m_data[2] = vA[2];
	m_data[3] = vB[0];
	m_data[4] = vB[1];
	m_data[5] = vB[2];
	m_data[6] = vC[0];
	m_data[7] = vC[1];
	m_data[8] = vC[2];
}

Mat3::Mat3( const float *data )
{
	m_data[0] = data[0];
	m_data[1] = data[1];
	m_data[2] = data[2];
	m_data[3] = data[3];
	m_data[4] = data[4];
	m_data[5] = data[5];
	m_data[6] = data[6];
	m_data[7] = data[7];
	m_data[8] = data[8];
}

Mat3::Mat3( const Mat3 &mat )
{
	m_data[0] = mat.m_data[0];
	m_data[1] = mat.m_data[1];
	m_data[2] = mat.m_data[2];
	m_data[3] = mat.m_data[3];
	m_data[4] = mat.m_data[4];
	m_data[5] = mat.m_data[5];
	m_data[6] = mat.m_data[6];
	m_data[7] = mat.m_data[7];
	m_data[8] = mat.m_data[8];
}

void Mat3::operator=( Mat3 other )
{
	m_data[0] = other.m_data[0];
	m_data[1] = other.m_data[1];
	m_data[2] = other.m_data[2];
	m_data[3] = other.m_data[3];
	m_data[4] = other.m_data[4];
	m_data[5] = other.m_data[5];
	m_data[6] = other.m_data[6];
	m_data[7] = other.m_data[7];
	m_data[8] = other.m_data[8];
}

bool Mat3::operator==( Mat3 other )
{
	for ( unsigned int i = 0; i < 9; i++ )
	{
		if ( m_data[i] != other.m_data[i] )
		{
			return false;
		}
	}
	return true;
}

Vec3 Mat3::GetRowVec( const unsigned int row ) const
{
	assert( row < 3 );
	Vec3 returnVec;
	returnVec[0] = m_data[row * 3];
	returnVec[1] = m_data[row * 3 + 1];
	returnVec[2] = m_data[row * 3 + 2];
	return returnVec;
}

const Vec3 &Mat3::RowVec( const unsigned int row ) const
{
	assert( row < 3 );
	return reinterpret_cast< const Vec3 & >( m_data[row * 3] );
}

void Mat3::SetRowVec( const unsigned int row, const Vec3 *vec )
{
	assert( row < 3 );
	m_data[row * 3] = ( *vec )[0];
	m_data[row * 3 + 1] = ( *vec )[1];
	m_data[row * 3 + 2] = ( *vec )[2];
}

Vec3 Mat3::GetColVec( const unsigned int col ) const
{
	assert( col < 3 );
	Vec3 returnVec;
	returnVec[0] = m_data[col * 3];
	returnVec[1] = m_data[col + 1 * 3];
	returnVec[2] = m_data[col + 1 * 3];
	return returnVec;
}

void Mat3::SetColVec( const unsigned int col, const Vec3 *vec )
{
	assert( col < 3 );
	m_data[col * 3] = ( *vec )[0];
	m_data[col + 1 * 3] = ( *vec )[1];
	m_data[col + 2 * 3] = ( *vec )[2];
}

float Mat3::GetComponent( unsigned int m, unsigned int n ) const
{
	assert( m < 3 && n < 3 );
	return m_data[m * 3 + n];
}

void Mat3::SetComponent( unsigned int m, unsigned int n, float val )
{
	assert( m < 3 && n < 3 );
	m_data[m * 3 + n] = val;
}

Mat3 Mat3::operator+( Mat3 other ) const
{
	Mat3 returnMat;
	returnMat.m_data[0] = m_data[0] + other.m_data[0];
	returnMat.m_data[1] = m_data[1] + other.m_data[1];
	returnMat.m_data[2] = m_data[2] + other.m_data[2];
	returnMat.m_data[3] = m_data[3] + other.m_data[3];
	returnMat.m_data[4] = m_data[4] + other.m_data[4];
	returnMat.m_data[5] = m_data[5] + other.m_data[5];
	returnMat.m_data[6] = m_data[6] + other.m_data[6];
	returnMat.m_data[7] = m_data[7] + other.m_data[7];
	returnMat.m_data[8] = m_data[8] + other.m_data[8];
	return returnMat;
}

void Mat3::operator+=( Mat3 other )
{
	m_data[0] += other.m_data[0];
	m_data[1] += other.m_data[1];
	m_data[2] += other.m_data[2];
	m_data[3] += other.m_data[3];
	m_data[4] += other.m_data[4];
	m_data[5] += other.m_data[5];
	m_data[6] += other.m_data[6];
	m_data[7] += other.m_data[7];
	m_data[8] += other.m_data[8];
}

Mat3 Mat3::operator-( Mat3 other ) const
{
	Mat3 returnMat;
	returnMat.m_data[0] = m_data[0] - other.m_data[0];
	returnMat.m_data[1] = m_data[1] - other.m_data[1];
	returnMat.m_data[2] = m_data[2] - other.m_data[2];
	returnMat.m_data[3] = m_data[3] - other.m_data[3];
	returnMat.m_data[4] = m_data[4] - other.m_data[4];
	returnMat.m_data[5] = m_data[5] - other.m_data[5];
	returnMat.m_data[6] = m_data[6] - other.m_data[6];
	returnMat.m_data[7] = m_data[7] - other.m_data[7];
	returnMat.m_data[8] = m_data[8] - other.m_data[8];
	return returnMat;
}

void Mat3::operator-=( Mat3 other )
{
	m_data[0] -= other.m_data[0];
	m_data[1] -= other.m_data[1];
	m_data[2] -= other.m_data[2];
	m_data[3] -= other.m_data[3];
	m_data[4] -= other.m_data[4];
	m_data[5] -= other.m_data[5];
	m_data[6] -= other.m_data[6];
	m_data[7] -= other.m_data[7];
	m_data[8] -= other.m_data[8];
}

Mat3 Mat3::operator*( float scalar ) const
{
	Mat3 returnMat;
	returnMat.m_data[0] = m_data[0] * scalar;
	returnMat.m_data[1] = m_data[1] * scalar;
	returnMat.m_data[2] = m_data[2] * scalar;
	returnMat.m_data[3] = m_data[3] * scalar;
	returnMat.m_data[4] = m_data[4] * scalar;
	returnMat.m_data[5] = m_data[5] * scalar;
	returnMat.m_data[6] = m_data[6] * scalar;
	returnMat.m_data[7] = m_data[7] * scalar;
	returnMat.m_data[8] = m_data[8] * scalar;
	return returnMat;
}

void Mat3::operator*=( float scalar )
{
	m_data[0] *= scalar;
	m_data[1] *= scalar;
	m_data[2] *= scalar;
	m_data[3] *= scalar;
	m_data[4] *= scalar;
	m_data[5] *= scalar;
	m_data[6] *= scalar;
	m_data[7] *= scalar;
	m_data[8] *= scalar;
}

Vec3 Mat3::operator*( Vec3 vec ) const
{
	Vec3 returnVec;
	returnVec[0] = m_data[0] * vec[0] + m_data[1] * vec[1] + m_data[2] * vec[2];
	returnVec[1] = m_data[3] * vec[0] + m_data[4] * vec[1] + m_data[5] * vec[2];
	returnVec[2] = m_data[6] * vec[0] + m_data[7] * vec[1] + m_data[8] * vec[2];
	return returnVec;
}

Mat3 Mat3::operator*( Mat3 other ) const
{
	Mat3 returnMat;
	returnMat[0] = m_data[0] * other.m_data[0] + m_data[1] * other.m_data[3] + m_data[2] * other.m_data[6];
	returnMat[1] = m_data[0] * other.m_data[1] + m_data[1] * other.m_data[4] + m_data[2] * other.m_data[7];
	returnMat[2] = m_data[0] * other.m_data[2] + m_data[1] * other.m_data[5] + m_data[2] * other.m_data[8];
	returnMat[3] = m_data[3] * other.m_data[0] + m_data[4] * other.m_data[3] + m_data[5] * other.m_data[6];
	returnMat[4] = m_data[3] * other.m_data[1] + m_data[4] * other.m_data[4] + m_data[5] * other.m_data[7];
	returnMat[5] = m_data[3] * other.m_data[2] + m_data[4] * other.m_data[5] + m_data[5] * other.m_data[8];
	returnMat[6] = m_data[6] * other.m_data[0] + m_data[7] * other.m_data[3] + m_data[8] * other.m_data[6];
	returnMat[7] = m_data[6] * other.m_data[1] + m_data[7] * other.m_data[4] + m_data[8] * other.m_data[7];
	returnMat[8] = m_data[6] * other.m_data[2] + m_data[7] * other.m_data[5] + m_data[8] * other.m_data[8];
	return returnMat;
}

void Mat3::operator*=( Mat3 other )
{
	const float a = m_data[0] * other.m_data[0] + m_data[1] * other.m_data[3] + m_data[2] * other.m_data[6];
	const float b = m_data[0] * other.m_data[1] + m_data[1] * other.m_data[4] + m_data[2] * other.m_data[7];
	const float c = m_data[0] * other.m_data[2] + m_data[1] * other.m_data[5] + m_data[2] * other.m_data[8];
	const float d = m_data[3] * other.m_data[0] + m_data[4] * other.m_data[3] + m_data[5] * other.m_data[6];
	const float e = m_data[3] * other.m_data[1] + m_data[4] * other.m_data[4] + m_data[5] * other.m_data[7];
	const float f = m_data[3] * other.m_data[2] + m_data[4] * other.m_data[5] + m_data[5] * other.m_data[8];
	const float g = m_data[6] * other.m_data[0] + m_data[7] * other.m_data[3] + m_data[8] * other.m_data[6];
	const float h = m_data[6] * other.m_data[1] + m_data[7] * other.m_data[4] + m_data[8] * other.m_data[7];
	const float i = m_data[6] * other.m_data[2] + m_data[7] * other.m_data[5] + m_data[8] * other.m_data[8];
	m_data[0] = a;
	m_data[1] = b;
	m_data[2] = c;
	m_data[3] = d;
	m_data[4] = e;
	m_data[5] = f;
	m_data[6] = g;
	m_data[7] = h;
	m_data[8] = i;
}

Mat3 Mat3::operator/( float scalar ) const
{
	Mat3 returnMat;
	returnMat.m_data[0] = m_data[0] / scalar;
	returnMat.m_data[1] = m_data[1] / scalar;
	returnMat.m_data[2] = m_data[2] / scalar;
	returnMat.m_data[3] = m_data[3] / scalar;
	returnMat.m_data[4] = m_data[4] / scalar;
	returnMat.m_data[5] = m_data[5] / scalar;
	returnMat.m_data[6] = m_data[6] / scalar;
	returnMat.m_data[7] = m_data[7] / scalar;
	returnMat.m_data[8] = m_data[8] / scalar;
	return returnMat;
}

void Mat3::operator/=( float scalar )
{
	m_data[0] /= scalar;
	m_data[1] /= scalar;
	m_data[2] /= scalar;
	m_data[3] /= scalar;
	m_data[4] /= scalar;
	m_data[5] /= scalar;
	m_data[6] /= scalar;
	m_data[7] /= scalar;
	m_data[8] /= scalar;
}

float Mat3::Determinant() const
{
	return m_data[0] * m_data[4] * m_data[8] +
		m_data[1] * m_data[5] * m_data[6] +
		m_data[2] * m_data[3] * m_data[7] -
		m_data[6] * m_data[4] * m_data[2] -
		m_data[7] * m_data[5] * m_data[0] -
		m_data[8] * m_data[3] * m_data[1];
}

bool Mat3::Inverse( Mat3 *inv )
{
	const float det = Determinant();

	if ( fabsf( det ) < FLT_EPSILON )
	{
		return false;
	}

	( *inv )[0] = m_data[4] * m_data[8] - m_data[5] * m_data[7];
	( *inv )[1] = m_data[2] * m_data[7] - m_data[1] * m_data[8];
	( *inv )[2] = m_data[1] * m_data[5] - m_data[2] * m_data[4];
	( *inv )[4] = m_data[0] * m_data[8] - m_data[2] * m_data[6];
	( *inv )[5] = m_data[2] * m_data[3] - m_data[0] * m_data[5];
	( *inv )[7] = m_data[1] * m_data[6] - m_data[0] * m_data[7];
	( *inv )[8] = m_data[0] * m_data[4] - m_data[1] * m_data[3];

	( *inv ) *= 1.0f / det;

	return true;
}

Mat3 Mat3::Inverse() const
{
	const float det = Determinant();
	assert( fabsf( det ) > FLT_EPSILON );

	Mat3 inv;

	inv[0] = m_data[4] * m_data[8] - m_data[5] * m_data[7];
	inv[1] = m_data[2] * m_data[7] - m_data[1] * m_data[8];
	inv[2] = m_data[1] * m_data[5] - m_data[2] * m_data[4];
	inv[4] = m_data[0] * m_data[8] - m_data[2] * m_data[6];
	inv[5] = m_data[2] * m_data[3] - m_data[0] * m_data[5];
	inv[7] = m_data[1] * m_data[6] - m_data[0] * m_data[7];
	inv[8] = m_data[0] * m_data[4] - m_data[1] * m_data[3];

	inv *= 1.0f / det;

	return inv;
}

Mat3 Mat3::Transpose() const
{
	Mat3 returnMat;
	returnMat.m_data[0] = m_data[0];
	returnMat.m_data[1] = m_data[3];
	returnMat.m_data[2] = m_data[6];
	returnMat.m_data[3] = m_data[1];
	returnMat.m_data[4] = m_data[4];
	returnMat.m_data[5] = m_data[7];
	returnMat.m_data[6] = m_data[2];
	returnMat.m_data[7] = m_data[5];
	returnMat.m_data[8] = m_data[8];
	return returnMat;
}

void Mat3::Transposed()
{
	for ( unsigned int i = 0; i < 3; i++ )
	{
		for ( unsigned int j = 0; j < 3; j++ )
		{
			float temp = m_data[j * 3 + i];
			m_data[j * 3 + i] = m_data[i * 3 + j];
			m_data[i * 3 + j] = temp;
		}
	}
}

MatN Mat3::as_MatN() const
{
	MatN returnMat( 3 );
	returnMat[0] = m_data[0];
	returnMat[1] = m_data[1];
	returnMat[2] = m_data[2];
	returnMat[3] = m_data[3];
	returnMat[4] = m_data[4];
	returnMat[5] = m_data[5];
	returnMat[6] = m_data[6];
	returnMat[7] = m_data[7];
	returnMat[8] = m_data[8];
	return returnMat;
}

Mat2 Mat3::as_Mat2() const
{
	Mat2 returnMat;
	returnMat[0] = m_data[0];
	returnMat[1] = m_data[1];
	returnMat[3] = m_data[3];
	returnMat[4] = m_data[4];
	return returnMat;
}

Mat4 Mat3::as_Mat4() const
{
	Mat4 returnMat;
	returnMat[0] = m_data[0];
	returnMat[1] = m_data[1];
	returnMat[2] = m_data[2];
	returnMat[3] = 0.0f;
	returnMat[4] = m_data[3];
	returnMat[5] = m_data[4];
	returnMat[6] = m_data[5];
	returnMat[7] = 0.0f;
	returnMat[8] = m_data[6];
	returnMat[9] = m_data[7];
	returnMat[10] = m_data[8];
	returnMat[11] = 0.0f;
	returnMat[12] = 0.0f;
	returnMat[13] = 0.0f;
	returnMat[14] = 0.0f;
	returnMat[15] = 1.0f;
	return returnMat;
}


/*
================================
Mat4
================================
*/
Mat4::Mat4()
{
	for ( unsigned int i = 0; i < 16; i++ )
	{
		m_data[i] = 0.0f;
	}
	m_data[0] = 1.0f;
	m_data[5] = 1.0f;
	m_data[10] = 1.0f;
	m_data[15] = 1.0f;
}

Mat4::Mat4( float v )
{
	for ( unsigned int i = 0; i < 16; i++ )
	{
		m_data[i] = v;
	}
}

Mat4::Mat4( const Vec4 &vA, const Vec4 &vB, const Vec4 &vC, const Vec4 &vD )
{
	m_data[0] = vA[0];
	m_data[1] = vA[1];
	m_data[2] = vA[2];
	m_data[3] = vA[3];
	m_data[4] = vB[0];
	m_data[5] = vB[1];
	m_data[6] = vB[2];
	m_data[7] = vB[3];
	m_data[8] = vC[0];
	m_data[9] = vC[1];
	m_data[10] = vC[2];
	m_data[11] = vC[3];
	m_data[12] = vD[0];
	m_data[13] = vD[1];
	m_data[14] = vD[2];
	m_data[15] = vD[3];
}

Mat4::Mat4( const float *data )
{
	for ( unsigned int i = 0; i < 16; i++ )
	{
		m_data[i] = data[i];
	}
}

Mat4::Mat4( const Mat4 &mat )
{
	for ( unsigned int i = 0; i < 16; i++ )
	{
		m_data[i] = mat.m_data[i];
	}
}

void Mat4::operator=( Mat4 other )
{
	for ( unsigned int i = 0; i < 16; i++ )
	{
		m_data[i] = other.m_data[i];
	}
}

bool Mat4::operator==( Mat4 other )
{
	for ( unsigned int i = 0; i < 16; i++ )
	{
		if ( m_data[i] != other.m_data[i] )
		{
			return false;
		}
	}
	return true;
}

Vec4 Mat4::GetRowVec( const unsigned int row ) const
{
	assert( row < 4 );
	Vec4 returnVec;
	returnVec[0] = m_data[row * 4];
	returnVec[1] = m_data[row * 4 + 1];
	returnVec[2] = m_data[row * 4 + 2];
	returnVec[3] = m_data[row * 4 + 3];
	return returnVec;
}

void Mat4::SetRowVec( const unsigned int row, const Vec4 *vec )
{
	assert( row < 4 );
	m_data[row * 4] = ( *vec )[0];
	m_data[row * 4 + 1] = ( *vec )[1];
	m_data[row * 4 + 2] = ( *vec )[2];
	m_data[row * 4 + 3] = ( *vec )[3];
}

Vec4 Mat4::GetColVec( const unsigned int col ) const
{
	assert( col < 4 );
	Vec4 returnVec;
	returnVec[0] = m_data[col * 4];
	returnVec[1] = m_data[col + 1 * 4];
	returnVec[2] = m_data[col + 2 * 4];
	returnVec[3] = m_data[col + 3 * 4];
	return returnVec;
}

void Mat4::SetColVec( const unsigned int col, const Vec4 *vec )
{
	assert( col < 4 );
	m_data[col * 4] = ( *vec )[0];
	m_data[col + 1 * 4] = ( *vec )[1];
	m_data[col + 2 * 4] = ( *vec )[2];
	m_data[col + 3 * 4] = ( *vec )[3];
}

float Mat4::GetComponent( unsigned int m, unsigned int n ) const
{
	assert( m < 4 && n < 4 );
	return m_data[m * 4 + n];
}

void Mat4::SetComponent( unsigned int m, unsigned int n, float val )
{
	assert( m < 4 && n < 4 );
	m_data[m * 4 + n] = val;
}

Mat4 Mat4::operator+( Mat4 other ) const
{
	Mat4 returnMat;
	for ( unsigned int i = 0; i < 16; i++ )
	{
		returnMat.m_data[i] = m_data[i] + other.m_data[i];
	}
	return returnMat;
}

void Mat4::operator+=( Mat4 other )
{
	for ( unsigned int i = 0; i < 16; i++ )
	{
		m_data[i] += other.m_data[i];
	}
}

Mat4 Mat4::operator-( Mat4 other ) const
{
	Mat4 returnMat;
	for ( unsigned int i = 0; i < 16; i++ )
	{
		returnMat.m_data[i] = m_data[i] - other.m_data[i];
	}
	return returnMat;
}

void Mat4::operator-=( Mat4 other )
{
	for ( unsigned int i = 0; i < 16; i++ )
	{
		m_data[i] -= other.m_data[i];
	}
}

Mat4 Mat4::operator*( float scalar ) const
{
	Mat4 returnMat;
	for ( unsigned int i = 0; i < 16; i++ )
	{
		returnMat.m_data[i] = m_data[i] * scalar;
	}
	return returnMat;
}

void Mat4::operator*=( float scalar )
{
	for ( unsigned int i = 0; i < 16; i++ )
	{
		m_data[i] *= scalar;
	}
}

Vec4 Mat4::operator*( Vec4 vec ) const
{
	Vec4 returnVec;
	returnVec[0] = m_data[0] * vec[0] + m_data[1] * vec[1] + m_data[2] * vec[2] + m_data[3] * vec[3];
	returnVec[1] = m_data[4] * vec[0] + m_data[5] * vec[1] + m_data[6] * vec[2] + m_data[7] * vec[3];
	returnVec[2] = m_data[8] * vec[0] + m_data[9] * vec[1] + m_data[10] * vec[2] + m_data[11] * vec[3];
	returnVec[3] = m_data[12] * vec[0] + m_data[13] * vec[1] + m_data[14] * vec[2] + m_data[15] * vec[3];
	return returnVec;
}

Mat4 Mat4::operator*( Mat4 other ) const
{
	Mat4 returnMat;

	returnMat[0] = m_data[0] * other.m_data[0] + m_data[1] * other.m_data[4] + m_data[2] * other.m_data[8] + m_data[3] * other.m_data[12];
	returnMat[1] = m_data[0] * other.m_data[1] + m_data[1] * other.m_data[5] + m_data[2] * other.m_data[9] + m_data[3] * other.m_data[13];
	returnMat[2] = m_data[0] * other.m_data[2] + m_data[1] * other.m_data[6] + m_data[2] * other.m_data[10] + m_data[3] * other.m_data[14];
	returnMat[3] = m_data[0] * other.m_data[3] + m_data[1] * other.m_data[7] + m_data[2] * other.m_data[11] + m_data[3] * other.m_data[15];

	returnMat[4] = m_data[4] * other.m_data[0] + m_data[5] * other.m_data[4] + m_data[6] * other.m_data[8] + m_data[7] * other.m_data[12];
	returnMat[5] = m_data[4] * other.m_data[1] + m_data[5] * other.m_data[5] + m_data[6] * other.m_data[9] + m_data[7] * other.m_data[13];
	returnMat[6] = m_data[4] * other.m_data[2] + m_data[5] * other.m_data[6] + m_data[6] * other.m_data[10] + m_data[7] * other.m_data[14];
	returnMat[7] = m_data[4] * other.m_data[3] + m_data[5] * other.m_data[7] + m_data[6] * other.m_data[11] + m_data[7] * other.m_data[15];

	returnMat[8] = m_data[8] * other.m_data[0] + m_data[9] * other.m_data[4] + m_data[10] * other.m_data[8] + m_data[11] * other.m_data[12];
	returnMat[9] = m_data[8] * other.m_data[1] + m_data[9] * other.m_data[5] + m_data[10] * other.m_data[9] + m_data[11] * other.m_data[13];
	returnMat[10] = m_data[8] * other.m_data[2] + m_data[9] * other.m_data[6] + m_data[10] * other.m_data[10] + m_data[11] * other.m_data[14];
	returnMat[11] = m_data[8] * other.m_data[3] + m_data[9] * other.m_data[7] + m_data[10] * other.m_data[11] + m_data[11] * other.m_data[15];

	returnMat[12] = m_data[12] * other.m_data[0] + m_data[13] * other.m_data[4] + m_data[14] * other.m_data[8] + m_data[15] * other.m_data[12];
	returnMat[13] = m_data[12] * other.m_data[1] + m_data[13] * other.m_data[5] + m_data[14] * other.m_data[9] + m_data[15] * other.m_data[13];
	returnMat[14] = m_data[12] * other.m_data[2] + m_data[13] * other.m_data[6] + m_data[14] * other.m_data[10] + m_data[15] * other.m_data[14];
	returnMat[15] = m_data[12] * other.m_data[3] + m_data[13] * other.m_data[7] + m_data[14] * other.m_data[11] + m_data[15] * other.m_data[15];

	return returnMat;
}

void Mat4::operator*=( Mat4 other )
{
	float a = m_data[0] * other.m_data[0] + m_data[1] * other.m_data[4] + m_data[2] * other.m_data[8] + m_data[3] * other.m_data[12];
	float b = m_data[0] * other.m_data[1] + m_data[1] * other.m_data[5] + m_data[2] * other.m_data[9] + m_data[3] * other.m_data[13];
	float c = m_data[0] * other.m_data[2] + m_data[1] * other.m_data[6] + m_data[2] * other.m_data[10] + m_data[3] * other.m_data[14];
	float d = m_data[0] * other.m_data[3] + m_data[1] * other.m_data[7] + m_data[2] * other.m_data[11] + m_data[3] * other.m_data[15];

	float e = m_data[4] * other.m_data[0] + m_data[5] * other.m_data[4] + m_data[6] * other.m_data[8] + m_data[7] * other.m_data[12];
	float f = m_data[4] * other.m_data[1] + m_data[5] * other.m_data[5] + m_data[6] * other.m_data[9] + m_data[7] * other.m_data[13];
	float g = m_data[4] * other.m_data[2] + m_data[5] * other.m_data[6] + m_data[6] * other.m_data[10] + m_data[7] * other.m_data[14];
	float h = m_data[4] * other.m_data[3] + m_data[5] * other.m_data[7] + m_data[6] * other.m_data[11] + m_data[7] * other.m_data[15];

	float i = m_data[8] * other.m_data[0] + m_data[9] * other.m_data[4] + m_data[10] * other.m_data[8] + m_data[11] * other.m_data[12];
	float j = m_data[8] * other.m_data[1] + m_data[9] * other.m_data[5] + m_data[10] * other.m_data[9] + m_data[11] * other.m_data[13];
	float k = m_data[8] * other.m_data[2] + m_data[9] * other.m_data[6] + m_data[10] * other.m_data[10] + m_data[11] * other.m_data[14];
	float l = m_data[8] * other.m_data[3] + m_data[9] * other.m_data[7] + m_data[10] * other.m_data[11] + m_data[11] * other.m_data[15];

	float m = m_data[12] * other.m_data[0] + m_data[13] * other.m_data[4] + m_data[14] * other.m_data[8] + m_data[15] * other.m_data[12];
	float n = m_data[12] * other.m_data[1] + m_data[13] * other.m_data[5] + m_data[14] * other.m_data[9] + m_data[15] * other.m_data[13];
	float o = m_data[12] * other.m_data[2] + m_data[13] * other.m_data[6] + m_data[14] * other.m_data[10] + m_data[15] * other.m_data[14];
	float p = m_data[12] * other.m_data[3] + m_data[13] * other.m_data[7] + m_data[14] * other.m_data[11] + m_data[15] * other.m_data[15];

	m_data[0] = a;
	m_data[1] = b;
	m_data[2] = c;
	m_data[3] = d;

	m_data[4] = e;
	m_data[5] = f;
	m_data[6] = g;
	m_data[7] = h;

	m_data[8] = i;
	m_data[9] = j;
	m_data[10] = k;
	m_data[11] = l;

	m_data[12] = m;
	m_data[13] = n;
	m_data[14] = o;
	m_data[15] = p;
}

Mat4 Mat4::operator/( float scalar ) const
{
	Mat4 returnMat;
	for ( unsigned int i = 0; i < 16; i++ )
	{
		returnMat.m_data[i] = m_data[i] / scalar;
	}
	return returnMat;
}

void Mat4::operator/=( float scalar )
{
	for ( unsigned int i = 0; i < 16; i++ )
	{
		m_data[i] *= scalar;
	}
}

float Mat4::Determinant()
{
	return m_data[3] * m_data[6] * m_data[9] * m_data[12] - m_data[2] * m_data[7] * m_data[9] * m_data[12] -
		m_data[3] * m_data[5] * m_data[10] * m_data[12] + m_data[1] * m_data[7] * m_data[10] * m_data[12] +
		m_data[2] * m_data[5] * m_data[11] * m_data[12] - m_data[1] * m_data[6] * m_data[11] * m_data[12] -
		m_data[3] * m_data[6] * m_data[8] * m_data[13] + m_data[2] * m_data[7] * m_data[8] * m_data[13] +
		m_data[3] * m_data[4] * m_data[10] * m_data[13] - m_data[0] * m_data[7] * m_data[10] * m_data[13] -
		m_data[2] * m_data[4] * m_data[11] * m_data[13] + m_data[0] * m_data[6] * m_data[11] * m_data[13] +
		m_data[3] * m_data[5] * m_data[8] * m_data[14] - m_data[1] * m_data[7] * m_data[8] * m_data[14] -
		m_data[3] * m_data[4] * m_data[9] * m_data[14] + m_data[0] * m_data[7] * m_data[9] * m_data[14] +
		m_data[1] * m_data[4] * m_data[11] * m_data[14] - m_data[0] * m_data[5] * m_data[11] * m_data[14] -
		m_data[2] * m_data[5] * m_data[8] * m_data[15] + m_data[1] * m_data[6] * m_data[8] * m_data[15] +
		m_data[2] * m_data[4] * m_data[9] * m_data[15] - m_data[0] * m_data[6] * m_data[9] * m_data[15] -
		m_data[1] * m_data[4] * m_data[10] * m_data[15] + m_data[0] * m_data[5] * m_data[10] * m_data[15];
}

bool Mat4::Inverse( Mat4 *inv )
{
//2x2 sub-determinants required to calculate 4x4 determinant
	float det2_01_01 = m_data[0] * m_data[5] - m_data[1] * m_data[4];
	float det2_01_02 = m_data[0] * m_data[6] - m_data[2] * m_data[4];
	float det2_01_03 = m_data[0] * m_data[8] - m_data[3] * m_data[4];
	float det2_01_12 = m_data[1] * m_data[6] - m_data[2] * m_data[5];
	float det2_01_13 = m_data[1] * m_data[7] - m_data[3] * m_data[5];
	float det2_01_23 = m_data[2] * m_data[7] - m_data[3] * m_data[6];

	//3x3 sub-determinants required to calculate 4x4 determinant
	float det3_201_012 = m_data[8] * det2_01_12 - m_data[9] * det2_01_02 + m_data[10] * det2_01_01;
	float det3_201_013 = m_data[8] * det2_01_13 - m_data[9] * det2_01_03 + m_data[11] * det2_01_01;
	float det3_201_023 = m_data[8] * det2_01_23 - m_data[10] * det2_01_03 + m_data[11] * det2_01_02;
	float det3_201_123 = m_data[9] * det2_01_23 - m_data[10] * det2_01_13 + m_data[11] * det2_01_12;

	const float det = ( -det3_201_123 * m_data[12] + det3_201_023 * m_data[13] - det3_201_013 * m_data[14] + det3_201_012 * m_data[15] );

	if ( ( fabsf( det ) < FLT_EPSILON ) )
	{
		return false;
	}

	//remaining 2x2 sub-determinants
	float det2_03_01 = m_data[0] * m_data[13] - m_data[1] * m_data[12];
	float det2_03_02 = m_data[0] * m_data[14] - m_data[2] * m_data[12];
	float det2_03_03 = m_data[0] * m_data[15] - m_data[3] * m_data[12];
	float det2_03_12 = m_data[1] * m_data[14] - m_data[2] * m_data[13];
	float det2_03_13 = m_data[1] * m_data[15] - m_data[3] * m_data[13];
	float det2_03_23 = m_data[2] * m_data[15] - m_data[3] * m_data[14];

	float det2_13_01 = m_data[4] * m_data[13] - m_data[5] * m_data[12];
	float det2_13_02 = m_data[4] * m_data[14] - m_data[6] * m_data[12];
	float det2_13_03 = m_data[4] * m_data[15] - m_data[7] * m_data[12];
	float det2_13_12 = m_data[5] * m_data[14] - m_data[6] * m_data[13];
	float det2_13_13 = m_data[5] * m_data[15] - m_data[7] * m_data[13];
	float det2_13_23 = m_data[6] * m_data[15] - m_data[7] * m_data[13];

	//remaining 3x3 sub-determinants
	float det3_203_012 = m_data[8] * det2_03_12 - m_data[9] * det2_03_02 + m_data[10] * det2_03_01;
	float det3_203_013 = m_data[8] * det2_03_13 - m_data[9] * det2_03_03 + m_data[11] * det2_03_01;
	float det3_203_023 = m_data[8] * det2_03_23 - m_data[10] * det2_03_03 + m_data[11] * det2_03_02;
	float det3_203_123 = m_data[9] * det2_03_23 - m_data[10] * det2_03_13 + m_data[11] * det2_03_12;

	float det3_213_012 = m_data[8] * det2_13_12 - m_data[9] * det2_13_02 + m_data[10] * det2_13_01;
	float det3_213_013 = m_data[8] * det2_13_13 - m_data[9] * det2_13_03 + m_data[11] * det2_13_01;
	float det3_213_023 = m_data[8] * det2_13_23 - m_data[10] * det2_13_03 + m_data[11] * det2_13_02;
	float det3_213_123 = m_data[9] * det2_13_23 - m_data[10] * det2_13_13 + m_data[11] * det2_13_12;

	float det3_301_012 = m_data[12] * det2_01_12 - m_data[13] * det2_01_02 + m_data[14] * det2_01_01;
	float det3_301_013 = m_data[12] * det2_01_13 - m_data[13] * det2_01_03 + m_data[15] * det2_01_01;
	float det3_301_023 = m_data[12] * det2_01_23 - m_data[13] * det2_01_03 + m_data[15] * det2_01_02;
	float det3_301_123 = m_data[13] * det2_01_23 - m_data[14] * det2_01_13 + m_data[15] * det2_01_12;

	( *inv )[0] = -det3_213_123 / det;
	( *inv )[1] = det3_213_023 / det;
	( *inv )[2] = -det3_213_013 / det;
	( *inv )[3] = det3_213_012 / det;

	( *inv )[4] = det3_203_123 / det;
	( *inv )[5] = -det3_203_023 / det;
	( *inv )[6] = det3_203_013 / det;
	( *inv )[7] = -det3_203_012 / det;

	( *inv )[8] = det3_301_123 / det;
	( *inv )[9] = -det3_301_023 / det;
	( *inv )[10] = det3_301_013 / det;
	( *inv )[11] = -det3_301_012 / det;

	( *inv )[12] = -det3_201_123 / det;
	( *inv )[13] = det3_201_023 / det;
	( *inv )[14] = -det3_201_013 / det;
	( *inv )[15] = det3_201_012 / det;

	return true;
}

Mat4 Mat4::Transpose() const
{
	Mat4 returnMat;

	returnMat.m_data[0] = m_data[0];
	returnMat.m_data[1] = m_data[4];
	returnMat.m_data[2] = m_data[8];
	returnMat.m_data[3] = m_data[12];

	returnMat.m_data[4] = m_data[1];
	returnMat.m_data[5] = m_data[5];
	returnMat.m_data[6] = m_data[9];
	returnMat.m_data[7] = m_data[13];

	returnMat.m_data[8] = m_data[2];
	returnMat.m_data[9] = m_data[6];
	returnMat.m_data[10] = m_data[10];
	returnMat.m_data[11] = m_data[14];

	returnMat.m_data[12] = m_data[3];
	returnMat.m_data[13] = m_data[7];
	returnMat.m_data[14] = m_data[11];
	returnMat.m_data[15] = m_data[15];

	return returnMat;
}

void Mat4::Transposed()
{
	for ( unsigned int i = 0; i < 4; i++ )
	{
		for ( unsigned int j = 0; j < 4; j++ )
		{
			float temp = m_data[j * 4 + i];
			m_data[j * 4 + i] = m_data[i * 4 + j];
			m_data[i * 4 + j] = temp;
		}
	}
}

void Mat4::Translate( const Vec3 &vPosition )
{
	m_data[3] = vPosition[0];
	m_data[7] = vPosition[1];
	m_data[11] = vPosition[2];

}

void Mat4::LookAt( const Vec3 look, const Vec3 up, const Vec3 pos )
{
	Vec3 r = look.Cross( up ).Normal();
	Vec3 d = look.Normal();
	Vec3 u = up.Normal();

	m_data[0] = r[0];
	m_data[1] = r[1];
	m_data[2] = r[2];
	m_data[3] = -r.Dot( pos );

	m_data[4] = u[0];
	m_data[5] = u[1];
	m_data[6] = u[2];
	m_data[7] = -u.Dot( pos );

	m_data[8] = -d[0];
	m_data[9] = -d[1];
	m_data[10] = -d[2];
	m_data[11] = d.Dot( pos );

	m_data[12] = 0.0f;
	m_data[13] = 0.0f;
	m_data[14] = 0.0f;
	m_data[15] = 1.0f;
}

void Mat4::Perspective( const float verticalFOV, const float aspect, const float near, const float far )
{
	const float tanHalfFOV = tanf( verticalFOV / 2.0f );

	m_data[0] = 1.0f / ( tanHalfFOV * aspect );
	m_data[1] = 0.0f;
	m_data[2] = 0.0f;
	m_data[3] = 0.0f;

	m_data[4] = 0.0f;
	m_data[5] = 1.0f / tanHalfFOV;
	m_data[6] = 0.0f;
	m_data[7] = 0.0f;

	m_data[8] = 0.0f;
	m_data[9] = 0.0f;
	m_data[10] = ( -far - near ) / ( far - near );
	m_data[11] = -1.0f;

	m_data[12] = 0.0f;
	m_data[13] = 0.0f;
	m_data[14] = ( -2.0f * far * near ) / ( far - near );
	m_data[15] = 0.0f;
}

void Mat4::Orthographic( const float left, const float right, const float bottom, const float top )
{
	const float near = -1.0f;
	const float far = 1.0f;
	Orthographic( left, right, bottom, top, near, far );
}

void Mat4::Orthographic( const float left, const float right, const float bottom, const float top, const float near, const float far )
{
	m_data[0] = 2.0f / ( right - left );
	m_data[1] = 0.0f;
	m_data[2] = 0.0f;
	m_data[3] = 0.0f;

	m_data[4] = 0.0f;
	m_data[5] = 2.0f / ( top - bottom );
	m_data[6] = 0.0f;
	m_data[7] = 0.0f;

	m_data[8] = 0.0f;
	m_data[9] = 0.0f;
	m_data[10] = -2.0f / ( far - near );
	m_data[11] = 0.0f;

	m_data[12] = ( right + left ) / ( right - left ) * -1.0f;
	m_data[13] = ( top + bottom ) / ( top - bottom ) * -1.0f;
	m_data[14] = ( far + near ) / ( far - near ) * -1.0f;
	m_data[15] = 1.0f;
}

MatN Mat4::as_MatN() const
{
	MatN returnMat( 4 );
	memcpy( returnMat.m_data, m_data, sizeof( float ) );
	return returnMat;
}

Mat2 Mat4::as_Mat2() const
{
	Mat2 returnMat;
	returnMat[0] = m_data[0];
	returnMat[1] = m_data[1];
	returnMat[3] = m_data[4];
	returnMat[4] = m_data[5];
	return returnMat;
}

Mat3 Mat4::as_Mat3() const
{
	Mat3 returnMat;
	returnMat[0] = m_data[0];
	returnMat[1] = m_data[1];
	returnMat[2] = m_data[2];
	returnMat[3] = m_data[4];
	returnMat[4] = m_data[5];
	returnMat[5] = m_data[6];
	returnMat[6] = m_data[8];
	returnMat[7] = m_data[9];
	returnMat[8] = m_data[10];
	return returnMat;
}