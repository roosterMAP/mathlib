#include "..\include\Vector.h"
#include "..\include\Matrix.h"
#include <math.h>
#include "assert.h"

float to_degrees( const float rads ) {
	return ( float )(rads * 180.0 / PI);
}

float to_radians( const float deg ) {
	return ( float )(deg * PI / 180.0);
}

/*
================================
VecN
================================
*/
VecN::VecN( unsigned int size ) {
	m_size = size;
	m_data = new float[m_size];
	for ( unsigned int i = 0; i < m_size; i++ ) {
		m_data[i] = 0.0f;
	}
}

VecN::VecN( unsigned int size, float n ) {
	m_size = size;
	m_data = new float[m_size];
	for ( unsigned int i = 0; i < m_size; i++ ) {
		m_data[i] = n;
	}
}

VecN::VecN( unsigned int size, const float * data ) {
	m_size = size;
	m_data = new float[m_size];
	for ( unsigned int i = 0; i < m_size; i++ ) {
		m_data[i] = data[i];
	}
}

VecN::VecN( const VecN &vec ) {
	m_size = vec.m_size;
	m_data = new float[m_size];
	for ( unsigned int i = 0; i < m_size; i++ ) {
		m_data[i] = vec.m_data[i];
	}
}

void VecN::operator=( VecN other ) {
	if ( m_size != other.m_size ) {
		delete[] m_data;
		m_data = new float[ other.m_size ];
	}

	for ( unsigned int i = 0; i < m_size; i++ ) {
		m_data[i] = other.m_data[i];
	}
}

bool VecN::operator==( VecN other ) {
	if ( m_size != other.m_size ) {
		return false;
	}

	for ( unsigned int i = 0; i < m_size; i++ ) {
		if ( m_data[i] != other.m_data[i] ) {
			return false;
		}
	}
	return true;
}

VecN VecN::operator+( VecN other ) const {
	assert( m_size == other.m_size );
	VecN returnVec( m_size );
	for ( unsigned int i = 0; i < m_size; i++ ) {
		returnVec[i] = m_data[i] + other.m_data[i];
	}
	return returnVec;
}

void VecN::operator+=( const VecN other ) {
	assert( m_size == other.m_size );
	for ( unsigned int i = 0; i < m_size; i++ ) {
		m_data[i] += other.m_data[i];
	}
}

VecN VecN::operator+( float n ) const {
	VecN returnVec( m_size );
	for ( unsigned int i = 0; i < m_size; i++ ) {
		returnVec[i] = m_data[i] + n;
	}
	return returnVec;
}

void VecN::operator+=( const float n ) {
	for ( unsigned int i = 0; i < m_size; i++ ) {
		m_data[i] += n;
	}
}

VecN VecN::operator-( VecN other ) const {
	assert( m_size == other.m_size );
	VecN returnVec( m_size );
	for ( unsigned int i = 0; i < m_size; i++ ) {
		returnVec[i] = m_data[i] - other.m_data[i];
	}
	return returnVec;
}

void VecN::operator-=( const VecN other ) {
	assert( m_size == other.m_size );
	for ( unsigned int i = 0; i < m_size; i++ ) {
		m_data[i] -= other.m_data[i];
	}
}

VecN VecN::operator-( float n ) const {
	VecN returnVec( m_size );
	for ( unsigned int i = 0; i < m_size; i++ ) {
		returnVec[i] = m_data[i] - n;
	}
	return returnVec;
}

void VecN::operator-=( const float n ) {
	for ( unsigned int i = 0; i < m_size; i++ ) {
		m_data[i] -= n;
	}
}

VecN VecN::operator*( const float n ) const {
	VecN returnVec( m_size );
	for ( unsigned int i = 0; i < m_size; i++ ) {
		returnVec[i] = m_data[i] * n;
	}
	return returnVec;
}

void VecN::operator*=( const float n ) {
	for ( unsigned int i = 0; i < m_size; i++ ) {
		m_data[i] *= n;
	}
}

VecN VecN::operator*( MatN mat ) const {
	assert( m_size == mat.RowCount() );
	const unsigned int colCount = mat.ColumnCount();
	VecN returnVec( colCount );
	for ( unsigned int i = 0; i < colCount; i++ ) {
		for ( unsigned int j = 0; j < m_size; j++ ) {
			float p1 = m_data[j];
			float p2 = mat[ i + j * colCount ];
			returnVec[i] += m_data[j] * mat[ i + j * colCount ];
		}
	}
	return returnVec;
}

void VecN::operator*=( const MatN mat ) {
	assert( m_size == mat.RowCount() );
	const unsigned int colCount = mat.ColumnCount();
	float * newData = new float[colCount];
	for ( unsigned int i = 0; i < colCount; i++ ) {
		for ( unsigned int j = 0; j < m_size; j++ ) {
			float p1 = m_data[j];
			float p2 = mat[ i + j * colCount ];
			newData[i] += m_data[j] * mat[ i + j * colCount ];
		}
	}
	delete[] m_data;
	m_data = newData;
}

VecN VecN::operator/( float n ) const {
	VecN returnVec( m_size );
	for ( unsigned int i = 0; i < m_size; i++ ) {
		returnVec[i] = m_data[i] / n;
	}
	return returnVec;
}

void VecN::operator/=( const float n ) {
	for ( unsigned int i = 0; i < m_size; i++ ) {
		m_data[i] /= n;
	}
}

float VecN::Dot( const VecN & other ) const {
	assert( m_size == other.m_size );
	float val = 0.0f;
	for ( unsigned int i = 0; i < m_size; i++ ) {
		val += m_data[i] * other[i];
	}
	return val;
}

float VecN::Length() const {
	return sqrt( Dot( *this ) );
}

VecN VecN::Normal() const {
	const float len = Length();
	VecN returnVec( m_size );
	for ( unsigned int i = 0; i < m_size; i++ ) {
		returnVec[i] = m_data[i] / len;
	}
	return returnVec;

}

void VecN::Normalize() {
	const float len = Length();
	VecN returnVec( m_size );
	for ( unsigned int i = 0; i < m_size; i++ ) {
		m_data[i] /= len;
	}
}

VecN VecN::Proj( const VecN & other ) {
	const float scalar = Dot( other.Normal() );
	return VecN( m_size, m_data ) * scalar;
}

float VecN::Angle( const VecN & other ) {
	const float d = Normal().Dot( other.Normal() );
	return acosf( d );
}

Vec2 VecN::as_Vec2() const {
	Vec2 returnVec;
	if( m_size < 2 ) {
		const int dif = 2 - m_size;
		for ( int i = 2; i > dif; i-- ) {
			returnVec[ i - 1 ] = 0.0f;
		}
		for ( unsigned int i = 0; i < m_size; i++ ) {
			returnVec[i] = m_data[i];
		}
	} else {
		returnVec[0] = m_data[0];
		returnVec[1] = m_data[1];
	}
	return returnVec;
}

Vec3 VecN::as_Vec3() const {
	Vec3 returnVec;
	if( m_size < 3 ) {
		const int dif = 3 - m_size;
		for ( int i = 3; i > dif; i-- ) {
			returnVec[ i - 1 ] = 0.0f;
		}
		for ( unsigned int i = 0; i < m_size; i++ ) {
			returnVec[i] = m_data[i];
		}
	} else {
		returnVec[0] = m_data[0];
		returnVec[1] = m_data[1];
		returnVec[2] = m_data[2];
	}
	return returnVec;
}

Vec4 VecN::as_Vec4() const {
	Vec4 returnVec;
	if( m_size < 4 ) {
		const int dif = 4 - m_size;
		for ( int i = 4; i > dif; i-- ) {
			returnVec[ i - 1 ] = 0.0f;
		}
		for ( int i = 0; i < m_size; i++ ) {
			returnVec[i] = m_data[i];
		}
	} else {
		returnVec[0] = m_data[0];
		returnVec[1] = m_data[1];
		returnVec[2] = m_data[2];
		returnVec[3] = m_data[3];
	}
	return returnVec;
}


/*
================================
Vec2
================================
*/
Vec2::Vec2( float n ) {
	m_size = 2;
	m_data = new float[2];
	m_data[0] = n;
	m_data[1] = n;
}

Vec2::Vec2( float x, float y ) {
	m_size = 2;
	m_data = new float[2];
	m_data[0] = x;
	m_data[1] = y;
}

Vec2::Vec2( const float * data ) {
	m_size = 2;
	m_data = new float[2];
	m_data[0] = data[0];
	m_data[1] = data[1];
}

Vec2::Vec2( const Vec2 &vec ) {
	m_size = vec.m_size;
	m_data = new float[m_size];
	m_data[0] = vec.m_data[0];
	m_data[1] = vec.m_data[1];
}

void Vec2::operator=( Vec2 other ) {
	m_data[0] = other.m_data[0];
	m_data[1] = other.m_data[1];
}

bool Vec2::operator==( Vec2 other ) {
	if ( m_data[0] != other.m_data[0] ||
		m_data[1] != other.m_data[1] ) {
		return false;
	}
	return true;
}

Vec2 Vec2::operator+( Vec2 other ) const {
	Vec2 returnVec;
	returnVec[0] = m_data[0] + other.m_data[0];
	returnVec[1] = m_data[1] + other.m_data[1];
	return returnVec;
}

void Vec2::operator+=( const Vec2 other ) {
	m_data[0] += other.m_data[0];
	m_data[1] += other.m_data[1];
}

Vec2 Vec2::operator+( float n ) const {
	Vec2 returnVec;
	returnVec[0] = m_data[0] + n;
	returnVec[1] = m_data[1] + n;
	return returnVec;
}

void Vec2::operator+=( const float n ) {
	m_data[0] += n;
	m_data[1] += n;
}

Vec2 Vec2::operator-( Vec2 other ) const {
	Vec2 returnVec;
	returnVec[0] = m_data[0] - other.m_data[0];
	returnVec[1] = m_data[1] - other.m_data[1];
	return returnVec;
}

void Vec2::operator-=( const Vec2 other ) {
	m_data[0] -= other.m_data[0];
	m_data[1] -= other.m_data[1];
}

Vec2 Vec2::operator-( float n ) const {
	Vec2 returnVec;
	returnVec[0] = m_data[0] - n;
	returnVec[1] = m_data[1] - n;
	return returnVec;
}

void Vec2::operator-=( const float n ) {
	m_data[0] -= n;
	m_data[1] -= n;
}

Vec2 Vec2::operator*( const float n ) const {
	Vec2 returnVec;
	returnVec[0] = m_data[0] * n;
	returnVec[1] = m_data[1] * n;
	return returnVec;
}

void Vec2::operator*=( const float n ) {
	m_data[0] *= n;
	m_data[1] *= n;
}

Vec2 Vec2::operator*( MatN mat ) const {
	assert( m_size == mat.RowCount() );
	Vec2 returnVec;
	returnVec[0] = m_data[0] * mat[0] + m_data[1] * mat[2];
	returnVec[1] = m_data[0] * mat[1] + m_data[1] * mat[3];
	return returnVec;
}


Vec2 Vec2::operator*( Mat2 mat ) const {
	Vec2 returnVec;
	returnVec[0] = m_data[0] * mat[0] + m_data[1] * mat[2];
	returnVec[1] = m_data[0] * mat[1] + m_data[1] * mat[3];
	return returnVec;
}

void Vec2::operator*=( const MatN mat ) {
	assert( m_size == mat.RowCount() );
	const float x = m_data[0] * mat[0] + m_data[1] * mat[2];
	const float y = m_data[0] * mat[1] + m_data[1] * mat[3];
	m_data[0] = x;
	m_data[1] = y;
}

void Vec2::operator*=( const Mat2 mat ) {
	const float x = m_data[0] * mat[0] + m_data[1] * mat[2];
	const float y = m_data[0] * mat[1] + m_data[1] * mat[3];
	m_data[0] = x;
	m_data[1] = y;
}

Vec2 Vec2::operator/( float n ) const {
	Vec2 returnVec;
	returnVec[0] = m_data[0] / n;
	returnVec[1] = m_data[1] / n;
	return returnVec;
}

void Vec2::operator/=( const float n ) {
	m_data[0] /= n;
	m_data[1] /= n;
}

float Vec2::Dot( const Vec2 & other ) const {
	return m_data[0] * other[0] + m_data[1] * other[1];
}

float Vec2::Length() const {
	return sqrt( Dot( *this ) );
}

Vec2 Vec2::Normal() const {
	const float len = Length();
	Vec2 returnVec;
	returnVec[0] = m_data[0] / len;
	returnVec[1] = m_data[1] / len;
	return returnVec;

}

void Vec2::Normalize() {
	const float len = Length();
	Vec2 returnVec;
	m_data[0] /= len;
	m_data[1] /= len;
}

Vec2 Vec2::Proj( const Vec2 & other ) {
	const float scalar = Dot( other.Normal() );
	return Vec2( m_data ) * scalar;
}

float Vec2::Angle( const Vec2 & other ) {
	const float d = Normal().Dot( other.Normal() );
	return acosf( d );
}

VecN Vec2::as_VecN() const {
	float * data = new float[2];
	data[0] = m_data[0];
	data[1] = m_data[1];
	return VecN( 2, data );
}

Vec3 Vec2::as_Vec3() const {
	return Vec3( m_data[0], m_data[1], 0.0f );
}

Vec4 Vec2::as_Vec4() const {
	return Vec4( m_data[0], m_data[1], 0.0f, 0.0f );
}


/*
================================
Vec3
================================
*/
Vec3::Vec3( float n ) {
	m_size = 3;
	m_data = new float[3];
	m_data[0] = n;
	m_data[1] = n;
	m_data[2] = n;
}

Vec3::Vec3( float x, float y, float z ) {
	m_size = 3;
	m_data = new float[3];
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
}

Vec3::Vec3( const float * data ) {
	m_size = 3;
	m_data = new float[3];
	m_data[0] = data[0];
	m_data[1] = data[1];
	m_data[2] = data[2];
}

Vec3::Vec3( const Vec3 &vec ) {
	m_size = vec.m_size;
	m_data = new float[m_size];
	m_data[0] = vec.m_data[0];
	m_data[1] = vec.m_data[1];
	m_data[2] = vec.m_data[2];
}

void Vec3::operator=( Vec3 other ) {
	m_data[0] = other.m_data[0];
	m_data[1] = other.m_data[1];
	m_data[2] = other.m_data[2];
}

bool Vec3::operator==( Vec3 other ) {
	if ( m_data[0] != other.m_data[0] ||
		m_data[1] != other.m_data[1] ||
		m_data[2] != other.m_data[2] ) {
		return false;
	}
	return true;
}

Vec3 Vec3::operator+( Vec3 other ) const {
	Vec3 returnVec;
	returnVec[0] = m_data[0] + other.m_data[0];
	returnVec[1] = m_data[1] + other.m_data[1];
	returnVec[2] = m_data[2] + other.m_data[2];
	return returnVec;
}

void Vec3::operator+=( const Vec3 other ) {
	m_data[0] += other.m_data[0];
	m_data[1] += other.m_data[1];
	m_data[2] += other.m_data[2];
}

Vec3 Vec3::operator+( float n ) const {
	Vec3 returnVec;
	returnVec[0] = m_data[0] + n;
	returnVec[1] = m_data[1] + n;
	returnVec[2] = m_data[2] + n;
	return returnVec;
}

void Vec3::operator+=( const float n ) {
	m_data[0] += n;
	m_data[1] += n;
	m_data[2] += n;
}

Vec3 Vec3::operator-( Vec3 other ) const {
	Vec3 returnVec;
	returnVec[0] = m_data[0] - other.m_data[0];
	returnVec[1] = m_data[1] - other.m_data[1];
	returnVec[2] = m_data[2] - other.m_data[2];
	return returnVec;
}

void Vec3::operator-=( const Vec3 other ) {
	m_data[0] -= other.m_data[0];
	m_data[1] -= other.m_data[1];
	m_data[2] -= other.m_data[2];
}

Vec3 Vec3::operator-( float n ) const {
	Vec3 returnVec;
	returnVec[0] = m_data[0] - n;
	returnVec[1] = m_data[1] - n;
	returnVec[2] = m_data[2] - n;
	return returnVec;
}

void Vec3::operator-=( const float n ) {
	m_data[0] -= n;
	m_data[1] -= n;
	m_data[2] -= n;
}

Vec3 Vec3::operator*( const float n ) const {
	Vec3 returnVec;
	returnVec[0] = m_data[0] * n;
	returnVec[1] = m_data[1] * n;
	returnVec[2] = m_data[2] * n;
	return returnVec;
}

void Vec3::operator*=( const float n ) {
	m_data[0] *= n;
	m_data[1] *= n;
	m_data[2] *= n;
}

Vec3 Vec3::operator*( MatN mat ) const {
	assert( m_size == mat.RowCount() );
	Vec3 returnVec;
	returnVec[0] = m_data[0] * mat[0] + m_data[1] * mat[3] + m_data[2] * mat[6];
	returnVec[1] = m_data[0] * mat[1] + m_data[1] * mat[4] + m_data[2] * mat[7];
	returnVec[2] = m_data[0] * mat[2] + m_data[1] * mat[5] + m_data[2] * mat[8];
	return returnVec;
}

Vec3 Vec3::operator*( Mat3 mat ) const {
	Vec3 returnVec;
	returnVec[0] = m_data[0] * mat[0] + m_data[1] * mat[3] + m_data[2] * mat[6];
	returnVec[1] = m_data[0] * mat[1] + m_data[1] * mat[4] + m_data[2] * mat[7];
	returnVec[2] = m_data[0] * mat[2] + m_data[1] * mat[5] + m_data[2] * mat[8];
	return returnVec;
}

Vec3 Vec3::operator*( Mat4 mat ) const {
	Vec3 returnVec;
	returnVec[0] = m_data[0] * mat[0] + m_data[1] * mat[4] + m_data[2] * mat[8];
	returnVec[1] = m_data[0] * mat[1] + m_data[1] * mat[5] + m_data[2] * mat[9];
	returnVec[2] = m_data[0] * mat[2] + m_data[1] * mat[6] + m_data[2] * mat[10];
	return returnVec;
}

void Vec3::operator*=( const MatN mat ) {
	assert( m_size == mat.RowCount() );
	const float x = m_data[0] * mat[0] + m_data[1] * mat[3] + m_data[2] * mat[6];
	const float y = m_data[0] * mat[1] + m_data[1] * mat[4] + m_data[2] * mat[7];
	const float z = m_data[0] * mat[2] + m_data[1] * mat[5] + m_data[2] * mat[8];
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
}

void Vec3::operator*=( const Mat3 mat ) {
	const float x = m_data[0] * mat[0] + m_data[1] * mat[3] + m_data[2] * mat[6];
	const float y = m_data[0] * mat[1] + m_data[1] * mat[4] + m_data[2] * mat[7];
	const float z = m_data[0] * mat[2] + m_data[1] * mat[5] + m_data[2] * mat[8];
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
}

void Vec3::operator*=( const Mat4 mat ) {
	const float x = m_data[0] * mat[0] + m_data[1] * mat[4] + m_data[2] * mat[8];
	const float y = m_data[0] * mat[1] + m_data[1] * mat[5] + m_data[2] * mat[9];
	const float z = m_data[0] * mat[2] + m_data[1] * mat[6] + m_data[2] * mat[10];
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
}

Vec3 Vec3::operator/( float n ) const {
	Vec3 returnVec;
	returnVec[0] = m_data[0] / n;
	returnVec[1] = m_data[1] / n;
	returnVec[2] = m_data[2] / n;
	return returnVec;
}

void Vec3::operator/=( const float n ) {
	m_data[0] /= n;
	m_data[1] /= n;
	m_data[2] /= n;
}

float Vec3::Dot( const Vec3 & other ) const {
	return m_data[0] * other[0] + m_data[1] * other[1] + m_data[2] * other[2];
}

float Vec3::Length() const {
	return sqrt( Dot( *this ) );
}

Vec3 Vec3::Cross( const Vec3 & other ) const {
	Vec3 returnVec;
	returnVec[0] = m_data[1] * other[2] - m_data[2] * other[1];
	returnVec[1] = m_data[2] * other[0] - m_data[0] * other[2];
	returnVec[2] = m_data[0] * other[1] - m_data[1] * other[0];
	return returnVec;
}

Vec3 Vec3::Normal() const {
	const float len = Length();
	Vec3 returnVec;
	returnVec[0] = m_data[0] / len;
	returnVec[1] = m_data[1] / len;
	returnVec[2] = m_data[2] / len;
	return returnVec;

}

void Vec3::Normalize() {
	const float len = Length();
	Vec3 returnVec;
	m_data[0] /= len;
	m_data[1] /= len;
	m_data[2] /= len;
}

Vec3 Vec3::Proj( const Vec3 & other ) {
	const float scalar = Dot( other.Normal() );
	return Vec3( m_data ) * scalar;
}

float Vec3::Angle( const Vec3 & other ) {
	const float d = Normal().Dot( other.Normal() );
	return acosf( d );
}

VecN Vec3::as_VecN() const {
	float * data = new float[3];
	data[0] = m_data[0];
	data[1] = m_data[1];
	data[2] = m_data[2];
	return VecN( 3, data );
}

Vec2 Vec3::as_Vec2() const {
	return Vec2( m_data[0], m_data[1] );
}

Vec4 Vec3::as_Vec4() const {
	return Vec4( m_data[0], m_data[1], 0.0f, 0.0f );
}


/*
================================
Vec4
================================
*/
Vec4::Vec4( float n ) {
	m_size = 4;
	m_data = new float[4];
	m_data[0] = n;
	m_data[1] = n;
	m_data[2] = n;
	m_data[3] = n;
}

Vec4::Vec4( float x, float y, float z, float w ) {
	m_size = 4;
	m_data = new float[4];
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
	m_data[3] = w;
}

Vec4::Vec4( const float * data ) {
	m_size = 4;
	m_data = new float[4];
	m_data[0] = data[0];
	m_data[1] = data[1];
	m_data[2] = data[2];
	m_data[3] = data[3];
}

Vec4::Vec4( const Vec4 &vec ) {
	m_size = vec.m_size;
	m_data = new float[m_size];
	m_data[0] = vec.m_data[0];
	m_data[1] = vec.m_data[1];
	m_data[2] = vec.m_data[2];
	m_data[3] = vec.m_data[3];
}

void Vec4::operator=( Vec4 other ) {
	m_data[0] = other.m_data[0];
	m_data[1] = other.m_data[1];
	m_data[2] = other.m_data[2];
	m_data[3] = other.m_data[3];
}

bool Vec4::operator==( Vec4 other ) {
	if ( m_data[0] != other.m_data[0] ||
		m_data[1] != other.m_data[1] ||
		m_data[2] != other.m_data[2] ||
		m_data[3] != other.m_data[3] ) {
		return false;
	}
	return true;
}

Vec4 Vec4::operator+( Vec4 other ) const {
	Vec4 returnVec;
	returnVec[0] = m_data[0] + other.m_data[0];
	returnVec[1] = m_data[1] + other.m_data[1];
	returnVec[2] = m_data[2] + other.m_data[2];
	returnVec[3] = m_data[3] + other.m_data[3];
	return returnVec;
}

void Vec4::operator+=( const Vec4 other ) {
	m_data[0] += other.m_data[0];
	m_data[1] += other.m_data[1];
	m_data[2] += other.m_data[2];
	m_data[3] += other.m_data[3];
}

Vec4 Vec4::operator+( float n ) const {
	Vec4 returnVec;
	returnVec[0] = m_data[0] + n;
	returnVec[1] = m_data[1] + n;
	returnVec[2] = m_data[2] + n;
	returnVec[3] = m_data[3] + n;
	return returnVec;
}

void Vec4::operator+=( const float n ) {
	m_data[0] += n;
	m_data[1] += n;
	m_data[2] += n;
	m_data[3] += n;
}

Vec4 Vec4::operator-( Vec4 other ) const {
	Vec4 returnVec;
	returnVec[0] = m_data[0] - other.m_data[0];
	returnVec[1] = m_data[1] - other.m_data[1];
	returnVec[2] = m_data[2] - other.m_data[2];
	returnVec[3] = m_data[3] - other.m_data[3];
	return returnVec;
}

void Vec4::operator-=( const Vec4 other ) {
	m_data[0] -= other.m_data[0];
	m_data[1] -= other.m_data[1];
	m_data[2] -= other.m_data[2];
	m_data[3] -= other.m_data[3];
}

Vec4 Vec4::operator-( float n ) const {
	Vec4 returnVec;
	returnVec[0] = m_data[0] - n;
	returnVec[1] = m_data[1] - n;
	returnVec[2] = m_data[2] - n;
	returnVec[3] = m_data[3] - n;
	return returnVec;
}

void Vec4::operator-=( const float n ) {
	m_data[0] -= n;
	m_data[1] -= n;
	m_data[2] -= n;
	m_data[3] -= n;
}

Vec4 Vec4::operator*( const float n ) const {
	Vec4 returnVec;
	returnVec[0] = m_data[0] * n;
	returnVec[1] = m_data[1] * n;
	returnVec[2] = m_data[2] * n;
	returnVec[3] = m_data[3] * n;
	return returnVec;
}

void Vec4::operator*=( const float n ) {
	m_data[0] *= n;
	m_data[1] *= n;
	m_data[2] *= n;
	m_data[3] *= n;
}

Vec4 Vec4::operator*( MatN mat ) const {
	assert( m_size == mat.RowCount() );
	Vec4 returnVec;
	returnVec[0] = m_data[0] * mat[0] + m_data[1] * mat[4] + m_data[2] * mat[8] + m_data[3] * mat[12];
	returnVec[1] = m_data[0] * mat[1] + m_data[1] * mat[5] + m_data[2] * mat[9] + m_data[3] * mat[13];
	returnVec[2] = m_data[0] * mat[2] + m_data[1] * mat[6] + m_data[2] * mat[10] + m_data[3] * mat[14];
	returnVec[3] = m_data[0] * mat[3] + m_data[1] * mat[7] + m_data[2] * mat[11] + m_data[3] * mat[15];
	return returnVec;
}

Vec4 Vec4::operator*( Mat4 mat ) const {
	Vec4 returnVec;
	returnVec[0] = m_data[0] * mat[0] + m_data[1] * mat[4] + m_data[2] * mat[8] + m_data[3] * mat[12];
	returnVec[1] = m_data[0] * mat[1] + m_data[1] * mat[5] + m_data[2] * mat[9] + m_data[3] * mat[13];
	returnVec[2] = m_data[0] * mat[2] + m_data[1] * mat[6] + m_data[2] * mat[10] + m_data[3] * mat[14];
	returnVec[3] = m_data[0] * mat[3] + m_data[1] * mat[7] + m_data[2] * mat[11] + m_data[3] * mat[15];
	return returnVec;
}

void Vec4::operator*=( const MatN mat ) {
	assert( m_size == mat.RowCount() );
	const float x = m_data[0] * mat[0] + m_data[1] * mat[4] + m_data[2] * mat[8] + m_data[3] * mat[12];
	const float y = m_data[0] * mat[1] + m_data[1] * mat[5] + m_data[2] * mat[9] + m_data[3] * mat[13];
	const float z = m_data[0] * mat[2] + m_data[1] * mat[6] + m_data[2] * mat[10] + m_data[3] * mat[14];
	const float w = m_data[0] * mat[3] + m_data[1] * mat[7] + m_data[2] * mat[11] + m_data[3] * mat[15];
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
	m_data[3] = w;
}

void Vec4::operator*=( const Mat4 mat ) {
	const float x = m_data[0] * mat[0] + m_data[1] * mat[4] + m_data[2] * mat[8] + m_data[3] * mat[12];
	const float y = m_data[0] * mat[1] + m_data[1] * mat[5] + m_data[2] * mat[9] + m_data[3] * mat[13];
	const float z = m_data[0] * mat[2] + m_data[1] * mat[6] + m_data[2] * mat[10] + m_data[3] * mat[14];
	const float w = m_data[0] * mat[3] + m_data[1] * mat[7] + m_data[2] * mat[11] + m_data[3] * mat[15];
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
	m_data[3] = w;
}

Vec4 Vec4::operator/( float n ) const {
	Vec4 returnVec;
	returnVec[0] = m_data[0] / n;
	returnVec[1] = m_data[1] / n;
	returnVec[2] = m_data[2] / n;
	return returnVec;
}

void Vec4::operator/=( const float n ) {
	m_data[0] /= n;
	m_data[1] /= n;
	m_data[2] /= n;
}

float Vec4::Dot( const Vec4 & other ) const {
	return m_data[0] * other[0] + m_data[1] * other[1] + m_data[2] * other[2] + m_data[3] * other[3];
}

float Vec4::Length() const {
	return sqrtf( Dot( *this ) );
}

Vec4 Vec4::Normal() const {
	const float len = Length();
	Vec4 returnVec;
	returnVec[0] = m_data[0] / len;
	returnVec[1] = m_data[1] / len;
	returnVec[2] = m_data[2] / len;
	returnVec[3] = m_data[3] / len;
	return returnVec;

}

void Vec4::Normalize() {
	const float len = Length();
	Vec4 returnVec;
	m_data[0] /= len;
	m_data[1] /= len;
	m_data[2] /= len;
	m_data[3] /= len;
}

Vec4 Vec4::Proj( const Vec4 & other ) {
	const float scalar = Dot( other.Normal() );
	return Vec4( m_data ) * scalar;
}

float Vec4::Angle( const Vec4 & other ) {
	const float d = Normal().Dot( other.Normal() );
	return acosf( d );
}

VecN Vec4::as_VecN() const {
	float * data = new float[4];
	data[0] = m_data[0];
	data[1] = m_data[1];
	data[2] = m_data[2];
	data[3] = m_data[3];
	return VecN( 4, data );
}

Vec2 Vec4::as_Vec2() const {
	return Vec2( m_data[0], m_data[1] );
}

Vec3 Vec4::as_Vec3() const {
	return Vec3( m_data[0], m_data[1], m_data[2] );
}