#include "..\include\Quaternion.h"
#include "Math.h"
#include <assert.h>

/*
================================
Quat::Quat
================================
*/
Quat::Quat() {
	//identity quaternion
	m_data[0] = 0.0f;
	m_data[1] = 0.0f;
	m_data[2] = 0.0f;
	m_data[3] = 1.0f;
}

/*
================================
Quat::Quat
================================
*/
Quat::Quat( const float n ) {
	m_data[0] = 0.0f;
	m_data[1] = 0.0f;
	m_data[2] = 0.0f;
	m_data[3] = n;
}

/*
================================
Quat::Quat
================================
*/
Quat::Quat( const float x, const float y, const float z, const float w ) {
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
	m_data[3] = w;
}

/*
================================
Quat::Quat
================================
*/
Quat::Quat( Vec3 pos ) {
	m_data[0] = pos[0];
	m_data[1] = pos[1];
	m_data[2] = pos[2];
	m_data[3] = 0.0f;
}

/*
================================
Quat::Quat
================================
*/
Quat::Quat( Vec3 axis, float radians ) {
	const float halfAngleRadians = 0.5f * radians;
	const float s = sinf( halfAngleRadians );
	axis.Normalize();
	m_data[0] = axis[0] * s;
	m_data[1] = axis[1] * s;
	m_data[2] = axis[2] * s;
	m_data[3] = cosf( halfAngleRadians );
}

/*
================================
Quat::Quat
================================
*/
Quat::Quat( const Quat &q ) {
	m_data[0] = q[0];
	m_data[1] = q[1];
	m_data[2] = q[2];
	m_data[3] = q[3];
}

/*
================================
Quat::operator=
================================
*/
void Quat::operator=( Quat other ) {
	m_data[0] = other[0];
	m_data[1] = other[1];
	m_data[2] = other[2];
	m_data[3] = other[3];
}

/*
================================
Quat::operator==
================================
*/
bool Quat::operator==( Quat other ) {
	if ( m_data[0] == other[0] && m_data[1] == other[1] && m_data[2] == other[2] && m_data[3] == other[3] ) {
		return true;
	}
	return false;
}

/*
================================
Quat::operator[]
================================
*/
float Quat::operator[]( int index ) const {
	assert( index >=0 && index < 4 );
	return m_data[index];
}

/*
================================
Quat::operator[]
================================
*/
float & Quat::operator[]( int index ) {
	assert( index >=0 && index < 4 );
	return m_data[index];
}

/*
================================
Quat::operator*
================================
*/
Quat Quat::operator*( const float n ) const {
	return Quat( m_data[0] * n, m_data[1] * n, m_data[2] * n, m_data[3] * n );
}

/*
================================
Quat::operator*
================================
*/
Quat Quat::operator*( const Quat rhs ) const {
	Quat temp;
	temp[0] = ( m_data[0] * rhs[3] ) + ( m_data[3] * rhs[0] ) + ( m_data[1] * rhs[2] ) - ( m_data[2] * rhs[1] ); // x
	temp[1] = ( m_data[1] * rhs[3] ) + ( m_data[3] * rhs[1] ) + ( m_data[2] * rhs[0] ) - ( m_data[0] * rhs[2] ); // y
	temp[2] = ( m_data[2] * rhs[3] ) + ( m_data[3] * rhs[2] ) + ( m_data[0] * rhs[1] ) - ( m_data[1] * rhs[0] ); // z
	temp[3] = ( m_data[3] * rhs[3] ) - ( m_data[0] * rhs[0] ) - ( m_data[1] * rhs[1] ) - ( m_data[2] * rhs[2] ); // w
	return temp;
}

/*
================================
Quat::operator*=
================================
*/
void Quat::operator*=( const float n ) {
	m_data[0] *= n;
	m_data[1] *= n;
	m_data[2] *= n;
	m_data[3] *= n;
}

/*
================================
Quat::operator*=
================================
*/
void Quat::operator*=( const Quat other ) {
	const float n0 = other[3] * m_data[3] - other[0] * m_data[0] - other[1] * m_data[1] - other[2] * m_data[2];	//w
	const float n1 = other[3] * m_data[0] + other[0] * m_data[3] - other[1] * m_data[2] + other[2] * m_data[1];	//x
	const float n2 = other[3] * m_data[1] + other[0] * m_data[2] + other[1] * m_data[3] - other[2] * m_data[0];	//y
	const float n3 = other[3] * m_data[2] - other[0] * m_data[1] + other[1] * m_data[0] + other[2] * m_data[3];	//z

	m_data[3] = n0;
	m_data[0] = n1;
	m_data[1] = n2;
	m_data[2] = n3;
}

/*
================================
Quat::operator/
================================
*/
Quat Quat::operator/( const float n ) const {
	return Quat( m_data[0] / n, m_data[1] / n, m_data[2] / n, m_data[3] / n );
}

/*
================================
Quat::operator/=
================================
*/
void Quat::operator/=( const float n ) {
	m_data[0] /= n;
	m_data[1] /= n;
	m_data[2] /= n;
	m_data[3] /= n;
}

/*
================================
Quat::Inverse
================================
*/
/*
Quat Quat::Inverse() const {
	Quat q( m_data[0] * -1.0f, m_data[1] * -1.0f, m_data[2] * -1.0f, m_data[3] * -1.0f );
	const float sum = q[0] + q[1] + q[2] + q[3];
	for ( unsigned int i = 0; i < 4; i++ ) {
		q[i] /= sum;
	}
	return q;
}
*/

void Quat::Invert() {
	*this *= 1.0f / LengthSquared();
	m_data[0] *= -1.0f;
	m_data[1] *= -1.0f;
	m_data[2] *= -1.0f;
}
Quat Quat::Inverse() const {
	Quat val( *this );
	val.Invert();
	return val;
}

/*
================================
Quat::Rotate
	-Rotate vector(not point) by quaternion.
================================
*/
void Quat::Rotate( Vec3 & v ) const {
	Quat q( v );
	q = *this * q * this->Conjugate();

	const float length = v.Length();
	v.Normalize();
	
	v[0] = q[0];
	v[1] = q[1];
	v[2] = q[2];

	v *= length;
}

/*
================================
Quat::Length
================================
*/
float Quat::Length() const {
	return sqrt( LengthSquared() );
}


/*
================================
Quat::LengthSquared
================================
*/
float Quat::LengthSquared() const {
	return m_data[0] * m_data[0] + m_data[1] * m_data[1] + m_data[2] * m_data[2] + m_data[3] * m_data[3];
}

/*
================================
Quat::Normalize
================================
*/
void Quat::Normalize() {
	const float len = Length();
	if ( len == 0.0 ) {
		m_data[0] = 0.0f;
		m_data[1] = 0.0f;
		m_data[2] = 0.0f;
		m_data[3] = 1.0f;
	} else {
		m_data[0] /= len;
		m_data[1] /= len;
		m_data[2] /= len;
		m_data[3] /= len;
	}
}

/*
================================
Quat::Normal
================================
*/
Quat Quat::Normal() const {
	Quat q = *this;

	const float len = Length();
	if ( len == 0.0 ) {
		q[0] = 0.0f;
		q[1] = 0.0f;
		q[2] = 0.0f;
		q[3] = 1.0f;
	} else {
		q[0] /= len;
		q[1] /= len;
		q[2] /= len;
		q[3] /= len;
	}

	return q;
}

/*
================================
Quat::Conjugate
================================
*/
Quat Quat::Conjugate() const {
	Quat q = *this;
	q[0] *= -1.0f;
	q[1] *= -1.0f;
	q[2] *= -1.0f;
	return q;
}

/*
================================
Quat::Slerp
================================
*/
void Quat::Slerp( const Quat & qa, const Quat & qb, double t, Quat & qm ) {
	//angle between qa and qb
	float cosHalfTheta = qa[0] * qb[0] + qa[1] * qb[1] + qa[2] * qb[2] + qa[3] * qb[3];

	// if qa=qb or qa=-qb then theta = 0 and we can return qa
	if ( fabs( cosHalfTheta ) >= 1.0f ){		
		qm[0] = qa[0];
		qm[1] = qa[1];
		qm[2] = qa[2];
		qm[3] = qa[3];
	}

	float halfTheta = acosf( cosHalfTheta );
	float sinHalfTheta = sqrt( 1.0f - cosHalfTheta * cosHalfTheta );

	// if theta = 180 degrees then result is not fully defined we could rotate around any axis normal to qa or qb
	if ( fabs( sinHalfTheta ) < 0.001f ) {		
		qm[0] = ( qa[0] * 0.5f + qb[0] * 0.5f );
		qm[1] = ( qa[1] * 0.5f + qb[1] * 0.5f );
		qm[2] = ( qa[2] * 0.5f + qb[2] * 0.5f );
		qm[3] = ( qa[3] * 0.5f + qb[3] * 0.5f );
	}

	float ratioA = sinf( ( 1.0f - t ) * halfTheta ) / sinHalfTheta;
	float ratioB = sinf( t * halfTheta ) / sinHalfTheta;

	qm[0] = ( qa[0] * ratioA + qb[0] * ratioB );
	qm[1] = ( qa[1] * ratioA + qb[1] * ratioB );
	qm[2] = ( qa[2] * ratioA + qb[2] * ratioB );
	qm[3] = ( qa[3] * ratioA + qb[3] * ratioB );
}

/*
================================
Quat::Rotate
	-Rotate point by theta radians about axis on center
================================
*/
void Quat::Rotate( const Vec3 & center, const Vec3 & axis, const float theta, Vec3 & point ) {
	Quat rotQuat( axis.Normal(), theta );
	Vec3 v = point - center;
	rotQuat.Rotate( v );
	point += center;
}

//GREG
Vec3 Quat::RotatePoint( const Vec3& rhs ) const {
	Quat vector( rhs[0], rhs[1], rhs[2], 0.0f);
	Quat final = *this * vector * Inverse() ;
	return Vec3( final[0], final[1], final[2]);
}


Mat3 Quat::as_Mat3() const
{
	float x = m_data[0];
	float y = m_data[1];
	float z = m_data[2];
	float w = m_data[3];

	float data[] = {
		1 - 2 * ( y * y + z * z ), 2 * ( x * y - w * z ), 2 * ( x * z + w * y ),
		2 * ( x * y + w * z ), 1 - 2 * ( x * x + z * z ), 2 * ( y * z - w * x ),
		2 * ( x * z - w * y ), 2 * ( y * z + w * x ), 1 - 2 * ( x * x + y * y )
	};

	// Calculate the matrix elements using the quaternion components
	return Mat3( data );
}


Mat4 Quat::as_Mat4( const Vec3 &vPosition ) const
{
	Mat4 mat = as_Mat3().as_Mat4();
	mat.Translate( vPosition );
	return mat;
}