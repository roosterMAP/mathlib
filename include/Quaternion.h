#pragma once
#ifndef __QUAT_H_INCLUDE__
#define __QUAT_H_INCLUDE__

#include "Vector.h"
#include "Matrix.h"

class Quat {
	public:
		Quat();
		Quat( const float n );
		Quat( const float x, const float y, const float z, const float w );
		Quat( Vec3 pos );
		Quat( Vec3 axis, float radians );
		Quat( const Quat &q );
		~Quat() {}

		void operator=( Quat other );
		bool  operator==( Quat other );
		bool  operator!=( Quat other ) { return !( *this == other ); }

		float operator[]( int index ) const;
		float & operator[]( int index );

		Quat operator*( const float n ) const;
		Quat operator*( const Quat other ) const;
		void operator*=( const float n );
		void operator*=( const Quat other );

		Quat operator/( const float n ) const;
		void operator/=( const float n );

		void Invert();
		Quat Inverse() const;
		void Rotate( Vec3 & v ) const;
		float Length() const;
		float LengthSquared() const;
		void Normalize();
		Quat Normal() const;
		Quat Conjugate() const;
		void Slerp( const Quat & qa, const Quat & qb, double t, Quat & qm );

		Mat3 as_Mat3() const;
		Mat4 as_Mat4( const Vec3 &vPosition ) const;

		static void Rotate( const Vec3 & center, const Vec3 & axis, const float theta, Vec3 & point );
		Vec3 RotatePoint( const Vec3 &rhs ) const;

	private:
		float m_data[4];
};

/*
//TestCode
Quat q1;
Quat q2( 5.0f );
Quat q3( 1.0f, 2.0f, 3.0f, 0.0f );
Vec3 dir{ 0.0f, 1.0f, 0.0f };
Quat q4( dir, 90.0f );
Quat q5( q4 );

if ( q4 == q5 ) {
	q5 = q3;
}
if ( q5 != q3 ) {
	return LXe_FAILED;
}

if ( q1[3] == 1.0f ) {
	q1[3] = 5.0f;
}
Quat q6 = q1 * 5.0f;
q1 *= 5.0f;
q6 = q6 / 5.0f;
q5 /= 5.0f;
q6 = q5.Inverse();
const float len = q6.Length();
if ( len != 1.0f ) {
	q6.Normalize();
}
q5 = q5.Normal();
Vec3 arm{ 1.0f, 0.0f, 0.0f };
q4.Rotate( arm ); //arm becomes [0.0, 0.0, -1.0]
*/

#endif