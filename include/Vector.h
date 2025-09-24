#pragma once

#ifndef __VECTOR_H_INCLUDE__
#define __VECTOR_H_INCLUDE__

class MatN;
class Mat2;
class Mat3;
class Mat4;
class Vec2;
class Vec3;
class Vec4;

#define EPSILON 0.00000001
#define PI 3.14159265359
#define PI2 6.283185307179586

float to_degrees( const float rads );
float to_radians( const float deg );
Vec3 BarycentricCoordinates2D( const Vec2 &vA, const Vec2 &vB, const Vec2 &vC, const Vec2 &vP );
Vec3 BarycentricCoordinates3D( const Vec3 &vA, const Vec3 &vB, const Vec3 &vC, const Vec3 &vP );
Vec4 BarycentricCoordinates3DTetra( const Vec3 &vA, const Vec3 &vB, const Vec3 &vC, const Vec3 &vD, const Vec3 &vP );

/*
================================
VecN::VecN
================================
*/
class VecN {
	public:
		VecN() { m_size = 0; }
		VecN( unsigned int size );
		VecN( unsigned int size, float n );
		VecN( unsigned int size, const float * data );
		~VecN() { delete[] m_data; }

		VecN( const VecN &vec );

		VecN& operator=( const VecN& rhs );
		bool operator==( const VecN &other ) const;
		bool operator!=( const VecN &other ) const { return !( *this == other ); }

		float operator[]( const int index ) const { return m_data[index]; }
		float & operator[]( const int index ) { return m_data[index]; }

		VecN operator+( const VecN &other ) const;
		VecN operator+( float n ) const;
		void operator+=( const VecN &other );
		void operator+=( const float n );

		VecN operator-( const VecN &other ) const;
		VecN operator-( float n ) const;
		void operator-=( const VecN &other );
		void operator-=( const float n );

		VecN operator*( const float n ) const;
		VecN operator*( const MatN &mat ) const;
		void operator*=( const float n );
		void operator*=( const MatN &mat );

		VecN operator/( float n ) const;
		void operator/=( const float n );

		unsigned int size() const { return m_size; }
		float Dot( const VecN & other ) const;
		float Length() const;
		VecN Normal() const;
		void Normalize();
		VecN Proj( const VecN & other );
		float Angle( const VecN & other );
		void Zero();

		Vec2 as_Vec2() const;
		Vec3 as_Vec3() const;
		Vec4 as_Vec4() const;

		const float * as_ptr() const { return m_data; }

	private:		
		unsigned int m_size;
		float * m_data = nullptr;

	friend class MatN;
};

/*
================================
Vec2::Vec2
================================
*/
class Vec2 {
	public:
		Vec2() { m_data[0] = 0.0f; m_data[1] = 0.0f; };
		Vec2( float n );
		Vec2( float x, float y );
		Vec2( const float * data );
		~Vec2() {};

		Vec2( const Vec2 &vec );

		void operator=( Vec2 other );
		bool operator==( const Vec2 &other ) const;
		bool operator!=( const Vec2 &other ) const { return !( *this == other ); }

		float operator[]( int index ) const { return m_data[index]; }
		float & operator[]( int index ) { return m_data[index]; }

		Vec2 operator+( Vec2 other ) const;
		Vec2 operator+( float n ) const;
		void operator+=( const Vec2 other );
		void operator+=( const float n );

		Vec2 operator-( Vec2 other ) const;
		Vec2 operator-( float n ) const;
		void operator-=( const Vec2 other );
		void operator-=( const float n );

		Vec2 operator*( const float n ) const;
		Vec2 operator*( Mat2 mat ) const;
		Vec2 operator*( MatN mat ) const;
		void operator*=( const float n );
		void operator*=( const Mat2 mat );
		void operator*=( const MatN mat );

		Vec2 operator/( float n ) const;
		void operator/=( const float n );

		unsigned int size() const { return 2; }
		float Dot( const Vec2 & other ) const;
		float Length() const;
		float LengthSquared() const { return m_data[0] * m_data[0] + m_data[1] * m_data[1]; }
		Vec2 Normal() const;
		void Normalize();
		Vec2 Proj( const Vec2 & other );
		float Angle( const Vec2 & other );
		void Zero() { m_data[0] = 0.0f; m_data[1] = 0.0f; }

		VecN as_VecN() const;
		Vec3 as_Vec3() const;
		Vec4 as_Vec4() const;

		static Vec2 X() { return Vec2( 1.0f, 0.0f ); }
		static Vec2 Y() { return Vec2( 0.0f, 1.0f ); }

		const float * as_ptr() const { return m_data; }

	private:
		float m_data[2];
};


class Vec3 {
	public:
		Vec3() { m_data[0] = 0.0f; m_data[1] = 0.0f; m_data[2] = 0.0f; };
		Vec3( float n );
		Vec3( float x, float y, float z );
		Vec3( const float f, const Vec2& vec );
		Vec3( const float * data );
		~Vec3() {};

		Vec3( const Vec3 &vec );

		void operator=( Vec3 other );
		bool operator==( const Vec3 &other ) const;
		bool operator!=( const Vec3 &other ) const { return !( *this == other ); }

		float operator[]( int index ) const { return m_data[index]; }
		float & operator[]( int index ) { return m_data[index]; }

		Vec3 operator+( const Vec3& other ) const;
		Vec3 operator+( float n ) const;
		void operator+=( const Vec3& other );
		void operator+=( const float n );

		Vec3 operator-( const Vec3 &other ) const;
		Vec3 operator-( const float n ) const;
		void operator-=( const Vec3 &other );
		void operator-=( const float n );

		Vec3 operator*( const float n ) const;
		Vec3 operator*( const MatN &mat ) const;
		Vec3 operator*( const Mat3 &mat ) const;
		Vec3 operator*( const Mat4 &mat ) const;
		void operator*=( const float n );
		void operator*=( const MatN &mat );
		void operator*=( const Mat3 &mat );
		void operator*=( const Mat4 &mat );

		Vec3 operator/( float n ) const;
		void operator/=( const float n );

		unsigned int size() const { return 3; }
		float Dot( const Vec3 & other ) const;
		float Length() const;
		float LengthSquared() const;
		Vec3 Cross( const Vec3 & other ) const;
		Vec3 Normal() const;
		void Normalize();
		Vec3 Proj( const Vec3 & other );
		float Angle( const Vec3 & other );
		void Zero() { m_data[0] = 0.0f; m_data[1] = 0.0f; m_data[2] = 0.0f; }
		bool AlmostEqual( const Vec3& other, float epsilon = EPSILON ) const;
		Vec3 GetOrthogonal() const;
		bool IsValid() const;

		VecN as_VecN() const;
		Vec2 as_Vec2() const;
		Vec4 as_Vec4() const;
		Vec2 YZ() const { return Vec2( m_data[1], m_data[2] ); }

		static Vec3 X() { return Vec3( 1.0f, 0.0f, 0.0f ); }
		static Vec3 Y() { return Vec3( 0.0f, 1.0f, 0.0f ); }
		static Vec3 Z() { return Vec3( 0.0f, 0.0f, 1.0f ); }
	
		const float * as_ptr() const { return m_data; }

	private:
		float m_data[3];
};


class Vec4 {
	public:
		Vec4()  { m_data[0] = 0.0f; m_data[1] = 0.0f; m_data[2] = 0.0f; m_data[3] = 0.0f; };
		Vec4( float n );
		Vec4( float x, float y, float z, float w );
		Vec4( const Vec3 &v, float w );
		Vec4( const float * data );
		Vec4( const float f, const Vec3& vec );
		~Vec4() {};

		Vec4( const Vec4 &vec );

		void operator=( Vec4 other );
		bool operator==( const Vec4 &other ) const;
		bool operator!=( const Vec4 &other ) const { return !( *this == other ); }

		float operator[]( int index ) const { return m_data[index]; }
		float & operator[]( int index ) { return m_data[index]; }

		Vec4 operator+( Vec4 other ) const;
		Vec4 operator+( float n ) const;
		void operator+=( const Vec4 other );
		void operator+=( const float n );

		Vec4 operator-( Vec4 other ) const;
		Vec4 operator-( float n ) const;
		void operator-=( const Vec4 other );
		void operator-=( const float n );

		Vec4 operator*( const float n ) const;
		Vec4 operator*( MatN mat ) const;
		Vec4 operator*( Mat4 mat ) const;
		void operator*=( const float n );
		void operator*=( const MatN mat );
		void operator*=( const Mat4 mat );

		Vec4 operator/( float n ) const;
		void operator/=( const float n );

		unsigned int size() const { return 4; }
		float Dot( const Vec4 & other ) const;
		float Length() const;
		Vec4 Normal() const;
		void Normalize();
		Vec4 Proj( const Vec4 & other );
		float Angle( const Vec4 & other );
		void Zero() { m_data[0] = 0.0f; m_data[1] = 0.0f; m_data[2] = 0.0f; m_data[3] = 0.0f; }

		VecN as_VecN() const;
		Vec2 as_Vec2() const;
		Vec3 as_Vec3() const;
		Vec3 YZW() const { return Vec3( m_data[1], m_data[2], m_data[3]); }

		const float * as_ptr() const { return m_data; }

	private:
		float m_data[4];
};

#endif