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

float to_degrees( const float rads );
float to_radians( const float deg );

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

		void operator=( VecN other );
		bool operator==( VecN other );
		bool operator!=( VecN other ) { return !( *this == other ); }

		float operator[]( int index ) const { return m_data[index]; }
		float & operator[]( int index ) { return m_data[index]; }

		VecN operator+( VecN other ) const;
		VecN operator+( float n ) const;
		void operator+=( const VecN other );
		void operator+=( const float n );

		VecN operator-( VecN other ) const;
		VecN operator-( float n ) const;
		void operator-=( const VecN other );
		void operator-=( const float n );

		VecN operator*( const float n ) const;
		VecN operator*( MatN mat ) const;
		void operator*=( const float n );
		void operator*=( const MatN mat );

		VecN operator/( float n ) const;
		void operator/=( const float n );

		unsigned int size() const { return m_size; }
		float Dot( const VecN & other ) const;
		float Length() const;
		VecN Normal() const;
		void Normalize();
		VecN Proj( const VecN & other );
		float Angle( const VecN & other );

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
		Vec2() { m_size = 2; m_data = new float[2]; }
		Vec2( float n );
		Vec2( float x, float y );
		Vec2( const float * data );
		~Vec2() { delete[] m_data; }

		Vec2( const Vec2 &vec );

		void operator=( Vec2 other );
		bool operator==( Vec2 other );
		bool operator!=( Vec2 other ) { return !( *this == other ); }

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

		unsigned int size() const { return m_size; }
		float Dot( const Vec2 & other ) const;
		float Length() const;
		Vec2 Normal() const;
		void Normalize();
		Vec2 Proj( const Vec2 & other );
		float Angle( const Vec2 & other );

		VecN as_VecN() const;
		Vec3 as_Vec3() const;
		Vec4 as_Vec4() const;

		const float * as_ptr() const { return m_data; }

	private:		
		unsigned int m_size;
		float * m_data = nullptr;
};


class Vec3 {
	public:
		Vec3() { m_size = 3; m_data = new float[3]; }
		Vec3( float n );
		Vec3( float x, float y, float z );
		Vec3( const float * data );
		~Vec3() { delete[] m_data; }

		Vec3( const Vec3 &vec );

		void operator=( Vec3 other );
		bool operator==( Vec3 other );
		bool operator!=( Vec3 other ) { return !( *this == other ); }

		float operator[]( int index ) const { return m_data[index]; }
		float & operator[]( int index ) { return m_data[index]; }

		Vec3 operator+( Vec3 other ) const;
		Vec3 operator+( float n ) const;
		void operator+=( const Vec3 other );
		void operator+=( const float n );

		Vec3 operator-( Vec3 other ) const;
		Vec3 operator-( float n ) const;
		void operator-=( const Vec3 other );
		void operator-=( const float n );

		Vec3 operator*( const float n ) const;
		Vec3 operator*( MatN mat ) const;
		Vec3 operator*( Mat3 mat ) const;
		Vec3 operator*( Mat4 mat ) const;
		void operator*=( const float n );
		void operator*=( const MatN mat );
		void operator*=( const Mat3 mat );
		void operator*=( const Mat4 mat );

		Vec3 operator/( float n ) const;
		void operator/=( const float n );

		unsigned int size() const { return m_size; }
		float Dot( const Vec3 & other ) const;
		float Length() const;
		Vec3 Cross( const Vec3 & other ) const;
		Vec3 Normal() const;
		void Normalize();
		Vec3 Proj( const Vec3 & other );
		float Angle( const Vec3 & other );

		VecN as_VecN() const;
		Vec2 as_Vec2() const;
		Vec4 as_Vec4() const;

		const float * as_ptr() const { return m_data; }

	private:
		unsigned int m_size;
		float * m_data = nullptr;
};


class Vec4 {
	public:
		Vec4() { m_size = 2; }
		Vec4( float n );
		Vec4( float x, float y, float z, float w );
		Vec4( const float * data );
		~Vec4() { delete[] m_data; }

		Vec4( const Vec4 &vec );

		void operator=( Vec4 other );
		bool operator==( Vec4 other );
		bool operator!=( Vec4 other ) { return !( *this == other ); }

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

		unsigned int size() const { return m_size; }
		float Dot( const Vec4 & other ) const;
		float Length() const;
		Vec4 Normal() const;
		void Normalize();
		Vec4 Proj( const Vec4 & other );
		float Angle( const Vec4 & other );

		VecN as_VecN() const;
		Vec2 as_Vec2() const;
		Vec3 as_Vec3() const;

		const float * as_ptr() const { return m_data; }

	private:
		unsigned int m_size;
		float * m_data = nullptr;
};

#endif