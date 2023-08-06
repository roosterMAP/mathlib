#pragma once

#ifndef __MATRIX_H_INCLUDE__
#define __MATRIX_H_INCLUDE__

#define EPSILON 0.00000001

class VecN;
class Vec2;
class Vec3;
class Vec4;
class Mat2;
class Mat3;
class Mat4;

/*
================================
MatN -row major matrix of arbitrary dimension
================================
*/
class MatN {
	public:		
		MatN() { m_row = 0; m_col = 0; }
		MatN( const unsigned int size );
		MatN( const unsigned int size, const float * data );
		MatN( const unsigned int row, const unsigned int col );
		MatN( const unsigned int row, const unsigned int col, const float * data );
		~MatN() { delete[] m_data; }

		MatN( const MatN &mat );

		void operator=( MatN other );
		bool operator==( MatN other );
		bool operator!=( MatN other ) { return !( *this == other ); }

		float operator[]( int i ) const { return m_data[i]; }
		float & operator[]( int i ) { return m_data[i]; }

		VecN GetRowVec( const unsigned int row ) const;
		void SetRowVec( const unsigned int row, const VecN * vec );

		VecN GetColVec( const unsigned int col ) const;
		void SetColVec( const unsigned int col, const VecN * vec );

		float GetComponent( unsigned int m, unsigned int n ) const;
		void SetComponent( unsigned int m, unsigned int n, float val );

		MatN operator+( MatN other ) const;
		void operator+=( MatN other );

		MatN operator-( MatN other ) const;
		void operator-=( MatN other );

		MatN operator*( float scalar ) const;
		void operator*=( float scalar );
		VecN operator*( VecN vec ) const;
		void operator*=( VecN vec );
		MatN operator*( MatN other ) const;
		void operator*=( MatN other );

		MatN operator/( float scalar ) const;
		void operator/=( float scalar );

		unsigned int ColumnCount() const { return m_col; }
		unsigned int RowCount() const { return m_row; }

		float Determinant();
		bool Inverse( MatN * inv );
		MatN Transpose() const;
		void Transposed();

		Mat2 as_Mat2() const;
		Mat3 as_Mat3() const;
		Mat4 as_Mat4() const;

		const float * as_ptr() const { return m_data; }

	private:
		float * m_data = nullptr;
		unsigned int m_col, m_row;

	friend class Mat4;
};


/*
================================
Mat2
================================
*/
class Mat2 {
	public:
		Mat2();
		Mat2( float v );
		Mat2( const float * data );
		~Mat2() { delete[] m_data; }

		Mat2( const Mat2 &mat );

		void operator=( Mat2 other );
		bool operator==( Mat2 other );
		bool operator!=( Mat2 other ) { return !( *this == other ); }

		float operator[]( int i ) const { return m_data[i]; }
		float & operator[]( int i ) { return m_data[i]; }

		Vec2 GetRowVec( const unsigned int row ) const;
		void SetRowVec( const unsigned int row, const Vec2 * vec );

		Vec2 GetColVec( const unsigned int col ) const;
		void SetColVec( const unsigned int col, const Vec2 * vec );

		float GetComponent( unsigned int m, unsigned int n ) const;
		void SetComponent( unsigned int m, unsigned int n, float val );

		Mat2 operator+( Mat2 other ) const;
		void operator+=( Mat2 other );

		Mat2 operator-( Mat2 other ) const;
		void operator-=( Mat2 other );

		Mat2 operator*( float scalar ) const;
		void operator*=( float scalar );
		Vec2 operator*( Vec2 vec ) const;
		Mat2 operator*( Mat2 other ) const;
		void operator*=( Mat2 other );

		Mat2 operator/( float scalar ) const;
		void operator/=( float scalar );

		float Determinant();
		bool Inverse( Mat2 * inv );
		Mat2 Transpose() const;
		void Transposed();

		MatN as_MatN() const;
		Mat3 as_Mat3() const;
		Mat4 as_Mat4() const;

		const float * as_ptr() const { return m_data; }

	private:
		float * m_data = nullptr;
		unsigned int m_size;
};


/*
================================
Mat3
================================
*/
class Mat3 {
	public:
		Mat3();
		Mat3( float v );
		Mat3( const float * data );
		~Mat3() { delete[] m_data; }

		Mat3( const Mat3 &mat );

		void operator=( Mat3 other );
		bool operator==( Mat3 other );
		bool operator!=( Mat3 other ) { return !( *this == other ); }

		float operator[]( int i ) const { return m_data[i]; }
		float & operator[]( int i ) { return m_data[i]; }

		Vec3 GetRowVec( const unsigned int row ) const;
		void SetRowVec( const unsigned int row, const Vec3 * vec );

		Vec3 GetColVec( const unsigned int col ) const;
		void SetColVec( const unsigned int col, const Vec3 * vec );

		float GetComponent( unsigned int m, unsigned int n ) const;
		void SetComponent( unsigned int m, unsigned int n, float val );

		Mat3 operator+( Mat3 other ) const;
		void operator+=( Mat3 other );

		Mat3 operator-( Mat3 other ) const;
		void operator-=( Mat3 other );

		Mat3 operator*( float scalar ) const;
		void operator*=( float scalar );
		Vec3 operator*( Vec3 vec ) const;
		Mat3 operator*( Mat3 other ) const;
		void operator*=( Mat3 other );

		Mat3 operator/( float scalar ) const;
		void operator/=( float scalar );

		float Determinant();
		bool Inverse( Mat3 * inv );
		Mat3 Transpose() const;
		void Transposed();

		MatN as_MatN() const;
		Mat2 as_Mat2() const;
		Mat4 as_Mat4() const;

		const float * as_ptr() const { return m_data; }

	private:
		float * m_data = nullptr;
		unsigned int m_size;
};


/*
================================
Mat4
================================
*/
class Mat4 {
	public:
		Mat4();
		Mat4( float v );
		Mat4( const float * data );
		~Mat4() { delete[] m_data; }

		Mat4( const Mat4 &mat );

		void operator=( Mat4 other );
		bool operator==( Mat4 other );
		bool operator!=( Mat4 other ) { return !( *this == other ); }

		float operator[]( int i ) const { return m_data[i]; }
		float & operator[]( int i ) { return m_data[i]; }

		Vec4 GetRowVec( const unsigned int row ) const;
		void SetRowVec( const unsigned int row, const Vec4 * vec );

		Vec4 GetColVec( const unsigned int col ) const;
		void SetColVec( const unsigned int col, const Vec4 * vec );

		float GetComponent( unsigned int m, unsigned int n ) const;
		void SetComponent( unsigned int m, unsigned int n, float val );

		Mat4 operator+( Mat4 other ) const;
		void operator+=( Mat4 other );

		Mat4 operator-( Mat4 other ) const;
		void operator-=( Mat4 other );

		Mat4 operator*( float scalar ) const;
		void operator*=( float scalar );
		Vec4 operator*( Vec4 vec ) const;
		Mat4 operator*( Mat4 other ) const;
		void operator*=( Mat4 other );

		Mat4 operator/( float scalar ) const;
		void operator/=( float scalar );

		float Determinant();
		bool Inverse( Mat4 * inv );
		Mat4 Transpose() const;
		void Transposed();

		void LookAt( const Vec3 look, const Vec3 up, const Vec3 pos );
		void Perspective( const float verticalFOV, const float aspect, const float near, const float far );
		void Orthographic( const float left, const float right, const float bottom, const float top );
		void Orthographic( const float left, const float right, const float bottom, const float top, const float near, const float far );

		MatN as_MatN() const;
		Mat2 as_Mat2() const;
		Mat3 as_Mat3() const;

		const float * as_ptr() const { return m_data; }

	private:
		float * m_data = nullptr;
		unsigned int m_size;
};
#endif