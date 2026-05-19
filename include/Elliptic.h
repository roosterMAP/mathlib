#pragma once

#ifndef __ELLIPTIC_H_INCLUDE__
#define __ELLIPTIC_H_INCLUDE__

#include <Math.h>

#define CARLSON_MAX_ITER 8

class Complex
{
public:
	Complex()
	{
		real = 0.0f;
		imag = 0.0f;
	}

	Complex( float r, float i )
	{
		real = r;
		imag = i;
	}

	float real;
	float imag;

	static Complex Zero()
	{
		return Complex( 0.0f, 0.0f );
	};

	static Complex One()
	{
		return Complex( 1.0f, 0.0f );
	};



	Complex operator+( const Complex& c ) const
	{
		return Complex( real + c.real, imag + c.imag );
	}

	void operator+=( const Complex& c )
	{
		real += c.real;
		imag += c.imag;
	}

	Complex operator-() const
	{
		return Complex(-real, -imag);
	}

	Complex operator-( const Complex& c ) const
	{
		return Complex( real - c.real, imag - c.imag );
	}

	Complex operator*( const float &s ) const
	{
		return Complex( real * s, imag * s );
	}

	Complex operator*( const Complex& c ) const
	{
		Complex r;
		r.real = real * c.real - imag * c.imag;
		r.imag = real * c.imag + imag * c.real;
		return r;
	}

	void operator*=( const Complex& c )
	{
		float r = real * c.real - imag * c.imag;
		float i = real * c.imag + imag * c.real;
		real = r;
		imag = i;
	}

	Complex operator/( const float& s ) const
	{
		return Complex( real / s, imag / s );
	}

	Complex operator/( const Complex& c ) const
	{
		Complex r( real, imag );
		r = r * c.Conjugate();
		Complex d = c * c.Conjugate();
		if ( d.real == 0.0f )
			d.real += 1e6f;
		r.real /= d.real;
		r.imag /= d.real;
		return r;
	}

	void operator/=( const float& s )
	{
		real /= s;
		imag /= s;
	}

	Complex Conjugate() const
	{
		return { real, -imag };
	}

	Complex Sqrt() const
	{
		float r = sqrtf( MagSqrd() );
		float theta = atan2f( imag, real );
		float isign = imag < 0.0f ? -1.0f : 1.0f;
		return Complex(
			sqrtf( 0.5f * (r + real) ),
			isign * sqrtf( 0.5f * (r - real) )
		);
	}

	Complex Sin() const
	{
		return Complex(
			sinf( real ) * cosh( imag ),
			cosf( real ) * sinh( imag )
		);
	}

	Complex Tanh() const
	{
		// For a complex number z = x + iy:
		// tanh(x + iy) = (sinh(2x) + i sin(2y)) / (cosh(2x) + cos(2y))

		float x2 = real * 2.0f;
		float y2 = imag * 2.0f;

		// Check for large Real part to avoid cosh/sinh overflow
		if (x2 > 20.0f) return Complex( 1.0f, 0.0f );
		if (x2 < -20.0f) return Complex( -1.0f, 0.0f );

		float den = coshf( x2 ) + cosf( y2 );

		// Guard against division by zero at poles
		if (fabsf( den ) < 1e-12f) {
			return Complex( 1e12f, 1e12f );
		}

		return Complex( sinhf( x2 ) / den, sinf( y2 ) / den );
	}

	bool Compare( const Complex& c ) const
	{
		return (*this - c ).MagSqrd() < 1e-6f;
	}

	float MagSqrd() const
	{
		return real * real + imag * imag;
	}

	float Mag() const
	{
		return sqrtf( MagSqrd() );
	}

	float isReal() const
	{
		return fabsf( imag ) < 1e-12;
	}

	static Complex FromSqrt( const float v )
	{
		if ( v > 0.0f )
			return Complex( sqrtf( v ), 0.0f );
		return Complex( 0.0f, sqrtf( -v ) );
	}
};


// Returns the roots of the given polynomial
void QuadraticRoots( float a, float b, float c, Complex* roots );
void CubicRoots( float a, float b, float c, float d, Complex* roots );
float DepressedCubicRoots( float p, float q, Complex* roots );
float QuarticRoots( float a, float b, float c, float d, float e, Complex* roots );
void SortRoots( Complex* pRoots, const unsigned int nRootCount, const bool bAscending );
bool ValidateRootFindingFunctions();

class Weierstrass
{
public:
	enum class DegenerateCase {
		TRIGONOMETRIC,
		HYPERBOLIC,
		NOT_DEGENERATE
	};

	Weierstrass() {};
	~Weierstrass() {};

	Weierstrass( float g2, float g3 )
	{
		m_g2 = g2;
		m_g3 = g3;

		float p = g2 / -4.0;
		float q = g3 / -4.0;
		DepressedCubicRoots( p, q, m_roots ); //update to complex
		SortRoots( m_roots, 3, false );
		TestDegenerateCases();

		m_omega1 = RealHalfPeriod();
	};

	void Initialize( float g2, float g3 )
	{
		m_g2 = g2;
		m_g3 = g3;

		float p = m_g2 / -4.0f;
		float q = m_g3 / -4.0f;
		DepressedCubicRoots( p, q, m_roots ); //update to complex
		SortRoots( m_roots, 3, false );
		TestDegenerateCases();

		m_omega1 = RealHalfPeriod();
	}

	const Complex &RootByIndex( unsigned int i ) const
	{
		return m_roots[i];
	}

	void TestDegenerateCases()
	{
		Complex e1 = m_roots[0];
		Complex e2 = m_roots[1];
		Complex e3 = m_roots[2];

		// degenerate cases
		const float epsilon = 1e-12f;
		bool e12 = (e1 - e2).MagSqrd() < epsilon;
		bool e23 = (e2 - e3).MagSqrd() < epsilon;

		if (e12)
			m_degencase = DegenerateCase::HYPERBOLIC;
		else if (e23)
			m_degencase = DegenerateCase::TRIGONOMETRIC;
		else
			m_degencase = DegenerateCase::NOT_DEGENERATE;
	}

	Complex Evaluate( const float z ) const;
	Complex EvaluateInv( const float z ) const;
	float RealHalfPeriod() const;
	float ImagHalfPeriod() const;
	float GetRealHalfPeriod() const { return m_omega1; }

	float m_g2;
	float m_g3;
	float m_omega1; //real half-period
	Complex m_roots[3];
	DegenerateCase m_degencase;
};


// Carlson symmetric forms
Complex Carlson_RC( const Complex& x, const Complex& y );
Complex Carlson_RD( const Complex& x, const Complex& y, const Complex& z );
Complex Carlson_RF( const Complex& x, const Complex& y, const Complex& z );
Complex Carlson_RJ( const Complex& x, const Complex& y, const Complex& z, const Complex& p );

float Carlson_RC( const float& x, const float& y );
float Carlson_RF( const float& x, const float& y, const float& z );
float Carlson_RD( const float& x, const float& y, const float& z );
float Carlson_RJ( const float& x, const float& y, const float& z, const float& p );

bool ValidateCarlsonFunctions();

class CarlsonForm
{
	 enum class FormType
	 {
		CUBIC,
		ONE_QUADRATIC_FACTOR,
		TWO_QUADRATIC_FACTORS
	 };

public:
	CarlsonForm() {};
	CarlsonForm( float e1, float e2, float e3, float e4, float y, float x );
	CarlsonForm( float e1, const Complex &e2, float e4, float y, float x );
	CarlsonForm( const Complex &e1, const Complex& e3, float e4, bool useB5ZeroCase, float y, float x );


	void Initialize( float e1, float e2, float e3, float e4, float y, float x );
	void Initialize( float e1, const Complex& e2, float e4, float y, float x );
	void Initialize( const Complex& e1, const Complex& e3, float e4, bool useB5ZeroCase, float y, float x );


	float I_1c() const;
	float I_3c() const;
	float K_2c() const;
	float H() const { return m_H; }
	float RF() const;
	float SIGMA() const;

	float a( int i ) const { return m_a[i-1]; }
	float c2( int i, int j ) const { return m_c2[i-1][j-1]; }
	float d( int i, int j ) const { return m_d[i-1][j-1]; }
	float r( int i, int j ) const { return m_d[i-1][j-1]; } //same as d because b_i=1
	float g() const { return m_g1; }
	float f() const { return m_f1; }
	float bet( int i ) const { return m_beti5[i-1]; }
	float del( int i, int j ) const { return m_del[i-1][j-1]; }

private:
	void W2();
	void Q2();
	void P2();
	void M2();
	void L2();

	float K_2c_quadratic() const;
	float K_2c_cubic() const;

	float m_a[4]; //cubic roots
	float m_f1, m_g1; //quadratic factor
	float m_f2, m_g2; //quadratic factor
	float m_y, m_x; //lower and upper limits of integration
	FormType m_case;

	float m_X1, m_X2, m_X3;
	float m_Y1, m_Y2, m_Y3;
	float m_U1_2, m_U2_2, m_U3_2;

	float m_X4_2, m_Y4_2;

	float m_W2, m_Q2, m_P2, m_M2, m_L2_n, m_L2_p;
	float m_U; //only used for casis with one quadratic factor

	float m_d[4][4];
	float m_c2[4][4];
	float m_e1, m_n1;
	float m_e2, m_n2;

	float m_H;
	float m_beti5[2];
	float m_del[2][2];
};

bool ValidateCarlsonSymmetricForms();


// Other elliptic functions
Complex Jacobi_Sn( const Complex& u, const Complex& k );
Complex Jacobi_InvSn( const Complex& u, const Complex& k );

#endif //__ELLIPTIC_H_INCLUDE__