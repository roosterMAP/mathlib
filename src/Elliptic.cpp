#include "..\include\Elliptic.h"
#include "..\include\Vector.h"
#include <math.h>
#include <vector>

Complex Weierstrass::Evaluate( const float z ) const
{
	Complex e1 = m_roots[0];
	Complex e2 = m_roots[1];
	Complex e3 = m_roots[2];

	//fold input z into the fundamental half-domain 0-omega1
	float new_z = fabsf( z );
	if ( new_z > m_omega1 )
		new_z = m_omega1 * 2.0f - new_z; //fold

	Complex cz( new_z, 0.0f );

	// Initial guess: leading term of the Laurent series
	Complex x = Complex::One() / (cz * cz);

	const int MAX_IT = 8;
	for (int i = 0; i < MAX_IT; i++)
	{
		Complex t1 = x - e1;
		Complex t2 = x - e2;
		Complex t3 = x - e3;

		Complex rf = Carlson_RF( t1, t2, t3 );

		// Residual
		Complex f = rf - cz;

		// 2. Analytical Derivative: df/dx = -0.5 / sqrt((x-e1)(x-e2)(x-e3))
		Complex df = Complex( -0.5f, 0.0f ) / (t1 * t2 * t3).Sqrt();

		Complex step = f / df;
		x = x - step;

		if (step.MagSqrd() < EPSILON )
			break;
	}

	return x;
}

Complex Weierstrass::EvaluateInv( const float z ) const
{
	Complex e1 = m_roots[0];
	Complex e2 = m_roots[1];
	Complex e3 = m_roots[2];

	// e2 = e3 → trigonometric degeneration
	if (m_degencase == DegenerateCase::HYPERBOLIC || m_degencase == DegenerateCase::TRIGONOMETRIC)
	{
		float a = e1.real;
		float b = e3.real;

		float s = sqrt( z - a );
		float denom = sqrt( a - b );

		return Complex( log( (s - denom) / (s + denom) ) / denom, 0.0f );
	}

	if ( m_roots[0].isReal() && m_roots[1].isReal() && m_roots[2].isReal() )
		return Complex( Carlson_RF( z - m_roots[0].real, z - m_roots[1].real, z - m_roots[2].real ), 0.0f );
	const Complex cz( z, 0.0f );
	return Carlson_RF( cz - m_roots[0], cz - m_roots[1], cz - m_roots[2] );
}

float Weierstrass::RealHalfPeriod() const
{
	Complex e1 = m_roots[0];
	Complex e2 = m_roots[1];
	Complex e3 = m_roots[2];

	if (m_degencase == DegenerateCase::HYPERBOLIC || m_degencase == DegenerateCase::TRIGONOMETRIC)
	{
		float x = (e1 - e3).real;
		return PI / sqrt( x );
	}

	return Carlson_RF( e1 - e3, e1 - e2, Complex::Zero() ).real;
}

float Weierstrass::ImagHalfPeriod() const
{
	Complex e1 = m_roots[0];
	Complex e2 = m_roots[1];
	Complex e3 = m_roots[2];

	if ( ( e2 - e3 ).MagSqrd() < EPSILON )
		return PI / sqrt( (e1 - e3).real );

	return Carlson_RF( e3 - e1, e3 - e2, Complex::Zero() ).imag;

}

void QuadraticRoots( float a, float b, float c, Complex* roots )
{
	float discriminant = b * b - 4.0f * a * c;
	if (fabsf( discriminant ) > EPSILON)
	{
		Complex sqrtD = Complex::FromSqrt( discriminant );
		Complex denom = Complex( 2.0f * a, 0.0f );

		roots[0] = (Complex( -b, 0.0f ) + sqrtD) / denom;
		roots[1] = (Complex( -b, 0.0f ) - sqrtD) / denom;
	}
	else
	{
		//double root
		roots[0] = Complex( -b / (2.0f * a), 0.0f );
		roots[1] = roots[0];
	}
}

void CubicRoots( float a, float b, float c, float d, Complex* roots )
{
	// Convert to depressed cubic t^3 + pt + q = 0
	float p = ( 3.0f * a * c - b * b ) / ( 3.0f * a * a );
	float q = ( 2.0f * b * b * b - 9.0f * a * b * c + 27.0f * a * a * d ) / ( 27.0f * a * a * a );
	float shift = b / (3.0f * a);

	DepressedCubicRoots( p, q, roots );

	for (int i = 0; i < 3; i++)
	{
		roots[i].real -= shift;
	}
}

float CardanosRoots( float p, float q, unsigned int k )
{
	float costheta = -q * 0.5f * sqrtf( 27.0f / -( p * p * p ) );
	float theta = acosf( fmaxf( -1.0f, fminf( 1.0f, costheta ) ) );
	return 2.0f * sqrtf( -p / 3.0f ) * cosf( (theta  + 2.0f * PI * (float)k ) / 3.0f );
}

float DepressedCubicRoots( float p, float q, Complex* roots )
{

	float qo2 = q * 0.5f;
	float po3 = p / 3.0f;

	float delta = qo2 * qo2 + po3 * po3 * po3; //descriminant
	const float epsilon = 1e-12f;
	if (delta > epsilon)
	{
		// 1 real root, two complex which are conjugate pairs
		float sqrt_disc = sqrtf( delta );
		float u = cbrtf( -qo2 + sqrt_disc );
		float v = cbrtf( -qo2 - sqrt_disc );

		roots[0] = { u + v, 0.0f }; //real
		float real_part = -(u + v) / 2.0f;
		float imag_part = (u - v) * sqrtf( 3.0f ) / 2.0f;
		roots[1] = { real_part, imag_part }; //complex
		roots[2] = { real_part, -imag_part }; //complex conjugate
	}
	else if (fabs( delta ) <= epsilon)
	{
		float u = cbrtf( -qo2 );

		if (fabs( u ) < epsilon)
		{
			// Triple root
			roots[0] = { 0.0f, 0.0f };
			roots[1] = { 0.0f, 0.0f };
			roots[2] = { 0.0f, 0.0f };
		} else
		{
			// One simple + one double root
			roots[0] = { 2.0f * u, 0.0f };
			roots[1] = { -u, 0.0f };
			roots[2] = { -u, 0.0f };
		}
	}
	else
	{
		// 3 real roots
		roots[0] = { CardanosRoots( p, q, 0 ), 0.0f };
		roots[1] = { CardanosRoots( p, q, 1 ), 0.0f };
		roots[2] = { CardanosRoots( p, q, 2 ), 0.0f };
	}

	return delta;
}

Complex HornersMethod( float* coeffs, unsigned int nCoeffs, Complex x )
{
	//evaluate polynomial
	Complex result = { coeffs[0], 0.0f};
	for (unsigned int i = 1; i < nCoeffs; ++i)
	{
		result = result * x;
		result.real += coeffs[i];
	}
	return result;
}

float QuarticRoots( float a, float b, float c, float d, float e, Complex* roots )
{
	// find roots using Durand-Kerner
	float coeffs[5] = {a, b, c, d, e};
	const unsigned int n = 4;

	if (fabsf( e ) < EPSILON)
	{
		roots[0] = Complex::Zero();
		CubicRoots( a, b, c, d, roots + 1 ); // solve cubic
		if ( fabsf( roots[1].imag ) < EPSILON && fabsf( roots[2].imag ) < EPSILON && fabsf( roots[3].imag ) < EPSILON)
			return 1.0f;
		return -1.0f;
	}

	// normalize the polynomial
	for (int i = 1; i < 5; i++)
		coeffs[i] /= a;
	coeffs[0] = 1.0f;

	//initial point distribution
	float radius = 1.0f; //circle radius
	float offset = 0.0f;
	for (unsigned int k = 1; k < n + 1; ++k)
	{
		radius = fmaxf( radius, powf( fabsf( coeffs[k] ), 1.0f / k ) );
	}
	radius += 1.0f;	
	offset = -b / ( n * a ); //circle center offset

	const float theta = 2.0f * PI / (float)n;
	const float theta_offset = theta * 0.5;
	Complex r = { offset + radius * cosf( theta * 0.0f + theta_offset ), radius * sinf( theta * 0.0f + theta_offset ) };
	Complex s = { offset + radius * cosf( theta * 1.0f + theta_offset ), radius * sinf( theta * 1.0f + theta_offset ) };
	Complex t = { offset + radius * cosf( theta * 2.0f + theta_offset ), radius * sinf( theta * 2.0f + theta_offset ) };
	Complex u = { offset + radius * cosf( theta * 3.0f + theta_offset ), radius * sinf( theta * 3.0f + theta_offset ) };

	unsigned int max_iter = 128;
	for (unsigned int i = 0; i < max_iter; ++i)
	{
		Complex r0 = r, s0 = s, t0 = t, u0 = u;
		Complex dr = HornersMethod( coeffs, 5, r0 ) / ((r0 - s0) * (r0 - t0) * (r0 - u0));
		Complex ds = HornersMethod( coeffs, 5, s0 ) / ((s0 - r0) * (s0 - t0) * (s0 - u0));
		Complex dt = HornersMethod( coeffs, 5, t0 ) / ((t0 - r0) * (t0 - s0) * (t0 - u0));
		Complex du = HornersMethod( coeffs, 5, u0 ) / ((u0 - r0) * (u0 - s0) * (u0 - t0));

		r = r0 - dr;
		s = s0 - ds;
		t = t0 - dt;
		u = u0 - du;

		// convergence check
		if ( dr.MagSqrd() < EPSILON && ds.MagSqrd() < EPSILON && dt.MagSqrd() < EPSILON && du.MagSqrd() < EPSILON )
			break;
	}

	roots[0] = r;
	roots[1] = s;
	roots[2] = t;
	roots[3] = u;

	if ( fabsf( r.imag ) < EPSILON && fabsf( s.imag ) < EPSILON && fabsf( t.imag ) < EPSILON && fabsf( u.imag ) < EPSILON )
		return 1.0f;
	return -1.0f;
}


Complex Carlson_RC( const Complex& x, const Complex& y )
{
	// Cauchy Principal Value for negative y
	if (y.real < 0.0f)
	{
		Complex prefactor = ( x / (x - y) ).Sqrt();
		return prefactor * Carlson_RC( x - y, y * -1.0f );
	}

	Complex A0 = (x + y * 2.0f) / 3.0f;
	Complex An = A0;
	Complex xn = x;
	Complex yn = y;
	unsigned int i;
	for (i = 0; i < CARLSON_MAX_ITER; ++i)
	{
		Complex lambda = xn.Sqrt() * yn.Sqrt() * 2.0f + yn;
		An = (An + lambda) * 0.25f;
		xn = (xn + lambda) * 0.25f;
		yn = (yn + lambda) * 0.25f;

		Complex dx = (An - xn) / An;
		Complex dy = (An - yn) / An;

		float err = fmaxf( dx.Mag(), dy.Mag() );
		if (err < 1e-6f) break;
	}

	Complex s = (yn - An) / An;

	Complex taylor = Complex::One() + s * s * 3.0f / 10.0f +
		s * s * s * 1.0f / 7.0f +
		s * s * s * s * 3.0f / 8.0f +
		s * s * s * s * s * 9.0f / 22.0f;

	return taylor / An.Sqrt();
}


float Carlson_RC( const float& x, const float& y )
{
	//Cauchy Principal Value for negative y
	if (y < 0.0f)
	{
		float prefactor = sqrtf( x / (x - y) );
		return prefactor * Carlson_RC( x - y, -y );
	}

	float nmz = 1.0f;
	float mt = fmaxf( x, y );
	if (mt > 1e6)
		nmz = mt;

	float xn = x / nmz;
	float yn = y / nmz;
	float A0 = (xn + 2.0f * yn) / 3.0f;
	float An = A0;
	unsigned int i;
	for (i = 0; i < CARLSON_MAX_ITER; ++i)
	{
		float lambda = 2.0f * sqrtf( xn ) * sqrtf( yn ) + yn;
		An = 0.25f * (An + lambda);
		xn = 0.25f * (xn + lambda);
		yn = 0.25f * (yn + lambda);

		float dx = (An - xn) / An;
		float dy = (An - yn) / An;

		float err = fmaxf( fabsf( dx ), fabsf( dy ) );
		if (err < 1e-5f) break;
	}

	float s = (yn - An) / An;

	float taylor = 1.0f + 3.0f / 10.0f * s * s +
					1.0f / 7.0f * s * s * s +
					3.0f / 8.0f * s * s * s * s +
					9.0f / 22.0f * s * s * s * s * s;

	return taylor / sqrtf( An ) / sqrtf( nmz );
}


Complex Carlson_RD( const Complex& x, const Complex& y, const Complex& z )
{
	Complex A0 = (x + y + z * 3.0f) / 5.0f;
	Complex An = A0;
	Complex xn = x;
	Complex yn = y;
	Complex zn = z;
	Complex sum( 0.0f, 0.0f );
	float fac = 1.0f;
	unsigned int i;
	for (i = 0; i < CARLSON_MAX_ITER; ++i)
	{
		Complex sx = xn.Sqrt();
		Complex sy = yn.Sqrt();
		Complex sz = zn.Sqrt();

		Complex lambda = sx * sy + sx * sz + sy * sz;

		sum += Complex( fac, 0.0f ) / ( sz * (zn + lambda));

		fac = 0.25f * fac;

		An = (An + lambda) * 0.25;
		xn = (xn + lambda) * 0.25;
		yn = (yn + lambda) * 0.25;
		zn = (zn + lambda) * 0.25;

		Complex dx = (An - xn) / An;
		Complex dy = (An - yn) / An;
		Complex dz = (An - zn) / An;

		float err = fmaxf( fmaxf( dx.Mag(), dy.Mag()), dz.Mag());
		if (err < 1e-5f) break;
	}

	Complex X = (xn - An) / An;
	Complex Y = (yn - An) / An;
	Complex Z = (X + Y) / -3.0f;

	Complex E2 = X * Y - Z * Z * 6.0f;
	Complex E3 = (X * Y * 3.0f - Z * Z * 8.0f) * Z;
	Complex E4 = (X * Y - Z * Z) * Z * Z * 3.0f;

	Complex taylor = Complex::One() - E2 * ( 3.0f / 14.0f ) +
		E3 * ( 1.0f / 6.0f ) +
		E2 * E2 * ( 9.0f / 88.0f ) -
		E4 * ( 3.0f / 22.0f );

	return sum * 3.0f + taylor / (An * An.Sqrt() ) * fac;
}


float Carlson_RD( const float& x, const float& y, const float& z )
{
	float nmz = 1.0f;
	float mt = fmaxf( x, fmaxf( y, z ) );
	if ( mt > 1e6 )
		nmz = mt;

	float xn = x / nmz;
	float yn = y / nmz;
	float zn = z / nmz;
	float A0 = (xn + yn + 3.0f * zn) / 5.0f;
	float An = A0;
	float sum = 0.0f;
	float fac = 1.0f;
	unsigned int i;
	for (i = 0; i < CARLSON_MAX_ITER; ++i)
	{
		float sx = sqrtf( xn );
		float sy = sqrtf( yn );
		float sz = sqrtf( zn );

		float lambda = sx * sy + sx * sz + sy * sz;

		sum += fac / (sqrtf( zn ) * (zn + lambda));

		fac = 0.25f * fac;

		An = 0.25f * (An + lambda);
		xn = 0.25f * (xn + lambda);
		yn = 0.25f * (yn + lambda);
		zn = 0.25f * (zn + lambda);		

		float dx = (An - xn) / An;
		float dy = (An - yn) / An;
		float dz = (An - zn) / An;

		float err = fmaxf( fmaxf( fabsf( dx ), fabsf( dy ) ), fabsf( dz ) );
		if (err < 1e-5f) break;
	}

	float X = (xn - An) / An;
	float Y = (yn - An) / An;
	float Z = -( X + Y ) / 3.0f;

	float E2 = X * Y - 6.0f * Z * Z;
	float E3 = ( X * Y * 3.0f - Z * Z * 8.0f ) * Z;
	float E4 = ( X * Y - Z * Z ) * Z * Z * 3.0f;

	float taylor = 1.0f - E2 * 3.0f / 14.0f +
					E3 * 1.0f / 6.0f +
					E2 * E2 * 9.0f / 88.0f -
					E4 * 3.0f / 22.0f;

	float result = 3.0f * sum + fac * taylor / (An * sqrtf( An ));
	return result / ( nmz * sqrtf( nmz ) );
}


// This must be Complex to handle the roots and the z-plane correctly
Complex Carlson_RF( const Complex& x, const Complex& y, const Complex& z )
{
	Complex xn = x;
	Complex yn = y;
	Complex zn = z;
	Complex lambda;

	Complex A0 = (xn + yn + zn) / 3.0f;
	Complex An = A0;

	int i = 0;
	for (i = 0; i < CARLSON_MAX_ITER; i++)
	{
		Complex sx = xn.Sqrt();
		Complex sy = yn.Sqrt();
		Complex sz = zn.Sqrt();

		lambda = sx * sy + sx * sz + sy * sz;

		An = (An + lambda) * 0.25f;
		xn = (xn + lambda) * 0.25f;
		yn = (yn + lambda) * 0.25f;
		zn = (zn + lambda) * 0.25f;

		// Check relative error against An's magnitude
		Complex dx = (An - xn) / An;
		Complex dy = (An - yn) / An;
		Complex dz = (An - zn) / An;
		float err = fmaxf( dx.MagSqrd(), fmaxf( dy.MagSqrd(), dz.MagSqrd() ) );
		if (err < 1e-5f) break;
	}

	Complex X = (An - xn) / An;
	Complex Y = (An - yn) / An;
	Complex Z = ( X + Y ) * -1.0f;


	Complex E2 = X * Y - Z * Z;
	Complex E3 = X * Y * Z;

	return ( Complex::One() - E2 / 10.0f + E3 / 14.0f + E2 * E2 / 24.0f - E2 * E3 * (3.0f / 44.0f)) / An.Sqrt();
}


float Carlson_RF( const float& x, const float& y, const float& z )
{
	float nmz = 1.0f;
	float mt = fmaxf( x, fmaxf( y, z ) );
	if (mt > 1e6)
		nmz = mt;

	float Xn = x / nmz;
	float Yn = y / nmz;
	float Zn = z / nmz;

	float lambda;
	float mu;
	float dx, dy, dz;
	unsigned int i = 0;
	float err = 1.0f;
	do
	{
		lambda = sqrtf(Xn * Yn) + sqrtf(Xn * Zn) + sqrtf(Yn * Zn);

		Xn = (Xn + lambda) * 0.25f;
		Yn = (Yn + lambda) * 0.25f;
		Zn = (Zn + lambda) * 0.25f;

		mu = (Xn + Yn + Zn) / 3.0f;

		dx = (mu - Xn) / mu;
		dy = (mu - Yn) / mu;
		dz = (mu - Zn) / mu;

		err = fmaxf( dx, fmaxf( dy, dz ) );

		if (i > CARLSON_MAX_ITER)
			break;
		i += 1;

	} while (sqrtf( err ) > 1e-5);

	float E2 = dx * dy + dy * dz + dz * dx;
	float E3 = dx * dy * dz;

	float taylor = ( 1.0f - E2 / 10.0f + E3 / 14.0f + E2 * E2 / 24.0f - E2 * E3 * 3.0f / 44.0f ) / sqrtf( mu );
	return taylor / sqrtf( nmz );
}


Complex Carlson_RJ( const Complex& x, const Complex& y, const Complex& z, const Complex& p )
{
	Complex xn = x;
	Complex yn = y;
	Complex zn = z;
	Complex pn = p;

	Complex A0 = (x + y + z + p * 2.0f) / 5.0f;
	Complex An = A0;
	Complex delta = (p - x) * (p - y) * (p - z);

	Complex sum( 0.0f, 0.0f );
	float fac = 1.0f;
	Complex dm( 0.0f, 0.0f );
	Complex em( 0.0f, 0.0f );
	Complex dx, dy, dz, dp;

	unsigned int i = 0;
	for (i = 0; i < CARLSON_MAX_ITER; ++i)
	{
		Complex sx = xn.Sqrt();
		Complex sy = yn.Sqrt();
		Complex sz = zn.Sqrt();
		Complex sp = pn.Sqrt();
		Complex lambda = sx * sy + sx * sz + sy * sz;

		dm = (sp + sx) * (sp + sy) * (sp + sz);
		float fac3 = fac * fac * fac;
		em = (delta * fac3) / (dm * dm);
		sum += Complex( fac, 0.0f ) / dm * Carlson_RC( Complex::One(), Complex::One() + em);
		fac *= 0.25;

		An = (An + lambda) * 0.25f;
		xn = (xn + lambda) * 0.25f;
		yn = (yn + lambda) * 0.25f;
		zn = (zn + lambda) * 0.25f;
		pn = (pn + lambda) * 0.25f;

		dx = (An - xn) / An;
		dy = (An - yn) / An;
		dz = (An - zn) / An;
		dp = (An - pn) / An;

		float err = fmaxf( fmaxf( fmaxf( dx.Mag(), dy.Mag()), dz.Mag()), dp.Mag());
		if (err < 1e-5f)
			break;
	}

	Complex X = dx;
	Complex Y = dy;
	Complex Z = dz;
	Complex P = dp;
	Complex XYZ = X * Y * Z;

	Complex E2 = X * Y + X * Z + Y * Z - P * P * 3.0f;
	Complex E3 = XYZ + E2 * P * 2.0f + P * P * P * 4.0f;
	Complex E4 = (XYZ * 2.0f + E2 * P + P * P * P * 3.0f) * P;

	Complex taylor = Complex::One()
		- E2 * (3.0f / 14.0f)
		+ E3 * (1.0f / 6.0f)
		+ (E2 * E2) * (9.0f / 88.0f)
		- E4 * (3.0f / 22.0f);

	return sum * 6.0f + (taylor * fac) / (An * An.Sqrt());
}


float Carlson_RJ( const float& x, const float& y, const float& z, const float& p )
{	
	if ( p < 0.0f )
	{
		// --- Cauchy Principal Value for p < 0 ---
		float xt = std::min( { x, y, z } );
		float zt = std::max( { x, y, z } );
		float yt = x + y + z - xt - zt;

		float a = 1.0f / (yt - p);
		float b = a * (zt - yt) * (yt - xt);
		float pt = yt + b;
		float rho = (xt * zt) / yt;
		float tau = (p * pt) / yt;

		float rcx;
		if (tau < 0.0f)
			rcx = sqrtf( rho / (rho - tau) ) * Carlson_RC( rho - tau, -tau );
		else
			rcx = Carlson_RC( rho, tau );
		float rfx = Carlson_RF( xt, yt, zt );
		float rjx = Carlson_RJ( xt, yt, zt, pt );

		return a * (b * rjx + 3.0f * (rcx - rfx));
	}

	float nmz = 1.0f;
	float mt = fmaxf( x, fmaxf( y, z ) );
	if (mt > 1e6)
		nmz = mt;

	float xn = x / nmz;
	float yn = y / nmz;
	float zn = z / nmz;
	float pn = p / nmz;

	float A0 = ( x + y + z + p * 2.0f ) / 5.0f;
	float An = A0;
	float delta = ( p - x ) * ( p - y ) * ( p - z );

	float sum = 0.0f;
	float fac = 1.0f;
	float dm = 0.0f;
	float em = 0.0f;
	float dx, dy, dz, dp;

	unsigned int i = 0;
	for (i = 0; i < CARLSON_MAX_ITER; ++i)
	{
		float sx = sqrtf( xn );
		float sy = sqrtf( yn );
		float sz = sqrtf( zn );
		float sp = sqrtf( pn );
		float lambda = sx * sy + sx * sz + sy * sz;

		dm = (sp + sx) * (sp + sy) * (sp + sz);
		float fac3 = fac * fac * fac;
		em = (fac3 * delta) / (dm * dm);
		sum += fac / dm * Carlson_RC( 1.0f, 1.0f + em );
		fac *= 0.25;		

		An = (An + lambda) * 0.25f;
		xn = (xn + lambda) * 0.25f;
		yn = (yn + lambda) * 0.25f;
		zn = (zn + lambda) * 0.25f;
		pn = (pn + lambda) * 0.25f;

		dx = (An - xn) / An;
		dy = (An - yn) / An;
		dz = (An - zn) / An;
		dp = (An - pn) / An;

		float err = fmaxf( fmaxf( fmaxf( fabsf( dx ), fabsf( dy ) ), fabsf( dz ) ), fabsf( dp ) );
		if (err < 1e-5f)
			break;
	}

	float X = dx;
	float Y = dy;
	float Z = dz;	
	float P = dp;
	float XYZ = X * Y * Z;	

	float E2 = X * Y + X * Z + Y * Z - P * P * 3.0f;
	float E3 = XYZ + E2 * P * 2.0f + P * P * P * 4.0f;
	float E4 = (XYZ * 2.0f + E2 * P + P * P * P * 3.0f) * P;

	float taylor = 1.0f
		- E2 * (3.0f / 14.0f)
		+ E3 * (1.0f / 6.0f)
		+ (E2 * E2) * (9.0f / 88.0f)
		- E4 * (3.0f / 22.0f);

	float result = sum * 6.0f + (taylor * fac) / (An * sqrtf( An));
	return result / (nmz * sqrtf( nmz ) );
}


Complex Jacobi_Sn( const Complex& u, const Complex& m ) //m = k * k
{
	float weight = fmaxf( 0.0f, fminf( m.real, 1.0f ) );
	Complex s = u.Sin() * (1.0f - weight) + u.Tanh() * weight;

	Complex o( 1.0f, 0.0f );
	for (int i = 0; i < 8; ++i)
	{
		Complex s2 = s * s;

		Complex x = o - s2;
		Complex y = o - m * s2;
		Complex z = o;

		Complex RF = Carlson_RF( x, y, z );

		Complex f = s * RF - u;

		Complex df = o / (x * y).Sqrt();

		s = s - f / df;

		if (f.MagSqrd() < EPSILON)
			break;
	}

	return s;
}


Complex Jacobi_InvSn( const Complex& u, const Complex& k )
{
	Complex x2 = u * u;
	Complex k2 = k * k;
	Complex one = Complex( 1.0f, 0.0f );
	return u * Carlson_RF( one - x2, one - k2 * x2, one );
}


void SortRoots( Complex* pRoots, const unsigned int nRootCount, bool bAscending )
{
	std::vector<Complex> realRoots;
	std::vector<Complex> complexRoots;

	for (int i = 0; i < nRootCount; ++i)
	{
		if (fabs( pRoots[i].imag ) < EPSILON)
			realRoots.emplace_back( pRoots[i].real, 0.0f );
		else
			complexRoots.push_back( pRoots[i] );
	}

	// Sort real roots
	const unsigned int nRealCount = realRoots.size();
	for (unsigned int i = 0; i < nRealCount; ++i)
	{
		for (unsigned int j = i + 1; j < nRealCount; ++j)
		{
			if ( (realRoots[j].real < realRoots[i].real ) == bAscending)
				std::swap( realRoots[i], realRoots[j] );
		}
	}

	// Ensure conjugate order
	if (complexRoots.size() == 2)
	{
		if (complexRoots[0].imag < 0)
			std::swap( complexRoots[0], complexRoots[1] );
	}

	// Reassemble
	int idx = 0;
	for (const Complex& r : realRoots)
		pRoots[idx++] = r;

	for (const Complex& c : complexRoots)
		pRoots[idx++] = c;
}


float test_weirstrass_zfunc( float z )
{
	return 1.0f / z * z;
}


bool ValidateRootFindingFunctions()
{
	bool success = true;

	// Test Cubic Roots
	Complex* roots = new Complex[3];

	CubicRoots( 1, 0, -3, 2, roots );
	success = roots[0].Compare( { -2.0f, 0.0f } ) &&
		roots[1].Compare( { 1.0f, 0.0f } );
	if (!success)
		return false;

	DepressedCubicRoots( -3, 2, roots );
	success = roots[0].Compare( { -2.0f, 0.0f } ) &&
		roots[1].Compare( { 1.0f, 0.0f } );
	if (!success)
		return false;

	CubicRoots( 1, -8, 9, 18, roots );
	success = roots[0].Compare( { 6.0f, 0.0f } ) &&
		roots[1].Compare( { -1.0f, 0.0f } ) &&
		roots[2].Compare( { 3.0f, 0.0f } );
	if (!success)
		return false;

	CubicRoots( 1, -1, 1, -1, roots );
	success = roots[0].Compare( { 1.0f, 0.0f } ) &&
		roots[1].Compare( { 0.0f, 1.0f } ) &&
		roots[2].Compare( { 0.0f, -1.0f } );
	if (!success)
		return false;

	CubicRoots( 1, 0, -3, 1, roots );
	float rr1 = 2.0f * cosf( 2.0f * PI / 9.0f );
	float rr2 = 2.0f * cosf( 4.0f * PI / 9.0f );
	float rr3 = 2.0f * cosf( 8.0f * PI / 9.0f );
	success = roots[0].Compare( { rr1, 0.0f } ) &&
		roots[1].Compare( { rr3, 0.0f } ) &&
		roots[2].Compare( { rr2, 0.0f } );
	if (!success)
		return false;

	delete[] roots;

	// Test Quartic Roots
	roots = new Complex[4];

	QuarticRoots( 1, -11, 41, -61, 30, roots );
	success = roots[0].Compare( { 5.0f, 0.0f } ) &&
		roots[1].Compare( { 1.0f, 0.0f } ) &&
		roots[2].Compare( { 2.0f, 0.0f } ) &&
		roots[3].Compare( { 3.0f, 0.0f } );
	if (!success)
		return false;

	/*
	* fails. DK struggles with repeated roots. Try Aberth method?
	QuarticRoots( 1, -2, -7, 8, 12, roots );
	success = roots[0].Compare( { 2.0f, 0.0f } ) &&
		roots[1].Compare( { 2.0f, 0.0f } ) &&
		roots[2].Compare( { -2.0f, 0.0f } ) &&
		roots[3].Compare( { 3.0f, 0.0f } );
	if (!success)
		return false;
	*/

	QuarticRoots( 1, 0, 0, 0, 1, roots );
	float a = 1.0f / sqrtf( 2.0f );
	success = roots[0].Compare( { a, a } ) &&
		roots[1].Compare( { -a, a } ) &&
		roots[2].Compare( { -a, -a } ) &&
		roots[3].Compare( { a, -a } );
	if (!success)
		return false;
}


bool ValidateCarlsonFunctions()
{
	bool ok = true;

	auto FloatEqual = []( float a, float b, float eps = 1e-5f ) -> bool
		{
			return fabsf( a - b ) <= eps;
		};

	// ============================================================
	// RF
	// ============================================================

	ok &= FloatEqual(
		Carlson_RF( 1.0f, 2.0f, 0.0f ),
		1.3110287771461f
	);

	ok &= Carlson_RF(
		Complex( 0.0f, 1.0f ),
		Complex( 0.0f, -1.0f ),
		Complex::Zero()
	).Compare(
		Complex( 1.8540746773014f, 0.0f )
	);

	ok &= Carlson_RF(
		Complex( -1.0f, 1.0f ),
		Complex( 0.0f, 1.0f ),
		Complex::Zero()
	).Compare(
		Complex( 0.79612586584234f, -1.2138566698365f )
	);

	ok &= FloatEqual(
		Carlson_RF( 2.0f, 3.0f, 4.0f ),
		0.58408284167715f
	);

	ok &= Carlson_RF(
		Complex( 0.0f, 1.0f ),
		Complex( 0.0f, -1.0f ),
		Complex( 2.0f, 0.0f )
	).Compare(
		Complex( 1.0441445654064f, 0.0f )
	);

	ok &= Carlson_RF(
		Complex( -1.0f, 1.0f ),
		Complex( 0.0f, 1.0f ),
		Complex( 1.0f, -1.0f )
	).Compare(
		Complex( 0.93912050218619f, -0.53296252018635f )
	);

	// ============================================================
	// RC
	// ============================================================

	ok &= FloatEqual(
		Carlson_RC( 0.0f, 0.25f ),
		3.1415926535898f
	);

	ok &= FloatEqual(
		Carlson_RC( 9.0f / 4.0f, 2.0f ),
		0.69314718055995f
	);

	ok &= Carlson_RC(
		Complex::Zero(),
		Complex( 0.0f, 1.0f )
	).Compare(
		Complex( 1.1107207345396f, -1.1107207345396f )
	);

	ok &= Carlson_RC(
		Complex( 0.0f, -1.0f ),
		Complex( 0.0f, 1.0f )
	).Compare(
		Complex( 1.2260849569072f, -0.34471136988768f )
	);

	ok &= FloatEqual(
		Carlson_RC( 0.25f, -2.0f ),
		0.23104906018665f
	);

	ok &= Carlson_RC(
		Complex( 0.0f, 1.0f ),
		Complex( -1.0f, 0.0f )
	).Compare(
		Complex( 0.77778596920447f, 0.19832484993429f )
	);

	// ============================================================
	// RD
	// ============================================================

	ok &= FloatEqual(
		Carlson_RD( 0.0f, 2.0f, 1.0f ),
		1.7972103521034f
	);

	ok &= FloatEqual(
		Carlson_RD( 2.0f, 3.0f, 4.0f ),
		0.16510527294261f
	);

	ok &= Carlson_RD(
		Complex( 0.0f, 1.0f ),
		Complex( 0.0f, -1.0f ),
		Complex( 2.0f, 0.0f )
	).Compare(
		Complex( 0.65933854154220f, 0.0f )
	);

	ok &= Carlson_RD(
		Complex::Zero(),
		Complex( 0.0f, 1.0f ),
		Complex( 0.0f, -1.0f )
	).Compare(
		Complex( 1.2708196271910f, 2.7811120159521f )
	);

	ok &= Carlson_RD(
		Complex::Zero(),
		Complex( -1.0f, 1.0f ),
		Complex( 0.0f, 1.0f )
	).Compare(
		Complex( -1.8577235439239f, -0.96193450888839f )
	);

	ok &= Carlson_RD(
		Complex( -2.0f, -1.0f ),
		Complex( 0.0f, -1.0f ),
		Complex( -1.0f, 1.0f )
	).Compare(
		Complex( 1.8249027393704f, -1.2218475784827f )
	);

	// ============================================================
	// RJ
	// ============================================================

	ok &= FloatEqual(
		Carlson_RJ( 0.0f, 1.0f, 2.0f, 3.0f ),
		0.77688623778582f
	);

	ok &= FloatEqual(
		Carlson_RJ( 2.0f, 3.0f, 4.0f, 5.0f ),
		0.14297579667157f
	);

	ok &= Carlson_RJ(
		Complex( 2.0f, 0.0f ),
		Complex( 3.0f, 0.0f ),
		Complex( 4.0f, 0.0f ),
		Complex( -1.0f, 1.0f )
	).Compare(
		Complex( 0.13613945827771f, -0.38207561624427f )
	);

	ok &= Carlson_RJ(
		Complex( 0.0f, 1.0f ),
		Complex( 0.0f, -1.0f ),
		Complex::Zero(),
		Complex( 2.0f, 0.0f )
	).Compare(
		Complex( 1.6490011662711f, 0.0f )
	);

	ok &= Carlson_RJ(
		Complex( -1.0f, 1.0f ),
		Complex( -1.0f, -1.0f ),
		Complex( 1.0f, 0.0f ),
		Complex( 2.0f, 0.0f )
	).Compare(
		Complex( 0.94148358841220f, 0.0f )
	);

	ok &= Carlson_RJ(
		Complex( 0.0f, 1.0f ),
		Complex( 0.0f, -1.0f ),
		Complex::Zero(),
		Complex( 1.0f, -1.0f )
	).Compare(
		Complex( 1.8260115229009f, 1.2290661908643f )
	);

	ok &= Carlson_RJ(
		Complex( -1.0f, 1.0f ),
		Complex( -1.0f, -1.0f ),
		Complex( 1.0f, 0.0f ),
		Complex( -3.0f, 1.0f )
	).Compare(
		Complex( -0.61127970812028f, -1.0684038390007f )
	);

	ok &= Carlson_RJ(
		Complex( -1.0f, 1.0f ),
		Complex( -2.0f, -1.0f ),
		Complex( 0.0f, -1.0f ),
		Complex( -1.0f, 1.0f )
	).Compare(
		Complex( 1.8249027393704f, -1.2218475784827f )
	);

	return ok;
}

CarlsonForm::CarlsonForm( float e1, float e2, float e3, float e4, float y, float x )
{
	Initialize( e1, e2, e3, e4, y, x );
}

void CarlsonForm::Initialize( float e1, float e2, float e3, float e4, float y, float x )
{
	m_a[0] = e1;
	m_a[1] = e2;
	m_a[2] = e3;
	m_a[3] = e4;

	m_x = x;
	m_y = y;

	m_case = FormType::CUBIC;

	//Carlson88(2.1)
	for (unsigned int i = 0; i < 4; ++i)
	{
		for (unsigned int j = 0; j < 4; ++j)
		{
			m_d[i][j] = m_a[i] - m_a[j];
		}
	}

	//Carlson88(2.2)
	m_Y1 = sqrtf( m_a[0] + m_y );
	m_Y2 = sqrtf( m_a[1] + m_y );
	m_Y3 = sqrtf( m_a[2] + m_y );
	m_Y4_2 = m_a[3] + m_y;

	// If x is infinite, then use Carlson88(2.10)
	m_U1_2 = m_Y1 * m_Y1;
	m_U2_2 = m_Y2 * m_Y2;
	m_U3_2 = m_Y3 * m_Y3;

	if (x != INFINITY)
	{
		//Carlson88(2.2)
		m_X1 = sqrtf( m_a[0] + m_x );
		m_X2 = sqrtf( m_a[1] + m_x );
		m_X3 = sqrtf( m_a[2] + m_x );
		m_X4_2 = m_a[3] + m_x;

		//Carlson88(2.3)
		const float limdel = m_x - m_y;
		m_U1_2 = (m_X1 * m_Y2 * m_Y3 + m_Y1 * m_X2 * m_X3) / limdel;
		m_U1_2 *= m_U1_2;
		m_U2_2 = m_U1_2 - d( 1, 2 ); //Carlson88(2.3)
		m_U3_2 = m_U1_2 - d( 1, 3 );
	}

	// Compute common terms
	W2();
	Q2();
	P2();
}


CarlsonForm::CarlsonForm( float e1, const Complex& e2, float e4, float y, float x )
{
	Initialize( e1, e2, e4, y, x );
}

void CarlsonForm::Initialize(float e1, const Complex & e2, float e4, float y, float x)
{
	m_a[0] = e1;
	m_a[3] = e4;
	m_g1 = -2.0f * e2.real;
	m_f1 = e2.real * e2.real + e2.imag * e2.imag;

	m_x = x;
	m_y = y;

	m_case = FormType::ONE_QUADRATIC_FACTOR;

	//Carlson91(2.4)
	m_e1 = sqrtf( m_f1 + m_g1 * x + x * x );
	m_n1 = sqrtf( m_f1 + m_g1 * y + y * y );

	//Carlson91(2.1)
	m_Y1 = sqrtf( m_a[0] + m_y);
	m_Y4_2 = m_a[3] + m_y;

	if (x != INFINITY)
	{
		//Carlson91(2.1)
		m_X1 = sqrtf( m_a[0] + m_x );
		m_X4_2 = m_a[3] + m_x;

		//Carlson91(3.3)
		m_U = (m_X1 * m_n1 + m_Y1 * m_e1) / (x - y);
	}
	else
	{
		//Carlson91(3.6)
		m_U = m_Y1;
	}

	//init d and c terms
	// use quadratic equation to solve for a2 and a3 in terms of g and f.
	// delta = sqrt( b^2 - 4ac )
	const float delta = sqrtf( m_g1 * m_g1 - 4.0f * m_f1 );
	m_a[1] = (-m_g1 + delta) / 2.0f;
	m_a[2] = (-m_g1 - delta) / 2.0f;

	//Carlson91(2.1, 2.3)
	for (unsigned int i = 0; i < 4; ++i)
	{
		for ( unsigned int j = 0; j < 4; ++j )
		{
			m_d[i][j] = m_a[i] - m_a[j];
			m_c[i][j] = sqrtf( 2.0f * m_f1 - m_g1 * ( m_a[i] + m_a[j] ) + 2.0f * m_a[i] * m_a[j]);
		}
	}

	// Compute common terms
	W2();
	Q2();
	P2();
	M2();
	L2();
}

CarlsonForm::CarlsonForm( const Complex& e1, const Complex& e3, float e4, bool useB5ZeroCase, float y, float x )
{
	Initialize( e1, e3, e4, useB5ZeroCase, y, x );
}


void CarlsonForm::Initialize( const Complex& e1, const Complex& e3, float e4, bool useB5ZeroCase, float y, float x )
{
	m_g1 = -2.0f * e1.real;
	m_f1 = e1.real * e1.real + e1.imag * e1.imag;

	m_g2 = -2.0f * e3.real;
	m_f2 = e3.real * e3.real + e3.imag * e3.imag;

	const float a5 = -e4;

	m_x = x;
	m_y = y;

	//Carlson92(2.1)
	m_e1 = sqrtf( m_f1 + m_g1 * m_x + m_x * m_x );
	m_e2 = sqrtf( m_f2 + m_g2 * m_x + m_x * m_x );
	m_n1 = sqrtf( m_f1 + m_g1 * m_y + m_y * m_y );
	m_n2 = sqrtf( m_f2 + m_g2 * m_y + m_y * m_y );

	//Carlson92(2.4)
	float xy = m_x - m_y;
	float xy2 = xy * xy;
	float theta1 = m_e1 * m_e1 + m_n1 * m_n1 - xy2;
	float theta2 = m_e2 * m_e2 + m_n2 * m_n2 - xy2;

	//Carlson92(2.5)
	float c1 = sqrtf( m_e1 * m_n1 * 2.0f + theta1 );
	float c2 = sqrtf( m_e2 * m_n2 * 2.0f + theta1 );

	//Carlson92(2.6)
	float u = (m_e1 * m_n2 + m_n1 * m_e2) / xy;
	float m = c1 * c2 / xy;
	float M2 = m * m;

	//Carlson92(2.7)
	float f[2] = { m_f1, m_f2 };
	float g[2] = { m_g1, m_g2 };
	for (unsigned int i = 0; i < 2; ++i)
	{
		for (unsigned int j = 0; j < 2; ++j)
		{
			m_del[i][j] = sqrtf( 2.0f * f[0] + 2.0f * f[1] - g[0] * g[1] );
		}
	}
	float del12_sqrd = del(1,2) * del( 1, 2 );

	//Carlson92(2.8)
	float delta = sqrtf( del12_sqrd * del12_sqrd - del(1,1) * del(1,1) * del(2,2) * del(2,2) );
	float delta_p = del12_sqrd + delta;
	float delta_n = del12_sqrd - delta;
	m_L2_p = M2 + delta_p;
	m_L2_n = M2 + delta_n;

	//Carlson92(2.11)
	float a15 = 2.0f * m_f1 - m_g1 * a5;
	float a25 = 2.0f * m_f2 - m_g2 * a5;
	float b15 = m_g1 - 2.0f * a5;
	float b25 = m_g2 - 2.0f * a5;
	m_beti5[0] = b15;
	m_beti5[1] = b25;

	//Carlson92(2.12)
	float y1 = useB5ZeroCase ? 1.0f : (a15 - b15 * a5) * 0.5f; //needed for H vs H_0 deffinition
	float y2 = (a25 - b25 * a5) * 0.5f;

	
	float A, w;
	if (useB5ZeroCase)
	{
		A = del(1,1) * del(1,1); //Carlson92(2.23)
		w = m_g1 - m_g2; //Carlson92(2.24)
	}
	else
	{
		A = del(1,1) * del(1,1) * y2 / y1; //Carlson92(2.13)
		w = (a15 * b25 - a25 * b15) / 2.0f; //Carlson92(2.14)
	}
	float O2 = M2 + A; //Carlson92(2.13)
	m_W2 = w * w;

	//Carlson92(2.15)
	float e5 = x + a5;
	float n5 = y + a5;

	if (useB5ZeroCase)
	{
		float ep1 = (m_g1 + 2.0f * m_x) / (m_e1 * 2.0f);
		float np1 = (m_g1 + 2.0f * m_y) / (m_n1 * 2.0f);
		m_X1 = (ep1 * m_e2 + np1 * m_n2) / -xy;
	}
	else
	{
		//Carlson92(2.17)
		m_X1 = m_n2 * e5 * (a15 + b15 * m_y) / m_n1;
		m_X1 += m_e2 * n5 * (a15 + b15 * m_x) / m_e1;
		m_X1 /= 2.0f * xy;
	}

	const float U2 = m_U * m_U;

	//Carlson92(2.18)
	float s = (M2 + del(1,2) * del(1,2)) * 0.5f - U2;
	float S2 = s * s;

	float T2, V2;
	float a2, b2;
	if ( useB5ZeroCase )
	{
		//Carlson92(2.26)
		float mu0 = 1.0f / (m_e1 * m_n1);
		float T0 = mu0 * s + 2.0f;
		T2 = T0 * T0;
		float V0 = mu0 * mu0 * (S2 + A * U2);
		V2 = V0 * V0;

		//Carlson92(2.27)
		a2 = s * O2 / U2 + A * m_U * 2.0f;
		a2 *= a2;
		b2 = (S2 / U2 + A) * (O2 * O2);
	}
	else
	{
		//Carlson92(2.19)
		float mu = (y1 * e5 * n5) / (m_e1 * m_n1);
		T2 = mu * s + (2.0f * y1 * y2);
		T2 *= T2;
		V2 = mu * mu * (S2 + A * U2);

		//Carlson92(2.20, 2.21)
		b2 = (S2 / U2 + A) * (O2 * O2);
		a2 = b2 + A * m_W2 / (y1 * y2);
	}

	//Carlson92(2.22)
	float RJ = Carlson_RJ( m_M2, m_L2_n, m_L2_p, O2 );
	float RC1 = Carlson_RC( a2, b2 );
	float RC2 = Carlson_RC( T2, V2 );
	m_H = del(1,1) * w * (RJ / 3.0f + RC1 / 2.0f) / (y1 * y1) - m_X1 * RC2;

}


float CarlsonForm::I_1c() const
{
	if (m_case == FormType::CUBIC)
	{
		//Carlson88(2.12)
		float RF = Carlson_RF( m_U3_2, m_U2_2, m_U1_2 );
		return RF * 2.0f;
	}
	else if (m_case == FormType::ONE_QUADRATIC_FACTOR)
	{
		// Carlson91(3.8)
		float RF = Carlson_RF( m_M2, m_L2_n, m_L2_p );
		return RF * 4.0f;
	}
}


float CarlsonForm::I_3c() const
{
	if (m_case == FormType::CUBIC)
	{
		float d12d13od14 = d( 1, 2 ) * d( 1, 3 ) / d( 1, 4 );

		// We know X4 and Y4 must always have the same sign.
		// If Y4 (and so X4) are complex, then the product X4*Y4 is negative.
		// This must be recovered to get the sign of RC
		bool m_YComplex = m_a[3] + m_y < 0.0f;
		const float RC_sign = m_YComplex ? -1.0f : 1.0f;

		//Carlson88(2.14)
		float RC = Carlson_RC( m_P2, m_Q2 ) * RC_sign;
		float RJ = Carlson_RJ( m_U3_2, m_U2_2, m_U1_2, m_W2 );
		return RJ * d12d13od14 * (-2.0f / 3.0f) + RC * 2.0f;
	}
	else if (m_case == FormType::ONE_QUADRATIC_FACTOR)
	{
		// Calrson91( 3.2 )
		float W2_p = m_M2 - (m_c[0][3] * m_c[0][3] + m_c[0][0] * m_c[3][3]) / d( 1, 4 );

		bool m_YComplex = m_a[3] + m_y < 0.0f;
		const float RC_sign = m_YComplex ? -1.0f : 1.0f;

		// Calrson91( 3.10 )
		float RC_UW = Carlson_RC( m_U * m_U, m_W2 );
		float RC_PQ = Carlson_RC( m_P2, m_Q2 ) * RC_sign;
		float RF = Carlson_RF( m_M2, m_L2_n, m_L2_p );
		float RJ = Carlson_RJ( m_M2, m_L2_n, m_L2_p, W2_p );
		float term1 = (2.0f * m_c[0][0]) / (3.0f * m_c[3][3]);
		float term2 = RJ * (-4.0f / d( 1, 4 ) * (m_c[0][3] * m_c[0][3] + m_c[0][0] * m_c[3][3])) - RF * 6.0F + RC_UW * 3.0f;
		float term3 = RC_PQ * 2.0f;
		float I_3c = term1 * term2 + term3;

		return I_3c;
	}
}

float CarlsonForm::K_2c_quadratic() const
{
	// Carlson91(3.1, 3.6)
	float A;
	float lastTerm_N2c = 0.0f;
	if (m_x == INFINITY)
	{
		// Carlson91(3.12) in K term we must solve for A(-1,1,1,-2)
		const float A_m1m1m10 = -1.0f / (m_Y1 * m_n1);
		const float A_m1m110 = -m_Y3 / ( m_Y1 * m_Y2 );
		const float A_m1m1m1m2 = -1.0f / (m_Y1 * m_n1 * m_Y4_2);
		const float A_m1m11m2 = m_d[2][3] * A_m1m1m1m2 + A_m1m1m10;
		const float A_m111m2 = m_d[1][3] * A_m1m11m2 + A_m1m110;
		A = A_m111m2;

	} else
	{
		//Carlson91(2.5) for A(-1,1,1,-2)
		A = m_e1 / (m_X1 * m_X4_2) - m_n1 / (m_Y1 * m_Y4_2);
		lastTerm_N2c = 2.0f / (m_X1 * m_Y1 * m_U);
	}

	//Calrson91( 3.5 )
	float bet1 = m_g1 - 2.0f * m_a[0];
	float rho = sqrtf( 2.0f ) * m_c[0][0] - bet1;

	// Carlson91( 3.11 )
	float term1 = sqrtf( 8.0f / (9.0f * m_c[0][0] * m_c[0][0]) );
	float RD = Carlson_RD( m_M2, m_L2_n, m_L2_p );
	float RF = Carlson_RF( m_M2, m_L2_n, m_L2_p );
	float term2 = RD * 4.0f * rho - RF * 6.0f + 3.0f / m_U;
	float N_2c = term1 * term2 + lastTerm_N2c;

	// Carlson91( 3.12 )
	float K_2c = N_2c * (m_c[0][0] * m_c[0][0] / 2.0f) - 2.0f * d( 1, 4 ) * A;

	return K_2c;
}


float CarlsonForm::K_2c_cubic() const
{
	float J_2c, A;
	if (m_x == INFINITY)
	{
		//Carlson88(2.17)
		float RD_term = (2.0f / 3.0f) * d( 1, 2 ) * d( 1, 3 ) * Carlson_RD( m_U3_2, m_U2_2, m_U1_2 );
		float last_term = (2.0f * d( 1, 3 ) * m_Y2) / (m_Y1 * m_Y3); //Carlson88(2.19)
		J_2c = RD_term + last_term;

		// A(1,1,-1,-2) is finite when x is inf when reducing using Recurance Relations.
		// Explained in another Carlson88-ATableOfEllipticInegralsOfTheThirdKind(4.8).
		const float A_m1m1m10 = -1.0f / (m_Y1 * m_Y2 * m_Y3);
		const float A_m1m1m1m2 = -1.0f / (m_Y1 * m_Y2 * m_Y3 * m_Y4_2);
		const float A_m11m1m2 = d( 2, 4 ) * A_m1m1m1m2 + A_m1m1m10;
		const float A_m11m10 = -m_Y2 / ( m_Y1 * m_Y3 );
		const float A_11m1m2 = d(1,4) * A_m11m1m2 + A_m11m10;
		A = A_11m1m2;
	}
	else
	{	
		//Carlson88(2.17)
		float RD_term = (2.0f / 3.0f) * d( 1, 2 ) * d( 1, 3 ) * Carlson_RD( m_U3_2, m_U2_2, m_U1_2 );
		float last_term = (2.0f * d( 1, 3 ) * m_X2 * m_Y2) / (m_X3 * m_Y3 * sqrtf( m_U1_2 ));
		J_2c = RD_term + last_term;

		// A[1,1,-1,-2] Carlson88(2.6)
		A = m_X1 * m_X2 / (m_X3 * m_X4_2) - m_Y1 * m_Y2 / (m_Y3 * m_Y4_2);
	}

	//Carlson88(2.59)
	return J_2c - 2.0f * d( 3, 4 ) * A;
}


float CarlsonForm::K_2c() const
{
	if (m_case == FormType::CUBIC)
		return K_2c_cubic();
	else if (m_case == FormType::ONE_QUADRATIC_FACTOR)
		return K_2c_quadratic();
}


void CarlsonForm::W2()
{
	if ( m_case == FormType::CUBIC )
	{	
		//Carlson88(2.4)
		m_W2 = m_U1_2 - d( 1, 2 ) * d( 1, 3 ) / d( 1, 4 );
		return;
	}
	else if (m_case == FormType::ONE_QUADRATIC_FACTOR)
	{
		//Carlson91(3.3)
		m_W2 = m_U * m_U - m_c[0][0] * m_c[0][0] / (2.0f * d( 1, 4 ));
		return;
	}
}


void CarlsonForm::Q2()
{
	if (m_x == INFINITY)
	{
		//Carlson88(2.11) and Carlson91(3.7)
		float q = m_Y4_2 / ( m_Y1 * m_Y1 );
		m_Q2 = q * m_W2;
		return;
	}

	//Carlson88(2.5) and Carlson91(3.4)
	// In case with one quadratic factor, we lose sign of Q since we only store the square.
	m_Q2 = ( m_X4_2 * m_Y4_2 * m_W2 ) / ( m_X1 * m_X1 * m_Y1 * m_Y1 );
	return;
}


void CarlsonForm::P2()
{
	if (m_case == FormType::CUBIC)
	{
		//Carlson88(2.5)
		m_P2 = m_Q2 + d( 2, 4 ) * d( 3, 4 ) / d( 1, 4 );
		return;
	}
	else if (m_case == FormType::ONE_QUADRATIC_FACTOR)
	{
		//Carlson91(3.4)
		m_P2 = m_Q2 + m_c[3][3] * m_c[3][3] / 2.0f / d( 1, 4 );
		return;
	}
}

void CarlsonForm::M2()
{
	if (m_case == FormType::ONE_QUADRATIC_FACTOR)
	{
		if ( m_x == INFINITY)
		{
			//Carlson91(3.6)
			m_M2 = 2.0f * m_n1 + m_g1 + 2.0f * m_y;
			return;
		}

		//Carlson91(3.1)
		float term = 2.0f * m_e1 * m_n1 + 2.0f * m_f1 + m_g1 * ( m_x + m_y ) + 2.0f * m_x * m_y;
		float m = (m_X1 + m_Y1) * sqrtf( term ) / (m_x - m_y);
		m_M2 = m * m;
	}
}

void CarlsonForm::L2()
{
	//Carlson91(3.2)
	float bet1 = m_g1 - 2.0f * m_a[0];
	m_L2_p = m_M2 - bet1 + sqrtf( 2.0f ) * m_c[0][0];
	m_L2_n = m_M2 - bet1 - sqrtf( 2.0f ) * m_c[0][0];
}

float CarlsonForm::RF() const
{
	if ( m_case != FormType::TWO_QUADRATIC_FACTORS )
		return 0.0f;

	return Carlson_RF( m_M2, m_L2_n, m_L2_p );
}

float CarlsonForm::SIGMA() const
{
	if (m_case != FormType::TWO_QUADRATIC_FACTORS)
		return 0.0f;
	
	float del11_sqrd = del( 1, 1 ) * del( 1, 1 );
	float del12_sqrd = del( 1, 2 ) * del( 1, 2 );

	//Carlson92(2.4)
	float xy = m_x - m_y;
	float xy2 = xy * xy;
	float theta1 = m_e1 * m_e1 + m_n1 * m_n1 - xy2;
	float theta2 = m_e2 * m_e2 + m_n2 * m_n2 - xy2;

	//Carlson92(2.8)
	float delta = sqrtf( del12_sqrd * del12_sqrd - del( 1, 1 ) * del( 1, 1 ) * del( 2, 2 ) * del( 2, 2 ) );
	float delta_p = del12_sqrd + delta;
	float delta_n = del12_sqrd - delta;

	//Carlson92(2.9)
	float RD = Carlson_RD( m_M2, m_L2_n, m_L2_p );
	float G = delta * delta_p * RD * (2.0f / 3.0f) +
		delta / m_U * 0.5f +
		(del12_sqrd * theta1 - del11_sqrd * theta2) / 4.0f * m_e1 * m_n1 * m_U;

	//Carlson92(2.2)
	float ep1 = (m_g1 + 2.0f * m_x) / (m_e1 * 2.0f);
	float np1 = (m_g1 + 2.0f * m_y) / (m_n1 * 2.0f);

	//Carlson92(2.3)
	float B = ep1 * m_e2 - np1 * m_n2;

	//Carlson92(2.10)
	return G - delta_p * RF() + B;
}


bool ValidateCarlsonSymmetricForms()
{
	bool ok = true;

	auto FloatEqual = []( float a, float b, float eps = 1e-5f ) -> bool
		{
			return fabsf( a - b ) <= eps;
		};

	CarlsonForm form;

	// ============================================================
	// [-1,-1,-1,-2]
	// ============================================================

	//Cubic
	//Integrate[Divide[1,\(40)4+t\(41)Sqrt[\(40)1+t\(41)\(40)2+t\(41)\(40)3+t\(41)]],{t,3.5,3.8}]
	form.Initialize( 1.0f, 2.0f, 3.0f, 4.0f, 3.5f, 3.8f );
	float fm1 = (form.I_3c() - form.I_1c()) / form.d( 1, 4 );
	ok &= FloatEqual(
		fm1,
		0.00296942f
	);

	//Cubic
	//Integrate[Divide[1,\(40)-4+t\(41)Sqrt[\(40)-1+t\(41)\(40)-2+t\(41)\(40)-3+t\(41)]],{t,3.5,3.8}]
	form.Initialize( -1.0f, -2.0f, -3.0f, -4.0f, 3.5f, 3.8f );
	float fm2 = (form.I_3c() - form.I_1c()) / form.d( 1, 4 );
	ok &= FloatEqual(
		fm2,
		-0.53436f
	);

	//One Quadratic Factor
	//Integrate[Divide[1,\(40)4+t\(41)Sqrt[\(40)1+t\(41)\(40)3 - 2t + Power[t,2]\(41)]],{t,3.5,3.8}]
	form.Initialize( 1.0f, Complex( 1.0f, sqrtf( 2.0f ) ), 4.0f, 3.5f, 3.8f );
	float fm3 = (form.I_3c() - form.I_1c()) / form.d( 1, 4 );
	ok &= FloatEqual(
		fm3,
		0.00606313f
	);

	//One Quadratic Factor
	//Integrate[Divide[1,\(40)-4+t\(41)Sqrt[\(40)-1+t\(41)\(40)3 + 2t + Power[t,2]\(41)]],{t,3.5,3.8}]
	form.Initialize( -1.0f, Complex( 1.0f, sqrtf( 2.0f ) ) * -1.0f, -4.0f, 3.5f, 3.8f );
	float fm4 = (form.I_3c() - form.I_1c()) / form.d( 1, 4 );
	ok &= FloatEqual(
		fm4,
		-0.114917f
	);


	// ============================================================
	// [-1,-1,-1,-4]
	// ============================================================

	//Cubic
	//Integrate[Divide[1,Square[\(40)4+t\(41)]Sqrt[\(40)1+t\(41)\(40)2+t\(41)\(40)3+t\(41)]],{t,0,8}]
	form.Initialize( 1.0f, 2.0f, 3.0f, 4.0f, 0.0f, 8.0f );
	float fm5;
	{
		const float I3c_term = form.I_3c() / (2.0f * form.d( 1, 4 )) * (1.0f / form.r( 1, 4 ) + 1.0f / form.r( 2, 4 ) + 1.0f / form.r( 3, 4 ));
		const float K2c_term = form.K_2c() / (2.0f * form.d( 1, 4 ) * form.d( 2, 4 ) * form.d( 3, 4 ));
		const float I1c_term = form.I_1c() * (1.0f - form.r( 1, 2 ) * form.r( 1, 3 ) / (2.0f * form.r( 2, 4 ) * form.r( 3, 4 ))) / (form.d( 1, 4 ) * form.d( 1, 4 ));
		fm5 = -I3c_term + K2c_term + I1c_term;
	}
	ok &= FloatEqual(
		fm5,
		0.0261882f
	);

	//Cubic
	//Integrate[Divide[1,Square[\(40)-4+t\(41)]Sqrt[\(40)-1+t\(41)\(40)-2+t\(41)\(40)-3+t\(41)]],{t,3.1,3.8}]
	form.Initialize( -1.0f, -2.0f, -3.0f, -4.0f, 3.1f, 3.8f );
	float fm6;
	{
		form.Initialize( 1.0f, 2.0f, 3.0f, 4.0, 3.0f, 8.0f );
		const float I3c_term = form.I_3c() / (2.0f * form.d( 1, 4 )) * (1.0f / form.r( 1, 4 ) + 1.0f / form.r( 2, 4 ) + 1.0f / form.r( 3, 4 ));
		const float K2c_term = form.K_2c() / (2.0f * form.d( 1, 4 ) * form.d( 2, 4 ) * form.d( 3, 4 ));
		const float I1c_term = form.I_1c() * (1.0f - form.r( 1, 2 ) * form.r( 1, 3 ) / (2.0f * form.r( 2, 4 ) * form.r( 3, 4 ))) / (form.d( 1, 4 ) * form.d( 1, 4 ));
		fm6 = -I3c_term + K2c_term + I1c_term;
	}
	ok &= FloatEqual(
		fm6,
		2.63298f
	);

	//One Quadratic Factor
	//Integrate[Divide[1,Square[\(40)4+t\(41)]Sqrt[\(40)1+t\(41)\(40)2+t\(41)\(40)3+t\(41)]],{t,0,8}]
	form.Initialize( 1.0f, Complex( 1.0f, sqrtf( 2.0f ) ), 4.0f, 0.0f, 8.0f );
	float fm7;
	{
		const float a1 = form.a( 1 );
		const float a4 = form.a( 4 );
		const float d14 = form.d( 1, 4 );
		float d24d34 = a4 * a4 + a4 * form.g() + form.f();
		float r123 = 1.0f / d14 - (form.g() + 2.0f * a4) / d24d34;
		float d12d13 = a1 * a1 + a1 * form.g() + form.f();

		const float I3c_term = form.I_3c() / (2.0f * d14) * r123;
		const float K2c_term = form.K_2c() / (2.0f * d14 * d24d34);
		const float I1c_term = form.I_1c() * (1.0f - d12d13 / (2.0f * d24d34)) / (d14 * d14);
		fm7 = -I3c_term + K2c_term + I1c_term;
	}
	ok &= FloatEqual(
		fm7,
		0.0261882f
	);

	return ok;
}