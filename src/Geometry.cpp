#include "..\include\Geometry.h"
#include <math.h>

unsigned int FurthestVertexFromPlane( const unsigned int iP1, const unsigned int iP2, const unsigned int iP3, bool *pInverted, const unsigned int m_nVCount, const Vec3 *pVertices )
{
	const Vec3 vNormal = ( pVertices[iP2] - pVertices[iP1] ).Cross( pVertices[iP3] - pVertices[iP1] ).Normal();
	const Vec3 &vPlanePos = pVertices[iP1];

	unsigned int iMaxIdx = 0;
	float fMaxDist = ( pVertices[0] - pVertices[iP1] ).Dot( vNormal );
	*pInverted = fMaxDist < 0.0f;
	fMaxDist = fabsf( fMaxDist );
	for ( unsigned int i = 1; i < m_nVCount; ++i )
	{
		const float fDist = ( pVertices[i] - pVertices[iP1] ).Dot( vNormal );
		if ( fabsf( fDist ) > fMaxDist )
		{
			*pInverted = fDist < 0.0f;
			fMaxDist = fabsf( fDist );
			iMaxIdx = i;
		}
	}
	return iMaxIdx;
}


float DistanceFromTriangle( const Vec3 &a, const Vec3 &b, const Vec3 &c, const Vec3 &pt )
{
	Vec3 ab = b - a;
	Vec3 ac = c - a;
	Vec3 normal = ab.Cross( ac );
	Vec3 ray = pt - a;
	return normal.Normal().Dot(ray);

	float dot = ray.Dot( normal );
	float normLenSq = normal.Dot( normal );
	if (normLenSq < 1e-12f) {
		// Degenerate triangle (normal is almost zero)
		return 1000000.0f;
	}

	return dot / sqrt( normLenSq );
}

float SignedTetrahedronVolume( const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d )
{
	const Vec3 vCross = ( b - a ).Cross( c - a );
	const float fDot = vCross.Dot( d - a );
	return fDot / 6.0f;
}

float DistanceFromLine( const Vec3& a, const Vec3& b, const Vec3& pt )
{
	Vec3 vAB = b - a;
	Vec3 vAP = pt - a;
	return sqrtf( vAB.Cross( vAP ).LengthSquared() / vAB.LengthSquared() );
}

unsigned int FurthestVertexFromLine( const unsigned int iEpt1, const unsigned int iEpt2, const unsigned int m_nVCount, const Vec3 *pVertices )
{
	const Vec3 vVec = pVertices[iEpt2] - pVertices[iEpt1];
	const float fVecLength = vVec.Length();

	unsigned int iMaxIdx = 0;
	float fMaxDist = -1.0f;
	for ( unsigned int i = 0; i < m_nVCount; ++i )
	{
		const Vec3 vPos = pVertices[i] - pVertices[iEpt1];
		const float fDist = fabsf( vVec[1] * vPos[0] - vVec[0] * vPos[1] ) / fVecLength;
		if ( fDist > fMaxDist )
		{
			fMaxDist = fDist;
			iMaxIdx = i;
		}
	}
	return iMaxIdx;
}

unsigned int FurthestVertexFromVertex( const unsigned int iSampleVidx, const unsigned int m_nVCount, const Vec3 *pVertices )
{
	const Vec3 &vSamplePos = pVertices[iSampleVidx];

	unsigned int iMaxIdx = 0;
	float fMaxDistSqr = ( vSamplePos - pVertices[0] ).LengthSquared();
	for ( unsigned int i = 1; i < m_nVCount; ++i )
	{
		const float fDistSqr = ( vSamplePos - pVertices[i] ).LengthSquared();
		if ( fDistSqr > fMaxDistSqr )
		{
			fMaxDistSqr = fDistSqr;
			iMaxIdx = i;
		}
	}
	return iMaxIdx;
}

Vec3 NormalFromThreeVIDx( unsigned int iVidx1, unsigned int iVidx2, unsigned int iVidx3, const Vec3 *pVertices )
{
	const Vec3 vVec1 = ( pVertices[iVidx2] - pVertices[iVidx1] ).Normal();
	const Vec3 vVec2 = ( pVertices[iVidx3] - pVertices[iVidx1] ).Normal();
	return ( vVec1.Cross( vVec2 ) ).Normal();
}


Vec3 ProjectPointToPlane( const Vec3& vP, const Vec3& vN, const Vec3& v )
{
	const Vec3 vN_normalized = vN.Normal();
	const float fDist = (v - vP).Dot( vN_normalized );
	return v - vN_normalized * fDist;
}