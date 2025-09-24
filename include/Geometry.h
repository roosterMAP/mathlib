#pragma once

#ifndef __GEOMETRY_H_INCLUDE__
#define __GEOMETRY_H_INCLUDE__

#include "Vector.h"

// returns the index of the vertex furhtest from the plane defined by the tri <iP1, iP2, iP3>. pInverted is true is the point is behind the plane.
unsigned int FurthestVertexFromPlane( const unsigned int iP1, const unsigned int iP2, const unsigned int iP3, bool *pInverted, const unsigned int m_nVCount, const Vec3 *pVertices );

// returns the signed distance from pt to the triangle <a, b, c>
float DistanceFromTriangle( const Vec3 &a, const Vec3 &b, const Vec3 &c, const Vec3 &pt );

// returns the signed volume of the tetrahedrion < a, b, c, d > with base triangle < a, b, c >
float SignedTetrahedronVolume( const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d );

// returns the shortest distace from pt to line <a, b>
float DistanceFromLine( const Vec3 &a, const Vec3 &b, const Vec3 &pt );

// returns the index of the vertex furthest from the line defined by < iEpt1, iEpt2 >
unsigned int FurthestVertexFromLine( const unsigned int iEpt1, const unsigned int iEpt2, const unsigned int m_nVCount, const Vec3 *pVertices );

// returns the index of the vertex furthest from position at index iSampleVidx in pVertices
unsigned int FurthestVertexFromVertex( const unsigned int iSampleVidx, const unsigned int m_nVCount, const Vec3 *pVertices );

// returns the normal of the triangle formed by triangle < iVidx1, iVidx2, iVidx3 >
Vec3 NormalFromThreeVIDx( unsigned int iVidx1, unsigned int iVidx2, unsigned int iVidx3, const Vec3 *pVertices );

Vec3 ProjectPointToPlane( const Vec3 &vP, const Vec3 &vN, const Vec3 &v );

#endif