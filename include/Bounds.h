#pragma once

#ifndef __BOUNDS_H_INCLUDE__
#define __BOUNDS_H_INCLUDE__


#include "Vector.h"

class Bounds
{
public:
	Bounds() { Clear(); }
	Bounds( const Bounds &rhs ) : m_vMin( rhs.m_vMin ), m_vMax( rhs.m_vMax ) {}
	const Bounds &operator=( const Bounds &rhs );
	~Bounds() {}
	void Clear() { m_vMin = Vec3( 1e6 ); m_vMax = Vec3( -1e6 ); }
	bool DoesIntersect( const Bounds &rhs ) const;
	void Expand( const Vec3 *pts, const int num );
	void Expand( const Vec3 &rhs );
	void Expand( const Bounds &rhs );
	float WidthX() const { return m_vMax[0] - m_vMin[0]; }
	float WidthY() const { return m_vMax[1] - m_vMin[1]; }
	float WidthZ() const { return m_vMax[2] - m_vMin[2]; }
public:
	Vec3 m_vMin;
	Vec3 m_vMax;
};

#endif //__BOUNDS_H_INCLUDE__