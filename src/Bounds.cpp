#include "..\include\Bounds.h"


const Bounds &Bounds::operator=( const Bounds &rhs )
{
	m_vMin = rhs.m_vMin;
	m_vMax = rhs.m_vMax;
	return *this;
}
bool Bounds::DoesIntersect( const Bounds &rhs ) const
{
	if ( m_vMax[0] < rhs.m_vMin[0] || m_vMax[1] < rhs.m_vMin[1] || m_vMax[2] < rhs.m_vMin[2] )
		return false;

	if ( rhs.m_vMax[0] < m_vMin[0] || rhs.m_vMax[1] < m_vMin[1] || rhs.m_vMax[2] < m_vMin[2] )
		return false;

	return true;
}

void Bounds::Expand( const Vec3 *pts, const int num )
{
	for ( int i = 0; i < num; i++ )
	{
		Expand( pts[i] );
	}
}

void Bounds::Expand( const Vec3& rhs )
{
	if ( rhs[0] < m_vMin[0] )
		m_vMin[0] = rhs[0];
	if ( rhs[1] < m_vMin[1] )
		m_vMin[1] = rhs[1];
	if ( rhs[2] < m_vMin[2] )
		m_vMin[2] = rhs[2];
	if ( rhs[0] > m_vMax[0] )
		m_vMax[0] = rhs[0];
	if ( rhs[1] > m_vMax[1] )
		m_vMax[1] = rhs[1];
	if ( rhs[2] > m_vMax[2] )
		m_vMax[2] = rhs[2];
}
void Bounds::Expand( const Bounds& rhs )
{
	Expand( rhs.m_vMin );
	Expand( rhs.m_vMax );
}