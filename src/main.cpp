#include "..\include\Vector.h"
#include "..\include\Matrix.h"
#include "..\include\Quaternion.h"

int main() {
	float a_data[15] = { 1, 5, 2, 1, 0, 2, 5, 8, 1, 1, 2, 2, 6, 5, 7 };
	float b_data[15] = { 1, 1, 5, 2, 9, 1, 2, 2, 6, 4 };
	MatN a( 3, 5, a_data );
	MatN b( 5, 2, b_data );

	VecN v1 = a.GetRowVec( 1 );
	for ( unsigned int i = 0; i < v1.size(); i++ ) {
		const float val = v1[i];
		v1[i] = 0.0f;
	}
	//a.SetRowVec( 1, &v1 );
	VecN v2 = a.GetColVec( 2 );
	for ( unsigned int i = 0; i < v2.size(); i++ ) {
		const float val = v2[i];
		v2[i] = 0.0f;
	}
	//a.SetColVec( 2, &v2 );
	
	MatN c = a * b;
	for ( unsigned int i = 0; i < c.ColumnCount() * c.RowCount(); i++ ) {
		float val = c[i];
		val = 0.0f;
	}

	float v3_data[15] = { 1, 5, 9, 2, 6 };
	VecN v3( 5, v3_data );
	VecN v4 = a * v3;
	for ( unsigned int i = 0; i < v4.size(); i++ ) {
		float val = v4[i];
		val = 0.0f;
	}

	float v5_data[15] = { 1, 5, 9, 2, 6 };
	VecN v5( 5, v5_data );
	float c_data[15] = { 1, 2, 2, 5, 5, 2, 2, 8, 6, 1, 1, 5, 0, 1, 7 };
	MatN g( 5, 3, c_data );
	VecN v6 = v5 * g;
	for ( unsigned int i = 0; i < v6.size(); i++ ) {
		float val = v6[i];
		val = 0.0f;
	}

	float d_data[15] = { 1, 5, 2, 2, 5, 8, 2, 2, 6 };
	float e_data[15] = { 1, 1, 5, 5, 2, 1, 9, 1, 1 };
	MatN d( 3, d_data );
	MatN e( 3, e_data );
	MatN f = d * e;
	float f_data[9] = { f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8] };

	float poopy[6][6];

	float h_data[36] = { 5, -2, 2, 7, 1, 3, 1, 0, 0, 3, 1, 1, -3, 1, 5, 0, 0, 2, 3, -1, -9, 4, 0, 9, 0, 0, 0, 0, 1, 4, 0, 0, 2, 2, 0, 1 };
	MatN h( 6, h_data );
	MatN h_inv( 6 );
	if ( !h.Inverse( &h_inv ) ) {
		return 0;
	}

	h.Transposed();

	return 1;

}

int n = main();