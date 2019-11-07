#include "vtk_math.h"

template<class T1, class T2>
inline void vtkQuaternionToMatrix3x3(const T1 quat[4], T2 A[3][3])
{
	T2 ww = quat[0] * quat[0];
	T2 wx = quat[0] * quat[1];
	T2 wy = quat[0] * quat[2];
	T2 wz = quat[0] * quat[3];

	T2 xx = quat[1] * quat[1];
	T2 yy = quat[2] * quat[2];
	T2 zz = quat[3] * quat[3];

	T2 xy = quat[1] * quat[2];
	T2 xz = quat[1] * quat[3];
	T2 yz = quat[2] * quat[3];

	T2 rr = xx + yy + zz;
	// normalization factor, just in case quaternion was not normalized
	T2 f = 1 / (ww + rr);
	T2 s = (ww - rr)*f;
	f *= 2;

	A[0][0] = xx * f + s;
	A[1][0] = (xy + wz)*f;
	A[2][0] = (xz - wy)*f;

	A[0][1] = (xy - wz)*f;
	A[1][1] = yy * f + s;
	A[2][1] = (yz + wx)*f;

	A[0][2] = (xz + wy)*f;
	A[1][2] = (yz - wx)*f;
	A[2][2] = zz * f + s;
}

//----------------------------------------------------------------------------
void vtkMath::QuaternionToMatrix3x3(const float quat[4], float A[3][3])
{
	vtkQuaternionToMatrix3x3(quat, A);
}

//----------------------------------------------------------------------------
void vtkMath::QuaternionToMatrix3x3(const double quat[4], double A[3][3])
{
	vtkQuaternionToMatrix3x3(quat, A);
}

void vtkMatrix4x4::Invert(const double inElements[16],
	double outElements[16])
{
	// inverse( original_matrix, inverse_matrix )
	// calculate the inverse of a 4x4 matrix
	//
	//     -1
	//     A  = ___1__ adjoint A
	//         det A
	//

	// calculate the 4x4 determinent
	// if the determinent is zero,
	// then the inverse matrix is not unique.

	double det = vtkMatrix4x4::Determinant(inElements);
	if (det == 0.0)
	{
		return;
	}

	// calculate the adjoint matrix
	vtkMatrix4x4::Adjoint(inElements, outElements);

	// scale the adjoint matrix to get the inverse
	for (int i = 0; i < 16; i++)
	{
		outElements[i] /= det;
	}
}

//----------------------------------------------------------------------------
double vtkMatrix4x4::Determinant(const double elem[16])
{
	double a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4;

	// assign to individual variable names to aid selecting
	//  correct elements

	a1 = elem[0];  b1 = elem[1];  c1 = elem[2];  d1 = elem[3];
	a2 = elem[4];  b2 = elem[5];  c2 = elem[6];  d2 = elem[7];
	a3 = elem[8];  b3 = elem[9];  c3 = elem[10]; d3 = elem[11];
	a4 = elem[12]; b4 = elem[13]; c4 = elem[14]; d4 = elem[15];

	return a1 * vtkMath::Determinant3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4)
		- b1 * vtkMath::Determinant3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4)
		+ c1 * vtkMath::Determinant3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4)
		- d1 * vtkMath::Determinant3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);
}

//----------------------------------------------------------------------------
void vtkMatrix4x4::Adjoint(const double elem[16], double outElem[16])
{
	//
	//   adjoint( original_matrix, inverse_matrix )
	//
	//     calculate the adjoint of a 4x4 matrix
	//
	//      Let  a   denote the minor determinant of matrix A obtained by
	//           ij
	//
	//      deleting the ith row and jth column from A.
	//
	//                    i+j
	//     Let  b   = (-1)    a
	//          ij            ji
	//
	//    The matrix B = (b  ) is the adjoint of A
	//                     ij
	//
	double a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4;

	// assign to individual variable names to aid
	// selecting correct values

	a1 = elem[0];  b1 = elem[1];  c1 = elem[2];  d1 = elem[3];
	a2 = elem[4];  b2 = elem[5];  c2 = elem[6];  d2 = elem[7];
	a3 = elem[8];  b3 = elem[9];  c3 = elem[10]; d3 = elem[11];
	a4 = elem[12]; b4 = elem[13]; c4 = elem[14]; d4 = elem[15];

	// row column labeling reversed since we transpose rows & columns

	outElem[0] =
		vtkMath::Determinant3x3(b2, b3, b4, c2, c3, c4, d2, d3, d4);
	outElem[4] =
		-vtkMath::Determinant3x3(a2, a3, a4, c2, c3, c4, d2, d3, d4);
	outElem[8] =
		vtkMath::Determinant3x3(a2, a3, a4, b2, b3, b4, d2, d3, d4);
	outElem[12] =
		-vtkMath::Determinant3x3(a2, a3, a4, b2, b3, b4, c2, c3, c4);

	outElem[1] =
		-vtkMath::Determinant3x3(b1, b3, b4, c1, c3, c4, d1, d3, d4);
	outElem[5] =
		vtkMath::Determinant3x3(a1, a3, a4, c1, c3, c4, d1, d3, d4);
	outElem[9] =
		-vtkMath::Determinant3x3(a1, a3, a4, b1, b3, b4, d1, d3, d4);
	outElem[13] =
		vtkMath::Determinant3x3(a1, a3, a4, b1, b3, b4, c1, c3, c4);

	outElem[2] =
		vtkMath::Determinant3x3(b1, b2, b4, c1, c2, c4, d1, d2, d4);
	outElem[6] =
		-vtkMath::Determinant3x3(a1, a2, a4, c1, c2, c4, d1, d2, d4);
	outElem[10] =
		vtkMath::Determinant3x3(a1, a2, a4, b1, b2, b4, d1, d2, d4);
	outElem[14] =
		-vtkMath::Determinant3x3(a1, a2, a4, b1, b2, b4, c1, c2, c4);

	outElem[3] =
		-vtkMath::Determinant3x3(b1, b2, b3, c1, c2, c3, d1, d2, d3);
	outElem[7] =
		vtkMath::Determinant3x3(a1, a2, a3, c1, c2, c3, d1, d2, d3);
	outElem[11] =
		-vtkMath::Determinant3x3(a1, a2, a3, b1, b2, b3, d1, d2, d3);
	outElem[15] =
		vtkMath::Determinant3x3(a1, a2, a3, b1, b2, b3, c1, c2, c3);
}

double vtkMath::Determinant3x3(double a1, double a2, double a3,
	double b1, double b2, double b3,
	double c1, double c2, double c3)
{
	return (a1 * vtkMath::Determinant2x2(b2, b3, c2, c3)
		- b1 * vtkMath::Determinant2x2(a2, a3, c2, c3)
		+ c1 * vtkMath::Determinant2x2(a2, a3, b2, b3));
}

double vtkMath::Determinant2x2(double a, double b, double c, double d) {
	return a * d - b * c;
}

void vtkMatrix4x4::DeepCopy(double destination[16], const double source[16])
{
	for (int i = 0; i < 16; i++)
	{
		destination[i] = source[i];
	}
}

template<class T1, class T2, class T3>
void vtkMatrix4x4MultiplyPoint(T1 elem[16], T2 in[4], T3 out[4])
{
	T3 v1 = in[0];
	T3 v2 = in[1];
	T3 v3 = in[2];
	T3 v4 = in[3];

	out[0] = v1 * elem[0] + v2 * elem[1] + v3 * elem[2] + v4 * elem[3];
	out[1] = v1 * elem[4] + v2 * elem[5] + v3 * elem[6] + v4 * elem[7];
	out[2] = v1 * elem[8] + v2 * elem[9] + v3 * elem[10] + v4 * elem[11];
	out[3] = v1 * elem[12] + v2 * elem[13] + v3 * elem[14] + v4 * elem[15];
}

void vtkMatrix4x4::MultiplyPoint(const double elements[16],
	const double in[4], double result[4])
{
	vtkMatrix4x4MultiplyPoint(elements, in, result);
}

void vtkMatrix4x4::Multiply4x4(const double a[16], const double b[16], double c[16])
{
	double tmp[16];

	for (int i = 0; i < 16; i += 4)
	{
		for (int j = 0; j < 4; j++)
		{
			tmp[i + j] = a[i + 0] * b[j + 0] +
				a[i + 1] * b[j + 4] +
				a[i + 2] * b[j + 8] +
				a[i + 3] * b[j + 12];
		}
	}

	for (int k = 0; k < 16; k++)
	{
		c[k] = tmp[k];
	}
}