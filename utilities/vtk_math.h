#pragma once

namespace vtkMath {

	void QuaternionToMatrix3x3(const float quat[4], float A[3][3]);
	void QuaternionToMatrix3x3(const double quat[4], double A[3][3]);
	double Determinant3x3(double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3);

	double Determinant2x2(double a, double b, double c, double d);
}

namespace vtkMatrix4x4 {
	void Invert(const double inElements[16], double outElements[16]);
	double Determinant(const double elem[16]);
	void Adjoint(const double elem[16], double outElem[16]);
	void DeepCopy(double destination[16], const double source[16]);

	void MultiplyPoint(const double elements[16], const double in[4], double out[4]);
	void Multiply4x4(const double a[16], const double b[16], double c[16]);
} 