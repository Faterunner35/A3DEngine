#include "Matrxi4x4.h"

namespace mathx {

	bool Inverse(const Matrix4x4& a, Matrix4x4& out, double tolerance)
	{
		const double m00 = a[0], m01 = a[1], m02 = a[2], m03 = a[3];
		const double m10 = a[4], m11 = a[5], m12 = a[6], m13 = a[7];
		const double m20 = a[8], m21 = a[9], m22 = a[10], m23 = a[11];
		const double m30 = a[12], m31 = a[13], m32 = a[14], m33 = a[15];

		// Шесть миноров 2x2 из строк 2,3 — по одному на каждую пару столбцов
		double v0 = m20 * m31 - m21 * m30;
		double v1 = m20 * m32 - m22 * m30;
		double v2 = m20 * m33 - m23 * m30;
		double v3 = m21 * m32 - m22 * m31;
		double v4 = m21 * m33 - m23 * m31;
		double v5 = m22 * m33 - m23 * m32;

		// Дополнения 3x3 собираются из готовых 2x2
		const double t00 = +(v5 * m11 - v4 * m12 + v3 * m13);
		const double t10 = -(v5 * m10 - v2 * m12 + v1 * m13);
		const double t20 = +(v4 * m10 - v2 * m11 + v0 * m13);
		const double t30 = -(v3 * m10 - v1 * m11 + v0 * m12);

		// det берётся из уже посчитанных дополнений — без отдельного прохода
		const double det = t00 * m00 + t10 * m01 + t20 * m02 + t30 * m03;
		if (std::abs(det) <= tolerance) return false;
		const double id = 1.0 / det;

		out[0] = t00 * id; out[4] = t10 * id; out[8] = t20 * id; out[12] = t30 * id;
		out[1] = -(v5 * m01 - v4 * m02 + v3 * m03) * id;
		out[5] = +(v5 * m00 - v2 * m02 + v1 * m03) * id;
		out[9] = -(v4 * m00 - v2 * m01 + v0 * m03) * id;
		out[13] = +(v3 * m00 - v1 * m01 + v0 * m02) * id;

		// Пересчёт миноров для строк 1,3
		v0 = m10 * m31 - m11 * m30;  v1 = m10 * m32 - m12 * m30;  v2 = m10 * m33 - m13 * m30;
		v3 = m11 * m32 - m12 * m31;  v4 = m11 * m33 - m13 * m31;  v5 = m12 * m33 - m13 * m32;
		out[2] = +(v5 * m01 - v4 * m02 + v3 * m03) * id;
		out[6] = -(v5 * m00 - v2 * m02 + v1 * m03) * id;
		out[10] = +(v4 * m00 - v2 * m01 + v0 * m03) * id;
		out[14] = -(v3 * m00 - v1 * m01 + v0 * m02) * id;

		// Пересчёт для строк 1,2
		v0 = m21 * m10 - m20 * m11;  v1 = m22 * m10 - m20 * m12;  v2 = m23 * m10 - m20 * m13;
		v3 = m22 * m11 - m21 * m12;  v4 = m23 * m11 - m21 * m13;  v5 = m23 * m12 - m22 * m13;
		out[3] = -(v5 * m01 - v4 * m02 + v3 * m03) * id;
		out[7] = +(v5 * m00 - v2 * m02 + v1 * m03) * id;
		out[11] = -(v4 * m00 - v2 * m01 + v0 * m03) * id;
		out[15] = +(v3 * m00 - v1 * m01 + v0 * m02) * id;
		return true;
	}

	double Determinant(const Matrix4x4& a)
	{
		const double v0 = a[8] * a[13] - a[9] * a[12];
		const double v1 = a[8] * a[14] - a[10] * a[12];
		const double v2 = a[8] * a[15] - a[11] * a[12];
		const double v3 = a[9] * a[14] - a[10] * a[13];
		const double v4 = a[9] * a[15] - a[11] * a[13];
		const double v5 = a[10] * a[15] - a[11] * a[14];
		return  (+(v5 * a[5] - v4 * a[6] + v3 * a[7])) * a[0]
			+ (-(v5 * a[4] - v2 * a[6] + v1 * a[7])) * a[1]
			+ (+(v4 * a[4] - v2 * a[5] + v0 * a[7])) * a[2]
			+ (-(v3 * a[4] - v1 * a[5] + v0 * a[6])) * a[3];
	}

	bool ToAffine(const Matrix4x4& m, Affine4x3& out)
	{
		if (!IsAffine(m)) return false;
		out = Affine4x3(m[0], m[1], m[2], m[3],
			m[4], m[5], m[6], m[7],
			m[8], m[9], m[10], m[11]);
		return true;
	}

	Matrix4x4 Perspective(double fovYRadians, double aspect,
		double zNear, double zFar, DepthRange range)
	{
		const double f = 1.0 / std::tan(0.5 * fovYRadians);
		Matrix4x4 m;
		for (int32_t i = 0; i < 16; ++i) m[i] = 0.0;
		m[M4Idx(0, 0)] = f / aspect;
		m[M4Idx(1, 1)] = f;
		m[M4Idx(3, 2)] = -1.0;               // w = -z, отсюда перспективное деление
		if (range == DepthRange::NegOneToOne) {
			m[M4Idx(2, 2)] = (zFar + zNear) / (zNear - zFar);
			m[M4Idx(2, 3)] = (2.0 * zFar * zNear) / (zNear - zFar);
		}
		else {
			m[M4Idx(2, 2)] = zFar / (zNear - zFar);
			m[M4Idx(2, 3)] = (zFar * zNear) / (zNear - zFar);
		}
		return m;
	}

	Matrix4x4 Orthographic(double l, double r, double b, double t,
		double zNear, double zFar, DepthRange range)
	{
		Matrix4x4 m;
		for (int32_t i = 0; i < 16; ++i) m[i] = 0.0;
		m[M4Idx(0, 0)] = 2.0 / (r - l);   m[M4Idx(0, 3)] = -(r + l) / (r - l);
		m[M4Idx(1, 1)] = 2.0 / (t - b);   m[M4Idx(1, 3)] = -(t + b) / (t - b);
		m[M4Idx(3, 3)] = 1.0;
		if (range == DepthRange::NegOneToOne) {
			m[M4Idx(2, 2)] = -2.0 / (zFar - zNear);  m[M4Idx(2, 3)] = -(zFar + zNear) / (zFar - zNear);
		}
		else {
			m[M4Idx(2, 2)] = -1.0 / (zFar - zNear);  m[M4Idx(2, 3)] = -zNear / (zFar - zNear);
		}
		return m;
	}

} // namespace mathx
