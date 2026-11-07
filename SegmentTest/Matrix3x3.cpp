#include "Matrix3x3.h"

namespace mathx
{

bool HasScale(const Matrix3x3& m, double tolerance) {
	for (int32_t c = 0; c < 3; ++c) {
		const double len2 = m(0, c) * m(0, c)
			+ m(1, c) * m(1, c)
			+ m(2, c) * m(2, c);
		if (std::abs(len2 - 1.0) > tolerance) return true;
	}
	return false;
}

bool IsOrthonormal(const Matrix3x3& m, double tolerance) {
	if (HasScale(m, tolerance)) return false;
	// Попарные скалярные произведения столбцов должны быть нулевыми
	for (int32_t a = 0; a < 3; ++a)
		for (int32_t b = a + 1; b < 3; ++b) {
			double d = 0.0;
			for (int32_t r = 0; r < 3; ++r) d += m(r, a) * m(r, b);
			if (std::abs(d) > tolerance) return false;
		}
	// det = +1, а не -1: иначе это отражение, а не поворот
	return mathx::Determinant(m) > 0.0;
}

bool Inverse(const Matrix3x3& m, Matrix3x3& out, double tolerance) {
	const Matrix3x3 adj = Adjugate(m);
	// det через уже посчитанные дополнения первой строки: adj транспонирована,
	// поэтому C00,C01,C02 лежат в adj[0], adj[3], adj[6].
	const double det = m[0] * adj[0] + m[1] * adj[3] + m[2] * adj[6];

	// Порог относительный (неравенство Адамара): |det| ≤ |c0||c1||c2|
	double n0 = 0, n1 = 0, n2 = 0;
	for (int r = 0; r < 3; ++r) { n0 += m[r * 3] * m[r * 3]; n1 += m[r * 3 + 1] * m[r * 3 + 1]; n2 += m[r * 3 + 2] * m[r * 3 + 2]; }
	if (std::abs(det) <= tolerance * std::sqrt(n0 * n1 * n2)) return false;

	const double invDet = 1.0 / det;
	for (uint16_t i = 0; i < 9; ++i) out[i] = adj[i] * invDet;
	return true;
}

void Orthonormalize(Matrix3x3& m) {
	// Столбец 0: нормируем
	double l0 = std::sqrt(m[0] * m[0] + m[3] * m[3] + m[6] * m[6]);
	m[0] /= l0; m[3] /= l0; m[6] /= l0;

	// Столбец 1: вычитаем проекцию на 0, нормируем
	double d01 = m[0] * m[1] + m[3] * m[4] + m[6] * m[7];
	m[1] -= d01 * m[0]; m[4] -= d01 * m[3]; m[7] -= d01 * m[6];
	double l1 = std::sqrt(m[1] * m[1] + m[4] * m[4] + m[7] * m[7]);
	m[1] /= l1; m[4] /= l1; m[7] /= l1;

	// Столбец 2 = столбец0 x столбец1 — сразу единичный и ортогональный обоим,
	// считать проекции не нужно.
	m[2] = m[3] * m[7] - m[6] * m[4];
	m[5] = m[6] * m[1] - m[0] * m[7];
	m[8] = m[0] * m[4] - m[3] * m[1];
}

Matrix3x3  CreateRotationFromAxis(const Vector<3, double>& vec, double angle)
{    
    Matrix3x3 mat;
	double s = std::sin(angle), c = std::cos(angle), d = 1.0 - c;
	mat[0] = c + d * (vec[0] * vec[0]);      mat[1] = d * (vec[0] * vec[1]) - s * vec[2];   mat[2] = d * (vec[0] * vec[2]) + s * vec[1];
	mat[3] = d * (vec[0] * vec[1]) + s * vec[2];   mat[4] = c + d * (vec[1] * vec[1]);      mat[5] = d * (vec[1] * vec[2]) - s * vec[0];
	mat[6] = d * (vec[0] * vec[2]) - s * vec[1];   mat[7] = d * (vec[1] * vec[2]) + s * vec[0];   mat[8] = c + d * (vec[2] * vec[2]);
	return mat;
}

} // namespace mathx
