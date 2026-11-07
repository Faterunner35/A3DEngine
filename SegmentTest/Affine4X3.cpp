#include "Affine4X3.h"

namespace mathx {

	bool Inverse(const Affine4x3& m, Affine4x3& out, double tolerance)
	{
		const Vector<3, double> a = m.Column(0);
		const Vector<3, double> b = m.Column(1);
		const Vector<3, double> c = m.Column(2);
		const Vector<3, double> t = m.Column(3);          // столбец переноса

		const Vector<3, double> r0 = Cross(b, c);
		const Vector<3, double> r1 = Cross(c, a);
		const Vector<3, double> r2 = Cross(a, b);

		const double det = DotProduct(r2, c);            // уже посчитан выше

		// Порог ОТНОСИТЕЛЬНЫЙ: |det| ≤ |a||b||c| (неравенство Адамара), и отношение
		// |det| / (|a||b||c|) не зависит от масштаба. Абсолютный порог отвергал бы
		// обратимые матрицы мелких объектов: det растёт как куб масштаба.
		const double bound = std::sqrt(DotProduct(a, a) * DotProduct(b, b) * DotProduct(c, c));
		if (std::abs(det) <= tolerance * bound) return false;
		const double inv = 1.0 / det;

		// Новый перенос: −R⁻¹·t — скалярные произведения строк R⁻¹ на t.
		out = Affine4x3(
			r0[0] * inv, r0[1] * inv, r0[2] * inv, -DotProduct(r0, t) * inv,
			r1[0] * inv, r1[1] * inv, r1[2] * inv, -DotProduct(r1, t) * inv,
			r2[0] * inv, r2[1] * inv, r2[2] * inv, -DotProduct(r2, t) * inv);
		return true;
	}

	std::optional<Affine4x3> TryInverse(const Affine4x3& m, double tolerance)
	{
		Affine4x3 out;
		if (!Inverse(m, out, tolerance)) return std::nullopt;
		return out;
	}

	bool MakeNormalMatrix(const Affine4x3& m, Matrix3x3& out, double tolerance)
	{
		const Vector<3, double> a = m.Column(0), b = m.Column(1), c = m.Column(2);
		const Vector<3, double> k0 = Cross(b, c), k1 = Cross(c, a), k2 = Cross(a, b);
		const double det = DotProduct(a, k0);
		const double bound = std::sqrt(DotProduct(a, a) * DotProduct(b, b) * DotProduct(c, c));
		if (std::abs(det) <= tolerance * bound) return false;
		const double s = det < 0.0 ? -1.0 : 1.0;
		out = Matrix3x3(s * k0[0], s * k1[0], s * k2[0],      // столбцы: k0, k1, k2
			s * k0[1], s * k1[1], s * k2[1],
			s * k0[2], s * k1[2], s * k2[2]);
		return true;
	}

	Affine4x3 InverseRigid(const Affine4x3& m)
	{
		assert(std::abs(Determinant(m) - 1.0) < 1e-6 &&
			"InverseRigid: есть масштаб или сдвиг — используй Inverse()");
		const double tx = m[3], ty = m[7], tz = m[11];
		// Rᵀ плюс перенос -Rᵀ·t
		return Affine4x3(
			m[0], m[4], m[8], -(m[0] * tx + m[4] * ty + m[8] * tz),
			m[1], m[5], m[9], -(m[1] * tx + m[5] * ty + m[9] * tz),
			m[2], m[6], m[10], -(m[2] * tx + m[6] * ty + m[10] * tz));
	}

	bool HasScale(const Affine4x3& m, double tolerance)
	{
		for (int32_t c = 0; c < 3; ++c) {
			const Vector<3, double> col = m.Column(c);
			if (std::abs(DotProduct(col, col) - 1.0) > tolerance) return true;
		}
		return false;
	}

} // namespace mathx
