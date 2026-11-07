#include "LinearTransform.h"

namespace mathx {

	double SpectralNorm(const Matrix3x3& m)
	{
		double A[3][3];
		for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
			double s = 0; for (int k = 0; k < 3; ++k) s += m(k, i) * m(k, j);
			A[i][j] = s;
		}
		const double p1 = A[0][1] * A[0][1] + A[0][2] * A[0][2] + A[1][2] * A[1][2];
		if (p1 == 0.0) return std::sqrt(std::max({ A[0][0], A[1][1], A[2][2] }));
		const double q = (A[0][0] + A[1][1] + A[2][2]) / 3.0;
		const double p2 = (A[0][0] - q) * (A[0][0] - q) + (A[1][1] - q) * (A[1][1] - q) + (A[2][2] - q) * (A[2][2] - q) + 2.0 * p1;
		const double p = std::sqrt(p2 / 6.0);
		double B[3][3];
		for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) B[i][j] = (A[i][j] - (i == j ? q : 0.0)) / p;
		double r = 0.5 * (B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1])
			- B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0])
			+ B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0]));
		r = std::clamp(r, -1.0, 1.0);
		return std::sqrt(q + 2.0 * p * std::cos(std::acos(r) / 3.0));
	}

	std::optional<AffineTransform> AffineTransform::TryCreate(const Affine4x3& m, double tolerance)
	{
		Affine4x3 inv;
		if (!mathx::Inverse(m, inv, tolerance)) return std::nullopt;
		return AffineTransform(m, inv);
	}

	std::optional<AffineTransform> AffineTransform::TryCreate(const Matrix3x3& L, const Vector<3, double>& t,
		double tolerance)
	{
		return TryCreate(Affine4x3(L(0, 0), L(0, 1), L(0, 2), t[0],
			L(1, 0), L(1, 1), L(1, 2), t[1],
			L(2, 0), L(2, 1), L(2, 2), t[2]), tolerance);
	}

	std::optional<AffineTransform> AffineTransform::TryCreate(const Matrix3x3& L, double tolerance)
	{
		return TryCreate(L, Vector<3, double>{0, 0, 0}, tolerance);
	}

	std::optional<AffineTransform> AffineTransform::FromTRS(const Vector<3, double>& position,
		const Quaternion& orientation,
		const Vector<3, double>& scale,
		double tolerance)
	{
		if (std::abs(scale[0]) <= tolerance || std::abs(scale[1]) <= tolerance ||
			std::abs(scale[2]) <= tolerance) return std::nullopt;
		return AffineTransform(mathx::MakeTransform(position, orientation, scale),
			MakeInverseTransform(position, orientation, scale));
	}

	bool AffineTransform::IsRigid(double tolerance) const
	{
		for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
			double s = 0; for (int k = 0; k < 3; ++k) s += m_[AIdx(k, i)] * m_[AIdx(k, j)];
			if (std::abs(s - (i == j ? 1.0 : 0.0)) > tolerance) return false;
		}
		return Determinant() > 0.0;
	}

	double AffineTransform::MaxStretch() const
	{
		return SpectralNorm(Matrix3x3(m_[0], m_[1], m_[2], m_[4], m_[5], m_[6], m_[8], m_[9], m_[10]));
	}

} // namespace mathx
