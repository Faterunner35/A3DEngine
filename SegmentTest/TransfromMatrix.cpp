#include "TransfromMatrix.h"

namespace mathx {

	Affine4x3 MakeInverseTransform(const Vector<3, double>& position,
		const Quaternion& orientation,
		const Vector<3, double>& scale)
	{
		assert(std::abs(scale[0]) > 1e-12 && std::abs(scale[1]) > 1e-12 &&
			std::abs(scale[2]) > 1e-12 && "MakeInverseTransform: нулевой масштаб");

		// Обратный поворот = сопряжение (для единичного кватерниона)
		const Quaternion inv = quat_conjugate(orientation);
		const double x = inv[0], y = inv[1], z = inv[2], w = inv[3];
		const double xx = x * x, yy = y * y, zz = z * z;
		const double xy = x * y, xz = x * z, yz = y * z;
		const double wx = w * x, wy = w * y, wz = w * z;

		const double isx = 1.0 / scale[0], isy = 1.0 / scale[1], isz = 1.0 / scale[2];

		// R⁻¹, затем масштабируем СТРОКИ на 1/s (S⁻¹ слева)
		const double r00 = (1.0 - 2.0 * (yy + zz)) * isx, r01 = (2.0 * (xy - wz)) * isx, r02 = (2.0 * (xz + wy)) * isx;
		const double r10 = (2.0 * (xy + wz)) * isy, r11 = (1.0 - 2.0 * (xx + zz)) * isy, r12 = (2.0 * (yz - wx)) * isy;
		const double r20 = (2.0 * (xz - wy)) * isz, r21 = (2.0 * (yz + wx)) * isz, r22 = (1.0 - 2.0 * (xx + yy)) * isz;

		const double px = position[0], py = position[1], pz = position[2];
		return Affine4x3(
			r00, r01, r02, -(r00 * px + r01 * py + r02 * pz),
			r10, r11, r12, -(r10 * px + r11 * py + r12 * pz),
			r20, r21, r22, -(r20 * px + r21 * py + r22 * pz));
	}

	bool DecomposeTransform(const Affine4x3& m,
		Vector<3, double>& outPosition,
		Quaternion& outOrientation,
		Vector<3, double>& outScale,
		double shearTolerance)
	{
		outPosition = GetTranslation(m);

		// Столбцы линейной части
		double c0[3] = { m[AIdx(0,0)], m[AIdx(1,0)], m[AIdx(2,0)] };
		double c1[3] = { m[AIdx(0,1)], m[AIdx(1,1)], m[AIdx(2,1)] };
		double c2[3] = { m[AIdx(0,2)], m[AIdx(1,2)], m[AIdx(2,2)] };

		double s0 = std::sqrt(c0[0] * c0[0] + c0[1] * c0[1] + c0[2] * c0[2]);
		double s1 = std::sqrt(c1[0] * c1[0] + c1[1] * c1[1] + c1[2] * c1[2]);
		double s2 = std::sqrt(c2[0] * c2[0] + c2[1] * c2[1] + c2[2] * c2[2]);

		if (s0 < 1e-12 || s1 < 1e-12 || s2 < 1e-12) return false;   // вырожденная

		// Нормируем столбцы -> должен получиться чистый поворот
		for (int i = 0; i < 3; ++i) { c0[i] /= s0; c1[i] /= s1; c2[i] /= s2; }

		// ПРОВЕРКА НА СДВИГ: столбцы обязаны быть попарно ортогональны.
		const double d01 = c0[0] * c1[0] + c0[1] * c1[1] + c0[2] * c1[2];
		const double d02 = c0[0] * c2[0] + c0[1] * c2[1] + c0[2] * c2[2];
		const double d12 = c1[0] * c2[0] + c1[1] * c2[1] + c1[2] * c2[2];
		if (std::abs(d01) > shearTolerance ||
			std::abs(d02) > shearTolerance ||
			std::abs(d12) > shearTolerance) {
			return false;   // есть сдвиг — раскладывается не в T·R·S
		}

		// Отражение: det < 0. Знак уходит в первый масштаб.
		const double det = c0[0] * (c1[1] * c2[2] - c1[2] * c2[1])
			- c1[0] * (c0[1] * c2[2] - c0[2] * c2[1])
			+ c2[0] * (c0[1] * c1[2] - c0[2] * c1[1]);
		if (det < 0.0) {
			s0 = -s0;
			for (int i = 0; i < 3; ++i) c0[i] = -c0[i];
		}

		outScale = Vector<3, double>{ s0, s1, s2 };

		// Шумейк по нормированным столбцам (R[r][c]: столбцы это c0,c1,c2)
		const double R00 = c0[0], R10 = c0[1], R20 = c0[2];
		const double R01 = c1[0], R11 = c1[1], R21 = c1[2];
		const double R02 = c2[0], R12 = c2[1], R22 = c2[2];

		Quaternion q;
		const double trace = R00 + R11 + R22;
		if (trace > 0.0) {
			double root = std::sqrt(trace + 1.0);
			q[3] = 0.5 * root; root = 0.5 / root;
			q[0] = (R21 - R12) * root; q[1] = (R02 - R20) * root; q[2] = (R10 - R01) * root;
		}
		else {
			const double d[3] = { R00, R11, R22 };
			const double Rm[9] = { R00,R01,R02, R10,R11,R12, R20,R21,R22 };
			static constexpr int32_t nxt[3] = { 1,2,0 };
			int32_t i = 0; if (d[1] > d[0]) i = 1; if (d[2] > d[i]) i = 2;
			const int32_t j = nxt[i], k = nxt[j];
			double root = std::sqrt(Rm[i * 3 + i] - Rm[j * 3 + j] - Rm[k * 3 + k] + 1.0);
			q[i] = 0.5 * root; root = 0.5 / root;
			q[3] = (Rm[k * 3 + j] - Rm[j * 3 + k]) * root;
			q[j] = (Rm[j * 3 + i] + Rm[i * 3 + j]) * root;
			q[k] = (Rm[k * 3 + i] + Rm[i * 3 + k]) * root;
		}
		outOrientation = q;
		return true;
	}

	std::optional<TransformParts> TryDecomposeTransform(
		const Affine4x3& m, double shearTolerance)
	{
		TransformParts p;
		if (!DecomposeTransform(m, p.position, p.orientation, p.scale, shearTolerance))
			return std::nullopt;
		return p;
	}

	void CreateRotorMatrix(Matrix3x3& mat, double ax, double ay, double az, double angle)
	{
		double s = std::sin(angle), c = std::cos(angle), d = 1.0 - c;
		mat[0] = c + d * (ax * ax);      mat[1] = d * (ax * ay) - s * az;   mat[2] = d * (ax * az) + s * ay;
		mat[3] = d * (ax * ay) + s * az;   mat[4] = c + d * (ay * ay);      mat[5] = d * (ay * az) - s * ax;
		mat[6] = d * (ax * az) - s * ay;   mat[7] = d * (ay * az) + s * ax;   mat[8] = c + d * (az * az);
	}

	Matrix3x3 CreateRotMatFromQuat(const Quaternion& quat)
	{


		Matrix3x3 mat_rot;

		// Все переменные строго double
		double x2 = quat[0] + quat[0];
		double y2 = quat[1] + quat[1];
		double z2 = quat[2] + quat[2];

		// Предвычисляем попарные произведения
		double xx = quat[0] * x2;   double xy = quat[0] * y2;   double xz = quat[0] * z2;
		double yy = quat[1] * y2;   double yz = quat[1] * z2;   double zz = quat[2] * z2;
		double wx = quat[3] * x2;   double wy = quat[3] * y2;   double wz = quat[3] * z2;

		// Главная диагональ
		mat_rot[0] = 1.0 - (yy + zz);
		mat_rot[4] = 1.0 - (xx + zz);
		mat_rot[8] = 1.0 - (xx + yy);

		// Внедиагональные элементы (Row-Major)
		mat_rot[1] = xy - wz;
		mat_rot[2] = xz + wy;
		mat_rot[3] = xy + wz;
		mat_rot[5] = yz - wx;
		mat_rot[6] = xz - wy;
		mat_rot[7] = yz + wx;

		return mat_rot;

	}

	Matrix3x3 MakeReflection(const Vector<3, double>& a)
	{
		assert(std::abs(a[0] * a[0] + a[1] * a[1] + a[2] * a[2] - 1.0) < 1e-9 &&
			"MakeReflection: нормаль должна быть единичной");
		// Форма Ленгьеля: -2a считается ОДИН раз, матрица симметрична,
		// поэтому каждый внедиагональный элемент вычисляется по разу.
		// 9 умножений вместо 12-15 при наивной записи формулы.
		const double x = a[0] * -2.0, y = a[1] * -2.0, z = a[2] * -2.0;
		const double axay = x * a[1], axaz = x * a[2], ayaz = y * a[2];
		return Matrix3x3(x * a[0] + 1.0, axay, axaz,
			axay, y * a[1] + 1.0, ayaz,
			axaz, ayaz, z * a[2] + 1.0);
	}

	Affine4x3 MakeReflectionPlane(const Vector<3, double>& a, double d)
	{
		const Matrix3x3 R = MakeReflection(a);
		return Affine4x3(R[0], R[1], R[2], -2.0 * d * a[0],
			R[3], R[4], R[5], -2.0 * d * a[1],
			R[6], R[7], R[8], -2.0 * d * a[2]);
	}

	Affine4x3 MakeReflectionThroughPoint(const Vector<3, double>& a,
		const Vector<3, double>& point)
	{
		return MakeReflectionPlane(a, -(a[0] * point[0] + a[1] * point[1] + a[2] * point[2]));
	}

	Matrix3x3 MakeInvolution(const Vector<3, double>& a)
	{
		assert(std::abs(a[0] * a[0] + a[1] * a[1] + a[2] * a[2] - 1.0) < 1e-9 &&
			"MakeInvolution: axis must be normalized");

		const double x2 = a[0] * 2.0;
		const double y2 = a[1] * 2.0;
		const double z2 = a[2] * 2.0;

		const double xy2 = x2 * a[1];
		const double xz2 = x2 * a[2];
		const double yz2 = y2 * a[2];

		return Matrix3x3(
			x2 * a[0] - 1.0, xy2, xz2,
			xy2, y2 * a[1] - 1.0, yz2,
			xz2, yz2, z2 * a[2] - 1.0
		);
	}

	Matrix3x3 MakeScaleAlong(const Vector<3, double>& a, double s)
	{
		assert(std::abs(DotProduct(a, a) - 1.0) < 1e-9 &&
			"MakeScaleAlong: направление должно быть единичным");
		const double k = s - 1.0;
		const double x = a[0] * k, y = a[1] * k, z = a[2] * k;   // (s−1)·a, один раз
		const double axay = x * a[1], axaz = x * a[2], ayaz = y * a[2];
		return Matrix3x3(x * a[0] + 1.0, axay, axaz,
			axay, y * a[1] + 1.0, ayaz,
			axaz, ayaz, z * a[2] + 1.0);
	}

	Matrix3x3 MakeScaleAlongAxes(const Vector<3, double>& u1, double s1,
		const Vector<3, double>& u2, double s2,
		const Vector<3, double>& u3, double s3)
	{
		assert(std::abs(DotProduct(u1, u2)) < 1e-9 &&
			std::abs(DotProduct(u1, u3)) < 1e-9 &&
			std::abs(DotProduct(u2, u3)) < 1e-9 &&
			"MakeScaleAlongAxes: оси должны быть взаимно ортогональны");
		Matrix3x3 m;
		for (int32_t r = 0; r < 3; ++r)
			for (int32_t c = 0; c < 3; ++c)
				m[r * 3 + c] = s1 * u1[r] * u1[c] + s2 * u2[r] * u2[c] + s3 * u3[r] * u3[c];
		return m;
	}

	Matrix3x3 MakeSkew(double angleRadians,
		const Vector<3, double>& a,
		const Vector<3, double>& b)
	{
		assert(std::abs(DotProduct(a, a) - 1.0) < 1e-9 &&
			std::abs(DotProduct(b, b) - 1.0) < 1e-9 &&
			"MakeSkew: направления должны быть единичными");
		assert(std::abs(DotProduct(a, b)) < 1e-9 &&
			"MakeSkew: a и b должны быть перпендикулярны");
		assert(std::abs(angleRadians) < 1.5707 &&
			"MakeSkew: угол должен быть строго меньше 90° — tan уходит в бесконечность");

		const double t = std::tan(angleRadians);
		const double ax = a[0] * t, ay = a[1] * t, az = a[2] * t;   // a·tanθ, один раз
		return Matrix3x3(ax * b[0] + 1.0, ax * b[1], ax * b[2],
			ay * b[0], ay * b[1] + 1.0, ay * b[2],
			az * b[0], az * b[1], az * b[2] + 1.0);
	}

	Matrix3x3 MakeSkewByTangent(double tangent,
		const Vector<3, double>& a,
		const Vector<3, double>& b)
	{
		assert(std::abs(DotProduct(a, b)) < 1e-9 && "MakeSkewByTangent: a ⊥ b");
		const double ax = a[0] * tangent, ay = a[1] * tangent, az = a[2] * tangent;
		return Matrix3x3(ax * b[0] + 1.0, ax * b[1], ax * b[2],
			ay * b[0], ay * b[1] + 1.0, ay * b[2],
			az * b[0], az * b[1], az * b[2] + 1.0);
	}

	bool HasNonOrthogonalColumns(const Matrix3x3& m, double tolerance)
	{
		Vector<3, double> c0{ m[0],m[3],m[6] }, c1{ m[1],m[4],m[7] }, c2{ m[2],m[5],m[8] };
		const double l0 = std::sqrt(DotProduct(c0, c0));
		const double l1 = std::sqrt(DotProduct(c1, c1));
		const double l2 = std::sqrt(DotProduct(c2, c2));
		if (l0 < 1e-12 || l1 < 1e-12 || l2 < 1e-12) return false;   // вырождена
		return std::abs(DotProduct(c0, c1)) / (l0 * l1) > tolerance
			|| std::abs(DotProduct(c0, c2)) / (l0 * l2) > tolerance
			|| std::abs(DotProduct(c1, c2)) / (l1 * l2) > tolerance;
	}

	bool HasSkew(const Matrix3x3& m, double tolerance)
	{
		return HasNonOrthogonalColumns(m, tolerance) && !IsSymmetric(m);
	}

} // namespace mathx
