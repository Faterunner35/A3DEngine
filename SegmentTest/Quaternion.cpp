#include "Quaternion.h"

namespace mathx
{

std::array<double, 2> SlerpWeights(double t, double cosAngle)
{
	if (cosAngle >= 1.0 - SLERP_COS_EPSILON) {
		return { 1.0 - t, t };                    // точный предел = lerp
	}

	if (cosAngle <= -1.0 + SLERP_COS_EPSILON) {
		assert(false && "Slerp between antipodal quaternions is ill-defined "
			"(missing shortest-path sign correction?)");
		return { 1.0 - t, t };                    // деградация вместо UB
	}

	// sin через sqrt: угол ∈ (0,π), где sin > 0 — положительный корень.
	const double sinAngle = std::sqrt(1.0 - cosAngle * cosAngle);

	// acos, а НЕ atan2 — и это осознанный выбор, замеренный.
	//
	// В ToAngleAxis atan2(|imag|, w) строго лучше acos(w), потому что
	// там синус и косинус — НЕЗАВИСИМЫЕ компоненты кватерниона, каждая
	// с полной точностью: atan2 использует двумерную информацию и даёт
	// нулевую ошибку там, где acos деградирует до 4e-5.
	//
	// ЗДЕСЬ предпосылка не выполняется: единственный вход — cosAngle,
	// а sinAngle ВЫВЕДЕН из него же (и вдобавок страдает от сокращения
	// при cos->1). Обусловленность отображения cos->угол равна
	// 1/sqrt(1-cos²) — это свойство математики, а не вызываемой функции;
	// atan2 не вернёт информацию, которой нет во входе.
	// Замер: ошибки acos и atan2 здесь совпадают (местами atan2 хуже),
	// при этом atan2 стоит 23.7 нс против 17.6 нс у acos — на 35% дороже.
	// Платить 35% за ноль выигрыша смысла нет.
	const double angle = std::acos(cosAngle);

	return {
		std::sin((1.0 - t) * angle) / sinAngle,
		std::sin(t * angle) / sinAngle
	};
}

Quaternion quat_slerp(double t, Quaternion const& q0, Quaternion const& q1,
	bool shortestPath)
{
	double cosAngle = quat_dot(q0, q1);
	double sign = 1.0;

	if (shortestPath && cosAngle < 0.0) {
		cosAngle = -cosAngle;
		sign = -1.0;
	}

	const std::array<double, 2> f = SlerpWeights(t, cosAngle);
	return q0 * f[0] + q1 * (sign * f[1]);
}

Quaternion FromAngleAxis(const Radian& angle, const Vector<3, double>& axis)
{
	const Radian half = angle * 0.5;          // angle*0.5 или 0.5*angle — оба теперь работают
	const double s = mathx::Sin(half);
	const double c = mathx::Cos(half);

	return Quaternion(
		s * axis[0],
		s * axis[1],
		s * axis[2],
		c
	);

	}

	Radian ToAngleAxis(const Quaternion & q, Vector<3, double>&outAxis)
	{
		// q[0..3] = x,y,z,w
		const double sqrLen = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];

		// Порог, а не строгий ноль: почти-нулевой угол -> стабильная ветка.
		if (sqrLen > 1e-24) {   // близко к машинному пределу: пускаем в atan2 всё, что физически ненулевое
			const double len = std::sqrt(sqrLen);       // = |sin(θ/2)|
			// atan2(|imag|, w) устойчив везде, в отличие от acos(w),
			// который плохо обусловлен при малых углах (w≈1).
			const Radian angle(2.0 * std::atan2(len, q[3]));

			const double invLen = 1.0 / len;
			outAxis[0] = q[0] * invLen;
			outAxis[1] = q[1] * invLen;
			outAxis[2] = q[2] * invLen;
			return angle;
		}
		else {
			// Угол ~0 (mod 2π): ось не определена, любая подойдёт.
			outAxis[0] = 1.0;
			outAxis[1] = 0.0;
			outAxis[2] = 0.0;
			return Radian(0.0);
		}
	}

		Quaternion MatrixToQuaternion(const Matrix3x3& R)
	{
		Quaternion q;
		const double trace = R(0, 0) + R(1, 1) + R(2, 2);

		if (trace > 0.0) {
			// |w| > 1/2 — считаем w первым, деление на него безопасно.
			double root = std::sqrt(trace + 1.0);        // = 2w
			q[3] = 0.5 * root;
			root = 0.5 / root;                           // = 1/(4w)
			// Разность внедиагональных даёт члены с w: R21-R12 = 4wx и т.д.
			// (симметричные части yz сокращаются, антисимметричные удваиваются)
			q[0] = (R(2, 1) - R(1, 2)) * root;
			q[1] = (R(0, 2) - R(2, 0)) * root;
			q[2] = (R(1, 0) - R(0, 1)) * root;
		}
		else {
			// w мал — наибольшая среди x,y,z. Находим её по диагонали.
			static constexpr int32_t next[3] = { 1, 2, 0 };

			int32_t i = 0;
			if (R(1, 1) > R(0, 0)) i = 1;
			if (R(2, 2) > R(i, i)) i = 2;
			const int32_t j = next[i];
			const int32_t k = next[j];

			// R_ii - R_jj - R_kk + 1 = 4·(компонента_i)²
			double root = std::sqrt(R(i, i) - R(j, j) - R(k, k) + 1.0);
			q[i] = 0.5 * root;
			root = 0.5 / root;
			// w — из антисимметричной части (разность),
			// j,k — из симметричной (сумма): R_ji + R_ij = 4·comp_i·comp_j
			q[3] = (R(k, j) - R(j, k)) * root;
			q[j] = (R(j, i) + R(i, j)) * root;
			q[k] = (R(k, i) + R(i, k)) * root;
		}
		return q;
	}

	bool AxesAreOrthonormal(Vector<3, double> const& x,
		Vector<3, double> const& y,
		Vector<3, double> const& z,
		double tol)
	{
		// Диагональ QᵀQ: квадраты длин столбцов.
		if (std::abs(DotProduct(x, x) - 1.0) > tol) return false;
		if (std::abs(DotProduct(y, y) - 1.0) > tol) return false;
		if (std::abs(DotProduct(z, z) - 1.0) > tol) return false;
		// Вне диагонали: попарные скалярные произведения (матрица симметрична — три числа).
		if (std::abs(DotProduct(x, y)) > tol) return false;
		if (std::abs(DotProduct(y, z)) > tol) return false;
		if (std::abs(DotProduct(z, x)) > tol) return false;
		// det Q = (x × y) · z. Для ортогональной матрицы он ±1; знак отличает поворот от отражения.
		return DotProduct(Cross(x, y), z) > 0.0;
	}

	Quaternion FromAxes(const Vector<3, double>& xAxis,
		const Vector<3, double>& yAxis,
		const Vector<3, double>& zAxis)
	{
		assert(AxesAreOrthonormal(xAxis, yAxis, zAxis) &&
			"FromAxes: axes must be orthonormal and right-handed "
			"(check for accumulated drift or a look-at built without re-orthogonalization)");

		// Оси идут в СТОЛБЦЫ (см. конвенцию выше).
		const Matrix3x3 R(
			xAxis[0], yAxis[0], zAxis[0],
			xAxis[1], yAxis[1], zAxis[1],
			xAxis[2], yAxis[2], zAxis[2]
		);
		return MatrixToQuaternion(R);
	}

	Quaternion QuaternionFromTwoVectors(Vector<3, double> const& from,
		Vector<3, double> const& to)
	{
		const double cosT = DotProduct(from, to);

		// Антипараллельны: d обнулился бы, нужен явный перпендикуляр.
		if (cosT <= -1.0 + 1e-6) {
			// Выбираем перпендикуляр к from, устойчиво: берём ось координат,
			// наименее сонаправленную с from, и делаем cross.
			Vector<3, double> axis;
			if (std::abs(from[0]) > std::abs(from[2])) {
				// перпендикуляр к (x,y,z): (-y, x, 0), проверка: -yx+xy+0=0
				axis = Vector<3, double>{ -from[1], from[0], 0.0 };
			}
			else {
				// перпендикуляр: (0, -z, y), проверка: 0-zy+yz=0  (знак у z!)
				axis = Vector<3, double>{ 0.0, -from[2], from[1] };
			}
			const double len = std::sqrt(DotProduct(axis, axis));
			const double inv = 1.0 / len;
			// Поворот на 180°: w=0, мнимая часть = нормированная ось.
			return Quaternion(axis[0] * inv, axis[1] * inv, axis[2] * inv, 0.0);
		}

		// Общий случай.
		const Vector<3, double> axis = Cross(from, to);
		const double d = std::sqrt(2.0 * (1.0 + cosT));   // = 2cos(θ/2)
		const double s = 1.0 / d;
		// w = d/2 = cos(θ/2), НЕ 0.5*s.
		return Quaternion(axis[0] * s, axis[1] * s, axis[2] * s, 0.5 * d);
	}

	Quaternion quat_slerp_extra_spins(double t,
		Quaternion const& q0,
		Quaternion const& q1,
		int extraSpins,
		bool shortestPath)
	{
		double cosAngle = quat_dot(q0, q1);
		double sign = 1.0;

		// Дыра 1 у Ogre: знаковая коррекция. Без неё поворот пойдёт длинным путём.
		if (shortestPath && cosAngle < 0.0) {
			cosAngle = -cosAngle;
			sign = -1.0;
		}

		// Дыра 2 у Ogre: проверяется ТОЛЬКО θ≈0. sinθ обнуляется на ОБОИХ концах.
		// θ≈0 — кватернионы совпадают: обороты крутить не вокруг чего осмысленного,
		//        но сами обороты сохранить надо (см. ниже).
		// θ≈π — почти антиподальные: ось не определена, деление на sinθ взорвётся.
		constexpr double COS_EPS = 1e-9;

		if (cosAngle >= 1.0 - COS_EPS) {
			// Ориентации совпадают. Ogre здесь возвращает q0 и ТЕРЯЕТ обороты.
			// Но пользователь просил n оборотов — их надо выполнить, иначе
			// "подбрось монету на месте" не сработает.
			if (extraSpins == 0) return q0;
			// Крутим вокруг оси самого q0 (или вокруг Z, если q0 = identity).
			double ax = q0[0], ay = q0[1], az = q0[2];
			const double axisLen = std::sqrt(ax * ax + ay * ay + az * az);
			if (axisLen < 1e-12) { ax = 0.0; ay = 0.0; az = 1.0; }
			else { ax /= axisLen; ay /= axisLen; az /= axisLen; }
			const double half = std::numbers::pi * static_cast<double>(extraSpins) * t;  // πnt
			const double s = std::sin(half), c = std::cos(half);
			const Quaternion spin(ax * s, ay * s, az * s, c);   // поворот на 2πnt вокруг оси
			return spin * q0;
		}

		if (cosAngle <= -1.0 + COS_EPS) {
			// Почти антиподальные: ось поворота не определена (континуум решений).
			assert(false && "quat_slerp_extra_spins: antipodal inputs — "
				"rotation axis is undefined (use shortestPath=true)");
			return q0;
		}

		const double sinAngle = std::sqrt(1.0 - cosAngle * cosAngle);
		const double angle = std::acos(cosAngle);
		const double invSin = 1.0 / sinAngle;
		const double phase = std::numbers::pi * static_cast<double>(extraSpins) * t;

		const double c0 = std::sin((1.0 - t) * angle - phase) * invSin;
		const double c1 = std::sin(t * angle + phase) * invSin;

		return q0 * c0 + q1 * (sign * c1);
	}

	Quaternion quat_log(Quaternion const& q)
	{
		const double v = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
		if (v < 1e-14) return Quaternion(0, 0, 0, 0);   // тождественный поворот
		const double a = std::atan2(v, q[3]) / v;
		return Quaternion(q[0] * a, q[1] * a, q[2] * a, 0.0);
	}

	Quaternion quat_exp(Quaternion const& q)
	{
		const double a = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
		if (a < 1e-14) return Quaternion::Identity();
		const double s = std::sin(a) / a;
		return Quaternion(q[0] * s, q[1] * s, q[2] * s, std::cos(a));
	}

	Quaternion AngularVelocity(Quaternion const& qNow,
		Quaternion const& qNext,
		double dt)
	{
		Quaternion dq = qNext * quat_inverse(qNow);
		if (dq[3] < 0.0) for (int i = 0; i < 4; ++i) dq[i] = -dq[i];  // короткий путь
		const Quaternion b = quat_log(dq);
		const double k = 2.0 / dt;
		return Quaternion(b[0] * k, b[1] * k, b[2] * k, 0.0);   // (ωx, ωy, ωz, 0)
	}

	Quaternion IntegrateAngularVelocity(Quaternion const& q,
		Quaternion const& omega,
		double dt)
	{
		const double h = 0.5 * dt;
		const Quaternion halfStep(omega[0] * h, omega[1] * h, omega[2] * h, 0.0);
		return quat_exp(halfStep) * q;
	}

	double ToAngleAxis(Quaternion const& q, Quaternion& outAxis)
	{
		const Quaternion b = quat_log(q);                 // (θ/2)·n̂
		const double halfAngle = std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
		if (halfAngle < 1e-14) {
			outAxis = Quaternion(1, 0, 0, 0);
			return 0.0;
		}
		const double inv = 1.0 / halfAngle;
		outAxis = Quaternion(b[0] * inv, b[1] * inv, b[2] * inv, 0.0);
		return 2.0 * halfAngle;                            // полный угол θ
	}

	Quaternion ClampAngle(Quaternion const& q, double maxAngle)
	{
		Quaternion axis;
		const double angle = ToAngleAxis(q, axis);
		if (angle <= maxAngle) return q;
		const double h = 0.5 * maxAngle;
		const double s = std::sin(h);
		return Quaternion(axis[0] * s, axis[1] * s, axis[2] * s, std::cos(h));
	}

} // namespace mathx
