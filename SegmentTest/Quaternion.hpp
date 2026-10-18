#pragma once

#include"Geomerty.h"
#include"Matrix3x3.hpp"
#include"Math.h"
#include <array>
#include <cstdint>
#include <cmath>
#include <cassert>
#include<bit>
#include<numbers>


// ============================================================================
//  Quaternion — хранение (x, y, z, w), w последний.
//  Конвенция GLM/Unity/GLSL vec4, совместима с GTE.
//
//  ПРИНЦИП РАЗДЕЛЕНИЯ:
//    в классе — только доступ к данным и особые значения;
//    операторы и геометрия — СВОБОДНЫЕ функции снаружи.
//  Причина не стилистическая: operator*(double, Quaternion) физически
//  не может быть методом (левый операнд не объект класса), а бинарные
//  операторы-члены принимают ОДИН параметр, не два.
// ============================================================================
namespace mathx
{ 
class Quaternion
{
public:
	Quaternion() : quaternion{ 0,0,0,1 } {}          // identity, не мусор

	Quaternion(double x, double y, double z, double w)
		: quaternion{ x, y, z, w } {
	}

	double const& operator[](int32_t i) const { return quaternion[i]; }
	double& operator[](int32_t i) { return quaternion[i]; }

	bool operator==(Quaternion const& q) const { return quaternion == q.quaternion; }
	bool operator!=(Quaternion const& q) const { return quaternion != q.quaternion; }
	// operator<, <=, >, >= намеренно отсутствуют: у кватернионов нет
	// естественного линейного порядка. Лексикографическое сравнение
	// array<double,4> скомпилируется, но семантически бессмысленно.

	// ---- Особые кватернионы ----
	static Quaternion Identity() { return Quaternion(0, 0, 0, 1); }
	// Нулевой — НЕ поворот. Только как аккумулятор перед суммированием.
	static Quaternion Zero() { return Quaternion(0, 0, 0, 0); }
	// Чисто мнимые базисные (w=0).
	static Quaternion I() { return Quaternion(1, 0, 0, 0); }
	static Quaternion J() { return Quaternion(0, 1, 0, 0); }
	static Quaternion K() { return Quaternion(0, 0, 1, 0); }

	double x() const { return quaternion[0]; }
	double y() const { return quaternion[1]; }
	double z() const { return quaternion[2]; }
	double w() const { return quaternion[3]; }

private:
	std::array<double, 4> quaternion;
};

// ============================================================================
//  Унарные операторы — ОДИН параметр, снаружи класса.
// ============================================================================
inline Quaternion operator+(Quaternion const& q) { return q; }

inline Quaternion operator-(Quaternion const& q) {
	return Quaternion(-q[0], -q[1], -q[2], -q[3]);
}

// ============================================================================
//  Составные присваивания — ДВА параметра, снаружи класса.
// ============================================================================
inline Quaternion& operator+=(Quaternion& q0, Quaternion const& q1) {
	for (int i = 0; i < 4; ++i) q0[i] += q1[i];
	return q0;
}
inline Quaternion& operator-=(Quaternion& q0, Quaternion const& q1) {
	for (int i = 0; i < 4; ++i) q0[i] -= q1[i];
	return q0;
}
inline Quaternion& operator*=(Quaternion& q, double scalar) {
	for (int i = 0; i < 4; ++i) q[i] *= scalar;
	return q;
}
// Деление на ноль -> обнуление, а не inf/NaN.
inline Quaternion& operator/=(Quaternion& q, double scalar) {
	if (scalar != 0.0) {
		for (int i = 0; i < 4; ++i) q[i] /= scalar;
	}
	else {
		for (int i = 0; i < 4; ++i) q[i] = 0.0;
	}
	return q;
}

// ============================================================================
//  Бинарные операторы
// ============================================================================
inline Quaternion operator+(Quaternion const& q0, Quaternion const& q1) {
	Quaternion r = q0; return r += q1;
}
inline Quaternion operator-(Quaternion const& q0, Quaternion const& q1) {
	Quaternion r = q0; return r -= q1;
}
inline Quaternion operator*(Quaternion const& q, double scalar) {
	Quaternion r = q; return r *= scalar;
}
// НЕ может быть методом: левый операнд — double, не объект класса.
inline Quaternion operator*(double scalar, Quaternion const& q) {
	Quaternion r = q; return r *= scalar;
}
inline Quaternion operator/(Quaternion const& q, double scalar) {
	Quaternion r = q; return r /= scalar;
}

// Hamilton product — композиция поворотов. Некоммутативен.
inline Quaternion operator*(Quaternion const& q0, Quaternion const& q1) {
	return Quaternion(
		+q0[0] * q1[3] + q0[1] * q1[2] - q0[2] * q1[1] + q0[3] * q1[0],
		-q0[0] * q1[2] + q0[1] * q1[3] + q0[2] * q1[0] + q0[3] * q1[1],
		+q0[0] * q1[1] - q0[1] * q1[0] + q0[2] * q1[3] + q0[3] * q1[2],
		-q0[0] * q1[0] - q0[1] * q1[1] - q0[2] * q1[2] + q0[3] * q1[3]
	);
}

// ============================================================================
//  Геометрические операции — по одному определению каждой.
// ============================================================================

// Dot как для Vec4, НЕ Hamilton product.
inline double quat_dot(Quaternion const& q0, Quaternion const& q1) {
	double dot = q0[0] * q1[0];
	for (int i = 1; i < 4; ++i) dot += q0[i] * q1[i];
	return dot;
}

inline double quat_length(Quaternion const& q) {
	return std::sqrt(quat_dot(q, q));
}

// In-place. Возвращает ПРЕДЫДУЩУЮ длину — диагностика дрейфа
// (далеко от 1 => накопилась ошибка). Защищена от деления на ноль.
inline double quat_normalize(Quaternion& q) {
	const double length = quat_length(q);
	if (length > 0.0) {
		q /= length;
	}
	else {
		for (int i = 0; i < 4; ++i) q[i] = 0.0;
	}
	return length;
}

// (x,y,z,w) -> (-x,-y,-z,w). Знак меняет только мнимая часть.
inline Quaternion quat_conjugate(Quaternion const& q) {
	return Quaternion(-q[0], -q[1], -q[2], +q[3]);
}

// conj(q)/|q|². Для ЕДИНИЧНОГО q достаточно quat_conjugate (без деления).
inline Quaternion quat_inverse(Quaternion const& q) {
	const double sqrLen = quat_dot(q, q);
	if (sqrLen > 0.0) return quat_conjugate(q) / sqrLen;
	return Quaternion::Zero();
}

// ============================================================================
//  Slerp
//
//  slerp(t) = [sin((1-t)θ)·q0 + sin(tθ)·q1] / sinθ,  θ = acos(q0·q1)
//
//  Доказанные свойства: концы; |slerp|=1; угол(slerp(t),q0)=t·θ;
//  геодезическая на S³ (∇_Ṙ Ṙ = 0, ускорение чисто радиальное).
//
//  Вырожденные случаи принципиально разные:
//    θ→0 : УСТРАНИМАЯ особенность — по Лопиталю slerp аналитически
//          продолжается в lerp (точный предел, не аппроксимация).
//    θ=π : НЕУСТРАНИМАЯ неоднозначность — антиподальные точки на S³
//          соединяет континуум геодезических. Лечится знаковой
//          коррекцией ДО вызова, не эпсилоном.
// ============================================================================

inline constexpr double SLERP_COS_EPSILON = 1e-12;

// Коэффициенты (a,b): slerp = a*q0 + b*q1.
// ВСЕ пути возвращают значение (в GTE ветка θ=π проваливается — UB).
inline std::array<double, 2> SlerpWeights(double t, double cosAngle)
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

// shortestPath: q и -q задают одно вращение (Ad: Spin(3)->SO(3) двулистно),
// из углов θ и π-θ берём меньший. Без коррекции — поворот "длинным путём".
inline Quaternion quat_slerp(double t, Quaternion const& q0, Quaternion const& q1,
	bool shortestPath = true)
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

// Дешевле slerp (нет тригонометрии), коммутативен, но скорость непостоянна.
// Годится для блендинга весов в skinning, где важна коммутативность.
inline Quaternion quat_nlerp(double t, Quaternion const& q0, Quaternion const& q1,
	bool shortestPath = true)
{
	double sign = 1.0;
	if (shortestPath && quat_dot(q0, q1) < 0.0) sign = -1.0;

	Quaternion result = q0 * (1.0 - t) + q1 * (sign * t);

	// Та же вырожденность: антиподальные входы при shortestPath=false дают
	// нулевой кватернион при t=0.5. quat_normalize его обнулит — не NaN,
	// но и не поворот. Возвращаем q0.
	if (quat_length(result) < 1e-12) {
		assert(false && "quat_nlerp: degenerate result (antipodal inputs "
			"with shortestPath=false)");
		return q0;
	}
	quat_normalize(result);
	return result;
}


inline Quaternion FromAngleAxis(const Radian& angle, const Vector<3, double>& axis)
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

	inline Radian ToAngleAxis(const Quaternion & q, Vector<3, double>&outAxis)
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

	inline Matrix3x3 QuaternionToMatrix(const Quaternion& q)
	{
		const double x = q[0], y = q[1], z = q[2], w = q[3];

		const double xx = x * x, yy = y * y, zz = z * z;
		const double xy = x * y, xz = x * z, yz = y * z;
		const double wx = w * x, wy = w * y, wz = w * z;

		return Matrix3x3(
			1.0 - 2.0 * (yy + zz), 2.0 * (xy - wz), 2.0 * (xz + wy),
			2.0 * (xy + wz), 1.0 - 2.0 * (xx + zz), 2.0 * (yz - wx),
			2.0 * (xz - wy), 2.0 * (yz + wx), 1.0 - 2.0 * (xx + yy)
		);
	}

	//Матрица->кватернион.Алгоритм Шумейка(Shoemake, SIGGRAPH 1987).
		//
		//  ЗАЧЕМ ВЕТВЛЕНИЕ. Наивно: tr(R) = 4w² - 1, откуда w = sqrt(tr+1)/2,
		//  затем x = (R21-R12)/(4w) и т.д. Но w = cos(θ/2) обнуляется при
		//  повороте на 180° — деление на ноль, а рядом с ним катастрофическая
		//  потеря точности. Условие tr > 0 равносильно |w| > 1/2, то есть
		//  отсекает опасную зону.
		//
		//  ИДЕЯ ШУМЕЙКА. Из x²+y²+z²+w² = 1 хотя бы одна компонента по модулю
		//  >= 1/2 (иначе сумма квадратов была бы < 1). Значит если w мал, велика
		//  какая-то из x,y,z — её и берём знаменателем.
		//
		//  КАК НАЙТИ НАИБОЛЬШУЮ, не зная кватерниона: из структуры матрицы
		//      R00 - R11 - R22 = 4x² - 1   (и циклически для y, z)
		//  то есть наибольший диагональный элемент указывает на наибольшую
		//  компоненту — достаточно сравнить диагональ.
		//
		//  ЦИКЛИЧЕСКАЯ ПЕРЕСТАНОВКА next = {1,2,0}: три случая (наибольший x/y/z)
		//  структурно одинаковы с точностью до сдвига i->j->k->i. Это следствие
		//  циклической симметрии кватернионной алгебры (ij=k, jk=i, ki=j),
		//  поэтому именно КРУГОВАЯ перестановка сохраняет знаки — обычная
		//  перестановка двух индексов дала бы отражение вместо поворота.
		//
		//  ПРЕДУСЛОВИЕ: R ортогональна, det = +1 (истинный поворот, не отражение).
		// ---------------------------------------------------------------------------
		inline Quaternion MatrixToQuaternion(const Matrix3x3& R)
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

	inline Vector<3, double> AxisX(const Quaternion& q)
	{
		const double x = q[0], y = q[1], z = q[2], w = q[3];
		const double ty = 2.0 * y, tz = 2.0 * z;
		const double twy = ty * w, twz = tz * w;
		const double txy = ty * x, txz = tz * x;
		const double tyy = ty * y, tzz = tz * z;
		return Vector<3, double>{ 1.0 - (tyy + tzz), txy + twz, txz - twy };
	}

	// Локальная ось Y = столбец 1.
	inline Vector<3, double> AxisY(const Quaternion& q)
	{
		const double x = q[0], y = q[1], z = q[2], w = q[3];
		const double tx = 2.0 * x, ty = 2.0 * y, tz = 2.0 * z;
		const double twx = tx * w, twz = tz * w;
		const double txx = tx * x, txy = ty * x;
		const double tyz = tz * y, tzz = tz * z;
		return Vector<3, double>{ txy - twz, 1.0 - (txx + tzz), tyz + twx };
	}

	// Локальная ось Z = столбец 2.
	inline Vector<3, double> AxisZ(const Quaternion& q)
	{
		const double x = q[0], y = q[1], z = q[2], w = q[3];
		const double tx = 2.0 * x, ty = 2.0 * y, tz = 2.0 * z;
		const double twx = tx * w, twy = ty * w;
		const double txx = tx * x, txz = tz * x;
		const double tyy = ty * y, tyz = tz * y;
		return Vector<3, double>{ txz + twy, tyz - twx, 1.0 - (txx + tyy) };
	}

	inline constexpr double kOrthonormalTolerance = 1e-9;
	
	inline bool AxesAreOrthonormal(Vector<3, double> const& x,
		Vector<3, double> const& y,
		Vector<3, double> const& z,
		double tol = kOrthonormalTolerance)
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

	inline Quaternion FromAxes(const Vector<3, double>& xAxis,
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


	inline Quaternion QuaternionFromTwoVectors(Vector<3, double> const& from,
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


	// Все три оси сразу — через матрицу, чтобы общие подвыражения считались один раз.
	inline std::array<Vector<3, double>, 3> ToAxes(const Quaternion& q)
	{
		const Matrix3x3 R=QuaternionToMatrix(q);
		return {
			Vector<3,double>{ R(0,0), R(1,0), R(2,0) },   // столбец 0
			Vector<3,double>{ R(0,1), R(1,1), R(2,1) },   // столбец 1
			Vector<3,double>{ R(0,2), R(1,2), R(2,2) }    // столбец 2
		};
	}


	inline Quaternion quat_slerp_extra_spins(double t,
		Quaternion const& q0,
		Quaternion const& q1,
		int extraSpins,
		bool shortestPath = true)
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


	// ============================================================================
//  Мост в алгебру Ли: log / exp и то, что на них строится.
//
//  Ротор = cos(θ/2) + sin(θ/2)·B̂  (скаляр в grade-0, бивектор в grade-2).
//  log извлекает бивекторную часть, масштабированную полууглом: (θ/2)·n̂.
//  exp — обратно. Это ФУНДАМЕНТ для физики (угловая скорость), SLAM
//  (разность поз на многообразии) и интерполяции — не анимационная деталь.
// ============================================================================

// log единичного кватерниона -> чисто мнимый (бивектор) (θ/2)·n̂, w=0.
// atan2(|imag|, w), НЕ acos(w): |imag| и w — независимые компоненты полной
// точности, поэтому atan2 строго лучше (нулевая ошибка на малых углах).
	inline Quaternion quat_log(Quaternion const& q)
	{
		const double v = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
		if (v < 1e-14) return Quaternion(0, 0, 0, 0);   // тождественный поворот
		const double a = std::atan2(v, q[3]) / v;
		return Quaternion(q[0] * a, q[1] * a, q[2] * a, 0.0);
	}

	// exp чисто мнимого кватерниона (бивектора) -> единичный. Обратно к quat_log.
	inline Quaternion quat_exp(Quaternion const& q)
	{
		const double a = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
		if (a < 1e-14) return Quaternion::Identity();
		const double s = std::sin(a) / a;
		return Quaternion(q[0] * s, q[1] * s, q[2] * s, std::cos(a));
	}

	// ---- Физика: угловая скорость ----

	// Угловая скорость (рад/с) как ВЕКТОР между двумя ориентациями за dt.
	// ω = 2·log(q_next · q_now⁻¹) / dt. Множитель 2 — потому что log даёт
	// полуугол (θ/2), а физическая скорость измеряется полным углом θ.
	// Направление вектора = ось вращения, длина = скорость. Готово для
	// гироскопических моментов, демпфирования, интегрирования инерции.
	inline Quaternion AngularVelocity(Quaternion const& qNow,
		Quaternion const& qNext,
		double dt)
	{
		Quaternion dq = qNext * quat_inverse(qNow);
		if (dq[3] < 0.0) for (int i = 0; i < 4; ++i) dq[i] = -dq[i];  // короткий путь
		const Quaternion b = quat_log(dq);
		const double k = 2.0 / dt;
		return Quaternion(b[0] * k, b[1] * k, b[2] * k, 0.0);   // (ωx, ωy, ωz, 0)
	}

	// Проинтегрировать ориентацию на шаг dt при угловой скорости omega
	// (вектор в мнимой части). q_next = exp(½·ω·dt) · q. Обратно к AngularVelocity.
	inline Quaternion IntegrateAngularVelocity(Quaternion const& q,
		Quaternion const& omega,
		double dt)
	{
		const double h = 0.5 * dt;
		const Quaternion halfStep(omega[0] * h, omega[1] * h, omega[2] * h, 0.0);
		return quat_exp(halfStep) * q;
	}

	// ---- Ось-угол как ОБЁРТКА над log ----
	//
	// Тот же (ось, угол), что log, но РАСПАКОВАННЫЙ: угол отдельно, ось отдельно.
	// Нужно, когда их обрабатывают ПОРОЗНЬ — например, clamp: обрезать угол,
	// ось оставить. Когда нужен единый вектор (физика) — бери quat_log напрямую.
	// outAxis пишется как единичная ось; при θ≈0 ось не определена -> (1,0,0).
	inline double ToAngleAxis(Quaternion const& q, Quaternion& outAxis)
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

	// Ограничить поворот максимальным углом (clamp). Демонстрирует, ЗАЧЕМ нужна
	// распакованная форма: угол трогаем, ось сохраняем. Сустав/камера/стабилизатор.
	inline Quaternion ClampAngle(Quaternion const& q, double maxAngle)
	{
		Quaternion axis;
		const double angle = ToAngleAxis(q, axis);
		if (angle <= maxAngle) return q;
		const double h = 0.5 * maxAngle;
		const double s = std::sin(h);
		return Quaternion(axis[0] * s, axis[1] * s, axis[2] * s, std::cos(h));
	}



}



