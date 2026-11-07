#pragma once

#include"Geomerty.h"
#include"Matrix3x3.h"
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
std::array<double, 2> SlerpWeights(double t, double cosAngle);

// shortestPath: q и -q задают одно вращение (Ad: Spin(3)->SO(3) двулистно),
// из углов θ и π-θ берём меньший. Без коррекции — поворот "длинным путём".
Quaternion quat_slerp(double t, Quaternion const& q0, Quaternion const& q1,
	bool shortestPath = true);

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


Quaternion FromAngleAxis(const Radian& angle, const Vector<3, double>& axis);

	Radian ToAngleAxis(const Quaternion & q, Vector<3, double>&outAxis);

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
		Quaternion MatrixToQuaternion(const Matrix3x3& R);

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
	
	bool AxesAreOrthonormal(Vector<3, double> const& x,
		Vector<3, double> const& y,
		Vector<3, double> const& z,
		double tol = kOrthonormalTolerance);

	Quaternion FromAxes(const Vector<3, double>& xAxis,
		const Vector<3, double>& yAxis,
		const Vector<3, double>& zAxis);


	Quaternion QuaternionFromTwoVectors(Vector<3, double> const& from,
		Vector<3, double> const& to);


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


	Quaternion quat_slerp_extra_spins(double t,
		Quaternion const& q0,
		Quaternion const& q1,
		int extraSpins,
		bool shortestPath = true);


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
	Quaternion quat_log(Quaternion const& q);

	// exp чисто мнимого кватерниона (бивектора) -> единичный. Обратно к quat_log.
	Quaternion quat_exp(Quaternion const& q);

	// ---- Физика: угловая скорость ----

	// Угловая скорость (рад/с) как ВЕКТОР между двумя ориентациями за dt.
	// ω = 2·log(q_next · q_now⁻¹) / dt. Множитель 2 — потому что log даёт
	// полуугол (θ/2), а физическая скорость измеряется полным углом θ.
	// Направление вектора = ось вращения, длина = скорость. Готово для
	// гироскопических моментов, демпфирования, интегрирования инерции.
	Quaternion AngularVelocity(Quaternion const& qNow,
		Quaternion const& qNext,
		double dt);

	// Проинтегрировать ориентацию на шаг dt при угловой скорости omega
	// (вектор в мнимой части). q_next = exp(½·ω·dt) · q. Обратно к AngularVelocity.
	Quaternion IntegrateAngularVelocity(Quaternion const& q,
		Quaternion const& omega,
		double dt);

	// ---- Ось-угол как ОБЁРТКА над log ----
	//
	// Тот же (ось, угол), что log, но РАСПАКОВАННЫЙ: угол отдельно, ось отдельно.
	// Нужно, когда их обрабатывают ПОРОЗНЬ — например, clamp: обрезать угол,
	// ось оставить. Когда нужен единый вектор (физика) — бери quat_log напрямую.
	// outAxis пишется как единичная ось; при θ≈0 ось не определена -> (1,0,0).
	double ToAngleAxis(Quaternion const& q, Quaternion& outAxis);

	// Ограничить поворот максимальным углом (clamp). Демонстрирует, ЗАЧЕМ нужна
	// распакованная форма: угол трогаем, ось сохраняем. Сустав/камера/стабилизатор.
	Quaternion ClampAngle(Quaternion const& q, double maxAngle);



}



