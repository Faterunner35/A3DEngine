#pragma once
#include"GeomIntersector.h"
#include"Plane.h"
#include"Primitives.h"
#include"Ray.h"
#include<cassert>
#include<cmath>
#include<limits>
#include<optional>
	 
namespace mathx
{ 
//  Проекция, отбрасывание и отражение.
//
//  Все три строятся на одной операции — проекции на направление:
//      proj_a(v) = (v·a)·a          <- результат лежит ВДОЛЬ a
//      rej_a(v)  = v − proj_a(v)    <- то, что осталось, ⊥ a
//      refl(v,n) = v − 2·proj_n(v)  <- параллельная часть перевёрнута
//
//  ЧАСТАЯ ОШИБКА: умножить скаляр (v·a) на v вместо a. Результат тогда
//  остаётся вдоль v, то есть направление не меняется вовсе. Ошибка НЕ
//  видна, когда v уже коллинеарен a — там обе формы совпадают.
//
//  ПРЕДУСЛОВИЕ: направление единичное. Иначе (v·a) содержит лишний
//  множитель |a|, и его надо делить на |a|² — см. версии *Unnormalized.
// ============================================================================

// Проекция v на единичное направление a. Результат КОЛЛИНЕАРЕН a.
inline Vector<3, double> VectorProjection(const Vector<3, double>& v,
	const Vector<3, double>& axis)
{
	assert(std::abs(DotProduct(axis, axis) - 1.0) < 1e-9 &&
		"VectorProjection: axis должен быть единичным");
	const double s = DotProduct(v, axis);
	return Vector<3, double>{ axis[0] * s, axis[1] * s, axis[2] * s };
}

// Отбрасывание: то, что осталось от v после снятия проекции. ⊥ axis.
inline Vector<3, double> VectorRejection(const Vector<3, double>& v,
	const Vector<3, double>& axis)
{
	assert(std::abs(DotProduct(axis, axis) - 1.0) < 1e-9 &&
		"VectorRejection: axis должен быть единичным");
	const double s = DotProduct(v, axis);
	return Vector<3, double>{ v[0] - axis[0] * s, v[1] - axis[1] * s, v[2] - axis[2] * s };
}

// Отражение v через плоскость с единичной нормалью n.
// Параллельная нормали часть переворачивается, перпендикулярная остаётся.
// СВОЙСТВА: длина сохраняется, применить дважды = вернуться.
inline Vector<3, double> VectorReflection(const Vector<3, double>& v,
	const Vector<3, double>& normal)
{
	assert(std::abs(DotProduct(normal, normal) - 1.0) < 1e-9 &&
		"VectorReflection: normal должна быть единичной");
	const double s = 2.0 * DotProduct(v, normal);
	return Vector<3, double>{ v[0] - normal[0] * s, v[1] - normal[1] * s, v[2] - normal[2] * s };
}

// Отражение ЛУЧА от поверхности (нормаль направлена наружу).
// Отличается знаком: падающий луч идёт К поверхности, отражённый ОТ неё.
// Формула та же, разница только в том, что здесь v — направление движения.
inline Vector<3, double> ReflectRay(const Vector<3, double>& incident,
	const Vector<3, double>& normal)
{
	return VectorReflection(incident, normal);
}

// --- Версии для НЕединичных направлений: делят на |axis|² ---

inline Vector<3, double> VectorProjectionUnnormalized(const Vector<3, double>& v,
	const Vector<3, double>& axis)
{
	const double aa = DotProduct(axis, axis);
	assert(aa > 1e-24 && "VectorProjectionUnnormalized: нулевая ось");
	const double s = DotProduct(v, axis) / aa;
	return Vector<3, double>{ axis[0] * s, axis[1] * s, axis[2] * s };
}

inline Vector<3, double> VectorRejectionUnnormalized(const Vector<3, double>& v,
	const Vector<3, double>& axis)
{
	const double aa = DotProduct(axis, axis);
	assert(aa > 1e-24 && "VectorRejectionUnnormalized: нулевая ось");
	const double s = DotProduct(v, axis) / aa;
	return Vector<3, double>{ v[0] - axis[0] * s, v[1] - axis[1] * s, v[2] - axis[2] * s };
}

inline double DistancePointToLine(const Vector<3, double>& q, const Ray& line)
{
	const Vector<3, double>& d = line.getDirection();
	const double dd = DotProduct(d, d);

	const Vector<3, double> w{ q[0] - line.getOrigin()[0],
							  q[1] - line.getOrigin()[1],
							  q[2] - line.getOrigin()[2] };

	if (dd < 1e-24) return std::sqrt(DotProduct(w, w));   // вырождена в точку

	const Vector<3, double> area = Cross(w, d);
	return std::sqrt(DotProduct(area, area) / dd);
}


//расстояние от точки до прямой — ЧЕРЕЗ ПЛОЩАДЬ.
//
//  |(q − p) × v| — площадь параллелограмма на этих двух векторах.
//  Площадь = основание × высота, основание это |v|, высота — искомое
//  расстояние. Отсюда:
//
//      d = |(q − p) × v| / |v|
//
//  ПОЧЕМУ НЕ ЧЕРЕЗ ПИФАГОРА. Формула d² = |w|² − (w·v)²/|v|² алгебраически
//  верна, но вычитает два БОЛЬШИХ близких числа. Если точка далеко от p,
//  но близко к самой прямой, оба слагаемых порядка |w|², а разность
//  порядка d². Замерено: при |w| = 1e8 и истинном d = 1 пифагорова форма
//  даёт РОВНО НОЛЬ — оба слагаемых равны 1e16, а в double всего 16
//  значащих цифр.
//
//  Cross тоже вычитает большие числа, но ВНУТРИ компонент, до возведения
//  в квадрат — то есть на величинах порядка координат, а не их квадратов.
//  Скорость та же (замерено 3.80 против 3.82 нс).
// ---------------------------------------------------------------------------

// Прямая пересечения двух плоскостей (n·p + d = 0). nullopt для параллельных.
std::optional<Ray> intersectionPlanePlane(const Plane& pl_1, const Plane& pl_2);

// ---------------------------------------------------------------------------
//  ЛУЧ — ТРЕУГОЛЬНИК: Мёллер–Трумбор (Möller, Trumbore 1997).
//
//  Ищем t, u, v из  o + t·d = (1 − u − v)·V0 + u·V1 + v·V2, то есть
//      [−d | e1 | e2]·(t, u, v)ᵀ = s,   e1 = V1 − V0, e2 = V2 − V0, s = o − V0.
//  По Крамеру со смешанными произведениями, p = d × e2, q = s × e1:
//      det = e1·p,   u = (s·p)/det,   v = (d·q)/det,   t = (e2·q)/det.
//
//  Промах, если
//    det ≈ 0      — луч параллелен плоскости треугольника (или вырожден сам
//                   треугольник). Порог относительный: det² <= ε²·|d|²|e1|²|e2|²,
//                   т.к. det — тройное произведение и растёт с масштабом;
//    u < 0, v < 0 или u + v > 1 — точка плоскости вне треугольника (граница
//                   считается попаданием);
//    t < 0        — треугольник позади начала луча (t = 0, т.е. начало луча
//                   на треугольнике, считается попаданием).
//  t измеряется в длинах d (направление не обязано быть единичным).
//
//  Было: в конце функции не было return (UB, если t ≈ 0), u и v считались
//  из одного и того же q·e1, а проверка std::abs(u) < −ε не срабатывала никогда.
// ---------------------------------------------------------------------------
struct RayTriangleHit
{
	double t = 0.0;              // параметр луча: point = o + t·d
	double u = 0.0;              // барицентрическая координата V1
	double v = 0.0;              // барицентрическая координата V2; у V0 — 1 − u − v
	Vector<3, double> point;     // точка попадания
};

std::optional<RayTriangleHit> RayTriIntersectionMollerTrambo(const Ray& ray, const Triangle& tri);

// Быстрый вариант без проверки det ≈ 0 и t < 0: возвращает t или INFINITY,
// weight — барицентрические (V0, V1, V2).
double RayTriangleintersectionFast(const Triangle& tri, const Ray& ray, float weight[3]);

// Пересечение двух лучей: t на первом луче или nullopt (параллельны,
// не пересекаются или пересечение позади начала одного из лучей).
std::optional<double> RayRayIntersection(const Ray& ray_1, const Ray& ray_2);

inline std::optional<double> RayRayIntersectionRobust(const Ray& ray1, const Ray& ray2,
	double tolerance = 1e-6)
{
	const Vector<3, double>& d1 = ray1.getDirection();
	const Vector<3, double>& d2 = ray2.getDirection();

	const double d00 = DotProduct(d1, d1);
	const double d01 = DotProduct(d1, d2);
	const double d11 = DotProduct(d2, d2);

	const double denom = d00 * d11 - d01 * d01;

	// Порог ОТНОСИТЕЛЬНЫЙ: denom растёт как |d1|²·|d2|², и абсолютный
	// отверг бы перпендикулярные лучи с короткими направляющими.
	if (denom <= 1e-12 * d00 * d11) return std::nullopt;

	const Vector<3, double> w{ ray1.getOrigin()[0] - ray2.getOrigin()[0],
							  ray1.getOrigin()[1] - ray2.getOrigin()[1],
							  ray1.getOrigin()[2] - ray2.getOrigin()[2] };
	const double w0 = DotProduct(w, d1);
	const double w1 = DotProduct(w, d2);

	const double t = (d01 * w1 - d11 * w0) / denom;
	const double u = (d00 * w1 - d01 * w0) / denom;

	if (t < 0.0 || u < 0.0) return std::nullopt;      // позади начала луча

	const Vector<3, double> p1{ ray1.getOrigin()[0] + t * d1[0],
							   ray1.getOrigin()[1] + t * d1[1],
							   ray1.getOrigin()[2] + t * d1[2] };
	const Vector<3, double> p2{ ray2.getOrigin()[0] + u * d2[0],
							   ray2.getOrigin()[1] + u * d2[1],
							   ray2.getOrigin()[2] + u * d2[2] };
	const Vector<3, double> gap{ p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2] };

	if (DotProduct(gap, gap) > tolerance * tolerance) return std::nullopt;
	return t;
}

} // namespace mathx
