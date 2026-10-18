#pragma once
#include"Geomerty.h"
#include"Matrix.h"
#include"Ray.h"
#include "Matrix3x3.hpp"
#include<cassert>
#include<optional>



//  Плоскость: n·p + d = 0, нормаль единичная.
//
//  d — расстояние от начала координат до плоскости со знаком МИНУС.
//  Подстановка точки f·p = n·p + d даёт её расстояние до плоскости со знаком.
//
//  Правило нуля: никаких явных копирования, присваивания и деструктора —
//  класс держит только два тривиальных поля, компилятор справится лучше.
namespace mathx
{
class Plane
{
public:
	Plane() = default;

	Plane(const Vector<3, double>& normal, double distance)
		: normal_(normal), distance_(distance) {
	}

	// Через три точки. Порядок обхода задаёт направление нормали.
	Plane(const Vector<3, double>& v0, const Vector<3, double>& v1,
		const Vector<3, double>& v2)
	{
		const Vector<3, double> e1{ v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2] };
		const Vector<3, double> e2{ v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2] };
		normal_ = Cross(e1, e2);
		const double len2 = DotProduct(normal_, normal_);
		assert(len2 > 1e-24 && "Plane: три точки коллинеарны — плоскость не определена");
		const double inv = 1.0 / std::sqrt(len2);
		normal_ = Vector<3, double>{ normal_[0] * inv, normal_[1] * inv, normal_[2] * inv };
		distance_ = -DotProduct(v0, normal_);
	}

	const Vector<3, double>& GetNormal() const { return normal_; }
	double GetDistance() const { return distance_; }

	// Расстояние от точки до плоскости СО ЗНАКОМ.
	// Положительное — точка с той стороны, куда смотрит нормаль.
	double SignedDistance(const Vector<3, double>& p) const {
		return DotProduct(normal_, p) + distance_;
	}

	// Проекция точки на плоскость: сдвиг на расстояние против нормали.
	Vector<3, double> ProjectPoint(const Vector<3, double>& p) const {
		const double h = SignedDistance(p);
		return Vector<3, double>{ p[0] - h * normal_[0], p[1] - h * normal_[1], p[2] - h * normal_[2] };
	}

private:
	Vector<3, double> normal_{ 0,0,1 };
	double distance_ = 0.0;
};

// ---------------------------------------------------------------------------
//  ПЕРЕСЕЧЕНИЕ ЛУЧА С ПЛОСКОСТЬЮ.
//
//  Подставляем P(t) = o + t·d в уравнение плоскости n·P + D = 0:
//      n·(o + t·d) + D = 0
//      n·o + t(n·d) + D = 0
//      t = −(n·o + D) / (n·d)
//
//  Числитель — это расстояние от начала луча до плоскости со знаком.
//  Знаменатель n·d — скорость приближения к плоскости вдоль луча.
//
//  ВЫРОЖДЕНИЕ: n·d = 0 значит луч ПАРАЛЛЕЛЕН плоскости. Тогда либо он
//  целиком в плоскости (числитель тоже ноль), либо не пересекает никогда.
//  Оба случая возвращают nullopt: пересечения в виде одной точки нет.
//
//  Порог ОТНОСИТЕЛЬНЫЙ: n·d зависит от длины направляющей, и абсолютный
//  отверг бы короткие направляющие как параллельные.
// ---------------------------------------------------------------------------
inline std::optional<double> RayPlaneIntersection(const Ray& ray, const Plane& plane,
	double parallelTolerance = 1e-12)
{
	const Vector<3, double>& d = ray.getDirection();
	const double nd = DotProduct(plane.GetNormal(), d);

	const double dd = DotProduct(d, d);
	if (nd * nd <= parallelTolerance * dd) return std::nullopt;   // параллелен

	const double t = -plane.SignedDistance(ray.getOrigin()) / nd;
	if (t < 0.0) return std::nullopt;                            // позади начала
	return t;
}

// Точка пересечения, если она есть.
inline std::optional<Vector<3, double>> RayPlanePoint(const Ray& ray, const Plane& plane,
	double parallelTolerance = 1e-12)
{
	auto t = RayPlaneIntersection(ray, plane, parallelTolerance);
	if (!t) return std::nullopt;
	const Vector<3, double>& d = ray.getDirection();
	return Vector<3, double>{ ray.getOrigin()[0] + (*t) * d[0],
		ray.getOrigin()[1] + (*t) * d[1],
		ray.getOrigin()[2] + (*t) * d[2] };
}

// ---------------------------------------------------------------------------
//  ПЕРЕСЕЧЕНИЕ ДВУХ ПЛОСКОСТЕЙ — это ПРЯМАЯ, а не точка.
//
//  Направление: n1 × n2 — перпендикулярно обеим нормалям, значит лежит
//  в обеих плоскостях.
//
//  Точка на прямой: прямая бесконечна, нужен однозначный представитель.
//  Берём БЛИЖАЙШУЮ К НАЧАЛУ координат — она лежит в плоскости, натянутой
//  на нормали, то есть перпендикулярна c.
//
//      c = n1 × n2,   a = c × n1,   b = n2 × c
//      p = (−d1·b − d2·a) / |c|²
//
//  ПОЧЕМУ РАБОТАЕТ. Базис {a, b} выбран так, что каждое слагаемое отвечает
//  ровно за одно уравнение и не мешает другому:
//      a·n1 = 0   (a = c × n1 перпендикулярен n1)
//      b·n2 = 0   (b = n2 × c перпендикулярен n2)
//      b·n1 = a·n2 = |c|²   (смешанные произведения при c = n1 × n2)
//
//  Подставляя в первое уравнение, член с a исчезает, остаётся
//  (−d1·|c|²)/|c|² = −d1 — ровно то, что требуется. Во втором симметрично
//  исчезает член с b.
//
//  Возвращает Ray: origin — точка на прямой, direction — её направление.
//  nullopt при параллельных плоскостях.
// ---------------------------------------------------------------------------
inline std::optional<Ray> PlanePlaneIntersection(const Plane& p1, const Plane& p2,
	double parallelTolerance = 1e-12)
{
	const Vector<3, double> c = Cross(p1.GetNormal(), p2.GetNormal());
	const double cc = DotProduct(c, c);
	if (cc <= parallelTolerance) return std::nullopt;   // нормали коллинеарны

	const Vector<3, double> a = Cross(c, p1.GetNormal());   // ⊥ c и ⊥ n1
	const Vector<3, double> b = Cross(p2.GetNormal(), c);   // ⊥ c и ⊥ n2
	const double inv = 1.0 / cc;

	// ВНИМАНИЕ на скобки: деление применяется ко ВСЕМУ числителю.
	const Vector<3, double> point{
		(-p1.GetDistance() * b[0] - p2.GetDistance() * a[0]) * inv,
		(-p1.GetDistance() * b[1] - p2.GetDistance() * a[1]) * inv,
		(-p1.GetDistance() * b[2] - p2.GetDistance() * a[2]) * inv };

	return Ray{ point, c };
}

// ---------------------------------------------------------------------------
//  Ортонормированный базис по нормали (Duff et al., "Building an Orthonormal
//  Basis, Revisited"). Без ветвлений и без деления на почти ноль:
//  знаменатель sign + n.z по модулю всегда >= 1.
//
//  Столбцы: t, b — касательные в плоскости, n — нормаль.
//  Это матрица перехода ИЗ локального базиса плоскости В мировой.
// ---------------------------------------------------------------------------
inline Matrix3x3 MakeTangentBasis(const Vector<3, double>& n)
{
	assert(std::abs(DotProduct(n, n) - 1.0) < 1e-9 && "MakeTangentBasis: нормаль не единична");
	const double sign = std::copysign(1.0, n[2]);
	const double a = -1.0 / (sign + n[2]);
	const double b = n[0] * n[1] * a;

	const Vector<3, double> t{ 1.0 + sign * n[0] * n[0] * a,  sign * b,               -sign * n[0] };
	const Vector<3, double> bt{ b,                       sign + n[1] * n[1] * a,   -n[1] };

	// столбцы: t, bt, n
	return Matrix3x3(t[0], bt[0], n[0],
		t[1], bt[1], n[1],
		t[2], bt[2], n[2]);
}


// ---------------------------------------------------------------------------
//  ПЕРЕСЕЧЕНИЕ ТРЁХ ПЛОСКОСТЕЙ — точка.
//
//  Каждая плоскость даёт уравнение n_i·p + d_i = 0. Три уравнения образуют
//  систему, где нормали стоят СТРОКАМИ:
//
//      | n1 |         | −d1 |
//      | n2 | · p  =  | −d2 |
//      | n3 |         | −d3 |
//
//  Обратная к матрице со строками n_i имеет векторные произведения в
//  СТОЛБЦАХ (следствие тождества (Mᵀ)⁻¹ = (M⁻¹)ᵀ), поэтому:
//
//      p = [ d1·(n3 × n2) + d2·(n1 × n3) + d3·(n2 × n1) ] / [n1,n2,n3]
//
//  Порядок в каждом cross перевёрнут, чтобы погасить минусы правой части.
//
//  СМЫСЛ СЛАГАЕМЫХ: n2 × n3 — направление прямой пересечения второй и
//  третьей плоскостей. Умножая на d1, задаём, как далеко идти вдоль этой
//  прямой до первой плоскости. Каждое слагаемое отвечает за одну плоскость.
//
//  ВЫРОЖДЕНИЕ. Знаменатель [n1,n2,n3] — смешанное произведение, то есть
//  ОБЪЁМ параллелепипеда на нормалях. Ноль означает компланарность нормалей,
//  и это покрывает ОБА случая сразу:
//    (а) две плоскости параллельны;
//    (б) ни одна пара не параллельна, но линии их попарного пересечения
//        параллельны — треугольная призма без вершины.
//  Геометрически: компланарные нормали означают общее направление, вдоль
//  которого все три плоскости бесконечны, поэтому единственной точки нет.
//
//  Порог ОТНОСИТЕЛЬНЫЙ: нормали единичны, так что объём сравнивается с 1.
// ---------------------------------------------------------------------------
inline std::optional<Vector<3, double>> ThreePlaneIntersection(
	const Plane& p1, const Plane& p2, const Plane& p3,
	double tolerance = 1e-9)
{
	const Vector<3, double>& n1 = p1.GetNormal();
	const Vector<3, double>& n2 = p2.GetNormal();
	const Vector<3, double>& n3 = p3.GetNormal();

	const Vector<3, double> c23 = Cross(n2, n3);
	const double det = DotProduct(n1, c23);          // = [n1,n2,n3]
	if (std::abs(det) <= tolerance) return std::nullopt;

	const Vector<3, double> c31 = Cross(n3, n1);
	const Vector<3, double> c12 = Cross(n1, n2);
	const double inv = 1.0 / det;

	// p = −(d1·c23 + d2·c31 + d3·c12) / det
	const double a = -p1.GetDistance() * inv;
	const double b = -p2.GetDistance() * inv;
	const double c = -p3.GetDistance() * inv;

	return Vector<3, double>{ a* c23[0] + b * c31[0] + c * c12[0],
		a* c23[1] + b * c31[1] + c * c12[1],
		a* c23[2] + b * c31[2] + c * c12[2] };
}

} // namespace mathx


