#include "Plane.h"

namespace mathx
{

Plane::Plane(const Vector<3, double>& v0, const Vector<3, double>& v1,
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

std::optional<Ray> PlanePlaneIntersection(const Plane& p1, const Plane& p2,
	double parallelTolerance)
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

std::optional<Vector<3, double>> ThreePlaneIntersection(
	const Plane& p1, const Plane& p2, const Plane& p3,
	double tolerance)
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
