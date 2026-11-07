#include "Primitives.h"
#include <algorithm>
#include <cassert>
#include <cmath>

namespace mathx
{

// =============================================================================
//  Triangle
// =============================================================================

// Ericson, RTCD 5.1.5. Треугольник ABC, точка P. Плоскость треугольника
// делится на 7 областей Вороного: 3 вершины, 3 ребра, грань. Область
// определяется знаками скалярных произведений, а не проекцией на плоскость:
//
//   d1 = AB·AP, d2 = AC·AP     вершина A:  d1 <= 0 и d2 <= 0
//   d3 = AB·BP, d4 = AC·BP     вершина B:  d3 >= 0 и d4 <= d3
//   d5 = AB·CP, d6 = AC·CP     вершина C:  d6 >= 0 и d5 <= d6
//
// Для рёбер нужны ещё «барицентрические площади» (с точностью до общего
// множителя |AB×AC|²):
//   vc = d1·d4 − d3·d2   ребро AB:  vc <= 0, d1 >= 0, d3 <= 0
//   vb = d5·d2 − d1·d6   ребро AC:  vb <= 0, d2 >= 0, d6 <= 0
//   va = d3·d6 − d5·d4   ребро BC:  va <= 0, d4 − d3 >= 0, d5 − d6 >= 0
// Иначе P проецируется внутрь грани, барицентрические (u, v, w) =
// (va, vb, vc) / (va + vb + vc), ответ A + v·AB + w·AC.
Vector<3, double> Triangle::ClosestPoint(const Vector<3, double>& point) const
{
	const Vector<3, double>& a = triangle_data[0];
	const Vector<3, double>& b = triangle_data[1];
	const Vector<3, double>& c = triangle_data[2];

	const Vector<3, double> ab = b - a;
	const Vector<3, double> ac = c - a;

	const Vector<3, double> ap = point - a;
	const double d1 = DotProduct(ab, ap);
	const double d2 = DotProduct(ac, ap);
	if (d1 <= 0.0 && d2 <= 0.0)
		return a;                                         // вершина A

	const Vector<3, double> bp = point - b;
	const double d3 = DotProduct(ab, bp);
	const double d4 = DotProduct(ac, bp);
	if (d3 >= 0.0 && d4 <= d3)
		return b;                                         // вершина B

	const double vc = d1 * d4 - d3 * d2;
	if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
	{
		const double v = d1 / (d1 - d3);
		return a + v * ab;                                // ребро AB
	}

	const Vector<3, double> cp = point - c;
	const double d5 = DotProduct(ab, cp);
	const double d6 = DotProduct(ac, cp);
	if (d6 >= 0.0 && d5 <= d6)
		return c;                                         // вершина C

	const double vb = d5 * d2 - d1 * d6;
	if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
	{
		const double w = d2 / (d2 - d6);
		return a + w * ac;                                // ребро AC
	}

	const double va = d3 * d6 - d5 * d4;
	if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
	{
		const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		return b + w * (c - b);                           // ребро BC
	}

	const double denom = 1.0 / (va + vb + vc);            // грань
	const double v = vb * denom;
	const double w = vc * denom;
	return a + v * ab + w * ac;
}

double Triangle::PointToTriangleDisctance(const Vector<3, double>& point) const
{
	return Length(point - ClosestPoint(point));
}

// =============================================================================
//  Polygon
// =============================================================================
namespace
{
	// Допуск совпадения вершин (покомпонентно), как positionEquals в Ogre.
	constexpr double kPolygonVertexTolerance = 1e-9;

	bool PositionEquals(const Vector<3, double>& a, const Vector<3, double>& b)
	{
		return std::abs(a[0] - b[0]) <= kPolygonVertexTolerance &&
			std::abs(a[1] - b[1]) <= kPolygonVertexTolerance &&
			std::abs(a[2] - b[2]) <= kPolygonVertexTolerance;
	}
}

Polygon::Polygon()
	: mVertexList(), mNormal(), mIsNormalSet(false)
{
}

// Нормаль по формуле Ньюэлла: сумма по рёбрам (x_i, y_i, z_i) -> (x_j, y_j, z_j)
//   N.x += (y_i − y_j)(z_i + z_j),  N.y += (z_i − z_j)(x_i + x_j),
//   N.z += (x_i − x_j)(y_i + y_j)
// Это удвоенная векторная площадь многоугольника. В отличие от одного
// Cross по трём вершинам, не ломается, если первые три вершины коллинеарны.
// Направление — по правилу правой руки для обхода вершин.
void Polygon::updateNormal() const
{
	assert(getVertexCount() >= 3 && "Polygon::updateNormal: нужно минимум 3 вершины");
	if (mIsNormalSet)
		return;

	Vector<3, double> n;
	const size_t count = mVertexList.size();
	for (size_t i = 0; i < count; ++i)
	{
		const Vector<3, double>& cur = mVertexList[i];
		const Vector<3, double>& nxt = mVertexList[(i + 1) % count];
		n[0] += (cur[1] - nxt[1]) * (cur[2] + nxt[2]);
		n[1] += (cur[2] - nxt[2]) * (cur[0] + nxt[0]);
		n[2] += (cur[0] - nxt[0]) * (cur[1] + nxt[1]);
	}
	Normalize(n);

	mNormal = n;
	mIsNormalSet = true;
}

void Polygon::insertVertex(const Vector<3, double>& vdata, size_t vertex)
{
	assert(vertex <= getVertexCount() && "Polygon::insertVertex: индекс вне диапазона");
	mVertexList.insert(mVertexList.begin() + static_cast<std::ptrdiff_t>(vertex), vdata);
	mIsNormalSet = false;
}

void Polygon::insertVertex(const Vector<3, double>& vdata)
{
	mVertexList.push_back(vdata);
	mIsNormalSet = false;
}

Vector<3, double> Polygon::CalculateNormalToTriangle(const std::array<Vector<3, double>, 3>& triangle)
{
	return UnitCross(triangle[1] - triangle[0], triangle[2] - triangle[0]);
}

const Vector<3, double>& Polygon::getVertex(size_t vertex) const
{
	assert(vertex < getVertexCount() && "Polygon::getVertex: индекс вне диапазона");
	return mVertexList[vertex];
}

void Polygon::setVertexListSize(size_t number_of_vertex)
{
	mVertexList.resize(number_of_vertex);
	mIsNormalSet = false;
}

void Polygon::setVertex(const Vector<3, double>& vdata, size_t vertexIndex)
{
	assert(vertexIndex < getVertexCount() && "Polygon::setVertex: индекс вне диапазона");
	mVertexList[vertexIndex] = vdata;
	mIsNormalSet = false;
}

// Удаляет вершину, совпадающую со следующей по обходу (включая пару
// последняя–первая). Как в Ogre::Polygon::removeDuplicates.
void Polygon::removeDuplicates()
{
	for (size_t i = 0; i < mVertexList.size() && mVertexList.size() > 1; )
	{
		const size_t next = (i + 1) % mVertexList.size();
		if (PositionEquals(mVertexList[i], mVertexList[next]))
			deleteVertex(i);
		else
			++i;
	}
}

size_t Polygon::getVertexCount() const
{
	return mVertexList.size();
}

const Vector<3, double>& Polygon::getNormal() const
{
	assert(getVertexCount() >= 3 && "Polygon::getNormal: нужно минимум 3 вершины");
	updateNormal();
	return mNormal;
}

void Polygon::deleteVertex(size_t vertex)
{
	assert(vertex < getVertexCount() && "Polygon::deleteVertex: индекс вне диапазона");
	mVertexList.erase(mVertexList.begin() + static_cast<std::ptrdiff_t>(vertex));
	mIsNormalSet = false;
}

// Выпуклый плоский многоугольник, вершины против часовой стрелки
// относительно n = getNormal(). Точка P внутри (граница включена), если
//   1) лежит в плоскости: |n·(P − V0)| мало;
//   2) для каждого ребра Vi -> Vj она не правее его: ((Vj − Vi) × (P − Vi))·n >= 0.
// Допуски относительные: обе величины масштабируются с размером фигуры.
bool Polygon::isPointInside(const Vector<3, double>& point) const
{
	const size_t count = getVertexCount();
	if (count < 3)
		return false;

	const Vector<3, double>& n = getNormal();
	if (DotProduct(n, n) == 0.0)
		return false;                                     // вырожденный многоугольник

	const double tolerance = 1e-9;
	const Vector<3, double> w0 = point - mVertexList[0];
	if (std::abs(DotProduct(n, w0)) > tolerance * std::max(1.0, Length(w0)))
		return false;                                     // не в плоскости

	for (size_t i = 0; i < count; ++i)
	{
		const Vector<3, double>& vi = mVertexList[i];
		const Vector<3, double> edge = mVertexList[(i + 1) % count] - vi;
		const Vector<3, double> w = point - vi;
		if (DotProduct(Cross(edge, w), n) < -tolerance * Length(edge) * std::max(1.0, Length(w)))
			return false;                                 // снаружи от ребра i
	}
	return true;
}

void Polygon::reset()
{
	mVertexList.clear();
	mNormal.MakeZero();
	mIsNormalSet = false;
}

// Равны, если вершины совпадают с точностью до циклического сдвига начала
// обхода (направление обхода — то же). Как в Ogre::Polygon::operator==.
bool Polygon::operator==(const Polygon& rhs) const
{
	const size_t count = getVertexCount();
	if (count != rhs.getVertexCount())
		return false;
	if (count == 0)
		return true;

	for (size_t start = 0; start < count; ++start)
	{
		if (!PositionEquals(mVertexList[0], rhs.mVertexList[start]))
			continue;
		bool equal = true;
		for (size_t i = 1; i < count && equal; ++i)
			equal = PositionEquals(mVertexList[i], rhs.mVertexList[(i + start) % count]);
		if (equal)
			return true;
	}
	return false;
}

std::ostream& operator<<(std::ostream& strm, const Polygon& poly)
{
	strm << "NUM VERTICES: " << poly.getVertexCount() << "\n";
	for (size_t j = 0; j < poly.getVertexCount(); ++j)
	{
		const Vector<3, double>& v = poly.getVertex(j);
		strm << "VERTEX " << j << ": (" << v[0] << ", " << v[1] << ", " << v[2] << ")\n";
	}
	return strm;
}

} // namespace mathx
