#pragma once
#include"GeomPPrimitives.h"
#include"Matrix.h"

// AABBBoundingVolume — не часть математической библиотеки, поэтому живёт вне
// mathx и квалифицирует математические типы явно (mathx::Vector, ...).

struct AABBBoundingVolume
{
	AABBBoundingVolume()
		:
		box{}
	{
	}

	// Ось разбиения (как GetSplittingAxis для AlignedBox в GTE): начало —
	// центр бокса, направление — ось с наибольшим полуразмером. Раньше
	// direction вычислялось в локальной переменной и терялось.
	void GetSplittingAxis(mathx::Vector<3, double>& origin, mathx::Vector<3, double>& direction) const
	{
		double const zero = 0.0;
		double const one = 1.0;
		double const half = 0.5;

		origin = half * (box.max + box.min);
		mathx::Vector<3, double> extents = half * (box.max - box.min);
		double projectionExtent = extents[0];
		direction = mathx::Vector<3, double>{ one, zero, zero };
		if (extents[1] > projectionExtent)
		{
			projectionExtent = extents[1];
			direction.SetNumbers({ zero, one, zero });
		}
		if (extents[2] > projectionExtent)
		{
			projectionExtent = extents[2];
			direction.SetNumbers({ zero, zero, one });
		}
	}


	inline void transformAffine(const mathx::Matrix<4, 4, double>& matrix);

	/// Returns whether or not this box intersects another.
	inline bool intersects(const AABBBoundingVolume& b2) const;

	/// Calculate the area of intersection of this box and another
	inline AABBBoundingVolume intersection(const AABBBoundingVolume& b2) const;

	/// Calculate the volume of this box
	inline double volume() const;

	/// Tests whether another box contained by this box.
	inline bool contains(const AABBBoundingVolume& other) const;

	/// Tests whether the given point contained by this box.
	inline bool contains(const mathx::Vector<3, double>& v) const;

	/// Returns the square of the minimum distance between a given point and any part of the box.
	inline double squaredDistance(const mathx::Vector<3, double>& v) const;

	/// Returns the minimum distance between a given point and any part of the box.
	inline double distance(const mathx::Vector<3, double>& v) const;

	/// Returns the radius of a sphere enclosing the aabb from the outside at center mCenter
	inline double getRadius() const;

	/// Returns the radius of a sphere enclosing the aabb from origin as center
	inline double getRadiusOrigin() const;


	mathx::AlignedBox<double> box;
};



