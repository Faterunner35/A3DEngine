#pragma once
#include"GeomPPrimitives.h"
#include"Matrix.h"


struct AABBBoundingVolume
{
	AABBBoundingVolume()
		:
		box{}
	{
	}

	void GetSplittingAxis(Vector<3, double>& origin) const
	{
		double zero = static_cast<float>(0);
		double one = static_cast<float>(1);
		double half = static_cast<float>(0.5);


		origin = half * (box.max + box.min);
		Vector<3, double> extents = half * (box.max - box.min);
		float projectionExtent = extents[0];
		Vector<3, double> direction({ one, zero, zero });
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


	inline void transformAffine(const Matrix<4, 4,double>& matrix);

	/// Returns whether or not this box intersects another.
	inline bool intersects(const AABBBoundingVolume& b2) const;

	/// Calculate the area of intersection of this box and another
	inline AABBBoundingVolume intersection(const AABBBoundingVolume& b2) const;

	/// Calculate the volume of this box
	inline double volume() const;

	/// Tests whether another box contained by this box.
	inline bool contains(const AABBBoundingVolume& other) const;

	/// Tests whether the given point contained by this box.
	inline bool contains(const Vector<3,double>& v) const;

	/// Returns the square of the minimum distance between a given point and any part of the box.
	inline double squaredDistance(const Vector<3,double>& v) const;

	/// Returns the minimum distance between a given point and any part of the box.
	inline double distance(const Vector<3,double>& v) const;

	/// Returns the radius of a sphere enclosing the aabb from the outside at center mCenter
	inline double getRadius() const;

	/// Returns the radius of a sphere enclosing the aabb from origin as center
	inline double getRadiusOrigin() const;


	AlignedBox<double>box;
};



