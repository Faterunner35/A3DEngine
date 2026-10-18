#pragma once
#include"Geomerty.h"


namespace mathx
{ 
class  Sphere
{

public:
    /** Standard constructor - creates a unit sphere around the origin.*/
    Sphere() : mRadius(1.0), mCenter({0,0,0}) {}
    /** Constructor allowing arbitrary spheres.
        @param center The center point of the sphere.
        @param radius The radius of the sphere.
    */
    Sphere(const Vector<3,double>& center, double radius) : mRadius(radius), mCenter(center) {}

    /** Returns the radius of the sphere. */
    double getRadius() const { return mRadius; }

    /** Sets the radius of the sphere. */
    void setRadius(double radius) { mRadius = radius; }

    /** Returns the center point of the sphere. */
    const Vector<3,double>& getCenter() const { return mCenter; }

    /** Sets the center point of the sphere. */
    void setCenter(const Vector<3,double>& center) { mCenter = center; }

    /** Returns whether or not this sphere intersects another sphere. */
    bool intersects(const Sphere& s) const
    {
        return mathx::Length(s.mCenter - mCenter) <=(s.mRadius + mRadius);
    }
    /** Returns whether or not this sphere intersects a box. */
    //bool intersects(const AABBBoundingVolume& box) const { return GeomCoreIntersector::intersects(*this, box); }
    ///** Returns whether or not this sphere intersects a plane. */
    //bool intersects(const Plane<3,double>& plane) const { return GeomCoreIntersector::intersects(*this, plane); }
    ///** Returns whether or not this sphere intersects a point. */
    //bool intersects(const Vector<3,double>& v) const
    //{
    //    return ((v - mCenter).squaredLength() <= Math::Sqr(mRadius));
    //}
    /** Merges another Sphere into the current sphere */
    void merge(const Sphere& oth)
    {
        Vector<3,double> diff = oth.getCenter() - mCenter;
        double    lengthSq = mathx::Length(diff);
        double    radiusDiff = oth.getRadius() - mRadius;

        // Early-out
        if (std::sqrt(radiusDiff) >= lengthSq)
        {
            // One fully contains the other
            if (radiusDiff <= 0.0f)
                return;  // no change
            else
            {
                mCenter = oth.getCenter();
                mRadius = oth.getRadius();
                return;
            }
        }

        double length = std::sqrt(lengthSq);
        double t = (length + radiusDiff) / (2.0f * length);
        mCenter = mCenter + diff * t;
        mRadius = 0.5f * (length + mRadius + oth.getRadius());
    }

    protected:
	double    mRadius;
	Vector<3, double> mCenter;
};

}

