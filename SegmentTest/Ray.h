#pragma once
#include"Geomerty.h"


namespace mathx
{



class  Ray
{
protected:
    Vector<3,double> mOrigin;
    Vector<3,double> mDirection;
    double direction;

public:
    Ray() : mOrigin({0,0,0}), mDirection({0,0,1}) {}
    Ray(const Vector<3,double>& origin, const Vector<3,double>& direction) :
        mOrigin(origin),
        mDirection(direction)
    {
    }

    /** Sets the origin of the ray. */
    void setOrigin(const Vector<3,double>& origin) { mOrigin = origin; }
    /** Gets the origin of the ray. */
    const Vector<3,double>& getOrigin() const { return mOrigin; }

    /** Sets the direction of the ray. */
    void setDirection(const Vector<3,double>& dir) { mDirection = dir; }
    /** Gets the direction of the ray. */
    const Vector<3,double>& getDirection() const { return mDirection; }

    /** Gets the position of a point t units along the ray. */
    Vector<3, double> getPoint(double t) const { return mOrigin + mDirection * t; }

    /** Gets the position of a point t units along the ray. */
    Vector<3,double> operator*(double t) const { return getPoint(t); }

    inline double  GetScalarDirection() const
    {
        return direction;
    }

    /** Tests whether this ray intersects the given plane.
    @return A pair structure where the first element indicates whether
        an intersection occurs, and if true, the second element will
        indicate the distance along the ray at which it intersects.
        This can be converted to a point in space by calling getPoint().
    */
    //std::pair<bool, double> intersects(const Plane& p) const { return Geom_Exp::GeomCoreIntersector::PlaneRayIntersection(*this, p); }
    ///** Tests whether this ray intersects the given plane bounded volume.
    //@return A pair structure where the first element indicates whether
    //an intersection occurs, and if true, the second element will
    //indicate the distance along the ray at which it intersects.
    //This can be converted to a point in space by calling getPoint().
    //*/
    //std::pair<bool, double> intersects(const AABBBoundingVolume& p) const
    //{
    //    return GeomCoreIntersector::intersects(*this, p.planes, p.outside == Plane<3,double>::POSITIVE_SIDE);
    //}
    ///** Tests whether this ray intersects the given sphere.
    //@return A pair structure where the first element indicates whether
    //    an intersection occurs, and if true, the second element will
    //    indicate the distance along the ray at which it intersects.
    //    This can be converted to a point in space by calling getPoint().
    //*/
    //std::pair<bool,double> intersects(const Sphere& s) const
    //{
    //    return Math::intersects(*this, s);
    //}
    ///** Tests whether this ray intersects the given box.
    //@return A pair structure where the first element indicates whether
    //    an intersection occurs, and if true, the second element will
    //    indicate the distance along the ray at which it intersects.
    //    This can be converted to a point in space by calling getPoint().
    //*/
    //std::pair<bool, double> intersects(const AxisAlignedBox& box) const
    //{
    //    return Math::intersects(*this, box);
    //}

};

}