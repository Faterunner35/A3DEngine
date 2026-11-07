#include "Sphere.h"

namespace mathx
{

void Sphere::merge(const Sphere& oth)
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

} // namespace mathx
