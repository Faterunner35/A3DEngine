#include "Sphere.h"

namespace mathx
{

// Наименьшая сфера, содержащая обе: эта (c1, r1) и oth (c2, r2).
//
//  d = |c2 − c1| — расстояние между центрами (длина, НЕ квадрат).
//
//  1) d + r2 <= r1 — oth целиком внутри этой сферы, ничего не меняем.
//  2) d + r1 <= r2 — эта сфера внутри oth, результат — oth.
//     Совпадающие центры (d = 0) всегда попадают в 1) или 2).
//  3) Иначе ответ касается обеих изнутри в двух дальних точках на прямой
//     центров: c1 − r1·u и c2 + r2·u, u = (c2 − c1)/d. Диаметр — расстояние
//     между ними, d + r1 + r2; центр — середина:
//         R = (d + r1 + r2) / 2
//         c = c1 + u·(R − r1)
//     Здесь d > 0: из d = 0 следует 1) или 2).
//
//  Было: sqrt бралась от разности радиусов (NaN при r2 < r1), а Length,
//  возвращающая длину, считалась квадратом длины.
void Sphere::merge(const Sphere& oth)
{
    const Vector<3,double> diff = oth.getCenter() - mCenter;
    const double dist = mathx::Length(diff);
    const double r1 = mRadius;
    const double r2 = oth.getRadius();

    if (dist + r2 <= r1)
        return;                       // oth внутри
    if (dist + r1 <= r2)
    {
        mCenter = oth.getCenter();    // эта сфера внутри oth
        mRadius = r2;
        return;
    }

    const double radius = 0.5 * (dist + r1 + r2);
    mCenter = mCenter + diff * ((radius - r1) / dist);
    mRadius = radius;
}

} // namespace mathx
