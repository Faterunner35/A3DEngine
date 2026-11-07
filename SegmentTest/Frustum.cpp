#include "Frustum.h"

double Frustum::CameraDepth(double zDevice) const
{
    const DepthParams dp = GetProjectionParamsAB();
    return mProjType == ProjectionType::Perspective ? dp.B / (zDevice - dp.A)
                                                    : (zDevice - dp.A) / dp.B;
}

std::array<mathx::Plane, Frustum::kPlaneCount> Frustum::ExtractPlanes(const mathx::Matrix4x4& A, ClipY clipY,
                                                                      DepthMode depth, bool infiniteFar)
{
    using V4 = mathx::Vector<4, double>;
    const double sy = ClipSignY(clipY);   // приводим строку y к «клип вниз»
    const V4 r0{ A(0, 0), A(0, 1), A(0, 2), A(0, 3) };
    const V4 rY{ sy * A(1, 0), sy * A(1, 1), sy * A(1, 2), sy * A(1, 3) };
    const V4 r2{ A(2, 0), A(2, 1), A(2, 2), A(2, 3) };
    const V4 r3{ A(3, 0), A(3, 1), A(3, 2), A(3, 3) };
    auto comb = [](const V4& u, const V4& v, double s) {
        return V4{ u[0] + s * v[0], u[1] + s * v[1], u[2] + s * v[2], u[3] + s * v[3] };
    };
    auto plane = [](const V4& v) {
        const double inv = 1.0 / std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        return mathx::Plane(mathx::Vector<3, double>{ v[0] * inv, v[1] * inv, v[2] * inv }, v[3] * inv);
    };
    std::array<mathx::Plane, kPlaneCount> pl;
    pl[kLeft]   = plane(comb(r3, r0, +1.0));   // x >= -w
    pl[kRight]  = plane(comb(r3, r0, -1.0));   // x <=  w
    pl[kTop]    = plane(comb(r3, rY, +1.0));   // y >= -w
    pl[kBottom] = plane(comb(r3, rY, -1.0));   // y <=  w
    const bool stdDepth = depth == DepthMode::Standard;
    const V4 zLow = r2, zHigh = comb(r3, r2, -1.0);   // z >= 0,  z <= w
    pl[kNear] = plane(stdDepth ? zLow : zHigh);
    pl[kFar]  = infiniteFar
        ? mathx::Plane(mathx::Vector<3, double>{ 0.0, 0.0, -1.0 }, std::numeric_limits<double>::infinity())
        : plane(stdDepth ? zHigh : zLow);
    return pl;
}

mathx::Matrix4x4 Frustum::MakeFrustumProjection(Radian fovY, double s, double n, double f,
                                                ClipY clipY)
{
    assert(s > 0.0 && n > 0.0 && f > n && "MakeFrustumProjection: нужны s > 0, 0 < n < f");
    const double g = FocalG(fovY), k = f / (f - n);
    return mathx::Matrix4x4(g / s, 0.0,              0.0, 0.0,
                     0.0,   ClipSignY(clipY) * g, 0.0, 0.0,
                     0.0,   0.0,              k,   -n * k,
                     0.0,   0.0,              1.0, 0.0);   // w = z
}

mathx::Matrix4x4 Frustum::MakeFrustumProjectionInverse(Radian fovY, double s, double n, double f,
                                                       ClipY clipY)
{
    assert(s > 0.0 && n > 0.0 && f > n);
    const double g = FocalG(fovY);
    return mathx::Matrix4x4(s / g, 0.0,                        0.0,               0.0,
                     0.0,   ClipSignY(clipY) / g,       0.0,               0.0,
                     0.0,   0.0,                        0.0,               1.0,
                     0.0,   0.0,                        (n - f) / (n * f), 1.0 / n);
}

mathx::Matrix4x4 Frustum::MakeReversedInfiniteProjection(Radian fovY, double s, double n,
                                                         ClipY clipY)
{
    assert(s > 0.0 && n > 0.0 && "MakeReversedInfiniteProjection: нужны s > 0, n > 0");
    const double g = FocalG(fovY);
    return mathx::Matrix4x4(g / s, 0.0,              0.0, 0.0,
                     0.0,   ClipSignY(clipY) * g, 0.0, 0.0,
                     0.0,   0.0,              0.0, n,
                     0.0,   0.0,              1.0, 0.0);
}

mathx::Matrix4x4 Frustum::MakeReversedInfiniteProjectionInverse(Radian fovY, double s, double n,
                                                                ClipY clipY)
{
    assert(s > 0.0 && n > 0.0);
    const double g = FocalG(fovY);
    return mathx::Matrix4x4(s / g, 0.0,                  0.0,     0.0,
                     0.0,   ClipSignY(clipY) / g, 0.0,     0.0,
                     0.0,   0.0,                  0.0,     1.0,
                     0.0,   0.0,                  1.0 / n, 0.0);
}

double Frustum::FocalG(Radian fovY)
{
    assert(fovY.ValueRadians() > 0.0 && fovY.ValueRadians() < kPi);
    return 1.0 / mathx::Tan(fovY * 0.5);
}

Frustum::DepthParams Frustum::GetProjectionParamsAB() const
{
    const bool   stdDepth = mDepthMode == DepthMode::Standard;
    const double n = mNearDist, f = mFarDist;
    if (mProjType == ProjectionType::Orthographic) {      // (6.48): линейная глубина
        const double k = 1.0 / (f - n);
        return stdDepth ? DepthParams{ -n * k, k } : DepthParams{ f * k, -k };
    }
    if (IsInfinite())                                      // (6.42) / (6.46)
        return stdDepth ? DepthParams{ 1.0 - kInfiniteEps, -n * (1.0 - kInfiniteEps) }
                        : DepthParams{ 0.0, n };
    const double k = 1.0 / (f - n);                        // (6.36) / (6.45)
    return stdDepth ? DepthParams{ f * k, -n * f * k } : DepthParams{ -n * k, n * f * k };
}

void Frustum::CalcProjectionParameters(double& left, double& right, double& top, double& bottom) const
{
    if (mProjType == ProjectionType::Perspective) {
        const double halfH = mathx::Tan(mFovy * 0.5);           // 1/g
        const double halfW = halfH * mAspect;                    // s/g
        // сдвиг задан в фокальной плоскости z = focal; на расстоянии 1 он в focal раз меньше
        const double offX = mFrustumOffset[0] / mFocalLength, offY = mFrustumOffset[1] / mFocalLength;
        left = -halfW + offX;  right  = halfW + offX;
        top  = -halfH + offY;  bottom = halfH + offY;
    } else {
        const double halfH = mOrthoHeight * 0.5, halfW = halfH * mAspect;
        left = -halfW;  right = halfW;  top = -halfH;  bottom = halfH;
    }
}

void Frustum::RebuildMatrices() const
{
    const bool persp = mProjType == ProjectionType::Perspective;
    if (persp)   // все лучи идут через камеру z = 0: ближняя плоскость обязана лежать перед ней
        assert(mNearDist > 0.0 && mFarDist > mNearDist && "перспектива: нужно 0 < n < f");
    else         // делим только на f - n; бесконечной ортографии не бывает (глубина линейна)
        assert(!IsInfinite() && mFarDist > mNearDist && "ортография: нужно n < f < inf");

    double l, r, t, b;
    CalcProjectionParameters(l, r, t, b);
    const DepthParams dp = GetProjectionParamsAB();
    const double sy = ClipSignY(mClipY);

    // Окно [l,r] x [t,b] -> [-1,1] x [-1,1]: одни и те же масштаб и сдвиг для обоих типов.
    // Перспектива: окно в тангенсах, x/z в [l,r]; симметричное -> a = g/s, by = ±g, c = dy = 0
    // (это (6.36)); сдвиг c, dy стоит в столбце z, потому что умножается на z, а потом /w = z.
    // Ортография: ящик (6.48); сдвиг стоит в столбце w = 1.
    const double iw = 1.0 / (r - l), ih = 1.0 / (b - t);
    const double a  = 2.0 * iw,      c  = -(r + l) * iw;
    const double by = sy * 2.0 * ih, dy = -sy * (b + t) * ih;

    if (persp) {
        mProjMatrix    = mathx::Matrix4x4(a,   0.0, c,    0.0,
                                   0.0, by,  dy,   0.0,
                                   0.0, 0.0, dp.A, dp.B,
                                   0.0, 0.0, 1.0,  0.0);           // w = z
        mInvProjMatrix = mathx::Matrix4x4(1.0 / a, 0.0,      0.0,        -c / a,
                                   0.0,     1.0 / by, 0.0,        -dy / by,
                                   0.0,     0.0,      0.0,        1.0,
                                   0.0,     0.0,      1.0 / dp.B, -dp.A / dp.B);
    } else {
        mProjMatrix    = mathx::Matrix4x4(a,   0.0, 0.0,  c,
                                   0.0, by,  0.0,  dy,
                                   0.0, 0.0, dp.B, dp.A,
                                   0.0, 0.0, 0.0,  1.0);           // w = 1
        mInvProjMatrix = mathx::Matrix4x4(1.0 / a, 0.0,      0.0,        -c / a,
                                   0.0,     1.0 / by, 0.0,        -dy / by,
                                   0.0,     0.0,      1.0 / dp.B, -dp.A / dp.B,
                                   0.0,     0.0,      0.0,        1.0);
    }
    mRecalcMatrices = false;
    mRecalcPlanes = true;
}

void Frustum::UpdatePlanes() const
{
    mFrustumPlanes = ExtractPlanes(mProjMatrix, mClipY, mDepthMode, IsInfinite());
    mRecalcPlanes = false;
}
