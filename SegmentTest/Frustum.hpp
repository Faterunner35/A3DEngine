#pragma once
#include "Math.h"          // Radian, mathx::Tan
#include "Projection.hpp"  // Matrix4x4, M4Idx, ClipY
#include "Plane.h"
#include <array>
#include <cassert>
#include <cmath>
#include <limits>

namespace mathx {

// Пирамида видимости в пространстве камеры: x вправо, y ВНИЗ, z ВПЕРЁД (FGED).
// Клип D3D12: -w <= x,y <= w, 0 <= z <= w. Положение камеры в мире — в Camera.
class Frustum
{
public:
    enum class ProjectionType { Perspective, Orthographic };
    enum class DepthMode { Standard,    // near -> 0, far -> 1
                           Reversed };  // near -> 1, far -> 0
    enum PlaneId { kLeft = 0, kRight, kTop, kBottom, kNear, kFar, kPlaneCount };
    struct DepthParams { double A, B; }; // перспектива: z_dev = A + B/z; ортография: A + B*z

    explicit Frustum(ClipY clipY = ClipY::Up, DepthMode depth = DepthMode::Reversed)
        : mClipY(clipY), mDepthMode(depth) {}

    // ---- форма: сеттеры только помечают «устарело» ----------------------------
    void SetFOVy(Radian fovY)              { mFovy = fovY; InvalidateFrustum(); }
    void SetAspectRatio(double aspect)     { mAspect = aspect; InvalidateFrustum(); }
    void SetNearClipDistance(double n)     { assert(n > 0.0); mNearDist = n; InvalidateFrustum(); }
    void SetFarClipDistance(double f)      { mFarDist = f; InvalidateFrustum(); }   // +inf — бесконечная
    void SetOrthoWindowHeight(double h)    { mOrthoHeight = h; InvalidateFrustum(); }
    void SetFrustumOffset(const Vector<2, double>& o) { mFrustumOffset = o; InvalidateFrustum(); }
    void SetFocalLength(double focal)      { assert(focal > 0.0); mFocalLength = focal; InvalidateFrustum(); }
    void SetProjectionType(ProjectionType t) { mProjType = t; InvalidateFrustum(); }
    void SetClipY(ClipY c)                 { mClipY = c; InvalidateFrustum(); }
    void SetDepthMode(DepthMode d)         { mDepthMode = d; InvalidateFrustum(); }

    Radian GetFOVy() const              { return mFovy; }
    double GetAspectRatio() const       { return mAspect; }
    double GetNearClipDistance() const  { return mNearDist; }
    double GetFarClipDistance() const   { return mFarDist; }
    bool   IsInfinite() const           { return std::isinf(mFarDist); }

    // ---- результаты: геттер досчитывает только устаревшее ---------------------
    const Matrix4x4& GetProjectionMatrix() const        { UpdateFrustum(); return mProjMatrix; }
    const Matrix4x4& GetInverseProjectionMatrix() const { UpdateFrustum(); return mInvProjMatrix; }
    const std::array<Plane, kPlaneCount>& GetFrustumPlanes() const
    {
        UpdateFrustum();
        if (mRecalcFrustumPlanes) UpdateFrustumPlanesImpl();
        return mFrustumPlanes;
    }

    DepthParams GetProjectionParamsAB() const
    {
        const bool stdDepth = mDepthMode == DepthMode::Standard;
        const double n = mNearDist, f = mFarDist;
        if (mProjType == ProjectionType::Orthographic) {
            const double k = 1.0 / (f - n);
            return stdDepth ? DepthParams{ -n * k, k } : DepthParams{ f * k, -k };
        }
        if (IsInfinite())
            return stdDepth ? DepthParams{ 1.0 - kInfiniteEps, -n * (1.0 - kInfiniteEps) } : DepthParams{ 0.0, n };
        const double k = 1.0 / (f - n);
        return stdDepth ? DepthParams{ f * k, -n * f * k } : DepthParams{ -n * k, n * f * k };
    }

    // Отсечение в пространстве камеры (в мире — через Camera).
    bool IsVisible(const Vector<3, double>& p) const
    {
        for (const Plane& pl : GetFrustumPlanes()) if (pl.SignedDistance(p) < 0.0) return false;
        return true;
    }
    bool IsVisible(const Vector<3, double>& center, double radius) const
    {
        for (const Plane& pl : GetFrustumPlanes()) if (pl.SignedDistance(center) < -radius) return false;
        return true;
    }

    // ---- твои функции: симметричная пирамида без сдвига -----------------------
    static Matrix4x4 CreatePerspectiveMatrix(Radian fovY, double aspect, double n, double f,
                                             ClipY clipY = ClipY::Up)
    {
        assert(fovY.ValueRadians() > 0.0 && aspect > 0.0 && n > 0.0 && f > n);
        const double g = 1.0 / mathx::Tan(fovY * 0.5);
        const double gy = (clipY == ClipY::Up) ? -g : g;
        const double k = f / (f - n);
        Matrix4x4 P;
        for (int32_t i = 0; i < 16; ++i) P[i] = 0.0;      // конструктор даёт единичную
        P[M4Idx(0, 0)] = g / aspect;
        P[M4Idx(1, 1)] = gy;
        P[M4Idx(2, 2)] = k;
        P[M4Idx(2, 3)] = -n * k;
        P[M4Idx(3, 2)] = 1.0;                              // w = z
        return P;
    }

    static Matrix4x4 MakeFrustumProjectionInverse(Radian fovY, double s, double n, double f,
                                                  ClipY clipY = ClipY::Up)
    {
        const double g = 1.0 / mathx::Tan(fovY * 0.5);
        const double gy = (clipY == ClipY::Up) ? -g : g;
        return Matrix4x4(s / g, 0.0, 0.0, 0.0,
                         0.0, 1.0 / gy, 0.0, 0.0,
                         0.0, 0.0, 0.0, 1.0,
                         0.0, 0.0, (n - f) / (n * f), 1.0 / n);
    }

    static Matrix4x4 MakeReversedInfiniteProjectionInverse(Radian fovY, double s, double n,
                                                           ClipY clipY = ClipY::Up)
    {
        const double g = 1.0 / mathx::Tan(fovY * 0.5);
        const double gy = (clipY == ClipY::Up) ? -g : g;
        return Matrix4x4(s / g, 0.0, 0.0, 0.0,
                         0.0, 1.0 / gy, 0.0, 0.0,
                         0.0, 0.0, 0.0, 1.0,
                         0.0, 0.0, 1.0 / n, 0.0);
    }

protected:
    // Твоя UpdateProJectionParam: окно на ближней плоскости. y камеры вниз,
    // поэтому top < bottom.
    void CalcProjectionParameters(double& left, double& right, double& top, double& bottom) const
    {
        if (mProjType == ProjectionType::Perspective) {
            const double tanY = mathx::Tan(mFovy * 0.5);
            const double tanX = tanY * mAspect;
            const double nearFocal = mNearDist / mFocalLength;
            const double offsetX = mFrustumOffset[0] * nearFocal;
            const double offsetY = mFrustumOffset[1] * nearFocal;
            const double halfW = tanX * mNearDist;
            const double halfH = tanY * mNearDist;
            left = -halfW + offsetX;   right = halfW + offsetX;
            top = -halfH + offsetY;    bottom = halfH + offsetY;
        } else {
            const double halfH = mOrthoHeight * 0.5, halfW = halfH * mAspect;
            left = -halfW; right = halfW; top = -halfH; bottom = halfH;
        }
    }

    void UpdateFrustum() const { if (mRecalcFrustum) UpdateFrustumImpl(); }

    void UpdateFrustumImpl() const
    {
        double l, r, t, b;
        CalcProjectionParameters(l, r, t, b);
        const DepthParams dp = GetProjectionParamsAB();
        const double sy = (mClipY == ClipY::Up) ? -1.0 : 1.0;   // D3D12: y клипа вверх

        if (mProjType == ProjectionType::Perspective) {
            const double n = mNearDist;
            const double a = 2.0 * n / (r - l), c = -(r + l) / (r - l);
            const double by = sy * 2.0 * n / (b - t), dy = -sy * (b + t) / (b - t);
            mProjMatrix = Matrix4x4(a, 0.0, c, 0.0,
                                    0.0, by, dy, 0.0,
                                    0.0, 0.0, dp.A, dp.B,
                                    0.0, 0.0, 1.0, 0.0);
            mInvProjMatrix = Matrix4x4(1.0 / a, 0.0, 0.0, -c / a,
                                       0.0, 1.0 / by, 0.0, -dy / by,
                                       0.0, 0.0, 0.0, 1.0,
                                       0.0, 0.0, 1.0 / dp.B, -dp.A / dp.B);
        } else {
            const double a = 2.0 / (r - l), c = -(r + l) / (r - l);
            const double by = sy * 2.0 / (b - t), dy = -sy * (b + t) / (b - t);
            mProjMatrix = Matrix4x4(a, 0.0, 0.0, c,
                                    0.0, by, 0.0, dy,
                                    0.0, 0.0, dp.B, dp.A,
                                    0.0, 0.0, 0.0, 1.0);
            mInvProjMatrix = Matrix4x4(1.0 / a, 0.0, 0.0, -c / a,
                                       0.0, 1.0 / by, 0.0, -dy / by,
                                       0.0, 0.0, 1.0 / dp.B, -dp.A / dp.B,
                                       0.0, 0.0, 0.0, 1.0);
        }
        mRecalcFrustum = false;
        mRecalcFrustumPlanes = true;
    }

    // Плоскости из строк P (камера), нормали внутрь.
    void UpdateFrustumPlanesImpl() const
    {
        auto row = [&](int32_t i, double k) {
            return Vector<4, double>{ k * mProjMatrix(i, 0), k * mProjMatrix(i, 1),
                                      k * mProjMatrix(i, 2), k * mProjMatrix(i, 3) };
        };
        auto add = [](const Vector<4, double>& u, const Vector<4, double>& v, double s) {
            return Vector<4, double>{ u[0] + s * v[0], u[1] + s * v[1], u[2] + s * v[2], u[3] + s * v[3] };
        };
        auto plane = [](const Vector<4, double>& v) {
            const double len = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            return Plane(Vector<3, double>{ v[0] / len, v[1] / len, v[2] / len }, v[3] / len);
        };
        const Vector<4, double> r0 = row(0, 1.0), r2 = row(2, 1.0), r3 = row(3, 1.0);
        const Vector<4, double> rY = row(1, mClipY == ClipY::Up ? -1.0 : 1.0);   // к «клип вниз»
        mFrustumPlanes[kLeft]   = plane(add(r3, r0, +1.0));
        mFrustumPlanes[kRight]  = plane(add(r3, r0, -1.0));
        mFrustumPlanes[kTop]    = plane(add(r3, rY, +1.0));
        mFrustumPlanes[kBottom] = plane(add(r3, rY, -1.0));
        const bool stdDepth = mDepthMode == DepthMode::Standard;
        const Vector<4, double> zLow = r2, zHigh = add(r3, r2, -1.0);   // z >= 0, z <= w
        mFrustumPlanes[kNear] = plane(stdDepth ? zLow : zHigh);
        mFrustumPlanes[kFar] = IsInfinite()
            ? Plane(Vector<3, double>{ 0.0, 0.0, -1.0 }, std::numeric_limits<double>::infinity())
            : plane(stdDepth ? zHigh : zLow);
        mRecalcFrustumPlanes = false;
    }

    void InvalidateFrustum() { mRecalcFrustum = true; mRecalcFrustumPlanes = true; }

    // ---- поля -------------------------------------------------------------------
    static constexpr double kInfiniteEps = 9.5367431640625e-07;   // 2^-20, Lengyel (6.42)

    ProjectionType     mProjType = ProjectionType::Perspective;
    Radian             mFovy{ ::detail::PI / 3.0 };
    double             mFarDist = 100000.0;
    double             mNearDist = 0.1;
    double             mAspect = 16.0 / 9.0;
    double             mOrthoHeight = 1000.0;
    Vector<2, double>  mFrustumOffset{ 0.0, 0.0 };   // в единицах камеры, y вниз
    double             mFocalLength = 1.0;
    ClipY              mClipY;
    DepthMode          mDepthMode;

    mutable Matrix4x4  mProjMatrix;
    mutable Matrix4x4  mInvProjMatrix;
    mutable std::array<Plane, kPlaneCount> mFrustumPlanes;
    mutable bool       mRecalcFrustum = true;
    mutable bool       mRecalcFrustumPlanes = true;
};

} // namespace mathx
