#pragma once
#include "Math.h"          // Radian, mathx::Tan, ::detail::PI
#include "Matrxi4x4.hpp"   // Matrix4x4
#include "Plane.h"
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>

namespace mathx {

// ============================================================================
//  Пирамида видимости — «объектив» камеры. Живёт только в пространстве камеры
//  FGED (Lengyel, т. 2, гл. 6): x вправо, y ВНИЗ, z ВПЕРЁД (правая тройка).
//  Положение камеры в мире и матрица вида V — в классе Camera.
//
//  Клип: -w <= x, y <= w,  0 <= z <= w  (глубина [0,1], D3D12 / Vulkan).
//  Матрицы хранятся по строкам, точка — столбец: p_clip = P * p_camera.
//
//  Кэш: сеттеры только помечают «устарело», геттеры досчитывают нужное.
//  P и P^-1 пересчитываются при смене параметров объектива, плоскости — лениво.
// ============================================================================
class Frustum
{
public:
    // ---- соглашения -----------------------------------------------------------
    enum class ClipY {
        Down,   // Vulkan, OpenGL + glClipControl(GL_UPPER_LEFT)
        Up      // D3D12, Metal: вторая строка P со знаком минус
    };
    enum class ProjectionType { Perspective, Orthographic };
    enum class DepthMode {
        Standard,   // near -> 0, far -> 1                     (6.36), тест LESS
        Reversed    // near -> 1, far -> 0  (6.45), (6.46): тест GREATER_EQUAL, очистка 0
    };
    enum PlaneId : int32_t { kLeft = 0, kRight, kTop, kBottom, kNear, kFar, kPlaneCount };
    struct DepthParams { double A, B; };   // перспектива: z_dev = A + B/z; ортография: z_dev = A + B*z

    explicit Frustum(ClipY clipY = ClipY::Up, DepthMode depth = DepthMode::Reversed) noexcept
        : mClipY(clipY), mDepthMode(depth) {}

    // ---- параметры объектива: сеттер только помечает «устарело» -------------------
    void SetFOVy(Radian fovY)                         { mFovy = fovY; Invalidate(); }
    void SetAspectRatio(double aspect)                { assert(aspect > 0.0); mAspect = aspect; Invalidate(); }
    // Перспектива требует 0 < n < f, ортография — любое n < f (ящик может начинаться
    // за камерой). Тип проекции можно сменить после n, поэтому проверка — в UpdateMatrices.
    void SetNearClipDistance(double n)                { mNearDist = n; Invalidate(); }
    void SetFarClipDistance(double f)                 { mFarDist = f; Invalidate(); }   // +inf — бесконечная
    void SetOrthoWindowHeight(double h)               { assert(h > 0.0); mOrthoHeight = h; Invalidate(); }
    void SetFrustumOffset(const Vector<2, double>& o) { mFrustumOffset = o; Invalidate(); }   // стерео, y вниз
    void SetFocalLength(double focal)                 { assert(focal > 0.0); mFocalLength = focal; Invalidate(); }
    void SetProjectionType(ProjectionType t)          { mProjType = t; Invalidate(); }
    void SetClipY(ClipY c)                            { mClipY = c; Invalidate(); }
    void SetDepthMode(DepthMode d)                    { mDepthMode = d; Invalidate(); }

    Radian         GetFOVy() const             { return mFovy; }
    double         GetAspectRatio() const      { return mAspect; }
    double         GetNearClipDistance() const { return mNearDist; }
    double         GetFarClipDistance() const  { return mFarDist; }
    double         GetOrthoWindowHeight() const{ return mOrthoHeight; }
    ProjectionType GetProjectionType() const   { return mProjType; }
    ClipY          GetClipY() const            { return mClipY; }
    DepthMode      GetDepthMode() const        { return mDepthMode; }
    bool           IsInfinite() const          { return std::isinf(mFarDist); }

    // ---- результаты: геттер досчитывает только устаревшее ----------------------
    const Matrix4x4& GetProjectionMatrix() const        { UpdateMatrices(); return mProjMatrix; }
    const Matrix4x4& GetInverseProjectionMatrix() const { UpdateMatrices(); return mInvProjMatrix; }
    const std::array<Plane, kPlaneCount>& GetFrustumPlanes() const   // нормали внутрь
    {
        UpdateMatrices();
        if (mRecalcPlanes) UpdatePlanes();
        return mFrustumPlanes;
    }

    // Строка глубины P для текущего режима.
    DepthParams GetProjectionParamsAB() const
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

    // ---- отсечение в пространстве камеры (в мире — через Camera) ----------------
    bool IsVisible(const Vector<3, double>& p) const
    {
        for (const Plane& pl : GetFrustumPlanes()) if (pl.SignedDistance(p) < 0.0) return false;
        return true;
    }
    bool IsVisible(const Vector<3, double>& center, double radius) const   // сфера
    {
        for (const Plane& pl : GetFrustumPlanes()) if (pl.SignedDistance(center) < -radius) return false;
        return true;
    }

    // ---- готовые матрицы без объекта (листинги FGED): симметричная пирамида ------
    // (6.36), листинг 6.1: near -> 0, far -> 1.   z_dev = A + B/z,  A = f/(f-n),  B = -nf/(f-n)
    static Matrix4x4 MakeFrustumProjection(Radian fovY, double s, double n, double f,
                                           ClipY clipY = ClipY::Up)
    {
        assert(s > 0.0 && n > 0.0 && f > n && "MakeFrustumProjection: нужны s > 0, 0 < n < f");
        const double g = FocalG(fovY), k = f / (f - n);
        return Matrix4x4(g / s, 0.0,              0.0, 0.0,
                         0.0,   ClipSignY(clipY) * g, 0.0, 0.0,
                         0.0,   0.0,              k,   -n * k,
                         0.0,   0.0,              1.0, 0.0);   // w = z
    }
    // (6.37): обратная к (6.36)
    static Matrix4x4 MakeFrustumProjectionInverse(Radian fovY, double s, double n, double f,
                                                  ClipY clipY = ClipY::Up)
    {
        assert(s > 0.0 && n > 0.0 && f > n);
        const double g = FocalG(fovY);
        return Matrix4x4(s / g, 0.0,                        0.0,               0.0,
                         0.0,   ClipSignY(clipY) / g,       0.0,               0.0,
                         0.0,   0.0,                        0.0,               1.0,
                         0.0,   0.0,                        (n - f) / (n * f), 1.0 / n);
    }
    // (6.46): обратная бесконечная, near -> 1, бесконечность -> 0.  z_dev = n/z.
    // Для направления (w = 0) строка (0,0,0,n) даёт z_clip = 0 ровно — epsilon не нужен
    // при тесте GREATER_EQUAL и очистке 0. Формат D32_FLOAT.
    static Matrix4x4 MakeReversedInfiniteProjection(Radian fovY, double s, double n,
                                                    ClipY clipY = ClipY::Up)
    {
        assert(s > 0.0 && n > 0.0 && "MakeReversedInfiniteProjection: нужны s > 0, n > 0");
        const double g = FocalG(fovY);
        return Matrix4x4(g / s, 0.0,              0.0, 0.0,
                         0.0,   ClipSignY(clipY) * g, 0.0, 0.0,
                         0.0,   0.0,              0.0, n,
                         0.0,   0.0,              1.0, 0.0);
    }
    static Matrix4x4 MakeReversedInfiniteProjectionInverse(Radian fovY, double s, double n,
                                                           ClipY clipY = ClipY::Up)
    {
        assert(s > 0.0 && n > 0.0);
        const double g = FocalG(fovY);
        return Matrix4x4(s / g, 0.0,                  0.0,     0.0,
                         0.0,   ClipSignY(clipY) / g, 0.0,     0.0,
                         0.0,   0.0,                  0.0,     1.0,
                         0.0,   0.0,                  1.0 / n, 0.0);
    }
    // (6.38): глубина камеры по значению буфера для (6.36)
    static double CameraDepthFromFrustum(double zDevice, double n, double f)
    {
        return n * f / ((n - f) * zDevice + f);
    }
    // Обращение z_dev = n/z для (6.46). zDevice = 0 — точка на бесконечности.
    static double CameraDepthFromReversedInfinite(double zDevice, double n)
    {
        return n / zDevice;
    }

private:
    static double FocalG(Radian fovY)                     // g = 1 / tan(fovY/2), вид сбоку
    {
        assert(fovY.ValueRadians() > 0.0 && fovY.ValueRadians() < ::detail::PI);
        return 1.0 / mathx::Tan(fovY * 0.5);
    }
    static double ClipSignY(ClipY c) { return c == ClipY::Up ? -1.0 : 1.0; }   // D3D12: y клипа вверх

    void Invalidate() { mRecalcMatrices = true; mRecalcPlanes = true; }

    // Окно на ближней плоскости (перспектива) или ящик по x, y (ортография). y вниз: top < bottom.
    void CalcProjectionParameters(double& left, double& right, double& top, double& bottom) const
    {
        if (mProjType == ProjectionType::Perspective) {
            const double halfH = mathx::Tan(mFovy * 0.5) * mNearDist;
            const double halfW = halfH * mAspect;
            const double nearFocal = mNearDist / mFocalLength;   // сдвиг задан в фокальной плоскости
            const double offX = mFrustumOffset[0] * nearFocal, offY = mFrustumOffset[1] * nearFocal;
            left = -halfW + offX;  right  = halfW + offX;
            top  = -halfH + offY;  bottom = halfH + offY;
        } else {
            const double halfH = mOrthoHeight * 0.5, halfW = halfH * mAspect;
            left = -halfW;  right = halfW;  top = -halfH;  bottom = halfH;
        }
    }

    void UpdateMatrices() const { if (mRecalcMatrices) RebuildMatrices(); }

    void RebuildMatrices() const
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

        // Общее для обоих типов: перспектива масштабирует окно на ближней плоскости
        // (множитель n, сдвиг уходит в столбец z), ортография — ящик (множитель 1,
        // сдвиг в столбец w).
        const double scale = persp ? mNearDist : 1.0;
        const double iw = 1.0 / (r - l), ih = 1.0 / (b - t);
        const double a  = 2.0 * scale * iw,      c  = -(r + l) * iw;
        const double by = sy * 2.0 * scale * ih, dy = -sy * (b + t) * ih;

        if (persp) {
            mProjMatrix    = Matrix4x4(a,   0.0, c,    0.0,
                                       0.0, by,  dy,   0.0,
                                       0.0, 0.0, dp.A, dp.B,
                                       0.0, 0.0, 1.0,  0.0);           // w = z
            mInvProjMatrix = Matrix4x4(1.0 / a, 0.0,      0.0,        -c / a,
                                       0.0,     1.0 / by, 0.0,        -dy / by,
                                       0.0,     0.0,      0.0,        1.0,
                                       0.0,     0.0,      1.0 / dp.B, -dp.A / dp.B);
        } else {
            mProjMatrix    = Matrix4x4(a,   0.0, 0.0,  c,
                                       0.0, by,  0.0,  dy,
                                       0.0, 0.0, dp.B, dp.A,
                                       0.0, 0.0, 0.0,  1.0);           // w = 1
            mInvProjMatrix = Matrix4x4(1.0 / a, 0.0,      0.0,        -c / a,
                                       0.0,     1.0 / by, 0.0,        -dy / by,
                                       0.0,     0.0,      1.0 / dp.B, -dp.A / dp.B,
                                       0.0,     0.0,      0.0,        1.0);
        }
        mRecalcMatrices = false;
        mRecalcPlanes = true;
    }

    // Плоскости из строк P (Грибб–Хартманн), в пространстве камеры, нормали внутрь.
    void UpdatePlanes() const
    {
        using V4 = Vector<4, double>;
        const Matrix4x4& P = mProjMatrix;
        const double sy = ClipSignY(mClipY);   // приводим строку y к «клип вниз»
        const V4 r0{ P(0, 0), P(0, 1), P(0, 2), P(0, 3) };
        const V4 rY{ sy * P(1, 0), sy * P(1, 1), sy * P(1, 2), sy * P(1, 3) };
        const V4 r2{ P(2, 0), P(2, 1), P(2, 2), P(2, 3) };
        const V4 r3{ P(3, 0), P(3, 1), P(3, 2), P(3, 3) };
        auto comb = [](const V4& u, const V4& v, double s) {
            return V4{ u[0] + s * v[0], u[1] + s * v[1], u[2] + s * v[2], u[3] + s * v[3] };
        };
        auto plane = [](const V4& v) {
            const double inv = 1.0 / std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            return Plane(Vector<3, double>{ v[0] * inv, v[1] * inv, v[2] * inv }, v[3] * inv);
        };
        mFrustumPlanes[kLeft]   = plane(comb(r3, r0, +1.0));   // x >= -w
        mFrustumPlanes[kRight]  = plane(comb(r3, r0, -1.0));   // x <=  w
        mFrustumPlanes[kTop]    = plane(comb(r3, rY, +1.0));   // y >= -w
        mFrustumPlanes[kBottom] = plane(comb(r3, rY, -1.0));   // y <=  w
        const bool stdDepth = mDepthMode == DepthMode::Standard;
        const V4 zLow = r2, zHigh = comb(r3, r2, -1.0);       // z >= 0,  z <= w
        mFrustumPlanes[kNear] = plane(stdDepth ? zLow : zHigh);
        mFrustumPlanes[kFar]  = IsInfinite()
            ? Plane(Vector<3, double>{ 0.0, 0.0, -1.0 }, std::numeric_limits<double>::infinity())
            : plane(stdDepth ? zHigh : zLow);
        mRecalcPlanes = false;
    }

    // ---- поля -------------------------------------------------------------------
    static constexpr double kInfiniteEps = 9.5367431640625e-07;   // 2^-20, Lengyel (6.42)

    ProjectionType     mProjType = ProjectionType::Perspective;
    Radian             mFovy{ ::detail::PI / 3.0 };
    double             mFarDist = 100000.0;
    double             mNearDist = 0.1;
    double             mAspect = 16.0 / 9.0;
    double             mOrthoHeight = 1000.0;
    Vector<2, double>  mFrustumOffset{ 0.0, 0.0 };   // в фокальной плоскости, y вниз
    double             mFocalLength = 1.0;
    ClipY              mClipY;
    DepthMode          mDepthMode;

    mutable Matrix4x4  mProjMatrix;
    mutable Matrix4x4  mInvProjMatrix;
    mutable std::array<Plane, kPlaneCount> mFrustumPlanes;
    mutable bool       mRecalcMatrices = true;
    mutable bool       mRecalcPlanes = true;
};

} // namespace mathx
