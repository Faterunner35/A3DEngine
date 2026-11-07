#pragma once
#include "Geomerty.h"      // Vector
#include "Math.h"          // Radian, Tan
#include "Matrxi4x4.h"   // Matrix4x4
#include "Plane.h"         // Plane
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>

// ============================================================================
//  Пирамида видимости — «объектив» камеры. Живёт в пространстве камеры
//  FGED (Lengyel, т. 2, гл. 6): x вправо, y ВНИЗ, z ВПЕРЁД (правая тройка).
//  Положение камеры в мире и матрица вида V — в классе Camera.
//  Плоскости в мире / объекте — из P·V / P·V·M (ExtractPlanes, FGED §6.3.5).
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
    enum PlaneId : std::int32_t { kLeft = 0, kRight, kTop, kBottom, kNear, kFar, kPlaneCount };

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
    void SetFrustumOffset(const mathx::Vector<2, double>& o) { mFrustumOffset = o; Invalidate(); }   // стерео, y вниз
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
    const mathx::Matrix4x4& GetProjectionMatrix() const        { UpdateMatrices(); return mProjMatrix; }
    const mathx::Matrix4x4& GetInverseProjectionMatrix() const { UpdateMatrices(); return mInvProjMatrix; }
    const std::array<mathx::Plane, kPlaneCount>& GetFrustumPlanes() const   // камера, нормали внутрь
    {
        UpdateMatrices();
        if (mRecalcPlanes) UpdatePlanes();
        return mFrustumPlanes;
    }
    // Плоскости в пространстве X, где toCamera: X -> камера (аффинная):
    //   toCamera = V -> мир,  V·M -> объект.  Не кэшируются.
    std::array<mathx::Plane, kPlaneCount> GetFrustumPlanes(const mathx::Matrix4x4& toCamera) const
    {
        return ExtractPlanes(GetProjectionMatrix() * toCamera, mClipY, mDepthMode, IsInfinite());
    }

    // Глубина камеры z по значению из буфера глубины — для любого режима
    // (Standard/Reversed, конечная/бесконечная, перспектива/ортография).
    //   перспектива z_dev = A + B/z  ->  z = B / (z_dev - A)          (6.38) в общем виде
    //   ортография  z_dev = A + B*z  ->  z = (z_dev - A) / B
    // Reversed + бесконечная: z_dev = 0 даёт +inf (точка на бесконечности).
    double CameraDepth(double zDevice) const
    {
        const DepthParams dp = GetProjectionParamsAB();
        return mProjType == ProjectionType::Perspective ? dp.B / (zDevice - dp.A)
                                                        : (zDevice - dp.A) / dp.B;
    }

    // ---- отсечение ----------------------------------------------------------------
    // Сфера (center, radius) хотя бы частично внутри; radius = 0 — точка.
    // center — в том же пространстве, что и плоскости. Плоскости нормированы:
    // для мира/объекта это верно и при V, M с масштабом (нормировка — после A).
    static bool IsSphereVisible(const std::array<mathx::Plane, kPlaneCount>& planes,
                                const mathx::Vector<3, double>& center, double radius = 0.0)
    {
        for (const mathx::Plane& pl : planes) if (pl.SignedDistance(center) < -radius) return false;
        return true;
    }
    // В пространстве камеры.
    bool IsVisible(const mathx::Vector<3, double>& center, double radius = 0.0) const
    {
        return IsSphereVisible(GetFrustumPlanes(), center, radius);
    }

    // ---- плоскости из ЛЮБОЙ матрицы X -> клип (Грибб–Хартманн, FGED §6.3.5, табл. 6.2) ----
    // Точка p в X внутри  <=>  c = A·p в клипе: -w <= x, y <= w, 0 <= z <= w.
    // Каждое неравенство линейно по строкам A: x >= -w  <=>  (r3 + r0)·p >= 0 —
    // это и есть плоскость в X. A = P -> камера, P·V -> мир, P·V·M -> объект.
    // Нормали внутрь, плоскости нормированы.
    //
    // infiniteFar: строка дальней плоскости вырождена — (0,0,0,n) у Reversed, почти
    // такая же у Standard (6.42). Это плоскость на бесконечности w = 0; аффинное
    // отображение её не двигает (последняя строка (0,0,0,1)), поэтому ставим
    // «всегда проходит» Plane((0,0,-1), +inf) в любом пространстве: расстояние +inf
    // от любой конечной точки, нормаль не важна. Для A = P·(аффинная) — верно.
    static std::array<mathx::Plane, kPlaneCount> ExtractPlanes(const mathx::Matrix4x4& A, ClipY clipY,
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

    // ---- готовые матрицы без объекта (листинги FGED): симметричная пирамида ------
    // (6.36), листинг 6.1: near -> 0, far -> 1.   z_dev = A + B/z,  A = f/(f-n),  B = -nf/(f-n)
    static mathx::Matrix4x4 MakeFrustumProjection(Radian fovY, double s, double n, double f,
                                           ClipY clipY = ClipY::Up)
    {
        assert(s > 0.0 && n > 0.0 && f > n && "MakeFrustumProjection: нужны s > 0, 0 < n < f");
        const double g = FocalG(fovY), k = f / (f - n);
        return mathx::Matrix4x4(g / s, 0.0,              0.0, 0.0,
                         0.0,   ClipSignY(clipY) * g, 0.0, 0.0,
                         0.0,   0.0,              k,   -n * k,
                         0.0,   0.0,              1.0, 0.0);   // w = z
    }
    // (6.37): обратная к (6.36)
    static mathx::Matrix4x4 MakeFrustumProjectionInverse(Radian fovY, double s, double n, double f,
                                                  ClipY clipY = ClipY::Up)
    {
        assert(s > 0.0 && n > 0.0 && f > n);
        const double g = FocalG(fovY);
        return mathx::Matrix4x4(s / g, 0.0,                        0.0,               0.0,
                         0.0,   ClipSignY(clipY) / g,       0.0,               0.0,
                         0.0,   0.0,                        0.0,               1.0,
                         0.0,   0.0,                        (n - f) / (n * f), 1.0 / n);
    }
    // (6.46): обратная бесконечная, near -> 1, бесконечность -> 0.  z_dev = n/z.
    // Для направления (w = 0) строка (0,0,0,n) даёт z_clip = 0 ровно — epsilon не нужен
    // при тесте GREATER_EQUAL и очистке 0. Формат D32_FLOAT.
    static mathx::Matrix4x4 MakeReversedInfiniteProjection(Radian fovY, double s, double n,
                                                    ClipY clipY = ClipY::Up)
    {
        assert(s > 0.0 && n > 0.0 && "MakeReversedInfiniteProjection: нужны s > 0, n > 0");
        const double g = FocalG(fovY);
        return mathx::Matrix4x4(g / s, 0.0,              0.0, 0.0,
                         0.0,   ClipSignY(clipY) * g, 0.0, 0.0,
                         0.0,   0.0,              0.0, n,
                         0.0,   0.0,              1.0, 0.0);
    }
    static mathx::Matrix4x4 MakeReversedInfiniteProjectionInverse(Radian fovY, double s, double n,
                                                           ClipY clipY = ClipY::Up)
    {
        assert(s > 0.0 && n > 0.0);
        const double g = FocalG(fovY);
        return mathx::Matrix4x4(s / g, 0.0,                  0.0,     0.0,
                         0.0,   ClipSignY(clipY) / g, 0.0,     0.0,
                         0.0,   0.0,                  0.0,     1.0,
                         0.0,   0.0,                  1.0 / n, 0.0);
    }

private:
    static double FocalG(Radian fovY)                     // g = 1 / tan(fovY/2), вид сбоку
    {
        assert(fovY.ValueRadians() > 0.0 && fovY.ValueRadians() < kPi);
        return 1.0 / mathx::Tan(fovY * 0.5);
    }
    static double ClipSignY(ClipY c) { return c == ClipY::Up ? -1.0 : 1.0; }   // D3D12: y клипа вверх

    void Invalidate() { mRecalcMatrices = true; mRecalcPlanes = true; }

    // Строка глубины P для текущего режима: перспектива z_dev = A + B/z, ортография z_dev = A + B*z.
    struct DepthParams { double A, B; };
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

    // Окно кадра, y вниз (top < bottom).
    //   перспектива: окно на расстоянии 1 от камеры, то есть ТАНГЕНСЫ углов до краёв
    //                (не зависит от n). Симметричное окно: halfH = 1/g, halfW = s/g.
    //   ортография:  ящик по x, y в единицах камеры.
    void CalcProjectionParameters(double& left, double& right, double& top, double& bottom) const
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

    // Плоскости в пространстве камеры: A = P.
    void UpdatePlanes() const
    {
        mFrustumPlanes = ExtractPlanes(mProjMatrix, mClipY, mDepthMode, IsInfinite());
        mRecalcPlanes = false;
    }

    // ---- поля -------------------------------------------------------------------
    static constexpr double kInfiniteEps = 9.5367431640625e-07;   // 2^-20, Lengyel (6.42)
    static constexpr double kPi = 3.14159265358979323846;

    ProjectionType     mProjType = ProjectionType::Perspective;
    Radian             mFovy{ kPi / 3.0 };
    double             mFarDist = 100000.0;
    double             mNearDist = 0.1;
    double             mAspect = 16.0 / 9.0;
    double             mOrthoHeight = 1000.0;
    mathx::Vector<2, double>  mFrustumOffset{ 0.0, 0.0 };   // в фокальной плоскости, y вниз
    double             mFocalLength = 1.0;
    ClipY              mClipY;
    DepthMode          mDepthMode;

    mutable mathx::Matrix4x4  mProjMatrix;
    mutable mathx::Matrix4x4  mInvProjMatrix;
    mutable std::array<mathx::Plane, kPlaneCount> mFrustumPlanes;
    mutable bool       mRecalcMatrices = true;
    mutable bool       mRecalcPlanes = true;
};

