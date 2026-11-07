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
    double CameraDepth(double zDevice) const;

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
                                                               DepthMode depth, bool infiniteFar);

    // ---- готовые матрицы без объекта (листинги FGED): симметричная пирамида ------
    // (6.36), листинг 6.1: near -> 0, far -> 1.   z_dev = A + B/z,  A = f/(f-n),  B = -nf/(f-n)
    static mathx::Matrix4x4 MakeFrustumProjection(Radian fovY, double s, double n, double f,
                                           ClipY clipY = ClipY::Up);
    // (6.37): обратная к (6.36)
    static mathx::Matrix4x4 MakeFrustumProjectionInverse(Radian fovY, double s, double n, double f,
                                                  ClipY clipY = ClipY::Up);
    // (6.46): обратная бесконечная, near -> 1, бесконечность -> 0.  z_dev = n/z.
    // Для направления (w = 0) строка (0,0,0,n) даёт z_clip = 0 ровно — epsilon не нужен
    // при тесте GREATER_EQUAL и очистке 0. Формат D32_FLOAT.
    static mathx::Matrix4x4 MakeReversedInfiniteProjection(Radian fovY, double s, double n,
                                                    ClipY clipY = ClipY::Up);
    static mathx::Matrix4x4 MakeReversedInfiniteProjectionInverse(Radian fovY, double s, double n,
                                                           ClipY clipY = ClipY::Up);

private:
    static double FocalG(Radian fovY);                    // g = 1 / tan(fovY/2), вид сбоку
    static double ClipSignY(ClipY c) { return c == ClipY::Up ? -1.0 : 1.0; }   // D3D12: y клипа вверх

    void Invalidate() { mRecalcMatrices = true; mRecalcPlanes = true; }

    // Строка глубины P для текущего режима: перспектива z_dev = A + B/z, ортография z_dev = A + B*z.
    struct DepthParams { double A, B; };
    DepthParams GetProjectionParamsAB() const;

    // Окно кадра, y вниз (top < bottom).
    //   перспектива: окно на расстоянии 1 от камеры, то есть ТАНГЕНСЫ углов до краёв
    //                (не зависит от n). Симметричное окно: halfH = 1/g, halfW = s/g.
    //   ортография:  ящик по x, y в единицах камеры.
    void CalcProjectionParameters(double& left, double& right, double& top, double& bottom) const;

    void UpdateMatrices() const { if (mRecalcMatrices) RebuildMatrices(); }

    void RebuildMatrices() const;

    // Плоскости в пространстве камеры: A = P.
    void UpdatePlanes() const;

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

