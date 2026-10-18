#pragma once
#include "Matrix3x3.hpp"
#include"Affine4X3.hpp"
#include"TransfromMatrix.hpp"
#include <optional>
#include <cmath>
#include <algorithm>

namespace mathx {

	// σ_max линейной части = √(наибольшее собственное число MᵀM), аналитически для 3×3.
	// Это максимальное растяжение, которое преобразование даёт какому-либо вектору.
	inline double SpectralNorm(const Matrix3x3& m)
	{
		double A[3][3];
		for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
			double s = 0; for (int k = 0; k < 3; ++k) s += m(k, i) * m(k, j);
			A[i][j] = s;
		}
		const double p1 = A[0][1] * A[0][1] + A[0][2] * A[0][2] + A[1][2] * A[1][2];
		if (p1 == 0.0) return std::sqrt(std::max({ A[0][0], A[1][1], A[2][2] }));
		const double q = (A[0][0] + A[1][1] + A[2][2]) / 3.0;
		const double p2 = (A[0][0] - q) * (A[0][0] - q) + (A[1][1] - q) * (A[1][1] - q) + (A[2][2] - q) * (A[2][2] - q) + 2.0 * p1;
		const double p = std::sqrt(p2 / 6.0);
		double B[3][3];
		for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) B[i][j] = (A[i][j] - (i == j ? q : 0.0)) / p;
		double r = 0.5 * (B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1])
			- B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0])
			+ B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0]));
		r = std::clamp(r, -1.0, 1.0);
		return std::sqrt(q + 2.0 * p * std::cos(std::acos(r) / 3.0));
	}

	// ============================================================================
	//  AffineTransform — аффинное преобразование F(p) = M·p + t вместе с обратным.
	//
	//  Обратное хранится, потому что нужно в двух местах:
	//    - перевод из мирового пространства в локальное (рейкаст по объекту);
	//    - нормали. Касательная переносится как t' = M·t, а нормаль должна
	//      сохранить свёртку с ней: n'·t' = n·t. Отсюда n' = (M⁻¹)ᵀ·n —
	//      транспонированная линейная часть ОБРАТНОГО.
	//
	//  Хранение обоих делает две операции бесплатными по обращению:
	//    Inverse()  — обмен матриц;
	//    A * B      — обратное произведения собирается как B⁻¹·A⁻¹.
	// ============================================================================
	class AffineTransform
	{
	public:
		// ---------------- создание ----------------

		static AffineTransform Identity() { return AffineTransform(Affine4x3::Identity(), Affine4x3::Identity()); }

		// Общая матрица: обращается численно. nullopt — линейная часть вырождена.
		static std::optional<AffineTransform> TryCreate(const Affine4x3& m, double tolerance = 1e-12)
		{
			Affine4x3 inv;
			if (!mathx::Inverse(m, inv, tolerance)) return std::nullopt;
			return AffineTransform(m, inv);
		}

		static std::optional<AffineTransform> TryCreate(const Matrix3x3& L, const Vector<3, double>& t,
			double tolerance = 1e-12)
		{
			return TryCreate(Affine4x3(L(0, 0), L(0, 1), L(0, 2), t[0],
				L(1, 0), L(1, 1), L(1, 2), t[1],
				L(2, 0), L(2, 1), L(2, 2), t[2]), tolerance);
		}

		// Линейное 3×3 — нулевой перенос
		static std::optional<AffineTransform> TryCreate(const Matrix3x3& L, double tolerance = 1e-12)
		{
			return TryCreate(L, Vector<3, double>{0, 0, 0}, tolerance);
		}

		// Из позиции, ориентации и масштаба. Обратное строится ТОЧНО, без
		// численного обращения: масштаб 1/s, сопряжённый кватернион, −позиция.
		static std::optional<AffineTransform> FromTRS(const Vector<3, double>& position,
			const Quaternion& orientation,
			const Vector<3, double>& scale,
			double tolerance = 1e-12)
		{
			if (std::abs(scale[0]) <= tolerance || std::abs(scale[1]) <= tolerance ||
				std::abs(scale[2]) <= tolerance) return std::nullopt;
			return AffineTransform(mathx::MakeTransform(position, orientation, scale),
				MakeInverseTransform(position, orientation, scale));
		}

		// ---------------- алгебра ----------------

		// Обратное: обмен матриц, обращение не считается
		AffineTransform Inverse() const { return AffineTransform(inv_, m_); }

		// Композиция: (A * B)(p) = A(B(p)). Обратное — B⁻¹·A⁻¹, тоже без обращения.
		friend AffineTransform operator*(const AffineTransform& a, const AffineTransform& b)
		{
			return AffineTransform(a.m_ * b.m_, b.inv_ * a.inv_);
		}

		// ---------------- вперёд: из локального в мировое ----------------

		Vector<3, double> ApplyPoint(const Vector<3, double>& p) const { return TransformPoint(m_, p); }
		Vector<3, double> ApplyDirection(const Vector<3, double>& d) const { return TransformDirection(m_, d); }

		// Нормаль: (M⁻¹)ᵀ·n — столбцы обратного вместо строк, затем нормировка
		Vector<3, double> ApplyNormal(const Vector<3, double>& n) const { return TransposedLinear(inv_, n); }

		// ---------------- назад: из мирового в локальное ----------------

		Vector<3, double> ApplyInversePoint(const Vector<3, double>& p) const { return TransformPoint(inv_, p); }
		Vector<3, double> ApplyInverseDirection(const Vector<3, double>& d) const { return TransformDirection(inv_, d); }

		// Нормаль обратно: ((M⁻¹)⁻¹)ᵀ·n = Mᵀ·n
		Vector<3, double> ApplyInverseNormal(const Vector<3, double>& n) const { return TransposedLinear(m_, n); }

		// ---------------- свойства ----------------

		const Affine4x3& Matrix()        const { return m_; }
		const Affine4x3& InverseMatrix() const { return inv_; }
		double Linear(int r, int c)      const { return m_[AIdx(r, c)]; }
		Vector<3, double> Translation()   const { return GetTranslation(m_); }
		double Determinant()             const { return mathx::Determinant(m_); }

		// det < 0: меняет ориентацию — важно для мешей
		bool Mirrors() const { return Determinant() < 0.0; }

		// Жёсткое: только поворот и перенос — столбцы ортонормированы и det = +1
		bool IsRigid(double tolerance = 1e-9) const
		{
			for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
				double s = 0; for (int k = 0; k < 3; ++k) s += m_[AIdx(k, i)] * m_[AIdx(k, j)];
				if (std::abs(s - (i == j ? 1.0 : 0.0)) > tolerance) return false;
			}
			return Determinant() > 0.0;
		}

		// Максимальное растяжение = спектральная норма линейной части (нужна сферам)
		double MaxStretch() const
		{
			return SpectralNorm(Matrix3x3(m_[0], m_[1], m_[2], m_[4], m_[5], m_[6], m_[8], m_[9], m_[10]));
		}

	private:
		AffineTransform(const Affine4x3& m, const Affine4x3& inv) : m_(m), inv_(inv) {}

		// Aᵀ·n по линейной части A, нормированное
		static Vector<3, double> TransposedLinear(const Affine4x3& a, const Vector<3, double>& n)
		{
			Vector<3, double> r;
			for (int i = 0; i < 3; ++i)
				r[i] = a[AIdx(0, i)] * n[0] + a[AIdx(1, i)] * n[1] + a[AIdx(2, i)] * n[2];
			const double k = 1.0 / std::sqrt(DotProduct(r, r));
			return Vector<3, double>{ r[0] * k, r[1] * k, r[2] * k };
		}

		Affine4x3 m_;      // F
		Affine4x3 inv_;    // F⁻¹
	};

}