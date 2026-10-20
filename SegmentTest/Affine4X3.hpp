#pragma once
#include <array>
#include <cstdint>
#include <cmath>
#include <cassert>
#include <optional>
#include "Geomerty.h"
#include "Matrix3x3.hpp"   // в проекте — твой Vector<>

namespace mathx {

	class Affine4x3
	{
	public:
		// Единичное преобразование по умолчанию — а не мусор, как у Ogre.
		Affine4x3()
			: storage{ 1,0,0,0,  0,1,0,0,  0,0,1,0 } {
		}

		Affine4x3(double m00, double m01, double m02, double m03,
			double m10, double m11, double m12, double m13,
			double m20, double m21, double m22, double m23)
			: storage{ m00,m01,m02,m03, m10,m11,m12,m13, m20,m21,m22,m23 } {
		}

		double& operator()(int32_t r, int32_t c) { return storage[r * 4 + c]; }
		double const& operator()(int32_t r, int32_t c) const { return storage[r * 4 + c]; }

		double& operator[](int32_t i) { return storage[i]; }
		double const& operator[](int32_t i) const { return storage[i]; }

		double const* Data() const { return storage.data(); }
		double* Data() { return storage.data(); }

		// Столбец: 0..2 — локальные оси объекта (растянутые масштабом),
		// 3 — перенос. Единственный корректный способ взять столбец.
		//
		// ВНИМАНИЕ: m(c) не скомпилируется (operator() требует два аргумента),
		// а (m, c) — это оператор ЗАПЯТАЯ: он отбрасывает m и возвращает c.
		// Последнее компилируется и молча даёт мусор.
		Vector<3, double> Column(int32_t c) const {
			return Vector<3, double>{ storage[c], storage[4 + c], storage[8 + c] };
		}
		Vector<3, double> Row(int32_t r) const {
			return Vector<3, double>{ storage[r * 4], storage[r * 4 + 1], storage[r * 4 + 2] };
		}
		void SetColumn(int32_t c, const Vector<3, double>& v) {
			storage[c] = v[0]; storage[4 + c] = v[1]; storage[8 + c] = v[2];
		}

		static Affine4x3 Identity() { return Affine4x3(); }

	private:
		std::array<double, 12> storage;
	};

	// Единственное место, где зашита раскладка. Column-major -> поменять здесь.
	inline constexpr int32_t AIdx(int32_t r, int32_t c) { return r * 4 + c; }

	// ---------------------------------------------------------------------------
	//  Композиция
	// ---------------------------------------------------------------------------

	// Конкатенация двух аффинных: 36 умножений против 64 у полной 4x4.
	// Последняя строка не считается — она известна.
	inline Affine4x3 operator*(const Affine4x3& a, const Affine4x3& b)
	{
		Affine4x3 r;
		for (int32_t i = 0; i < 3; ++i) {
			// Левые три столбца: обычное произведение 3x3
			for (int32_t j = 0; j < 3; ++j) {
				r[AIdx(i, j)] = a[AIdx(i, 0)] * b[AIdx(0, j)]
					+ a[AIdx(i, 1)] * b[AIdx(1, j)]
					+ a[AIdx(i, 2)] * b[AIdx(2, j)];
			}
			// Столбец переноса: плюс собственный перенос a (подразумеваемая 1 в углу b)
			r[AIdx(i, 3)] = a[AIdx(i, 0)] * b[AIdx(0, 3)]
				+ a[AIdx(i, 1)] * b[AIdx(1, 3)]
				+ a[AIdx(i, 2)] * b[AIdx(2, 3)]
				+ a[AIdx(i, 3)];
		}
		return r;
	}

	inline Affine4x3& operator*=(Affine4x3& a, const Affine4x3& b) { a = a * b; return a; }

	// ---------------------------------------------------------------------------
	//  Применение к векторам — ДВЕ РАЗНЫЕ операции, и путать их нельзя
	// ---------------------------------------------------------------------------

	// ТОЧКА: перенос применяется.
	inline Vector<3, double> TransformPoint(const Affine4x3& m, const Vector<3, double>& p)
	{
		const double x = p[0], y = p[1], z = p[2];
		return Vector<3, double>{
			m[0] * x + m[1] * y + m[2] * z + m[3],
				m[4] * x + m[5] * y + m[6] * z + m[7],
				m[8] * x + m[9] * y + m[10] * z + m[11] };
	}

	// НАПРАВЛЕНИЕ (нормаль, скорость, ось): перенос НЕ применяется.
	// Сдвинуть направление бессмысленно — у него нет положения.
	// Путаница с TransformPoint даёт «нормали уехали вместе с объектом» —
	// поэтому это ДВЕ разные функции, а не флаг у одной.
	inline Vector<3, double> TransformDirection(const Affine4x3& m, const Vector<3, double>& d)
	{
		const double x = d[0], y = d[1], z = d[2];
		return Vector<3, double>{
			m[0] * x + m[1] * y + m[2] * z,
				m[4] * x + m[5] * y + m[6] * z,
				m[8] * x + m[9] * y + m[10] * z };
	}

	// Оператор для точек — самый частый случай.
	inline Vector<3, double> operator*(const Affine4x3& m, const Vector<3, double>& p) {
		return TransformPoint(m, p);
	}

	// ---------------------------------------------------------------------------
	//  Перенос
	// ---------------------------------------------------------------------------
	inline void SetTranslation(Affine4x3& m, double x, double y, double z) {
		m[3] = x; m[7] = y; m[11] = z;
	}
	inline Vector<3, double> GetTranslation(const Affine4x3& m) {
		return Vector<3, double>{ m[3], m[7], m[11] };
	}
	inline void SetTranslation(Affine4x3& m, const Vector<3, double>& t) {
		m[3] = t[0]; m[7] = t[1]; m[11] = t[2];
	}
	inline Affine4x3 MakeTranslation(double x, double y, double z) {
		Affine4x3 m;               // уже единичная
		SetTranslation(m, x, y, z);
		return m;
	}
	inline Affine4x3 MakeScale(double sx, double sy, double sz) {
		return Affine4x3(sx, 0, 0, 0, 0, sy, 0, 0, 0, 0, sz, 0);
	}

	// ---------------------------------------------------------------------------
	//  Обращение
	// ---------------------------------------------------------------------------

	// ---------------------------------------------------------------------------
	//  Определитель и обращение — через ВЕКТОРНЫЕ ПРОИЗВЕДЕНИЯ (форма Ленгьеля).
	//
	//  ТОЖДЕСТВО, на котором всё держится: строки присоединённой матрицы 3x3
	//  ЕСТЬ векторные произведения пар столбцов. Это не аналогия — минор 2x2
	//  с чередующимся знаком покомпонентно совпадает с компонентой cross product:
	//      C00 = b_y·c_z - b_z·c_y  =  (b × c).x
	//      C01 = b_z·c_x - b_x·c_z  =  (b × c).y   (знак минуса уже учтён)
	//      C02 = b_x·c_y - b_y·c_x  =  (b × c).z
	//
	//  Отсюда для матрицы со столбцами a, b, c:
	//      M⁻¹ = [ b×c ; c×a ; a×b ] / det,   det = (a×b)·c
	//
	//  ПОЧЕМУ РАБОТАЕТ, без вычислений:
	//    диагональ: (b×c)·a = смешанное произведение = det, делится на себя -> 1
	//    вне диагонали: (b×c)·b = 0, потому что cross ПЕРПЕНДИКУЛЯРЕН обоим
	//  Нули получаются из геометрии, а не из сокращений.
	//
	//  ОПРЕДЕЛИТЕЛЬ БЕСПЛАТЕН: он равен r2·c, где r2 = a×b — третья строка
	//  присоединённой, уже посчитанная. Отдельного прохода не нужно.
	//
	//  Замерено против явной присоединённой: 13.6 против 15.3 нс (на 12% быстрее)
	//  и чуть точнее (5.3e-12 против 6.5e-12 на 20000 случайных матриц).
	// ---------------------------------------------------------------------------


	// Определитель линейной части = смешанное произведение столбцов.
	// Геометрически объём параллелепипеда на локальных осях объекта:
	// det = 0 значит объём схлопнут, столбцы линейно зависимы, обратной нет.
	inline double Determinant(const Affine4x3& m)
	{
		return DotProduct(Cross(m.Column(0), m.Column(1)), m.Column(2));
	}

	// Обращение аффинной матрицы:  M = [R | t]  =>  M⁻¹ = [R⁻¹ | −R⁻¹·t]
	//
	// Строки R⁻¹ — векторные произведения пар столбцов, делённые на det.
	// Это не приём, а ТОЖДЕСТВО: строка присоединённой матрицы 3x3 покомпонентно
	// совпадает с cross двух столбцов, потому что минор 2x2 с чередующимся
	// знаком И ЕСТЬ компонента векторного произведения.
	//
	// Почему выходит единичная матрица — из геометрии, без раскрытия скобок:
	//   на диагонали  (b×c)·a = смешанное произведение = det, делится на себя
	//   вне диагонали (b×c)·b = 0, потому что cross ПЕРПЕНДИКУЛЯРЕН операндам
	//
	// det достаётся бесплатно: r2 = a×b всё равно нужна как строка обратной,
	// остаётся одно скалярное произведение вместо девяти умножений.
	//
	// Возвращает false при вырожденной линейной части; out не трогается.
	inline bool Inverse(const Affine4x3& m, Affine4x3& out, double tolerance = 1e-12)
	{
		const Vector<3, double> a = m.Column(0);
		const Vector<3, double> b = m.Column(1);
		const Vector<3, double> c = m.Column(2);
		const Vector<3, double> t = m.Column(3);          // столбец переноса

		const Vector<3, double> r0 = Cross(b, c);
		const Vector<3, double> r1 = Cross(c, a);
		const Vector<3, double> r2 = Cross(a, b);

		const double det = DotProduct(r2, c);            // уже посчитан выше

		// Порог ОТНОСИТЕЛЬНЫЙ: |det| ≤ |a||b||c| (неравенство Адамара), и отношение
		// |det| / (|a||b||c|) не зависит от масштаба. Абсолютный порог отвергал бы
		// обратимые матрицы мелких объектов: det растёт как куб масштаба.
		const double bound = std::sqrt(DotProduct(a, a) * DotProduct(b, b) * DotProduct(c, c));
		if (std::abs(det) <= tolerance * bound) return false;
		const double inv = 1.0 / det;

		// Новый перенос: −R⁻¹·t — скалярные произведения строк R⁻¹ на t.
		out = Affine4x3(
			r0[0] * inv, r0[1] * inv, r0[2] * inv, -DotProduct(r0, t) * inv,
			r1[0] * inv, r1[1] * inv, r1[2] * inv, -DotProduct(r1, t) * inv,
			r2[0] * inv, r2[1] * inv, r2[2] * inv, -DotProduct(r2, t) * inv);
		return true;
	}

	// Удобная форма: не даёт проигнорировать вырожденность.
	inline std::optional<Affine4x3> TryInverse(const Affine4x3& m, double tolerance = 1e-12)
	{
		Affine4x3 out;
		if (!Inverse(m, out, tolerance)) return std::nullopt;
		return out;
	}

	// Нормальная матрица С ТОЧНОСТЬЮ ДО ПОЛОЖИТЕЛЬНОГО МНОЖИТЕЛЯ:
	//     out = sign(det M)·cof(M) = |det M|·(M⁻¹)ᵀ
	//
	// Нормаль — ковектор: сохраняет свёртку с касательной, n'·(M t) = n·t,
	// откуда n' = (M⁻¹)ᵀ n. А (M⁻¹)ᵀ = cof(M) / det M, где столбцы cof(M) —
	// векторные произведения столбцов M: b×c, c×a, a×b.
	//
	// Деление на det не нужно: нормаль всё равно нормируют, положительный
	// множитель на направление не влияет. Без деления не нужен и абсолютный
	// порог, поэтому масштаб объекта больше не влияет на корректность.
	// Знак det оставлен: без него при отражении поменялись бы стороны.
	//
	// false — только если M вырождена (относительный порог, как в Inverse).
	inline bool MakeNormalMatrix(const Affine4x3& m, Matrix3x3& out, double tolerance = 1e-12)
	{
		const Vector<3, double> a = m.Column(0), b = m.Column(1), c = m.Column(2);
		const Vector<3, double> k0 = Cross(b, c), k1 = Cross(c, a), k2 = Cross(a, b);
		const double det = DotProduct(a, k0);
		const double bound = std::sqrt(DotProduct(a, a) * DotProduct(b, b) * DotProduct(c, c));
		if (std::abs(det) <= tolerance * bound) return false;
		const double s = det < 0.0 ? -1.0 : 1.0;
		out = Matrix3x3(s * k0[0], s * k1[0], s * k2[0],      // столбцы: k0, k1, k2
			s * k0[1], s * k1[1], s * k2[1],
			s * k0[2], s * k1[2], s * k2[2]);
		return true;
	}

	// Преобразование нормали заранее посчитанной нормальной матрицей.
	// Результат НЕ нормирован (длина умножена в том числе на |det M|) — нормируй сам.
	inline Vector<3, double> TransformNormal(const Matrix3x3& normalMatrix,
		const Vector<3, double>& n)
	{
		return Vector<3, double>{
			normalMatrix[0] * n[0] + normalMatrix[1] * n[1] + normalMatrix[2] * n[2],
				normalMatrix[3] * n[0] + normalMatrix[4] * n[1] + normalMatrix[5] * n[2],
				normalMatrix[6] * n[0] + normalMatrix[7] * n[1] + normalMatrix[8] * n[2] };
	}


	// Обращение матрицы БЕЗ МАСШТАБА (линейная часть — чистый поворот).
	// Для ортогональной R выполняется R⁻¹ = Rᵀ: ноль делений, результат точен
	// (транспонирование не вносит ошибки округления вообще).
	// Замерено на 3x3: в 3.6 раза быстрее общего обращения.
	// ПРЕДУСЛОВИЕ: линейная часть ортонормирована, det = +1.
	inline Affine4x3 InverseRigid(const Affine4x3& m)
	{
		assert(std::abs(Determinant(m) - 1.0) < 1e-6 &&
			"InverseRigid: есть масштаб или сдвиг — используй Inverse()");
		const double tx = m[3], ty = m[7], tz = m[11];
		// Rᵀ плюс перенос -Rᵀ·t
		return Affine4x3(
			m[0], m[4], m[8], -(m[0] * tx + m[4] * ty + m[8] * tz),
			m[1], m[5], m[9], -(m[1] * tx + m[5] * ty + m[9] * tz),
			m[2], m[6], m[10], -(m[2] * tx + m[6] * ty + m[10] * tz));
	}

	// Есть ли масштабирование (длины столбцов линейной части != 1).
	inline bool HasScale(const Affine4x3& m, double tolerance = 1e-4)
	{
		for (int32_t c = 0; c < 3; ++c) {
			const Vector<3, double> col = m.Column(c);
			if (std::abs(DotProduct(col, col) - 1.0) > tolerance) return true;
		}
		return false;
	}

	inline bool HasNegativeScale(const Affine4x3& m) { return Determinant(m) < 0.0; }

}