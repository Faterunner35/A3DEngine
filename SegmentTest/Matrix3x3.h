#pragma once

#include"Geomerty.h"
#include<array>
#include<iostream>
#include<cassert>

namespace mathx
{ 

class Matrix3x3 {
public:
	Matrix3x3() : storage{} {}   // явно нулевая инициализация, а не мусор

	Matrix3x3(double e00, double e01, double e02,
		double e10, double e11, double e12,
		double e20, double e21, double e22)
		: storage{ e00,e01,e02, e10,e11,e12, e20,e21,e22 } {
	}
	// explicit убран: конструктору от 9 аргументов нечего случайно
	// конвертировать, explicit здесь бесполезен.

	double& operator()(uint16_t r, uint16_t c) { return storage[r * 3 + c]; }
	double const& operator()(uint16_t r, uint16_t c) const { return storage[r * 3 + c]; }

	// Линейный доступ mat[i] — то, чего не хватало для CreateRotorMatrix.
	double& operator[](uint16_t i) { return storage[i]; }
	double const& operator[](uint16_t i) const { return storage[i]; }

	void SetRow(uint16_t r, const Vector<3, double>& v) {
		for (int i = 0; i < 3; ++i) storage[r * 3 + i] = v[i];
	}
	void SetColumn(uint16_t c, const Vector<3, double>& v) {
		for (int i = 0; i < 3; ++i) storage[i * 3 + c] = v[i];
	}

	// ---- Утилиты ----

	static Matrix3x3 Identity3() {
		return Matrix3x3(1, 0, 0, 0, 1, 0, 0, 0, 1);
	}

	

	static Vector<3, double> GetRow(const Matrix3x3& m, int32_t r) {
		return Vector<3, double>{ m(r, 0), m(r, 1), m(r, 2) };
	}
	static Vector<3, double> GetColumn(const Matrix3x3& m, int32_t c) {
		return Vector<3, double>{ m(0, c), m(1, c), m(2, c) };
	}

	


	// ---- Определитель и обращение ----
	//
	// Обращение через присоединённую матрицу: A⁻¹ = adj(A)/det(A).
	// adj(A) — транспонированная матрица алгебраических дополнений; для 3x3
	// каждое дополнение это определитель 2x2, то есть пара произведений.
	//
	// ТРЮК: определитель раскладывается по первой СТРОКЕ через уже посчитанные
	// дополнения: det = a00*C00 + a01*C01 + a02*C02. Поэтому считаем девять
	// дополнений один раз, а det достаём из трёх из них — без пересчёта.

	// Матрица алгебраических дополнений, СРАЗУ транспонированная (то есть adj).
	

	// Ортогональна ли (с точностью до допуска) и правая ли тройка.
	

	// Ортонормировка Грама-Шмидта по столбцам. Лечит накопленный дрейф после
	// многих перемножений: столбцы перестают быть единичными и ортогональными.
	
	
	// Свой operator= не пишем — компилятор сам генерирует корректное
	// почленное копирование std::array. Ручная реализация тут не нужна
	// и в исходном виде была невалидна (см. разбор).

private:
	std::array<double, 9> storage;   // no alignas unless you commit to 3x4 SIMD layout
};


inline Matrix3x3 operator+(const Matrix3x3& a, const Matrix3x3& b) {
	Matrix3x3 r;
	for (uint16_t i = 0; i < 9; ++i) r[i] = a[i] + b[i];
	return r;
}
inline Matrix3x3 operator-(const Matrix3x3& a, const Matrix3x3& b) {
	Matrix3x3 r;
	for (uint16_t i = 0; i < 9; ++i) r[i] = a[i] - b[i];
	return r;
}
inline Matrix3x3 operator-(const Matrix3x3& a) {
	Matrix3x3 r;
	for (uint16_t i = 0; i < 9; ++i) r[i] = -a[i];
	return r;
}
inline Matrix3x3 operator*(const Matrix3x3& m, double s) {
	Matrix3x3 r;
	for (uint16_t i = 0; i < 9; ++i) r[i] = m[i] * s;
	return r;
}
// Свободная функция обязательна: левый операнд double не может быть методом.
inline Matrix3x3 operator*(double s, const Matrix3x3& m) { return m * s; }

// Композиция поворотов. Некоммутативна: A*B != B*A.
// Порядок задаёт, в ЧЬИХ осях довернули — мировых или локальных.
inline Matrix3x3 operator*(const Matrix3x3& a, const Matrix3x3& b) {
	Matrix3x3 r;
	for (int32_t i = 0; i < 3; ++i)
		for (int32_t j = 0; j < 3; ++j) {
			double s = 0.0;
			for (int32_t k = 0; k < 3; ++k) s += a(i, k) * b(k, j);
			r(i, j) = s;
		}
	return r;
}

// M * v (вектор-столбец) — согласовано с конвенцией v' = M*v.
inline Vector<3, double> operator*(const Matrix3x3& m, const Vector<3, double>& v) {
	return Vector<3, double>{
		m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
			m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
			m[6] * v[0] + m[7] * v[1] + m[8] * v[2] };
}

inline bool operator==(const Matrix3x3& a, const Matrix3x3& b) {
	for (uint16_t i = 0; i < 9; ++i) if (a[i] != b[i]) return false;
	return true;
}
inline bool operator!=(const Matrix3x3& a, const Matrix3x3& b) { return !(a == b); }


inline Matrix3x3 Adjugate(const Matrix3x3& m) {
	return Matrix3x3(
		m[4] * m[8] - m[5] * m[7], m[2] * m[7] - m[1] * m[8], m[1] * m[5] - m[2] * m[4],
		m[5] * m[6] - m[3] * m[8], m[0] * m[8] - m[2] * m[6], m[2] * m[3] - m[0] * m[5],
		m[3] * m[7] - m[4] * m[6], m[1] * m[6] - m[0] * m[7], m[0] * m[4] - m[1] * m[3]);
}

inline double Determinant(const Matrix3x3& m) {
	// Разложение по первой строке; элементы adj[0], adj[3], adj[6] — это
	// дополнения C00, C01, C02 (adj уже транспонирована, отсюда шаг 3).
	return m[0] * (m[4] * m[8] - m[5] * m[7])
		- m[1] * (m[3] * m[8] - m[5] * m[6])
		+ m[2] * (m[3] * m[7] - m[4] * m[6]);
}

inline Matrix3x3 Transpose(const Matrix3x3& m) {
	return Matrix3x3(m[0], m[3], m[6],
		m[1], m[4], m[7],
		m[2], m[5], m[8]);
}


inline bool HasScale(const Matrix3x3& m, double tolerance = 1e-4) {
	for (int32_t c = 0; c < 3; ++c) {
		const double len2 = m(0, c) * m(0, c)
			+ m(1, c) * m(1, c)
			+ m(2, c) * m(2, c);
		if (std::abs(len2 - 1.0) > tolerance) return true;
	}
	return false;
}


inline bool IsOrthonormal(const Matrix3x3& m, double tolerance = 1e-6) {
	if (HasScale(m, tolerance)) return false;
	// Попарные скалярные произведения столбцов должны быть нулевыми
	for (int32_t a = 0; a < 3; ++a)
		for (int32_t b = a + 1; b < 3; ++b) {
			double d = 0.0;
			for (int32_t r = 0; r < 3; ++r) d += m(r, a) * m(r, b);
			if (std::abs(d) > tolerance) return false;
		}
	// det = +1, а не -1: иначе это отражение, а не поворот
	return mathx::Determinant(m) > 0.0;
}

inline Matrix3x3 InverseOrthogonal(const Matrix3x3& m) {
	assert(std::abs(Determinant(m) - 1.0) < 1e-6 &&
		"InverseOrthogonal: матрица не ортогональна (есть масштаб/сдвиг) — "
		"используй Inverse()");
	return Transpose(m);
}

// Обращение общей матрицы. Возвращает false, если |det| <= tolerance
// (матрица вырождена — обратной не существует), out при этом не трогается.
inline bool Inverse(const Matrix3x3& m, Matrix3x3& out, double tolerance = 1e-12) {
	const Matrix3x3 adj = Adjugate(m);
	// det через уже посчитанные дополнения первой строки: adj транспонирована,
	// поэтому C00,C01,C02 лежат в adj[0], adj[3], adj[6].
	const double det = m[0] * adj[0] + m[1] * adj[3] + m[2] * adj[6];

	// Порог относительный (неравенство Адамара): |det| ≤ |c0||c1||c2|
	double n0 = 0, n1 = 0, n2 = 0;
	for (int r = 0; r < 3; ++r) { n0 += m[r * 3] * m[r * 3]; n1 += m[r * 3 + 1] * m[r * 3 + 1]; n2 += m[r * 3 + 2] * m[r * 3 + 2]; }
	if (std::abs(det) <= tolerance * std::sqrt(n0 * n1 * n2)) return false;

	const double invDet = 1.0 / det;
	for (uint16_t i = 0; i < 9; ++i) out[i] = adj[i] * invDet;
	return true;
}

// Обращение ОРТОГОНАЛЬНОЙ матрицы (чистого поворота).
//
// Для ортогональной M выполняется M⁻¹ = Mᵀ — это следствие MᵀM = I.
// Ноль умножений и делений против ~30 у общего обращения, плюс результат
// точен: транспонирование не вносит ошибки округления вообще.
// ПРЕДУСЛОВИЕ: матрица ортогональна с det = +1. Если есть масштаб или
// сдвиг — результат будет молча неверным, поэтому assert в debug.


// ---- Ортогональность и масштаб ----

// Есть ли масштабирование: длины столбцов (локальных осей) должны быть 1.



inline void Orthonormalize(Matrix3x3& m) {
	// Столбец 0: нормируем
	double l0 = std::sqrt(m[0] * m[0] + m[3] * m[3] + m[6] * m[6]);
	m[0] /= l0; m[3] /= l0; m[6] /= l0;

	// Столбец 1: вычитаем проекцию на 0, нормируем
	double d01 = m[0] * m[1] + m[3] * m[4] + m[6] * m[7];
	m[1] -= d01 * m[0]; m[4] -= d01 * m[3]; m[7] -= d01 * m[6];
	double l1 = std::sqrt(m[1] * m[1] + m[4] * m[4] + m[7] * m[7]);
	m[1] /= l1; m[4] /= l1; m[7] /= l1;

	// Столбец 2 = столбец0 x столбец1 — сразу единичный и ортогональный обоим,
	// считать проекции не нужно.
	m[2] = m[3] * m[7] - m[6] * m[4];
	m[5] = m[6] * m[1] - m[0] * m[7];
	m[8] = m[0] * m[4] - m[3] * m[1];
}

// Тензорное (внешнее) произведение: P = u * vᵀ. Нужно для построения
// проекторов и матриц отражения.
inline Matrix3x3 TensorProduct(const Vector<3, double>& u, const Vector<3, double>& v) {
	return Matrix3x3(u[0] * v[0], u[0] * v[1], u[0] * v[2],
		u[1] * v[0], u[1] * v[1], u[1] * v[2],
		u[2] * v[0], u[2] * v[1], u[2] * v[2]);
}


inline Matrix3x3  CreateRotationFromAxis(const Vector<3, double>& vec, double angle)
{    
    Matrix3x3 mat;
	double s = std::sin(angle), c = std::cos(angle), d = 1.0 - c;
	mat[0] = c + d * (vec[0] * vec[0]);      mat[1] = d * (vec[0] * vec[1]) - s * vec[2];   mat[2] = d * (vec[0] * vec[2]) + s * vec[1];
	mat[3] = d * (vec[0] * vec[1]) + s * vec[2];   mat[4] = c + d * (vec[1] * vec[1]);      mat[5] = d * (vec[1] * vec[2]) - s * vec[0];
	mat[6] = d * (vec[0] * vec[2]) - s * vec[1];   mat[7] = d * (vec[1] * vec[2]) + s * vec[0];   mat[8] = c + d * (vec[2] * vec[2]);
	return mat;
}




}
