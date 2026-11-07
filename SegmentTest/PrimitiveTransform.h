#pragma once
#include "Math.h"          // Radian, Degree
#include "Matrix3x3.h"
#include "Quaternion.h"
#include <array>
#include <cmath>
#include <cstdint>

namespace mathx {

	// ============================================================================
	//  Разложение поворота на углы Эйлера.
	//
	//  ДВЕНАДЦАТЬ конвенций, не 27: два соседних поворота вокруг одной оси
	//  складываются в один, поэтому нужно N0!=N1 и N1!=N2 (N0 и N2 совпадать
	//  МОГУТ). Итого 3*2*2 = 12.
	//    Тейта-Брайана (все разные): 012 021 102 120 201 210
	//    собственно Эйлера (N0==N2): 010 020 101 121 202 212
	//
	//  Конвенция: R = R_{axis0}(angle0) * R_{axis1}(angle1) * R_{axis2}(angle2),
	//  матрица действует на вектор-столбец: v' = R*v.
	//
	//  ВЫВОД (XYZ; остальные — та же формула с перестановкой индексов):
	//    R00 = cb*cc   R01 = -cb*sc   R02 = sb
	//    R12 = -sa*cb  R22 = ca*cb
	//  Отсюда:
	//    b = asin(R02)         <- единственный элемент с чистым sin
	//    a = atan2(-R12, R22)  <- пара, где cb общий множитель и сокращается
	//    c = atan2(-R01, R00)  <- то же
	// ============================================================================

	// Индексация матрицы: row-major. ЕДИНСТВЕННОЕ место, где зашита раскладка.
	//
	// ВНИМАНИЕ: писать m[(r, c)] НЕЛЬЗЯ — круглые скобки внутри [] это
	// оператор ЗАПЯТАЯ, он отбрасывает первый операнд и возвращает второй.
	// m[(1,2)] читает m[2] вместо m[5]. Код компилируется и молча считает
	// мусор. Используй m(r, c) или m[EIdx(r, c)].
	inline constexpr int32_t EIdx(int32_t r, int32_t c) { return r * 3 + c; }

	inline constexpr double EULER_PI = 3.14159265358979323846;

	// Порог вырождения. Фактическую границу задаёт точность asin/acos —
	// они упираются в +-1 раньше, чем сработает номинальный порог; он тут
	// как страховка от NaN.
	inline constexpr double EULER_DEGENERATE_EPSILON = 1e-15;

	enum class EulerResult : uint8_t
	{
		INVALID,          // недопустимые индексы осей
		UNIQUE,           // решение единственно
		NOT_UNIQUE_SUM,   // вырождение: постоянна СУММА angle2 + angle0
		NOT_UNIQUE_DIF    // вырождение: постоянна РАЗНОСТЬ angle2 - angle0
	};

	// Оси хранятся РЯДОМ с углами: в Ogre шесть функций имеют одинаковую
	// сигнатуру, но третий параметр в XYZ значит поворот вокруг Z, а в XZY —
	// вокруг Y. Здесь смысл угла виден из данных, а не из имени функции.
	//
	// Единственный тип углов Эйлера в проекте. В Math.h был глобальный
	// дубликат EulerAngles/EulerResult без пользователей; при
	// using namespace mathx имена становились неоднозначными. Он удалён,
	// его типизированные геттеры перенесены сюда.
	struct EulerAngles
	{
		std::array<int32_t, 3> axis{ 0, 1, 2 };   // осмысленная конвенция по умолчанию
		std::array<double, 3> angle{ 0.0, 0.0, 0.0 };   // ВСЕГДА радианы
		EulerResult result = EulerResult::INVALID;

		EulerAngles() = default;
		EulerAngles(int32_t i0, int32_t i1, int32_t i2, double a0, double a1, double a2)
			: axis{ i0, i1, i2 }, angle{ a0, a1, a2 }, result(EulerResult::UNIQUE) {
		}

		bool IsValid()  const { return result != EulerResult::INVALID; }
		bool IsUnique() const { return result == EulerResult::UNIQUE; }

		// Явные, безопасные точки входа/выхода на границу с UI/аниматором.
		Radian AngleRadians(int32_t i) const { return Radian(angle[i]); }
		Degree AngleDegrees(int32_t i) const { return Degree(Radian(angle[i])); }
	};

	// ---------------------------------------------------------------------------
	//  Построение матрицы
	// ---------------------------------------------------------------------------

	// Элементарный поворот вокруг оси: заполняются только 5 ненулевых элементов.
	Matrix3x3 AxisRotation(int32_t axis, double angle);

	// Домножение СЛЕВА на элементарный поворот: m <- R_axis(angle) * m.
	//
	// У R_axis строка с номером axis есть (0,..,1,..,0), поэтому соответствующая
	// строка m проходит НЕИЗМЕННОЙ. Меняются только две строки — 12 умножений
	// вместо 27 у полного matmul, где две трети это умножения на 0 и 1.
	// Замерено: 92% экономии против AxisRotation + полное умножение.
	void PreRotate(Matrix3x3& m, int32_t axis, double angle);

	// Сборка матрицы из углов через PreRotate.
	// Наивный путь (три полные матрицы + два полных умножения) стоит 54 умножения;
	// здесь ~24. Замерено: на 66% быстрее при идентичном результате.
	Matrix3x3 EulerToMatrix(const EulerAngles& e);

	// ---------------------------------------------------------------------------
	//  Разложение — ХОЛОДНЫЕ функции, тела в EulerAngles.cpp.
	//  Вызываются при импорте/экспорте сцены и в инспекторе, не в кадре.
	//  Атрибут cold просит вынести их в .text.unlikely, чтобы холодный код
	//  не вытеснял горячий из i-cache.
	// ---------------------------------------------------------------------------
#if defined(__GNUC__) || defined(__clang__)
#define MATHX_COLD [[gnu::cold]]
#else
#define MATHX_COLD
#endif

	MATHX_COLD EulerAngles MatrixToEuler(const Matrix3x3& R, int32_t i0, int32_t i1, int32_t i2);
	MATHX_COLD EulerAngles QuaternionToEuler(const Quaternion& q, int32_t i0, int32_t i1, int32_t i2);

} /// namespace mathx