#pragma once

#include "Quaternion.hpp"
#include"PrimitiveTransform.hpp"
// подстрой под имя своего файла конверсий

namespace mathx {
	namespace {

		// Циклическая ли перестановка: 012 -> 120 -> 201.
		// Циклическая сохраняет ориентацию, транспозиция двух индексов даёт
		// отражение — отсюда общий знак sigma во всех формулах.
		inline bool IsCyclic(int32_t i0, int32_t i1, int32_t i2)
		{
			return (i1 == (i0 + 1) % 3) && (i2 == (i1 + 1) % 3);
		}

		// Защита от NaN: накопленная ошибка даёт 1.0000000001, и asin/acos
		// вернули бы NaN без единого признака ошибки.
		inline double Clamp11(double v)
		{
			return v < -1.0 ? -1.0 : (v > 1.0 ? 1.0 : v);
		}

	} // anonymous namespace

	EulerAngles MatrixToEuler(const Matrix3x3& R, int32_t i0, int32_t i1, int32_t i2)
	{
		EulerAngles e;
		e.axis = { i0, i1, i2 };

		const bool inRange = (0 <= i0 && i0 < 3) && (0 <= i1 && i1 < 3) && (0 <= i2 && i2 < 3);
		if (!inRange || i0 == i1 || i1 == i2) {
			e.result = EulerResult::INVALID;
			return e;
		}

		if (i0 != i2) {
			// ---- Тейт-Брайан: три разные оси ----
			const double sigma = IsCyclic(i0, i1, i2) ? 1.0 : -1.0;
			const double sinA1 = Clamp11(sigma * R[EIdx(i0, i2)]);

			if (std::abs(sinA1) < 1.0 - EULER_DEGENERATE_EPSILON) {
				e.angle[1] = std::asin(sinA1);
				e.angle[0] = std::atan2(-sigma * R[EIdx(i1, i2)], R[EIdx(i2, i2)]);
				e.angle[2] = std::atan2(-sigma * R[EIdx(i0, i1)], R[EIdx(i0, i0)]);
				e.result = EulerResult::UNIQUE;
			}
			else {
				// Оси i0 и i2 коллинеарны: определена лишь комбинация крайних углов.
				// Соглашение — angle0 = 0, всё сваливаем в angle2.
				e.angle[0] = 0.0;
				e.angle[1] = (sinA1 > 0.0) ? (0.5 * EULER_PI) : (-0.5 * EULER_PI);
				e.angle[2] = std::atan2(sigma * R[EIdx(i1, i0)], R[EIdx(i1, i1)]);
				// +pi/2 -> оси СОнаправлены -> сумма; -pi/2 -> встречны -> разность
				e.result = (sinA1 > 0.0) ? EulerResult::NOT_UNIQUE_SUM
					: EulerResult::NOT_UNIQUE_DIF;
			}
		}
		else {
			// ---- Собственно Эйлер: i0 == i2 ----
			// Средний угол через acos от диагонального элемента,
			// вырождение при angle1 = 0 или pi (а не +-pi/2).
			const int32_t k = 3 - i0 - i1;                     // свободная ось
			const double sigma = (i1 == (i0 + 1) % 3) ? 1.0 : -1.0;
			const double cosA1 = Clamp11(R[EIdx(i0, i0)]);

			if (std::abs(cosA1) < 1.0 - EULER_DEGENERATE_EPSILON) {
				e.angle[1] = std::acos(cosA1);
				e.angle[0] = std::atan2(R[EIdx(i1, i0)], -sigma * R[EIdx(k, i0)]);
				e.angle[2] = std::atan2(R[EIdx(i0, i1)], sigma * R[EIdx(i0, k)]);
				e.result = EulerResult::UNIQUE;
			}
			else {
				e.angle[0] = 0.0;
				e.angle[1] = (cosA1 > 0.0) ? 0.0 : EULER_PI;
				// При angle1 = pi средний поворот разворачивает ось, и третий
				// поворот идёт в обратную сторону — отсюда смена знака.
				const double sgn = (cosA1 > 0.0) ? 1.0 : -1.0;
				e.angle[2] = std::atan2(sgn * sigma * R[EIdx(k, i1)], R[EIdx(i1, i1)]);
				e.result = (cosA1 > 0.0) ? EulerResult::NOT_UNIQUE_SUM
					: EulerResult::NOT_UNIQUE_DIF;
			}
		}
		return e;
	}

	EulerAngles QuaternionToEuler(const Quaternion& q, int32_t i0, int32_t i1, int32_t i2)
	{
		return MatrixToEuler(QuaternionToMatrix(q), i0, i1, i2);
	}

}


