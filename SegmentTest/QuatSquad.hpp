#pragma once
#include "Quaternion.hpp"
#include <cmath>
#include <cassert>

namespace mathx {

	// ============================================================================
	//  Squad — сплайновая интерполяция ориентаций.
	//
	//  ЗАЧЕМ. Цепочка slerp между ключевыми кадрами — это ломаная: внутри
	//  сегмента движение ровное, но в узле скорость меняется СКАЧКОМ и объект
	//  дёргается. Squad проводит через те же ключи гладкую кривую.
	//  Нужен при РЕДКИХ ключах (позы на 1, 12, 24 кадре при 60 fps).
	//  При плотных ключах разница не видна, а slerp втрое дешевле.
	// ============================================================================

	// quat_log / quat_exp — в quaternion.hpp (нужны физике независимо).

	// Контрольные точки для узла q1 по трём подряд идущим ключам.
	//
	//   P0 = q0⁻¹·q1  — «откуда пришли»,  P1 = q1⁻¹·q2 — «куда идём»
	//   arg = ¼(log P0 − log P1)  — оценка касательной в узле; это в точности
	//         Катмулл-Ром, только в касательном пространстве вместо прямой
	//   A = q1·exp(arg)   — исходящая контрольная точка
	//   B = q1·exp(−arg)  — входящая, симметричная
	//
	// Коэффициент ¼ подобран так, чтобы производная squad-кривой в узле
	// совпала с этой касательной — это и даёт непрерывность скорости.
	//
	// ПРЕДУСЛОВИЕ: q0, q1, q2 единичные.
	inline void quat_intermediate(Quaternion const& q0,
		Quaternion const& q1,
		Quaternion const& q2,
		Quaternion& outA,
		Quaternion& outB)
	{
		assert(std::fabs(quat_length(q0) - 1.0) < 1e-6 &&
			std::fabs(quat_length(q1) - 1.0) < 1e-6 &&
			std::fabs(quat_length(q2) - 1.0) < 1e-6 &&
			"quat_intermediate: inputs must be unit quaternions");

		const Quaternion P0 = quat_conjugate(q0) * q1;   // для единичного conj == inverse
		const Quaternion P1 = quat_conjugate(q1) * q2;

		const Quaternion l0 = quat_log(P0);
		const Quaternion l1 = quat_log(P1);
		const Quaternion arg = (l0 - l1) * 0.25;

		outA = q1 * quat_exp(arg);
		outB = q1 * quat_exp(-arg);
	}

	// Сферическая кубическая интерполяция: кривая от p до q, форму задают
	// контрольные точки a (исходящая из p) и b (входящая в q).
	//
	// Множитель 2t(1−t) равен нулю на концах => внешний slerp там возвращает
	// первый аргумент, то есть кривая проходит через p и q ТОЧНО, контрольные
	// точки на концы не влияют. Максимум 0.5 при t=0.5 — там они максимально
	// оттягивают кривую. Это сферический аналог де Кастельжо.
	inline Quaternion quat_squad(double t,
		Quaternion const& p,
		Quaternion const& a,
		Quaternion const& b,
		Quaternion const& q,
		bool shortestPath = true)
	{
		const double outerT = 2.0 * t * (1.0 - t);
		const Quaternion s1 = quat_slerp(t, p, q, shortestPath);   // прямая дуга
		const Quaternion s2 = quat_slerp(t, a, b, shortestPath);   // дуга контрольных точек
		return quat_slerp(outerT, s1, s2, shortestPath);
	}

} // namespace mathx