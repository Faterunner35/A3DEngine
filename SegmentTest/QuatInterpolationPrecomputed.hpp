#pragma once
#include "Quaternion.hpp"
#include <vector>
#include <span>
#include <cassert>

namespace mathx {

	// ============================================================================
	//  ПРОТОТИП: slerp с предвычислением для анимации.
	//
	//  ЦЕНА, которую платим (осознанно): базовый quat_slerp самодостаточен —
	//  дал два кватерниона, получил результат. Здесь появляется КОНТРАКТ между
	//  двумя фазами: подготовка клипа и покадровое проигрывание. Следствия:
	//    - предвычисленные данные надо где-то держать и инвалидировать;
	//    - изменил ключ -> кэш протух, кто-то обязан пересчитать;
	//    - функция принимает cosA НА ВЕРУ, проверить дёшево нельзя.
	//  Поэтому это НЕ замена quat_slerp, а отдельный путь для случая, когда
	//  профиль показал, что интерполяция — узкое место.
	//
	//  ЧТО ВЫНОСИТСЯ ИЗ ГОРЯЧЕГО ЦИКЛА:
	//    - знаковая коррекция (Dot < 0 -> инвертировать ключ) — делается один
	//      раз для всей последовательности при загрузке;
	//    - сам Dot(q0,q1) — пары ключей не меняются между кадрами;
	//    - опционально acos(cosA) -> угол, если считать коэффициенты напрямую.
	//  В покадровом пути остаётся только вычисление весов.
	// ============================================================================

	// ---------------------------------------------------------------------------
	//  Фаза 1: подготовка клипа (один раз при загрузке)
	// ---------------------------------------------------------------------------

	// Предвычисленные данные для ОДНОГО сегмента между соседними ключами.
	struct SlerpSegment
	{
		double cosAngle = 1.0;   // Dot(k[i], k[i+1]) после знаковой коррекции, >= 0
		double angle = 0.0;   // acos(cosAngle) — вынесен из покадрового пути
		double invSin = 0.0;   // 1/sin(angle); 0 если сегмент вырожден
		bool   degenerate = true; // угол ~0 -> использовать lerp вместо slerp
	};

	// Подготовка последовательности ключей.
	// ВАЖНО: keys модифицируется на месте — знаковая коррекция переворачивает
	// ключи так, чтобы соседние образовывали острый угол. После этого в рантайме
	// проверка Dot < 0 не нужна вообще.
	// Возвращает segments.size() == keys.size() - 1.
	inline std::vector<SlerpSegment> PrepareSlerpSequence(std::span<Quaternion> keys)
	{
		std::vector<SlerpSegment> segments;
		if (keys.size() < 2) return segments;
		segments.reserve(keys.size() - 1);

		for (size_t i = 0; i + 1 < keys.size(); ++i) {
			// Знаковая коррекция: разворачиваем СЛЕДУЮЩИЙ ключ, чтобы угол был острым.
			// Делается цепочкой, поэтому порядок обхода важен (слева направо).
			double c = quat_dot(keys[i], keys[i + 1]);
			if (c < 0.0) {
				for (int j = 0; j < 4; ++j) keys[i + 1][j] = -keys[i + 1][j];
				c = -c;
			}

			SlerpSegment s;
			s.cosAngle = c;
			if (c >= 1.0 - SLERP_COS_EPSILON) {
				s.degenerate = true;               // ключи совпали -> lerp
			}
			else {
				s.degenerate = false;
				s.angle = std::acos(c);
				s.invSin = 1.0 / std::sqrt(1.0 - c * c);
			}
			segments.push_back(s);
		}
		return segments;
	}

	// ---------------------------------------------------------------------------
	//  Фаза 2: покадровое проигрывание (горячий путь)
	// ---------------------------------------------------------------------------

	// Интерполяция внутри подготовленного сегмента.
	// ПРЕДУСЛОВИЕ: q0, q1 — соседние ключи ПОСЛЕ PrepareSlerpSequence
	// (знаковая коррекция уже применена, поэтому shortestPath не нужен),
	// seg — соответствующий им сегмент.
	// В горячем пути: два sin, никаких acos/dot/ветвлений по знаку.
	inline Quaternion quat_slerp_precomputed(double t,
		Quaternion const& q0,
		Quaternion const& q1,
		SlerpSegment const& seg)
	{
		if (seg.degenerate) {
			// Точный предел slerp при theta->0 (правило Лопиталя), не аппроксимация.
			Quaternion r = q0 * (1.0 - t) + q1 * t;
			quat_normalize(r);
			return r;
		}
		const double f0 = std::sin((1.0 - t) * seg.angle) * seg.invSin;
		const double f1 = std::sin(t * seg.angle) * seg.invSin;
		return q0 * f0 + q1 * f1;
	}

	// ---------------------------------------------------------------------------
	//  Батчевый путь — под Ogre-Next-стиль (SoA).
	//
	//  В Ogre-Next ArrayQuaternion держит компоненты раздельно (x[], y[], z[], w[])
	//  и обрабатывает 4 кости за раз в SIMD-регистре. Здесь интерфейс без SIMD,
	//  но с ТОЙ ЖЕ раскладкой, чтобы позже подставить интринсики не меняя API.
	//  Раскладка важнее самого SIMD: она определяет, сколько полезных данных
	//  приходит за одну cache line.
	// ---------------------------------------------------------------------------

	// Кватернионы в SoA-раскладке: 4 отдельных непрерывных массива.
	struct QuaternionSoA
	{
		std::vector<double> x, y, z, w;

		size_t size() const { return w.size(); }
		void resize(size_t n) { x.resize(n); y.resize(n); z.resize(n); w.resize(n); }

		Quaternion Get(size_t i) const { return Quaternion(x[i], y[i], z[i], w[i]); }
		void Set(size_t i, Quaternion const& q) { x[i] = q[0]; y[i] = q[1]; z[i] = q[2]; w[i] = q[3]; }
	};

	// Интерполяция целой позы скелета: N костей, у каждой свой сегмент.
	// Один t на всю позу (обычная ситуация: клип проигрывается синхронно).
	// Внутренний цикл трогает только непрерывные массивы -> префетчер угадывает.
	inline void SlerpBatch(double t,
		QuaternionSoA const& from,
		QuaternionSoA const& to,
		std::span<const SlerpSegment> segments,
		QuaternionSoA& out)
	{
		const size_t n = from.size();
		assert(to.size() == n && segments.size() == n && "SlerpBatch: size mismatch");
		out.resize(n);

		for (size_t i = 0; i < n; ++i) {
			const SlerpSegment& s = segments[i];
			double f0, f1;
			if (s.degenerate) {
				f0 = 1.0 - t;  f1 = t;
			}
			else {
				f0 = std::sin((1.0 - t) * s.angle) * s.invSin;
				f1 = std::sin(t * s.angle) * s.invSin;
			}
			// Раздельные массивы: каждая строка — линейный проход
			double rx = from.x[i] * f0 + to.x[i] * f1;
			double ry = from.y[i] * f0 + to.y[i] * f1;
			double rz = from.z[i] * f0 + to.z[i] * f1;
			double rw = from.w[i] * f0 + to.w[i] * f1;
			if (s.degenerate) {
				const double inv = 1.0 / std::sqrt(rx * rx + ry * ry + rz * rz + rw * rw);
				rx *= inv; ry *= inv; rz *= inv; rw *= inv;
			}
			out.x[i] = rx; out.y[i] = ry; out.z[i] = rz; out.w[i] = rw;
		}
	}

} // namespace mat