#pragma once
#include <random>
#include <cstdint>

// ============================================================================
//  math_random.hpp — ГСЧ как объект, не как глобальная статика.
//
//  Почему не как в Ogre (static UnitRandom/RangeRandom + RandomValueProvider):
//    - глобальное состояние => не потокобезопасно без блокировок;
//    - невозможно воспроизвести прогон в тестах (нет управления seed);
//    - две подсистемы, дёргающие один генератор, влияют друг на друга.
//  Объект решает всё три: у каждого потока/подсистемы свой, seed явный.
// ============================================================================

namespace mathx {

	class RandomGenerator {
	public:
		// Детерминированный seed — для тестов и воспроизводимых прогонов.
		explicit RandomGenerator(uint64_t seed) : engine(seed) {}

		// Недетерминированный seed — для рантайма.
		RandomGenerator() : engine(std::random_device{}()) {}

		// [0, 1]
		double Unit() { return unitDist(engine); }

		// [low, high]
		double Range(double low, double high) {
			return low + (high - low) * Unit();
		}

		// [-1, 1]
		double Symmetric() { return 2.0 * Unit() - 1.0; }

		// Целое в [low, high] включительно.
		int32_t RangeInt(int32_t low, int32_t high) {
			std::uniform_int_distribution<int32_t> d(low, high);
			return d(engine);
		}

		// Нормальное (гауссово) распределение — заменяет Ogre::gaussianDistribution.
		double Gaussian(double mean = 0.0, double stddev = 1.0) {
			std::normal_distribution<double> d(mean, stddev);
			return d(engine);
		}

		void Reseed(uint64_t seed) { engine.seed(seed); }

	private:
		std::mt19937_64 engine;
		std::uniform_real_distribution<double> unitDist{ 0.0, 1.0 };
	};

} // namespace mathx
