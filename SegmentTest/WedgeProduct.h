
//     ^   wedge (внешнее)            a ^ b
//     &   antiwedge (regressive)     a & b   (meet/join)
//     *   geometric product          a * b   (использует метрику)
//     +, -                           линейные
//     *  (float)                     скаляр
//
//  Named-методы: .wedge() .antiwedge() .geometric() .dual() .grade() ...
//  Фабрики (static): Multivector::point/line/plane/e(i)/scalar.
//
//  Сигнатура метрики: METRIC[i] = e_i^2 in {+1,-1,0}.
//     PGA 3D  Cl(3,0,1) : {+1,+1,+1,0}  (e4^2 = 0)  — по умолчанию.
//     Евклид Cl(3,0)   : {+1,+1,+1,+1}.
//
//  Оси (n=4): e1=0b0001, e2=0b0010, e3=0b0100, e4=0b1000.
//  Точка: x e1 + y e2 + z e3 + w e4.  Плоскость: R ^ S ^ T.
//  Требует C++20.
// ============================================================================
#pragma once
#include <cstdint>
#include <array>
#include <bit>
#include <string>
#include <cmath>
#include <algorithm>
#include <ostream>

class Multivector {
public:
	static constexpr int DIM = 4;
	static constexpr int COUNT = 1 << DIM;      // 16 блейдов

	// коэффициент при каждом блейде, индекс = битовая маска
	std::array<float, COUNT> c{};

	// ======================= конструкторы / фабрики =========================
	Multivector() = default;

	static Multivector e(int i) {                 // базисный вектор e_i (1..DIM)
		Multivector m; m.c[1u << (i - 1)] = 1.0f; return m;
	}
	static Multivector scalar(float s) {
		Multivector m; m.c[0] = s; return m;
	}
	static Multivector point(float x, float y, float z, float w = 1.0f) {
		Multivector p;
		p.c[0b0001] = x; p.c[0b0010] = y; p.c[0b0100] = z; p.c[0b1000] = w;
		return p;
	}
	static Multivector line(const Multivector& R, const Multivector& S) {
		return R ^ S;                              // grade-2
	}
	static Multivector plane(const Multivector& R, const Multivector& S,
		const Multivector& T) {
		return R ^ S ^ T;                          // grade-3
	}

	// ============================ метрика ===================================
	static std::array<int, DIM>& metric() {
		static std::array<int, DIM> m = { +1, +1, +1, 0 };   // PGA 3D
		return m;
	}
	static void set_metric(std::array<int, DIM> m) { metric() = m; }

	// ======================== доступ по маске ================================
	float& operator[](uint32_t mask) { return c[mask]; }
	float  operator[](uint32_t mask) const { return c[mask]; }

	// ======================== линейные операции =============================
	Multivector operator+(const Multivector& o) const {
		Multivector r; for (uint32_t m = 0; m < COUNT; ++m) r.c[m] = c[m] + o.c[m]; return r;
	}
	Multivector operator-(const Multivector& o) const {
		Multivector r; for (uint32_t m = 0; m < COUNT; ++m) r.c[m] = c[m] - o.c[m]; return r;
	}
	Multivector operator*(float s) const {
		Multivector r; for (uint32_t m = 0; m < COUNT; ++m) r.c[m] = c[m] * s; return r;
	}

	// ========================= grade-проекция ===============================
	Multivector grade(int g) const {
		Multivector r;
		for (uint32_t m = 0; m < COUNT; ++m)
			if (std::popcount(m) == g) r.c[m] = c[m];
		return r;
	}
	int top_grade() const {
		int g = -1;
		for (uint32_t m = 0; m < COUNT; ++m)
			if (std::fabs(c[m]) > 1e-7f) g = std::max<int>(g, std::popcount(m));
		return g;
	}

	// ===================== ГЛАВНЫЕ ПРОИЗВЕДЕНИЯ =============================

	Multivector wedge(const Multivector& o) const {          // ^  внешнее
		return bilinear(*this, o, [](uint32_t a, uint32_t b) -> BladeProduct {
			return { a | b, wedge_sign(a, b) };
			});
	}
	Multivector operator^(const Multivector& o) const { return wedge(o); }

	Multivector geometric(const Multivector& o) const {      // *  клиффордово
		const auto& met = metric();
		return bilinear(*this, o, [&](uint32_t a, uint32_t b) {
			return geo_blade(a, b, met);
			});
	}
	Multivector operator*(const Multivector& o) const { return geometric(o); }

	// Reverse (антиавтоморфизм ~): обращение порядка векторных сомножителей.
	// На грейде k действует знаком (-1)^{k(k-1)/2}: скаляр/вектор +, бивектор/тривектор -.
	Multivector reverse() const {
		Multivector r;
		for (uint32_t m = 0; m < COUNT; ++m) {
			int k = std::popcount(m);
			int sign = ((k * (k - 1) / 2) & 1) ? -1 : +1;
			r.c[m] = sign * c[m];
		}
		return r;
	}

	Multivector dual() const {                                // right complement
		constexpr uint32_t I = COUNT - 1;                     // 0b1111
		Multivector r;
		for (uint32_t m = 0; m < COUNT; ++m) {
			if (c[m] == 0.0f) continue;
			uint32_t comp = I ^ m;
			r.c[comp] = wedge_sign(m, comp) * c[m];
		}
		return r;
	}

	Multivector antiwedge(const Multivector& o) const {       // &  regressive
		return dual().wedge(o.dual()).dual();
	}
	Multivector operator&(const Multivector& o) const { return antiwedge(o); }

	// ============ CONTRACTIONS (свёртки) как grade-проекции geometric ========
	// Для grade-a * grade-b полное произведение содержит грейды |a-b|..a+b.
	// Contraction выбирает нужную часть.

	// Scalar product: <u * v>_0  — только скалярная часть.
	// Для двух grade-1 векторов даёт евклидов dot по метрике.
	Multivector scalar_product(const Multivector& o) const {
		return geometric(o).grade(0);
	}

	// Left contraction u ⌋ v: проекция на grade (b - a).
	// "u contracted INTO v" — удаляет u-часть из v, понижает грейд.
	// Ноль, если grade(u) > grade(v).
	Multivector left_contract(const Multivector& o) const {
		Multivector r;
		for (int a = 0; a <= DIM; ++a) {
			for (int b = 0; b <= DIM; ++b) {
				int target = b - a;
				if (target < 0) continue;                          // left: b>=a
				r = r + grade(a).geometric(o.grade(b)).grade(target);
			}
		}
		return r;
	}

	// Right contraction u ⌊ v: проекция на grade (a - b). Симметрично.
	Multivector right_contract(const Multivector& o) const {
		Multivector r;
		for (int a = 0; a <= DIM; ++a) {
			for (int b = 0; b <= DIM; ++b) {
				int target = a - b;
				if (target < 0) continue;
				r = r + grade(a).geometric(o.grade(b)).grade(target);
			}
		}
		return r;
	}

	// Симметричный dot Hestenes'а: |<uv>_(|a-b|)|. Часто именно это имеют в виду.
	// Для векторов совпадает со scalar_product, для высших — самая младшая grade.
	Multivector dot(const Multivector& o) const {
		Multivector r;
		for (int a = 0; a <= DIM; ++a) {
			for (int b = 0; b <= DIM; ++b) {
				int target = std::abs(a - b);
				r = r + grade(a).geometric(o.grade(b)).grade(target);
			}
		}
		return r;
	}

	// ======================== утилиты =======================================
	bool approx_equal(const Multivector& o, float eps = 1e-5f) const {
		for (uint32_t m = 0; m < COUNT; ++m)
			if (std::fabs(c[m] - o.c[m]) > eps) return false;
		return true;
	}
	std::string str() const {
		std::string s; bool first = true;
		for (uint32_t m = 0; m < COUNT; ++m) {
			if (std::fabs(c[m]) < 1e-7f) continue;
			if (!first) s += " + ";
			s += std::to_string(c[m]) + "*" + blade_name(m);
			first = false;
		}
		return first ? "0" : s;
	}

private:
	struct BladeProduct { uint32_t blade; int sign; };

	static int wedge_sign(uint32_t a, uint32_t b) {           // знак перестановки
		if (a & b) return 0;
		int swaps = 0; a >>= 1;
		while (a) { swaps += std::popcount(a & b); a >>= 1; }
		return (swaps & 1) ? -1 : +1;
	}

	static BladeProduct geo_blade(uint32_t a, uint32_t b,     // инверсии + метрика
		const std::array<int, DIM>& met) {
		int sign = 1;
		{
			uint32_t aa = a >> 1; int sw = 0;
			while (aa) { sw += std::popcount(aa & b); aa >>= 1; }
			if (sw & 1) sign = -sign;
		}
		uint32_t common = a & b;
		while (common) {
			int i = std::countr_zero(common);
			int m = met[i];
			if (m == 0) return { 0u, 0 };
			if (m < 0)  sign = -sign;
			common &= common - 1;
		}
		return { a ^ b, sign };
	}

	template <class Combine>
	static Multivector bilinear(const Multivector& x, const Multivector& y,
		Combine combine) {
		Multivector r;
		for (uint32_t a = 0; a < COUNT; ++a) {
			if (x.c[a] == 0.0f) continue;
			for (uint32_t b = 0; b < COUNT; ++b) {
				if (y.c[b] == 0.0f) continue;
				auto pr = combine(a, b);
				if (pr.sign == 0) continue;
				r.c[pr.blade] += pr.sign * x.c[a] * y.c[b];
			}
		}
		return r;
	}

	static std::string blade_name(uint32_t m) {
		if (m == 0) return "1";
		std::string s;
		for (int i = 0; i < DIM; ++i)
			if (m & (1u << i)) s += "e" + std::to_string(i + 1);
		return s;
	}
};

inline Multivector operator*(float s, const Multivector& m) { return m * s; }
inline std::ostream& operator<<(std::ostream& os, const Multivector& mv) {
	return os << mv.str();
}