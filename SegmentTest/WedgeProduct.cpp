#include "WedgeProduct.h"

Multivector Multivector::grade(int g) const {
	Multivector r;
	for (uint32_t m = 0; m < COUNT; ++m)
		if (std::popcount(m) == g) r.c[m] = c[m];
	return r;
}

int Multivector::top_grade() const {
	int g = -1;
	for (uint32_t m = 0; m < COUNT; ++m)
		if (std::fabs(c[m]) > 1e-7f) g = std::max<int>(g, std::popcount(m));
	return g;
}

Multivector Multivector::reverse() const {
	Multivector r;
	for (uint32_t m = 0; m < COUNT; ++m) {
		int k = std::popcount(m);
		int sign = ((k * (k - 1) / 2) & 1) ? -1 : +1;
		r.c[m] = sign * c[m];
	}
	return r;
}

Multivector Multivector::dual() const {
	constexpr uint32_t I = COUNT - 1;                     // 0b1111
	Multivector r;
	for (uint32_t m = 0; m < COUNT; ++m) {
		if (c[m] == 0.0f) continue;
		uint32_t comp = I ^ m;
		r.c[comp] = wedge_sign(m, comp) * c[m];
	}
	return r;
}

Multivector Multivector::scalar_product(const Multivector& o) const {
	return geometric(o).grade(0);
}

Multivector Multivector::left_contract(const Multivector& o) const {
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

Multivector Multivector::right_contract(const Multivector& o) const {
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

Multivector Multivector::dot(const Multivector& o) const {
	Multivector r;
	for (int a = 0; a <= DIM; ++a) {
		for (int b = 0; b <= DIM; ++b) {
			int target = std::abs(a - b);
			r = r + grade(a).geometric(o.grade(b)).grade(target);
		}
	}
	return r;
}

bool Multivector::approx_equal(const Multivector& o, float eps) const {
	for (uint32_t m = 0; m < COUNT; ++m)
		if (std::fabs(c[m] - o.c[m]) > eps) return false;
	return true;
}

std::string Multivector::str() const {
	std::string s; bool first = true;
	for (uint32_t m = 0; m < COUNT; ++m) {
		if (std::fabs(c[m]) < 1e-7f) continue;
		if (!first) s += " + ";
		s += std::to_string(c[m]) + "*" + blade_name(m);
		first = false;
	}
	return first ? "0" : s;
}

std::string Multivector::blade_name(uint32_t m) {
	if (m == 0) return "1";
	std::string s;
	for (int i = 0; i < DIM; ++i)
		if (m & (1u << i)) s += "e" + std::to_string(i + 1);
	return s;
}

std::ostream& operator<<(std::ostream& os, const Multivector& mv) {
	return os << mv.str();
}
