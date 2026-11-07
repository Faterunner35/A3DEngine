#include "Math.h"

// ---- Реализации, требующие обоих классов сразу (как в оригинале) ----

Radian::Radian(const Degree& d) : radian(d.ValueDegrees()* detail::ANGLE_DEG_TO_RAD) {}

Radian& Radian::operator=(const Degree& d) {
	radian = d.ValueDegrees() * detail::ANGLE_DEG_TO_RAD;
	return *this;
}
Radian Radian::operator+(const Degree& d) const {
	return Radian(radian + d.ValueDegrees() * detail::ANGLE_DEG_TO_RAD);
}
Radian& Radian::operator+=(const Degree& d) {
	radian += d.ValueDegrees() * detail::ANGLE_DEG_TO_RAD;
	return *this;
}
Radian Radian::operator-(const Degree& d) const {
	return Radian(radian - d.ValueDegrees() * detail::ANGLE_DEG_TO_RAD);
}
Radian& Radian::operator-=(const Degree& d) {
	radian -= d.ValueDegrees() * detail::ANGLE_DEG_TO_RAD;
	return *this;
}

double Radian::ValueDegrees() const { return radian * detail::ANGLE_RAD_TO_DEG; }
double Degree::ValueRadians() const { return degree * detail::ANGLE_DEG_TO_RAD; }

namespace mathx {

	double WrapAngle(double radians) {
		double r = std::fmod(radians + PI, TWO_PI);
		if (r < 0.0) r += TWO_PI;
		return r - PI;
	}

	double WrapUnit(double v) {
		double intPart;
		double frac = std::modf(v, &intPart);
		if (frac < 0.0) frac += 1.0;
		return frac;
	}

} // namespace mathx
