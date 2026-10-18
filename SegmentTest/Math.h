#pragma once
#include"Geomerty.h"
#include <ostream>
#include<cassert>
#include<bit>


	// Константы конверсии — определены здесь, чтобы angle.hpp был
	// самодостаточным (math.hpp инклудит его, не наоборот).
	namespace detail {
		inline constexpr double ANGLE_PI = 3.14159265358979323846;
		inline constexpr double ANGLE_DEG_TO_RAD = ANGLE_PI / 180.0;
		inline constexpr double ANGLE_RAD_TO_DEG = 180.0 / ANGLE_PI;
		inline constexpr double PI = 3.14159265358979323846;
		
	}

	class Degree;

	class Radian
	{
	public:
		explicit Radian(double r = 0) : radian(r) {}
		Radian(const Degree& d);

		Radian& operator=(double f) { radian = f; return *this; }
		Radian& operator=(const Degree& d);

		double ValueDegrees() const;
		double ValueRadians() const { return radian; }

		const Radian& operator+() const { return *this; }
		Radian        operator+(const Radian& r) const { return Radian(radian + r.radian); }
		Radian        operator+(const Degree& d) const;
		Radian& operator+=(const Radian& r) { radian += r.radian; return *this; }
		Radian& operator+=(const Degree& d);

		Radian  operator-() const { return Radian(-radian); }
		Radian  operator-(const Radian& r) const { return Radian(radian - r.radian); }
		Radian  operator-(const Degree& d) const;
		Radian& operator-=(const Radian& r) { radian -= r.radian; return *this; }
		Radian& operator-=(const Degree& d);

		Radian  operator*(double f) const { return Radian(radian * f); }
		Radian  operator*(const Radian& r) const { return Radian(radian * r.radian); }
		Radian& operator*=(double f) { radian *= f; return *this; }

		Radian  operator/(double f) const { return Radian(radian / f); }
		Radian& operator/=(double f) { radian /= f; return *this; }

		bool operator<(const Radian& r)  const { return radian < r.radian; }
		bool operator<=(const Radian& r) const { return radian <= r.radian; }
		bool operator==(const Radian& r) const { return radian == r.radian; }
		bool operator!=(const Radian& r) const { return radian != r.radian; }
		bool operator>=(const Radian& r) const { return radian >= r.radian; }
		bool operator>(const Radian& r)  const { return radian > r.radian; }

		friend std::ostream& operator<<(std::ostream& o, const Radian& v) {
			return o << v.radian << " rad";
		}

	private:
		double radian;
	};

	class Degree
	{
	public:
		explicit Degree(double d = 0) : degree(d) {}
		Degree(const Radian& r) : degree(r.ValueDegrees()) {}

		Degree& operator=(double f) { degree = f; return *this; }
		Degree& operator=(const Radian& r) { degree = r.ValueDegrees(); return *this; }

		double ValueDegrees() const { return degree; }
		double ValueRadians() const;

		const Degree& operator+() const { return *this; }
		Degree        operator+(const Degree& d) const { return Degree(degree + d.degree); }
		Degree        operator+(const Radian& r) const { return Degree(degree + r.ValueDegrees()); }
		Degree& operator+=(const Degree& d) { degree += d.degree; return *this; }
		Degree& operator+=(const Radian& r) { degree += r.ValueDegrees(); return *this; }

		Degree  operator-() const { return Degree(-degree); }
		Degree  operator-(const Degree& d) const { return Degree(degree - d.degree); }
		Degree  operator-(const Radian& r) const { return Degree(degree - r.ValueDegrees()); }
		Degree& operator-=(const Degree& d) { degree -= d.degree; return *this; }
		Degree& operator-=(const Radian& r) { degree -= r.ValueDegrees(); return *this; }

		Degree  operator*(double f) const { return Degree(degree * f); }
		Degree  operator*(const Degree& d) const { return Degree(degree * d.degree); }
		Degree& operator*=(double f) { degree *= f; return *this; }

		Degree  operator/(double f) const { return Degree(degree / f); }
		Degree& operator/=(double f) { degree /= f; return *this; }

		bool operator<(const Degree& d)  const { return degree < d.degree; }
		bool operator<=(const Degree& d) const { return degree <= d.degree; }
		bool operator==(const Degree& d) const { return degree == d.degree; }
		bool operator!=(const Degree& d) const { return degree != d.degree; }
		bool operator>=(const Degree& d) const { return degree >= d.degree; }
		bool operator>(const Degree& d)  const { return degree > d.degree; }

		friend std::ostream& operator<<(std::ostream& o, const Degree& v) {
			return o << v.degree << " deg";
		}

	private:
		double degree;
	};

	// ---- Реализации, требующие обоих классов сразу (как в оригинале) ----


	inline Radian::Radian(const Degree& d) : radian(d.ValueDegrees()* detail::ANGLE_DEG_TO_RAD) {}

	inline Radian& Radian::operator=(const Degree& d) {
		radian = d.ValueDegrees() * detail::ANGLE_DEG_TO_RAD;
		return *this;
	}
	inline Radian Radian::operator+(const Degree& d) const {
		return Radian(radian + d.ValueDegrees() * detail::ANGLE_DEG_TO_RAD);
	}
	inline Radian& Radian::operator+=(const Degree& d) {
		radian += d.ValueDegrees() * detail::ANGLE_DEG_TO_RAD;
		return *this;
	}
	inline Radian Radian::operator-(const Degree& d) const {
		return Radian(radian - d.ValueDegrees() * detail::ANGLE_DEG_TO_RAD);
	}
	inline Radian& Radian::operator-=(const Degree& d) {
		radian -= d.ValueDegrees() * detail::ANGLE_DEG_TO_RAD;
		return *this;
	}

	inline double Radian::ValueDegrees() const { return radian * detail::ANGLE_RAD_TO_DEG; }
	inline double Degree::ValueRadians() const { return degree * detail::ANGLE_DEG_TO_RAD; }



enum class EulerResult
{
	// The solution is invalid (incorrect axis indices).
	INVALID,

	// The solution is unique.
	UNIQUE,

	// The solution is not unique.  A sum of angles is constant.
	NOT_UNIQUE_SUM,

	// The solution is not unique.  A difference of angles is constant.
	NOT_UNIQUE_DIF
};

class EulerAngles {
public:
	EulerAngles() : axis{ 0,0,0 }, angle{ 0,0,0 }, result(EulerResult::INVALID) {}
	EulerAngles(int32_t i0, int32_t i1, int32_t i2, double a0, double a1, double a2)
		: axis{ i0,i1,i2 }, angle{ a0,a1,a2 }, result(EulerResult::UNIQUE) {
	}

	std::array<int32_t, 3> axis;
	std::array<double, 3> angle;   // ВСЕГДА радианы — внутреннее соглашение, как в AxisAngle
	EulerResult result;

	// Явные, безопасные точки входа/выхода на границу с UI/аниматором:
	Radian AngleRadians(int i) const { return Radian(angle[i]); }
	Degree AngleDegrees(int i) const { return Degree(angle[i] * 180.0 / detail::PI); }
};




namespace mathx {

	// ---------------------------------------------------------------------------
	//  Константы
	// ---------------------------------------------------------------------------
	inline constexpr double PI = 3.14159265358979323846;
	inline constexpr double TWO_PI = 2.0 * PI;
	inline constexpr double HALF_PI = 0.5 * PI;
	inline constexpr double INV_PI = 1.0 / PI;
	inline constexpr double DEG_TO_RAD = PI / 180.0;
	inline constexpr double RAD_TO_DEG = 180.0 / PI;

	inline constexpr double POS_INFINITY = std::numeric_limits<double>::infinity();
	inline constexpr double NEG_INFINITY = -std::numeric_limits<double>::infinity();
	inline constexpr double EPSILON = std::numeric_limits<double>::epsilon();

	// ---------------------------------------------------------------------------
	//  Скалярные операции
	//  std::abs/sqrt/pow/exp/log НЕ оборачиваем — вызывай напрямую.
	// ---------------------------------------------------------------------------
	template <typename T>
	inline constexpr T Sqr(T v) { return v * v; }

	template <typename T>
	inline constexpr T Sign(T v) {
		return (v > T(0)) ? T(1) : ((v < T(0)) ? T(-1) : T(0));
	}

	template <typename T>
	inline constexpr T Clamp(T val, T minVal, T maxVal) {
		assert(minVal <= maxVal && "Invalid clamp range");
		return std::max(std::min(val, maxVal), minVal);
	}

	// Зажим в [0,1] — saturate из шейдеров.
	template <typename T>
	inline constexpr T Saturate(T t) { return Clamp(t, T(0), T(1)); }

	// Lerp(a,b,0)=a, Lerp(a,b,1)=b.
	template <typename T, typename S>
	inline constexpr T Lerp(const T& a, const T& b, const S& w) {
		return a + w * (b - a);
	}

	// Гладкая интерполяция Эрмита между edge0 и edge1.
	template <typename T>
	inline constexpr T Smoothstep(T edge0, T edge1, T x) {
		const T t = Saturate((x - edge0) / (edge1 - edge0));
		return t * t * (T(3) - T(2) * t);
	}

	// Используется в нормализации векторов достаточно часто, чтобы иметь имя.
	inline double InvSqrt(double v) { return 1.0 / std::sqrt(v); }

	// ---------------------------------------------------------------------------
	//  Проверки и сравнения
	// ---------------------------------------------------------------------------
	inline bool IsNaN(double f) { return std::isnan(f); }
	inline bool IsFinite(double f) { return std::isfinite(f); }

	// Абсолютный допуск — годится для величин порядка 1.
	inline bool NearEqual(double a, double b, double tolerance = 1e-9) {
		return std::abs(b - a) <= tolerance;
	}

	// Относительный — корректен для любых масштабов.
	inline bool NearEqualRelative(double a, double b, double relTolerance = 1e-9) {
		const double diff = std::abs(a - b);
		const double scale = std::max(std::abs(a), std::abs(b));
		return diff <= relTolerance * std::max(scale, 1.0);
	}

	// ---------------------------------------------------------------------------
	//  Тригонометрия — принимает Radian, чтобы нельзя было передать градусы.
	//  Прямые обёртки над std::, без lookup-таблиц.
	// ---------------------------------------------------------------------------
	inline double Sin(Radian a) { return std::sin(a.ValueRadians()); }
	inline double Cos(Radian a) { return std::cos(a.ValueRadians()); }
	inline double Tan(Radian a) { return std::tan(a.ValueRadians()); }

	// Обратные возвращают Radian, не голый double.
	// ASin/ACos зажимают домен: накопленная ошибка даёт 1.0000001,
	// от чего std::asin вернул бы NaN.
	inline Radian ASin(double v) { return Radian(std::asin(Clamp(v, -1.0, 1.0))); }
	inline Radian ACos(double v) { return Radian(std::acos(Clamp(v, -1.0, 1.0))); }
	inline Radian ATan(double v) { return Radian(std::atan(v)); }
	inline Radian ATan2(double y, double x) { return Radian(std::atan2(y, x)); }

	// Оба сразу — частый паттерн (построение ротора, матрицы поворота).
	inline void SinCos(Radian a, double& outSin, double& outCos) {
		const double r = a.ValueRadians();
		outSin = std::sin(r);
		outCos = std::cos(r);
	}

	// Abs/Sign сохраняют тип угла.
	// Sqrt(Radian) намеренно отсутствует: корень из угла не имеет
	// геометрического смысла (в Ogre это артефакт слепого перегружания).
	inline Radian Abs(Radian a) { return Radian(std::abs(a.ValueRadians())); }
	inline Degree Abs(Degree d) { return Degree(std::abs(d.ValueDegrees())); }
	inline Radian Sign(Radian a) { return Radian(Sign(a.ValueRadians())); }
	inline Degree Sign(Degree d) { return Degree(Sign(d.ValueDegrees())); }

	// ---------------------------------------------------------------------------
	//  Заворачивание диапазонов
	// ---------------------------------------------------------------------------
	// Угол в [-PI, PI).
	inline double WrapAngle(double radians) {
		double r = std::fmod(radians + PI, TWO_PI);
		if (r < 0.0) r += TWO_PI;
		return r - PI;
	}

	// Значение в [0,1).
	inline double WrapUnit(double v) {
		double intPart;
		double frac = std::modf(v, &intPart);
		if (frac < 0.0) frac += 1.0;
		return frac;
	}

	// ---------------------------------------------------------------------------
	//  Целочисленные и битовые утилиты
	// ---------------------------------------------------------------------------
	inline constexpr uint32_t NextPow2(uint32_t x) {
		return (x <= 1u) ? 1u : std::bit_ceil(x);
	}

	inline constexpr bool IsPow2(uint32_t x) {
		return x != 0u && (x & (x - 1u)) == 0u;
	}

	// Итеративный, не рекурсивный — рекурсия здесь тратит стек зря.
	inline constexpr uint64_t Gcd(uint64_t a, uint64_t b) {
		while (b != 0u) { uint64_t t = b; b = a % b; a = t; }
		return a;
	}

	// (a/gcd)*b, а не (a*b)/gcd — защита от переполнения.
	inline constexpr uint64_t Lcm(uint64_t a, uint64_t b) {
		return (a == 0u || b == 0u) ? 0u : (a / Gcd(a, b)) * b;
	}

	// Выравнивание вверх до кратного alignment (степень двойки).
	inline constexpr uint64_t AlignUp(uint64_t value, uint64_t alignment) {
		return (value + alignment - 1u) & ~(alignment - 1u);
	}

} //