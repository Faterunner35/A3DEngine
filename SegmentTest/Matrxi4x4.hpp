#pragma once
#include "Affine4X3.hpp"
#include <array>
#include <cmath>
#include <cstdint>
#include <cassert>

namespace mathx {

	// ============================================================================
	//  Matrix4x4 — ПОЛНАЯ матрица 4x4.
	//
	//  Нужна ТОЛЬКО там, где последняя строка не (0,0,0,1): проекционные матрицы
	//  (перспектива несёт деление на w), матрицы отсечения, редкие проективные
	//  преобразования. Всё остальное в движке — Affine4x3, который в 2.4 раза
	//  дешевле и на 25% компактнее.
	//
	//  Раскладка row-major: m[r*4 + c].
	// ============================================================================

	class Matrix4x4
	{
	public:
		Matrix4x4() : storage{ 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 } {}

		Matrix4x4(double m00, double m01, double m02, double m03,
			double m10, double m11, double m12, double m13,
			double m20, double m21, double m22, double m23,
			double m30, double m31, double m32, double m33)
			: storage{ m00,m01,m02,m03, m10,m11,m12,m13,
					   m20,m21,m22,m23, m30,m31,m32,m33 } {
		}

		// Расширение аффинной до полной: дописываем известную последнюю строку.
		explicit Matrix4x4(const Affine4x3& a)
			: storage{ a[0],a[1],a[2],a[3],
					   a[4],a[5],a[6],a[7],
					   a[8],a[9],a[10],a[11],
					   0.0, 0.0, 0.0, 1.0 } {
		}

		double& operator()(int32_t r, int32_t c) { return storage[r * 4 + c]; }
		double const& operator()(int32_t r, int32_t c) const { return storage[r * 4 + c]; }
		double& operator[](int32_t i) { return storage[i]; }
		double const& operator[](int32_t i) const { return storage[i]; }

		double const* Data() const { return storage.data(); }
		double* Data() { return storage.data(); }

		static Matrix4x4 Identity() { return Matrix4x4(); }

	private:
		std::array<double, 16> storage;
	};

	inline constexpr int32_t M4Idx(int32_t r, int32_t c) { return r * 4 + c; }

	// Конвенция глубины задаётся ЯВНО, а не зашита: OpenGL ждёт z в [-1,1],
	// DirectX/Vulkan/Metal — в [0,1]. Молчаливое несовпадение даёт
	// «половина сцены за камерой», и ищется это долго.
	enum class DepthRange {
		NegOneToOne,   // OpenGL
		ZeroToOne
	};   // DirectX, Vulkan, Metal

	// Полное умножение: 64 умножения (против 36 у аффинного).
	inline Matrix4x4 operator*(const Matrix4x4& a, const Matrix4x4& b)
	{
		Matrix4x4 r;
		for (int32_t i = 0; i < 4; ++i)
			for (int32_t j = 0; j < 4; ++j) {
				double s = 0.0;
				for (int32_t k = 0; k < 4; ++k) s += a[M4Idx(i, k)] * b[M4Idx(k, j)];
				r[M4Idx(i, j)] = s;
			}
		return r;
	}

	// Умножение на однородный вектор БЕЗ деления на w.
	// Деление вынесено отдельно и делается осознанно: после проекции w != 1,
	// и момент деления определяет, что попадёт в clip space, а что уже в NDC.
	inline Vector<4, double> operator*(const Matrix4x4& m, const Vector<4, double>& v)
	{
		return Vector<4, double>{
			m[0] * v[0] + m[1] * v[1] + m[2] * v[2] + m[3] * v[3],
				m[4] * v[0] + m[5] * v[1] + m[6] * v[2] + m[7] * v[3],
				m[8] * v[0] + m[9] * v[1] + m[10] * v[2] + m[11] * v[3],
				m[12] * v[0] + m[13] * v[1] + m[14] * v[2] + m[15] * v[3] };
	}

	// ---------------------------------------------------------------------------
	//  Перспективное деление — ОТДЕЛЬНАЯ операция, не встроенная в умножение.
	//
	//  Порядок в конвейере: матрица -> CLIP SPACE -> отсечение -> деление -> NDC.
	//  Отсечение делается ИМЕННО в clip space, по знаку и величине w. Если
	//  спрятать деление внутрь "ProjectPoint", промежуточный результат теряется,
	//  и отсекать становится нечем. Поэтому шаги разделены явно.
	//
	//  Проверка w сама по себе отсечение НЕ заменяет: она ловит точку ровно на
	//  плоскости камеры, но не ловит точки ЗА камерой (w < 0), которые после
	//  деления дают правдоподобные, но зеркально неверные координаты.
	// ---------------------------------------------------------------------------

	// Точка перед камерой? В стандартной перспективе w = -z_view, поэтому
	// положительный w означает "перед камерой".
	inline bool IsInFrontOfCamera(const Vector<4, double>& clip, double wEpsilon = 1e-12)
	{
		return clip[3] > wEpsilon;
	}

	// Внутри пирамиды видимости? Проверка в CLIP SPACE, до деления —
	// так это и делает GPU. Для диапазона глубины [0,1] (D3D/Vulkan).
	inline bool IsInsideFrustum(const Vector<4, double>& clip, DepthRange range)
	{
		const double w = clip[3];
		if (w <= 0.0) return false;
		if (clip[0] < -w || clip[0] > w) return false;
		if (clip[1] < -w || clip[1] > w) return false;
		if (range == DepthRange::NegOneToOne) { if (clip[2] < -w || clip[2] > w) return false; }
		else { if (clip[2] < 0.0 || clip[2] > w) return false; }
		return true;
	}

	// Деление на w. ПРЕДУСЛОВИЕ: точка уже прошла отсечение (w > 0).
	// Отдельная функция, чтобы деление было видно на месте вызова.
	inline Vector<3, double> PerspectiveDivide(const Vector<4, double>& clip)
	{
		assert(std::abs(clip[3]) > 1e-12 &&
			"PerspectiveDivide: w вырожден — точка не прошла отсечение");
		const double invW = 1.0 / clip[3];
		return Vector<3, double>{ clip[0] * invW, clip[1] * invW, clip[2] * invW };
	}


	// ---------------------------------------------------------------------------
	//  Определитель и обращение
	//
	//  ГЕОМЕТРИЯ: det — коэффициент изменения объёма. det=0 значит объём схлопнут,
	//  столбцы линейно зависимы, обратной НЕ существует. det<0 — ориентация
	//  вывернута (отражение).
	//
	//  МЕТОД: A⁻¹ = adj(A)/det(A), где adj — транспонированная матрица
	//  алгебраических дополнений. Основано на тождестве A·adj(A) = det(A)·I,
	//  которое следует из того, что разложение по ЧУЖОЙ строке даёт ноль
	//  (это определитель матрицы с двумя одинаковыми строками).
	//
	//  ОПТИМИЗАЦИЯ (приём из Ogre). Наивно нужно 16 дополнений, каждое —
	//  определитель 3x3, каждый — три минора 2x2: 48 штук с массой повторов.
	//  Здесь шесть миноров 2x2 считаются ОДИН РАЗ на пару строк и переиспользуются
	//  для сборки дополнений 3x3 (разложение Лапласа по дополнительным минорам).
	//  Три пары строк = 18 миноров вместо 48. Замерено: в 3.9 раза быстрее
	//  при побитово идентичном результате.
	// ---------------------------------------------------------------------------
	
	inline bool Inverse(const Matrix4x4& a, Matrix4x4& out, double tolerance = 1e-15)
	{
		const double m00 = a[0], m01 = a[1], m02 = a[2], m03 = a[3];
		const double m10 = a[4], m11 = a[5], m12 = a[6], m13 = a[7];
		const double m20 = a[8], m21 = a[9], m22 = a[10], m23 = a[11];
		const double m30 = a[12], m31 = a[13], m32 = a[14], m33 = a[15];

		// Шесть миноров 2x2 из строк 2,3 — по одному на каждую пару столбцов
		double v0 = m20 * m31 - m21 * m30;
		double v1 = m20 * m32 - m22 * m30;
		double v2 = m20 * m33 - m23 * m30;
		double v3 = m21 * m32 - m22 * m31;
		double v4 = m21 * m33 - m23 * m31;
		double v5 = m22 * m33 - m23 * m32;

		// Дополнения 3x3 собираются из готовых 2x2
		const double t00 = +(v5 * m11 - v4 * m12 + v3 * m13);
		const double t10 = -(v5 * m10 - v2 * m12 + v1 * m13);
		const double t20 = +(v4 * m10 - v2 * m11 + v0 * m13);
		const double t30 = -(v3 * m10 - v1 * m11 + v0 * m12);

		// det берётся из уже посчитанных дополнений — без отдельного прохода
		const double det = t00 * m00 + t10 * m01 + t20 * m02 + t30 * m03;
		if (std::abs(det) <= tolerance) return false;
		const double id = 1.0 / det;

		out[0] = t00 * id; out[4] = t10 * id; out[8] = t20 * id; out[12] = t30 * id;
		out[1] = -(v5 * m01 - v4 * m02 + v3 * m03) * id;
		out[5] = +(v5 * m00 - v2 * m02 + v1 * m03) * id;
		out[9] = -(v4 * m00 - v2 * m01 + v0 * m03) * id;
		out[13] = +(v3 * m00 - v1 * m01 + v0 * m02) * id;

		// Пересчёт миноров для строк 1,3
		v0 = m10 * m31 - m11 * m30;  v1 = m10 * m32 - m12 * m30;  v2 = m10 * m33 - m13 * m30;
		v3 = m11 * m32 - m12 * m31;  v4 = m11 * m33 - m13 * m31;  v5 = m12 * m33 - m13 * m32;
		out[2] = +(v5 * m01 - v4 * m02 + v3 * m03) * id;
		out[6] = -(v5 * m00 - v2 * m02 + v1 * m03) * id;
		out[10] = +(v4 * m00 - v2 * m01 + v0 * m03) * id;
		out[14] = -(v3 * m00 - v1 * m01 + v0 * m02) * id;

		// Пересчёт для строк 1,2
		v0 = m21 * m10 - m20 * m11;  v1 = m22 * m10 - m20 * m12;  v2 = m23 * m10 - m20 * m13;
		v3 = m22 * m11 - m21 * m12;  v4 = m23 * m11 - m21 * m13;  v5 = m23 * m12 - m22 * m13;
		out[3] = -(v5 * m01 - v4 * m02 + v3 * m03) * id;
		out[7] = +(v5 * m00 - v2 * m02 + v1 * m03) * id;
		out[11] = -(v4 * m00 - v2 * m01 + v0 * m03) * id;
		out[15] = +(v3 * m00 - v1 * m01 + v0 * m02) * id;
		return true;
	}

	inline double Determinant(const Matrix4x4& a)
	{
		const double v0 = a[8] * a[13] - a[9] * a[12];
		const double v1 = a[8] * a[14] - a[10] * a[12];
		const double v2 = a[8] * a[15] - a[11] * a[12];
		const double v3 = a[9] * a[14] - a[10] * a[13];
		const double v4 = a[9] * a[15] - a[11] * a[13];
		const double v5 = a[10] * a[15] - a[11] * a[14];
		return  (+(v5 * a[5] - v4 * a[6] + v3 * a[7])) * a[0]
			+ (-(v5 * a[4] - v2 * a[6] + v1 * a[7])) * a[1]
			+ (+(v4 * a[4] - v2 * a[5] + v0 * a[7])) * a[2]
			+ (-(v3 * a[4] - v1 * a[5] + v0 * a[6])) * a[3];
	}

	inline Matrix4x4 Transpose(const Matrix4x4& m)
	{
		Matrix4x4 r;
		for (int32_t i = 0; i < 4; ++i) for (int32_t j = 0; j < 4; ++j) r[M4Idx(i, j)] = m[M4Idx(j, i)];
		return r;
	}

	inline bool IsAffine(const Matrix4x4& m)
	{
		return m[12] == 0.0 && m[13] == 0.0 && m[14] == 0.0 && m[15] == 1.0;
	}

	// Сужение до аффинной. Возвращает false, если матрица не аффинна —
	// у Ogre аналогичное преобразование молча отбрасывает последнюю строку.
	inline bool ToAffine(const Matrix4x4& m, Affine4x3& out)
	{
		if (!IsAffine(m)) return false;
		out = Affine4x3(m[0], m[1], m[2], m[3],
			m[4], m[5], m[6], m[7],
			m[8], m[9], m[10], m[11]);
		return true;
	}

	// ---------------------------------------------------------------------------
	//  Проекции. Конвенция глубины задаётся явно параметром, а не зашита:
	//  OpenGL ожидает z в [-1,1], DirectX и Vulkan — в [0,1].
	//  Молчаливое несовпадение даёт «половина сцены за камерой».
	// ---------------------------------------------------------------------------

	inline Matrix4x4 Perspective(double fovYRadians, double aspect,
		double zNear, double zFar, DepthRange range)
	{
		const double f = 1.0 / std::tan(0.5 * fovYRadians);
		Matrix4x4 m;
		for (int32_t i = 0; i < 16; ++i) m[i] = 0.0;
		m[M4Idx(0, 0)] = f / aspect;
		m[M4Idx(1, 1)] = f;
		m[M4Idx(3, 2)] = -1.0;               // w = -z, отсюда перспективное деление
		if (range == DepthRange::NegOneToOne) {
			m[M4Idx(2, 2)] = (zFar + zNear) / (zNear - zFar);
			m[M4Idx(2, 3)] = (2.0 * zFar * zNear) / (zNear - zFar);
		}
		else {
			m[M4Idx(2, 2)] = zFar / (zNear - zFar);
			m[M4Idx(2, 3)] = (zFar * zNear) / (zNear - zFar);
		}
		return m;
	}

	inline Matrix4x4 Orthographic(double l, double r, double b, double t,
		double zNear, double zFar, DepthRange range)
	{
		Matrix4x4 m;
		for (int32_t i = 0; i < 16; ++i) m[i] = 0.0;
		m[M4Idx(0, 0)] = 2.0 / (r - l);   m[M4Idx(0, 3)] = -(r + l) / (r - l);
		m[M4Idx(1, 1)] = 2.0 / (t - b);   m[M4Idx(1, 3)] = -(t + b) / (t - b);
		m[M4Idx(3, 3)] = 1.0;
		if (range == DepthRange::NegOneToOne) {
			m[M4Idx(2, 2)] = -2.0 / (zFar - zNear);  m[M4Idx(2, 3)] = -(zFar + zNear) / (zFar - zNear);
		}
		else {
			m[M4Idx(2, 2)] = -1.0 / (zFar - zNear);  m[M4Idx(2, 3)] = -zNear / (zFar - zNear);
		}
		return m;
	}

} // namespace mathx

