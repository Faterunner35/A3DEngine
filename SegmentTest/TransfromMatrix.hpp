#pragma once
#include "affine4x3.hpp"
#include "quaternion.hpp"
#include"Matrix3x3.hpp"
#include <cmath>
#include <cassert>
#include <optional>

namespace mathx {

	// ============================================================================
	//  transform_ops — сборка и разбор аффинного преобразования.
	//
	//  ПОЧЕМУ ОТДЕЛЬНЫЙ ФАЙЛ. Эти функции тянут Quaternion, а affine4x3.hpp
	//  должен оставаться дешёвым для всех, кому нужна просто матрица (физика,
	//  рендер, обход сцены). Замерено: включение конверсий в матричный заголовок
	//  удваивает время компиляции каждой единицы трансляции.
	//
	//  ЭТО МАТЕМАТИКА, НЕ СЦЕНА. Компонент Transform в слое сцены (родитель,
	//  dirty-флаги, кэш мировой матрицы) — другая сущность; он ИСПОЛЬЗУЕТ
	//  эти функции, но живёт выше.
	//
	//  ПОРЯДОК ОПЕРАЦИЙ: масштаб -> поворот -> перенос, то есть M = T·R·S.
	//  Следствия: перенос не зависит от осей поворота, масштаб не влияет на
	//  величину переноса, поворот и масштаб центрированы в начале координат.
	// ============================================================================

	// Сборка: M = T · R · S.
	// Линейная часть = R·S, а умножение справа на диагональную S масштабирует
	// СТОЛБЦЫ матрицы поворота (столбцы — это локальные оси объекта).
	inline Affine4x3 MakeTransform(const Vector<3, double>& position,
		const Quaternion& orientation,
		const Vector<3, double>& scale)
	{
		const double x = orientation[0], y = orientation[1], z = orientation[2], w = orientation[3];
		const double xx = x * x, yy = y * y, zz = z * z;
		const double xy = x * y, xz = x * z, yz = y * z;
		const double wx = w * x, wy = w * y, wz = w * z;

		// Матрица поворота, каждый столбец сразу умножен на свой масштаб
		const double sx = scale[0], sy = scale[1], sz = scale[2];
		return Affine4x3(
			(1.0 - 2.0 * (yy + zz)) * sx, (2.0 * (xy - wz)) * sy, (2.0 * (xz + wy)) * sz, position[0],
			(2.0 * (xy + wz)) * sx, (1.0 - 2.0 * (xx + zz)) * sy, (2.0 * (yz - wx)) * sz, position[1],
			(2.0 * (xz - wy)) * sx, (2.0 * (yz + wx)) * sy, (1.0 - 2.0 * (xx + yy)) * sz, position[2]);
	}

	// Сборка ОБРАТНОГО преобразования из тех же данных, без обращения матрицы.
	// M⁻¹ = S⁻¹·R⁻¹·T⁻¹: масштаб 1/s, поворот сопряжённым кватернионом,
	// перенос -position, применённый ПОСЛЕ поворота и масштаба.
	// Дешевле и точнее, чем построить M и обратить: нет деления на det.
	inline Affine4x3 MakeInverseTransform(const Vector<3, double>& position,
		const Quaternion& orientation,
		const Vector<3, double>& scale)
	{
		assert(std::abs(scale[0]) > 1e-12 && std::abs(scale[1]) > 1e-12 &&
			std::abs(scale[2]) > 1e-12 && "MakeInverseTransform: нулевой масштаб");

		// Обратный поворот = сопряжение (для единичного кватерниона)
		const Quaternion inv = quat_conjugate(orientation);
		const double x = inv[0], y = inv[1], z = inv[2], w = inv[3];
		const double xx = x * x, yy = y * y, zz = z * z;
		const double xy = x * y, xz = x * z, yz = y * z;
		const double wx = w * x, wy = w * y, wz = w * z;

		const double isx = 1.0 / scale[0], isy = 1.0 / scale[1], isz = 1.0 / scale[2];

		// R⁻¹, затем масштабируем СТРОКИ на 1/s (S⁻¹ слева)
		const double r00 = (1.0 - 2.0 * (yy + zz)) * isx, r01 = (2.0 * (xy - wz)) * isx, r02 = (2.0 * (xz + wy)) * isx;
		const double r10 = (2.0 * (xy + wz)) * isy, r11 = (1.0 - 2.0 * (xx + zz)) * isy, r12 = (2.0 * (yz - wx)) * isy;
		const double r20 = (2.0 * (xz - wy)) * isz, r21 = (2.0 * (yz + wx)) * isz, r22 = (1.0 - 2.0 * (xx + yy)) * isz;

		const double px = position[0], py = position[1], pz = position[2];
		return Affine4x3(
			r00, r01, r02, -(r00 * px + r01 * py + r02 * pz),
			r10, r11, r12, -(r10 * px + r11 * py + r12 * pz),
			r20, r21, r22, -(r20 * px + r21 * py + r22 * pz));
	}

	// Разбор: M -> позиция, ориентация, масштаб.
	//
	// Позиция — прямо из столбца переноса.
	// Масштаб — длины СТОЛБЦОВ линейной части (столбцы = локальные оси).
	// Поворот — те же столбцы, нормированные, через Шумейка в кватернион.
	//
	// ОГРАНИЧЕНИЕ. Разложение неоднозначно, если в матрице есть СДВИГ (shear):
	// тогда столбцы не ортогональны, и представить её как T·R·S нельзя вовсе.
	// Функция это обнаруживает и возвращает false, вместо того чтобы молча
	// выдать бессмысленный результат (Ogre в decomposition() такой проверки
	// не делает). Для матриц со сдвигом нужна QDU или полярная декомпозиция.
	//
	// ОТРИЦАТЕЛЬНЫЙ МАСШТАБ. При det < 0 матрица содержит отражение, которое
	// в тройку (позиция, кватернион, масштаб) с положительным масштабом не
	// раскладывается. Договорённость: знак уходит в X-масштаб.
	struct TransformParts
	{
		Vector<3, double> position{};
		Quaternion       orientation{};
		Vector<3, double> scale{};
	};

	inline bool DecomposeTransform(const Affine4x3& m,
		Vector<3, double>& outPosition,
		Quaternion& outOrientation,
		Vector<3, double>& outScale,
		double shearTolerance = 1e-6)
	{
		outPosition = GetTranslation(m);

		// Столбцы линейной части
		double c0[3] = { m[AIdx(0,0)], m[AIdx(1,0)], m[AIdx(2,0)] };
		double c1[3] = { m[AIdx(0,1)], m[AIdx(1,1)], m[AIdx(2,1)] };
		double c2[3] = { m[AIdx(0,2)], m[AIdx(1,2)], m[AIdx(2,2)] };

		double s0 = std::sqrt(c0[0] * c0[0] + c0[1] * c0[1] + c0[2] * c0[2]);
		double s1 = std::sqrt(c1[0] * c1[0] + c1[1] * c1[1] + c1[2] * c1[2]);
		double s2 = std::sqrt(c2[0] * c2[0] + c2[1] * c2[1] + c2[2] * c2[2]);

		if (s0 < 1e-12 || s1 < 1e-12 || s2 < 1e-12) return false;   // вырожденная

		// Нормируем столбцы -> должен получиться чистый поворот
		for (int i = 0; i < 3; ++i) { c0[i] /= s0; c1[i] /= s1; c2[i] /= s2; }

		// ПРОВЕРКА НА СДВИГ: столбцы обязаны быть попарно ортогональны.
		const double d01 = c0[0] * c1[0] + c0[1] * c1[1] + c0[2] * c1[2];
		const double d02 = c0[0] * c2[0] + c0[1] * c2[1] + c0[2] * c2[2];
		const double d12 = c1[0] * c2[0] + c1[1] * c2[1] + c1[2] * c2[2];
		if (std::abs(d01) > shearTolerance ||
			std::abs(d02) > shearTolerance ||
			std::abs(d12) > shearTolerance) {
			return false;   // есть сдвиг — раскладывается не в T·R·S
		}

		// Отражение: det < 0. Знак уходит в первый масштаб.
		const double det = c0[0] * (c1[1] * c2[2] - c1[2] * c2[1])
			- c1[0] * (c0[1] * c2[2] - c0[2] * c2[1])
			+ c2[0] * (c0[1] * c1[2] - c0[2] * c1[1]);
		if (det < 0.0) {
			s0 = -s0;
			for (int i = 0; i < 3; ++i) c0[i] = -c0[i];
		}

		outScale = Vector<3, double>{ s0, s1, s2 };

		// Шумейк по нормированным столбцам (R[r][c]: столбцы это c0,c1,c2)
		const double R00 = c0[0], R10 = c0[1], R20 = c0[2];
		const double R01 = c1[0], R11 = c1[1], R21 = c1[2];
		const double R02 = c2[0], R12 = c2[1], R22 = c2[2];

		Quaternion q;
		const double trace = R00 + R11 + R22;
		if (trace > 0.0) {
			double root = std::sqrt(trace + 1.0);
			q[3] = 0.5 * root; root = 0.5 / root;
			q[0] = (R21 - R12) * root; q[1] = (R02 - R20) * root; q[2] = (R10 - R01) * root;
		}
		else {
			const double d[3] = { R00, R11, R22 };
			const double Rm[9] = { R00,R01,R02, R10,R11,R12, R20,R21,R22 };
			static constexpr int32_t nxt[3] = { 1,2,0 };
			int32_t i = 0; if (d[1] > d[0]) i = 1; if (d[2] > d[i]) i = 2;
			const int32_t j = nxt[i], k = nxt[j];
			double root = std::sqrt(Rm[i * 3 + i] - Rm[j * 3 + j] - Rm[k * 3 + k] + 1.0);
			q[i] = 0.5 * root; root = 0.5 / root;
			q[3] = (Rm[k * 3 + j] - Rm[j * 3 + k]) * root;
			q[j] = (Rm[j * 3 + i] + Rm[i * 3 + j]) * root;
			q[k] = (Rm[k * 3 + i] + Rm[i * 3 + k]) * root;
		}
		outOrientation = q;
		return true;
	}


	// Удобная форма: не даёт проигнорировать неудачу разбора (сдвиг/вырождение).
	inline std::optional<TransformParts> TryDecomposeTransform(
		const Affine4x3& m, double shearTolerance = 1e-6)
	{
		TransformParts p;
		if (!DecomposeTransform(m, p.position, p.orientation, p.scale, shearTolerance))
			return std::nullopt;
		return p;
	}


	void CreateRotorMatrix(Matrix3x3& mat, double ax, double ay, double az, double angle)
	{
		double s = std::sin(angle), c = std::cos(angle), d = 1.0 - c;
		mat[0] = c + d * (ax * ax);      mat[1] = d * (ax * ay) - s * az;   mat[2] = d * (ax * az) + s * ay;
		mat[3] = d * (ax * ay) + s * az;   mat[4] = c + d * (ay * ay);      mat[5] = d * (ay * az) - s * ax;
		mat[6] = d * (ax * az) - s * ay;   mat[7] = d * (ay * az) + s * ax;   mat[8] = c + d * (az * az);
	}

	Matrix3x3 CreateRotMatFromQuat(const Quaternion& quat)
	{


		Matrix3x3 mat_rot;

		// Все переменные строго double
		double x2 = quat[0] + quat[0];
		double y2 = quat[1] + quat[1];
		double z2 = quat[2] + quat[2];

		// Предвычисляем попарные произведения
		double xx = quat[0] * x2;   double xy = quat[0] * y2;   double xz = quat[0] * z2;
		double yy = quat[1] * y2;   double yz = quat[1] * z2;   double zz = quat[2] * z2;
		double wx = quat[3] * x2;   double wy = quat[3] * y2;   double wz = quat[3] * z2;

		// Главная диагональ
		mat_rot[0] = 1.0 - (yy + zz);
		mat_rot[4] = 1.0 - (xx + zz);
		mat_rot[8] = 1.0 - (xx + yy);

		// Внедиагональные элементы (Row-Major)
		mat_rot[1] = xy - wz;
		mat_rot[2] = xz + wy;
		mat_rot[3] = xy + wz;
		mat_rot[5] = yz - wx;
		mat_rot[6] = xz - wy;
		mat_rot[7] = yz + wx;

		return mat_rot;

	}
	

	inline Matrix3x3 OuterProduct(const Vector<3, double>& a, const Vector<3, double>& b)
	{
		return Matrix3x3(a[0] * b[0], a[0] * b[1], a[0] * b[2],
			a[1] * b[0], a[1] * b[1], a[1] * b[2],
			a[2] * b[0], a[2] * b[1], a[2] * b[2]);
	}

	// ---------------------------------------------------------------------------
	//  Отражение через плоскость, проходящую ЧЕРЕЗ НАЧАЛО КООРДИНАТ,
	//  перпендикулярную единичному вектору a.
	//
	//  ВЫВОД: раскладываем v = v_перп + v_парал относительно a и меняем знак
	//  у параллельной части:  v' = v_перп − v_парал.
	//    v_парал = (a·v)a = (a ⊗ a)v
	//    v_перп  = v − (a ⊗ a)v = (I − a ⊗ a)v
	//    v' = (I − a⊗a)v − (a⊗a)v = (I − 2 a⊗a) v
	//
	//  det = −1 всегда: по тождеству det(I + a⊗b) = 1 + a·b при b = −2a
	//  получаем 1 + a·(−2a) = 1 − 2 = −1.
	//
	//  ПРЕДУСЛОВИЕ: a единичный.
	// ---------------------------------------------------------------------------
	inline Matrix3x3 MakeReflection(const Vector<3, double>& a)
	{
		assert(std::abs(a[0] * a[0] + a[1] * a[1] + a[2] * a[2] - 1.0) < 1e-9 &&
			"MakeReflection: нормаль должна быть единичной");
		// Форма Ленгьеля: -2a считается ОДИН раз, матрица симметрична,
		// поэтому каждый внедиагональный элемент вычисляется по разу.
		// 9 умножений вместо 12-15 при наивной записи формулы.
		const double x = a[0] * -2.0, y = a[1] * -2.0, z = a[2] * -2.0;
		const double axay = x * a[1], axaz = x * a[2], ayaz = y * a[2];
		return Matrix3x3(x * a[0] + 1.0, axay, axaz,
			axay, y * a[1] + 1.0, ayaz,
			axaz, ayaz, z * a[2] + 1.0);
	}

	// ---------------------------------------------------------------------------
	//  Отражение через плоскость на ПРОИЗВОЛЬНОМ расстоянии от начала координат.
	//  Плоскость задана уравнением a·p + d = 0 (a единичный).
	//
	//  Схема: перенести плоскость в начало, отразить, вернуть обратно.
	//    M = T(−d·a) · Reflect(a) · T(d·a)
	//  Раскрывая, столбец переноса получается как −2d·a: точка на расстоянии
	//  d от плоскости уходит на −d, то есть сдвигается на 2d против нормали.
	// ---------------------------------------------------------------------------
	inline Affine4x3 MakeReflectionPlane(const Vector<3, double>& a, double d)
	{
		const Matrix3x3 R = MakeReflection(a);
		return Affine4x3(R[0], R[1], R[2], -2.0 * d * a[0],
			R[3], R[4], R[5], -2.0 * d * a[1],
			R[6], R[7], R[8], -2.0 * d * a[2]);
	}

	// Плоскость через точку p с нормалью a: d = −(a·p).
	inline Affine4x3 MakeReflectionThroughPoint(const Vector<3, double>& a,
		const Vector<3, double>& point)
	{
		return MakeReflectionPlane(a, -(a[0] * point[0] + a[1] * point[1] + a[2] * point[2]));
	}

	// ---------------------------------------------------------------------------
	//  Involution — обратный случай: гасится ПЕРПЕНДИКУЛЯРНАЯ часть,
	//  параллельная остаётся:  v' = v_парал − v_перп = (2 a⊗a − I) v.
	//
	//  Это в точности ОТРИЦАНИЕ матрицы отражения: M_invol = −M_reflect.
	//
	//  И отрицание меняет знак определителя: в 3D оно умножает каждую из трёх
	//  строк на −1, значит det умножается на (−1)³ = −1. Поэтому отражение
	//  (det = −1) превращается в ПОВОРОТ (det = +1), а именно — в поворот на
	//  180° вокруг оси a (проверено сравнением с кватернионом: расхождение 9e-16).
	//
	//  Название «involution» — из свойства M² = I: применить дважды значит
	//  вернуться. Отражения и повороты на 180° — её примеры.
	// ---------------------------------------------------------------------------
	inline Matrix3x3 MakeInvolution(const Vector<3, double>& a)
	{
		assert(std::abs(a[0] * a[0] + a[1] * a[1] + a[2] * a[2] - 1.0) < 1e-9 &&
			"MakeInvolution: axis must be normalized");

		const double x2 = a[0] * 2.0;
		const double y2 = a[1] * 2.0;
		const double z2 = a[2] * 2.0;

		const double xy2 = x2 * a[1];
		const double xz2 = x2 * a[2];
		const double yz2 = y2 * a[2];

		return Matrix3x3(
			x2 * a[0] - 1.0, xy2, xz2,
			xy2, y2 * a[1] - 1.0, yz2,
			xz2, yz2, z2 * a[2] - 1.0
		);
	}

	// ---------------------------------------------------------------------------
	//  Масштабирование ВДОЛЬ произвольного направления: параллельная часть
	//  умножается на s, перпендикулярная не трогается.
	//    v' = v_перп + s·v_парал = (I + (s−1) a⊗a) v
	//  Отражение — частный случай при s = −1, involution — при s = ... нет,
	//  involution требует масштаба −1 у ПЕРПЕНДИКУЛЯРНОЙ части.
	// ---------------------------------------------------------------------------
	

	inline Matrix3x3 MakeScaleAlong(const Vector<3, double>& a, double s)
	{
		assert(std::abs(DotProduct(a, a) - 1.0) < 1e-9 &&
			"MakeScaleAlong: направление должно быть единичным");
		const double k = s - 1.0;
		const double x = a[0] * k, y = a[1] * k, z = a[2] * k;   // (s−1)·a, один раз
		const double axay = x * a[1], axaz = x * a[2], ayaz = y * a[2];
		return Matrix3x3(x * a[0] + 1.0, axay, axaz,
			axay, y * a[1] + 1.0, ayaz,
			axaz, ayaz, z * a[2] + 1.0);
	}

	// Масштаб вдоль ТРЁХ взаимно ортогональных направлений.
	// M = s1·(u1⊗u1) + s2·(u2⊗u2) + s3·(u3⊗u3)
	//
	// Это спектральное разложение: проекторы на собственные направления,
	// взвешенные собственными значениями. det = s1·s2·s3.
	// Обычный осевой масштаб — частный случай при u_i = базисные орты.
	inline Matrix3x3 MakeScaleAlongAxes(const Vector<3, double>& u1, double s1,
		const Vector<3, double>& u2, double s2,
		const Vector<3, double>& u3, double s3)
	{
		assert(std::abs(DotProduct(u1, u2)) < 1e-9 &&
			std::abs(DotProduct(u1, u3)) < 1e-9 &&
			std::abs(DotProduct(u2, u3)) < 1e-9 &&
			"MakeScaleAlongAxes: оси должны быть взаимно ортогональны");
		Matrix3x3 m;
		for (int32_t r = 0; r < 3; ++r)
			for (int32_t c = 0; c < 3; ++c)
				m[r * 3 + c] = s1 * u1[r] * u1[c] + s2 * u2[r] * u2[c] + s3 * u3[r] * u3[c];
		return m;
	}

	// Осевой масштаб уже есть в matrix3x3_ops.hpp как MakeScale(sx,sy,sz).
	// Он ДЕШЕВЛЕ обеих формул выше: применение к вектору стоит 3 умножения
	// вместо 9, и хранить надо 3 числа вместо 9. Общие формулы нужны только
	// для НАКЛОННЫХ осей.

	// ---------------------------------------------------------------------------
	//  СДВИГ (skew / shear).
	//
	//      v' = v + a·(b·v)·tanθ
	//
	//  Читается так: ИЗМЕРИТЬ проекцию v на направление b, СДВИНУТЬ на
	//  пропорциональную величину вдоль направления a. Величина сдвига линейна
	//  по расстоянию вдоль b — отсюда наклон, а не параллельный перенос.
	//
	//  Почему tanθ, а не sinθ: сдвиг это ОТНОШЕНИЕ смещения вдоль a к
	//  расстоянию вдоль b, то есть противолежащий катет к прилежащему.
	//  При θ → 90° уходит в бесконечность и фигура вырождается.
	//
	//  СВОЙСТВА:
	//    det = 1 — объём сохраняется точно
	//    углы НЕ сохраняются — квадрат становится параллелограммом
	//    a⊗b при a ⊥ b НЕсимметрична, в отличие от всех формул выше
	//
	//  Отсюда главное следствие: сдвиг НЕ раскладывается на поворот и масштаб.
	//  В движке он появляется сам — от неравномерного масштаба родителя,
	//  применённого к повёрнутому ребёнку.
	//
	//  ПРЕДУСЛОВИЕ: a и b единичные и взаимно перпендикулярные.
	// ---------------------------------------------------------------------------
	inline Matrix3x3 MakeSkew(double angleRadians,
		const Vector<3, double>& a,
		const Vector<3, double>& b)
	{
		assert(std::abs(DotProduct(a, a) - 1.0) < 1e-9 &&
			std::abs(DotProduct(b, b) - 1.0) < 1e-9 &&
			"MakeSkew: направления должны быть единичными");
		assert(std::abs(DotProduct(a, b)) < 1e-9 &&
			"MakeSkew: a и b должны быть перпендикулярны");
		assert(std::abs(angleRadians) < 1.5707 &&
			"MakeSkew: угол должен быть строго меньше 90° — tan уходит в бесконечность");

		const double t = std::tan(angleRadians);
		const double ax = a[0] * t, ay = a[1] * t, az = a[2] * t;   // a·tanθ, один раз
		return Matrix3x3(ax * b[0] + 1.0, ax * b[1], ax * b[2],
			ay * b[0], ay * b[1] + 1.0, ay * b[2],
			az * b[0], az * b[1], az * b[2] + 1.0);
	}

	// Сдвиг, заданный тангенсом напрямую — когда угол не нужен.
	// Избавляет от вызова tan: в анимации и редакторе величину сдвига
	// часто задают отношением, а не углом.
	inline Matrix3x3 MakeSkewByTangent(double tangent,
		const Vector<3, double>& a,
		const Vector<3, double>& b)
	{
		assert(std::abs(DotProduct(a, b)) < 1e-9 && "MakeSkewByTangent: a ⊥ b");
		const double ax = a[0] * tangent, ay = a[1] * tangent, az = a[2] * tangent;
		return Matrix3x3(ax * b[0] + 1.0, ax * b[1], ax * b[2],
			ay * b[0], ay * b[1] + 1.0, ay * b[2],
			az * b[0], az * b[1], az * b[2] + 1.0);
	}

	// Ортогональны ли столбцы матрицы.
	//
	// ЧТО ЭТО ЗНАЧИТ НА САМОМ ДЕЛЕ: столбцы ортогональны тогда и только тогда,
	// когда матрица раскладывается на ПОВОРОТ и ОСЕВОЙ масштаб (R·S с
	// диагональной S). Именно это проверяет DecomposeTransform.
	//
	// ВАЖНО, что это НЕ «есть ли сдвиг»:
	//   сдвиг                     -> столбцы неортогональны  (не раскладывается)
	//   масштаб вдоль НАКЛОННОЙ оси -> столбцы тоже неортогональны!
	//
	// Второй случай не сдвиг: матрица симметрична и раскладывается как
	// R·S·Rᵀ (спектрально), просто НЕ в виде R·S с диагональной S. Отличить
	// одно от другого можно по симметрии: масштаб вдоль любой оси симметричен,
	// сдвиг — нет.
	inline bool HasNonOrthogonalColumns(const Matrix3x3& m, double tolerance = 1e-6)
	{
		Vector<3, double> c0{ m[0],m[3],m[6] }, c1{ m[1],m[4],m[7] }, c2{ m[2],m[5],m[8] };
		const double l0 = std::sqrt(DotProduct(c0, c0));
		const double l1 = std::sqrt(DotProduct(c1, c1));
		const double l2 = std::sqrt(DotProduct(c2, c2));
		if (l0 < 1e-12 || l1 < 1e-12 || l2 < 1e-12) return false;   // вырождена
		return std::abs(DotProduct(c0, c1)) / (l0 * l1) > tolerance
			|| std::abs(DotProduct(c0, c2)) / (l0 * l2) > tolerance
			|| std::abs(DotProduct(c1, c2)) / (l1 * l2) > tolerance;
	}

	// Симметрична ли матрица. Все преобразования на a⊗a (масштаб, отражение,
	// проекция) симметричны; сдвиг на a⊗b при a ⊥ b — нет.
	// Вместе с проверкой выше это отличает сдвиг от наклонного масштаба.
	inline bool IsSymmetric(const Matrix3x3& m, double tolerance = 1e-9)
	{
		return std::abs(m[1] - m[3]) <= tolerance
			&& std::abs(m[2] - m[6]) <= tolerance
			&& std::abs(m[5] - m[7]) <= tolerance;
	}

	// Настоящая проверка на сдвиг: столбцы неортогональны И матрица несимметрична.
	inline bool HasSkew(const Matrix3x3& m, double tolerance = 1e-6)
	{
		return HasNonOrthogonalColumns(m, tolerance) && !IsSymmetric(m);
	}

}