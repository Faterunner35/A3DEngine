#pragma once
#include "Matrxi4x4.hpp"
#include <cassert>
#include <cmath>

namespace mathx {

	// ============================================================================
	//  Перспективные проекции в соглашениях FGED (Lengyel, т. 2, гл. 6).
	//
	//  Пространство камеры: x вправо, y ВНИЗ, z ВПЕРЁД (правая тройка).
	//  Вьюпорт — прямоугольник y в [-1,1], x в [-s,s] на плоскости z = g, где
	//      g = 1 / tan(fovY/2)          (вид сбоку: катеты g и 1)
	//      s = ширина / высота
	//  Хранение по строкам, точка — столбец: p_clip = P * p_camera, w_clip = z.
	//
	//  Клип-пространство по y у API разное:
	//      ClipY::Down — Vulkan; OpenGL с glClipControl(GL_UPPER_LEFT, ...)
	//      ClipY::Up   — D3D12, Metal: вторая строка берётся со знаком минус (g -> -g)
	//  Глубина в обоих случаях [0,1] (для OpenGL — glClipControl(..., GL_ZERO_TO_ONE)).
	//
	//  Старая Perspective() в Matrxi4x4.hpp — в соглашениях OpenGL (камера смотрит
	//  в -z, y вверх) и с этими функциями не смешивается.
	// ============================================================================

	enum class ClipY {
		Down,   // Vulkan, OpenGL + glClipControl(GL_UPPER_LEFT)
		Up      // D3D12, Metal
	};

	// (6.36), листинг 6.1 FGED: ближняя плоскость -> глубина 0, дальняя -> 1.
	//     z_dev = A + B/z,   A = f/(f-n),   B = -n f/(f-n)
	inline Matrix4x4 MakeFrustumProjection(double fovY, double s, double n, double f,
		ClipY clipY = ClipY::Up)
	{
		assert(fovY > 0.0 && fovY < 3.14159265358979 && s > 0.0 && n > 0.0 && f > n &&
			"MakeFrustumProjection: нужны 0 < fovY < pi, s > 0, 0 < n < f");
		const double g = 1.0 / std::tan(fovY * 0.5);
		const double gy = (clipY == ClipY::Up) ? -g : g;
		const double k = f / (f - n);
		return Matrix4x4(g / s, 0.0, 0.0, 0.0,
			0.0, gy, 0.0, 0.0,
			0.0, 0.0, k, -n * k,
			0.0, 0.0, 1.0, 0.0);
	}

	// Обратная бесконечная проекция: ближняя плоскость -> 1, бесконечность -> 0.
	//     z_dev = n / z
	// Даёт на порядки лучшую точность float-глубины, чем (6.36), и не требует
	// поправки epsilon из (6.42): для направления (w = 0) строка (0,0,0,n) даёт
	// z_clip = 0 ровно. Требует: формат D32_FLOAT, тест GREATER_EQUAL, очистку 0.
	inline Matrix4x4 MakeReversedInfiniteProjection(double fovY, double s, double n,
		ClipY clipY = ClipY::Up)
	{
		assert(fovY > 0.0 && fovY < 3.14159265358979 && s > 0.0 && n > 0.0 &&
			"MakeReversedInfiniteProjection: нужны 0 < fovY < pi, s > 0, n > 0");
		const double g = 1.0 / std::tan(fovY * 0.5);
		const double gy = (clipY == ClipY::Up) ? -g : g;
		return Matrix4x4(g / s, 0.0, 0.0, 0.0,
			0.0, gy, 0.0, 0.0,
			0.0, 0.0, 0.0, n,
			0.0, 0.0, 1.0, 0.0);
	}

	// ---------------------------------------------------------------------------
	//  Обратные матрицы в явном виде — (6.37) FGED. Общее обращение 4x4 не нужно:
	//  матрица блочно-диагональная, блоки обращаются по отдельности.
	//
	//      P = | diag(g/s, gy)   0 |        P^-1 = | diag(s/g, 1/gy)   0     |
	//          |      0          M |               |      0          M^-1    |
	//
	//  Применение: экранная точка (X, Y, z_dev, 1) -> P^-1 -> деление на w
	//  даёт точку в пространстве камеры (x, y, z). Так по буферу глубины
	//  восстанавливают позицию в deferred-рендере.
	// ---------------------------------------------------------------------------

	// Обратная к MakeFrustumProjection.
	//     M = | A  B |,   M^-1 = | 0          1   |
	//         | 1  0 |           | (n-f)/(nf) 1/n |
	inline Matrix4x4 MakeFrustumProjectionInverse(double fovY, double s, double n, double f,
		ClipY clipY = ClipY::Up)
	{
		assert(fovY > 0.0 && fovY < 3.14159265358979 && s > 0.0 && n > 0.0 && f > n &&
			"MakeFrustumProjectionInverse: нужны 0 < fovY < pi, s > 0, 0 < n < f");
		const double g = 1.0 / std::tan(fovY * 0.5);
		const double gy = (clipY == ClipY::Up) ? -g : g;
		return Matrix4x4(s / g, 0.0, 0.0, 0.0,
			0.0, 1.0 / gy, 0.0, 0.0,
			0.0, 0.0, 0.0, 1.0,
			0.0, 0.0, (n - f) / (n * f), 1.0 / n);
	}

	// Обратная к MakeReversedInfiniteProjection.
	//     M = | 0  n |,   M^-1 = | 0    1 |
	//         | 1  0 |           | 1/n  0 |
	inline Matrix4x4 MakeReversedInfiniteProjectionInverse(double fovY, double s, double n,
		ClipY clipY = ClipY::Up)
	{
		assert(fovY > 0.0 && fovY < 3.14159265358979 && s > 0.0 && n > 0.0 &&
			"MakeReversedInfiniteProjectionInverse: нужны 0 < fovY < pi, s > 0, n > 0");
		const double g = 1.0 / std::tan(fovY * 0.5);
		const double gy = (clipY == ClipY::Up) ? -g : g;
		return Matrix4x4(s / g, 0.0, 0.0, 0.0,
			0.0, 1.0 / gy, 0.0, 0.0,
			0.0, 0.0, 0.0, 1.0,
			0.0, 0.0, 1.0 / n, 0.0);
	}

	// Глубина в пространстве камеры по значению из буфера глубины.
	// (6.38): обращение z_dev = A + B/z для MakeFrustumProjection.
	inline double CameraDepthFromFrustum(double zDevice, double n, double f)
	{
		return n * f / ((n - f) * zDevice + f);
	}

	// Обращение z_dev = n/z для MakeReversedInfiniteProjection.
	// zDevice = 0 — точка на бесконечности.
	inline double CameraDepthFromReversedInfinite(double zDevice, double n)
	{
		return n / zDevice;
	}

} // namespace mathx
