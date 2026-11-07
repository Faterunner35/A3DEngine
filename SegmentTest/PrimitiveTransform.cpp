#include "PrimitiveTransform.h"

namespace mathx {

	Matrix3x3 AxisRotation(int32_t axis, double angle)
	{
		const double c = std::cos(angle), s = std::sin(angle);
		const int32_t j = (axis + 1) % 3, k = (axis + 2) % 3;

		Matrix3x3 m;                        // конструктор по умолчанию обнуляет
		m(axis, axis) = 1.0;
		m(j, j) = c;   m(j, k) = -s;
		m(k, j) = s;   m(k, k) = c;
		return m;
	}

	void PreRotate(Matrix3x3& m, int32_t axis, double angle)
	{
		const double c = std::cos(angle), s = std::sin(angle);
		const int32_t j = (axis + 1) % 3, k = (axis + 2) % 3;

		for (int32_t col = 0; col < 3; ++col) {
			const double vj = m(j, col), vk = m(k, col);
			m(j, col) = c * vj - s * vk;
			m(k, col) = s * vj + c * vk;
		}
	}

	Matrix3x3 EulerToMatrix(const EulerAngles& e)
	{
		Matrix3x3 m = AxisRotation(e.axis[2], e.angle[2]);
		PreRotate(m, e.axis[1], e.angle[1]);
		PreRotate(m, e.axis[0], e.angle[0]);
		return m;
	}

} // namespace mathx
