#include "QuatInterpolationPrecomputed.h"

namespace mathx {

	std::vector<SlerpSegment> PrepareSlerpSequence(std::span<Quaternion> keys)
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

	void SlerpBatch(double t,
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

} // namespace mathx
