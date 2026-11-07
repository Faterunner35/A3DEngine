#include "QuatSquad.h"

namespace mathx {

	void quat_intermediate(Quaternion const& q0,
		Quaternion const& q1,
		Quaternion const& q2,
		Quaternion& outA,
		Quaternion& outB)
	{
		assert(std::fabs(quat_length(q0) - 1.0) < 1e-6 &&
			std::fabs(quat_length(q1) - 1.0) < 1e-6 &&
			std::fabs(quat_length(q2) - 1.0) < 1e-6 &&
			"quat_intermediate: inputs must be unit quaternions");

		const Quaternion P0 = quat_conjugate(q0) * q1;   // для единичного conj == inverse
		const Quaternion P1 = quat_conjugate(q1) * q2;

		const Quaternion l0 = quat_log(P0);
		const Quaternion l1 = quat_log(P1);
		const Quaternion arg = (l0 - l1) * 0.25;

		outA = q1 * quat_exp(arg);
		outB = q1 * quat_exp(-arg);
	}

	Quaternion quat_squad(double t,
		Quaternion const& p,
		Quaternion const& a,
		Quaternion const& b,
		Quaternion const& q,
		bool shortestPath)
	{
		const double outerT = 2.0 * t * (1.0 - t);
		const Quaternion s1 = quat_slerp(t, p, q, shortestPath);   // прямая дуга
		const Quaternion s2 = quat_slerp(t, a, b, shortestPath);   // дуга контрольных точек
		return quat_slerp(outerT, s1, s2, shortestPath);
	}

} // namespace mathx
