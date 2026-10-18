#pragma once
#include"Primitives.h"
#include"Plane.h"
#include"Ray.h"
#include"Sphere.h"
#include<tuple>



namespace mathx
{ 

class GeomCoreIntersector
{
public:

	GeomCoreIntersector() = default;

	struct Result
	{
		Result()
			:
			distance(static_cast<double>(0)),
			sqrDistance(static_cast<double>(0)),
			parameter{ static_cast<double>(0), static_cast<double>(0) },
			closest{ Vector<3,double>::Zero(), Vector<3,double>::Zero() }
		{
		}

		double distance, sqrDistance;
		std::array<double, 2> parameter;
		std::array<Vector<3, double>, 2> closest;
	};

	inline std::pair<Vector<3, double>, bool> PlaneRayIntersection(const Plane& plane_, const Segment<3, double>& seg)
	{
		Vector<3, double> point_intersect;
		Vector<3,double> normal_pl=plane_.GetNormal();
		double scalar_dot_1 = DotProduct(normal_pl, seg.p[1]);

		if (scalar_dot_1 >= std::numeric_limits<double>::epsilon())
		{
			point_intersect = seg.p[0] - seg.p[1] * DotProduct(seg.p[0], normal_pl) / scalar_dot_1;

			return std::make_pair(point_intersect, true);

		}

		else
		{
			return std::make_pair(point_intersect, false);
		}


	}

	double RayTriangleIntersect(const Ray& ray, const Triangle& tri);

	//void static convex_hull_creation(std::vector<Vector<3, double>> points_vector, Triangle triangle, std::vector<Triangle>& _array_triangles);


	std::pair<bool, Vector<3,double>> RaySphereIntersection(const Ray& ray, const Sphere& s);

	double ProjectPointOnSegment(const Segment<3, double>& segment, const Vector<3, double>& point)
	{

		Vector<3, double> dist_segment = segment.p[1] - segment.p[0];
		Vector<3, double> point_segment_dist = point - segment.p[0];

		Vector<3, double> a = Cross(point_segment_dist, segment.p[1]);
		return (std::sqrt(DotProduct(a, a) / DotProduct(segment.p[1], segment.p[1])));

	}


	Vector<3, double> PointSegmentDistance(Vector<3, double>& point, Segment<3, double>& segment)
	{
		Vector<3, double>closest_point;
		Vector<3, double> dist_segment = segment.p[1] - segment.p[0];
		double t = DotProduct(point - segment.p[1], dist_segment);
		if (t > 0.0f)
			t = 1.0f;
		else
			t = 0.0f;

	    closest_point = segment.p[0] + t * dist_segment;
		return closest_point;

	}

	double PriciseSegmentSegmentDistance(Segment<3,double>& segment_1, Segment<3,double>& segment_2)
	{
		auto dot_1 = DotProduct(segment_1.p[1], segment_1.p[1]);
		auto dot_2 = DotProduct(segment_2.p[1], segment_2.p[1]);
		auto dot_3 = DotProduct(segment_1.p[1], segment_2.p[1]);
		auto segnet_diff = segment_2.p[0] - segment_1.p[0];

		auto determinant = 1 / std::pow(DotProduct(segment_1.p[1], segment_2.p[1]), 2) - dot_1 * dot_2;

		auto t_1 = determinant * (((-dot_2 * dot_2) * DotProduct(segnet_diff, segment_1.p[1])) + (-dot_3 * DotProduct(segnet_diff, segment_1.p[1])));

		auto t_2 = determinant * ((dot_3 * DotProduct(segment_1.p[1], segnet_diff) + (dot_1 * dot_1) * DotProduct(segnet_diff, segment_2.p[1])));

		return std::sqrt(DotProduct((segnet_diff + (segment_2.p[1] * t_2) - (segment_1.p[1] * t_1)), (segnet_diff + (segment_2.p[1] * t_2) - (segment_1.p[1] * t_1))));

	}

	double DistanceSegment(Segment<3, double>& segment_1, Segment<3, double>& segment_2)
	{

		auto dp = segment_2.p[0] - segment_1.p[0];
		auto v12 = DotProduct(segment_1.p[1], segment_1.p[1]);
		auto v22 = DotProduct(segment_2.p[1], segment_2.p[1]);
		auto v1v2 = DotProduct(segment_1.p[1], segment_2.p[1]);
		auto det = v1v2 * v1v2 - v12 * v22;

		if (std::fabs(det) > FLT_MIN)
			det = 1.0F / det;
		double dpvl = DotProduct(dp, segment_1.p[1]);
		double dpv2 = DotProduct(dp, segment_2.p[1]);
		double tl = (v1v2 * dpv2 - v22 * dpvl) * det;
		double t2 = (v12 * dpv2 - v1v2 * dpvl) * det;
		return std::sqrt(DotProduct((dp + (segment_2.p[1] * t2) - (segment_1.p[1] * tl)), (dp + (segment_2.p[1] * t2) - (segment_1.p[1] * tl))));

	}



	template<typename T>
	bool Intersects(std::array<Vector<3, double>, 3> U,
		std::array<Vector<3, double>, 3>V, Vector<3, T> segment[2]) const
	{
		// Compute the plane normal for triangle U.
		Vector<3, T> edge1 = U[1] - U[0];
		Vector<3, T> edge2 = U[2] - U[0];
		Vector<3, T> normal = UnitCross(edge1, edge2);

		// Test whether the edges of triangle V transversely intersect the
		// plane of triangle U.
		double d[3];
		int32_t positive = 0, negative = 0, zero = 0;
		for (int32_t i = 0; i < 3; ++i)
		{
			d[i] = DotProduct(normal, V[i] - U[0]);
			if (d[i] > 0.0f)
			{
				++positive;
			}
			else if (d[i] < 0.0f)
			{
				++negative;
			}
			else
			{
				++zero;
			}
		}
		// positive + negative + zero == 3
		//We will use linear interpolation to find segmets c=a+t*b,t-interpolation factor,d will be scalar factor here
		if (positive > 0 && negative > 0)
		{
			if (positive == 2)  // and negative == 1
			{
				if (d[0] < 0.0f)
				{
					segment[0] = (d[1] * V[0] - d[0] * V[1]) / (d[1] - d[0]);
					segment[1] = (d[2] * V[0] - d[0] * V[2]) / (d[2] - d[0]);
				}
				else if (d[1] < 0.0f)
				{
					segment[0] = (d[0] * V[1] - d[1] * V[0]) / (d[0] - d[1]);
					segment[1] = (d[2] * V[1] - d[1] * V[2]) / (d[2] - d[1]);
				}
				else  // d[2] < 0.0f
				{
					segment[0] = (d[0] * V[2] - d[2] * V[0]) / (d[0] - d[2]);
					segment[1] = (d[1] * V[2] - d[2] * V[1]) / (d[1] - d[2]);
				}
			}
			else if (negative == 2)  // and positive == 1
			{
				if (d[0] > 0.0f)
				{
					segment[0] = (d[1] * V[0] - d[0] * V[1]) / (d[1] - d[0]);
					segment[1] = (d[2] * V[0] - d[0] * V[2]) / (d[2] - d[0]);
				}
				else if (d[1] > 0.0f)
				{
					segment[0] = (d[0] * V[1] - d[1] * V[0]) / (d[0] - d[1]);
					segment[1] = (d[2] * V[1] - d[1] * V[2]) / (d[2] - d[1]);
				}
				else  // d[2] > 0.0f
				{
					segment[0] = (d[0] * V[2] - d[2] * V[0]) / (d[0] - d[2]);
					segment[1] = (d[1] * V[2] - d[2] * V[1]) / (d[1] - d[2]);
				}
			}
			else  // positive == 1, negative == 1, zero == 1
			{
				if (d[0] == 0.0f)
				{
					segment[0] = V[0];
					segment[1] = (d[2] * V[1] - d[1] * V[2]) / (d[2] - d[1]);
				}
				else if (d[1] == 0.0f)
				{
					segment[0] = V[1];
					segment[1] = (d[0] * V[2] - d[2] * V[0]) / (d[0] - d[2]);
				}
				else  // d[2] == 0.0f
				{
					segment[0] = V[2];
					segment[1] = (d[1] * V[0] - d[0] * V[1]) / (d[1] - d[0]);
				}
			}
			return true;
		}

		// Triangle V does not transversely intersect triangle U, although it is
		// possible a vertex or edge of V is just touching U.  In this case, we
		// do not call this an intersection.
		return false;
	}

	bool TriangleIntersection(std::array<Vector<3, double>, 3> tri_1, std::array<Vector<3, double>, 3> tri_2) const
	{
		Vector<3, double> S0[2], S1[2];
		if (Intersects(tri_1, tri_2, S0) && Intersects(tri_1, tri_2, S1))
		{
			// Theoretically, the segments lie on the same line.  A direction D
			// of the line is the Cross(NormalOf(U),NormalOf(V)).  We choose the
			// average A of the segment endpoints as the line origin.
			Vector<3, double> uNormal = Cross(tri_1[1] - tri_1[0], tri_1[2] - tri_1[0]);
			Vector<3, double> vNormal = Cross(tri_2[1] - tri_2[0], tri_2[2] - tri_2[0]);
			Vector<3, double> D = UnitCross(uNormal, vNormal);
			Vector<3, double> A = 0.25 * (S0[0] + S0[1] + S1[0] + S1[1]);

			// Each segment endpoint is of the form A + t*D.  Compute the
			// t-values to obtain I0 = [t0min,t0max] for S0 and I1 = [t1min,t1max]
			// for S1.  The segments intersect when I0 overlaps I1.  Although this
			// application acts as a "test intersection" query, in fact the
			// construction here is a "find intersection" query.
			double t00 = DotProduct(D, S0[0] - A); double t01 = DotProduct(D, S0[1] - A);
			double t10 = DotProduct(D, S1[0] - A); double t11 = DotProduct(D, S1[1] - A);
			auto I0 = std::minmax(t00, t01);
			auto I1 = std::minmax(t10, t11);
			return (I0.second > I1.first && I0.first < I1.second);
		}

		return false;
	}


	bool overlap(float min1, float max1, float min2, float max2) {
		return !(max1 < min2 || max2 < min1);
	}




	bool isIntersecting( OOB& box1, OOB& box2)
	{
		// All potential separating axes
		std::array<Vector<3,double>,15> axes;

		// Add the face normals of both boxes
		for (int i = 0; i < 3; ++i)
		{
			axes[i] = box1.GetOOBAxis()[i];
			axes[i + 3] = box2.GetOOBAxis()[i];
		}

		// Add the cross products of the edges of the two boxes
		int index = 6;
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j)
			{
				axes[index++] = Cross(box1.GetOOBAxis()[i],box2.GetOOBAxis()[j]);
			}
		}

		// Test all axes for separation
		for (const auto& axis : axes)
		{
			// Skip zero-length axes (can happen with parallel edges)
			if (std::abs(axis[0]) < 1e-6 && std::abs(axis[1]) < 1e-6 && std::abs(axis[2]) < 1e-6) {
				continue;
			}

			// Project both boxes onto the axis
			float min1, max1, min2, max2;
		    box1.projectOntoAxis(box1, axis, min1, max1);
		    box2.projectOntoAxis(box2, axis, min2, max2);

			// Check for overlap
			if (!overlap(min1, max1, min2, max2)) {
				return false; // Found a separating axis, no intersection
			}

		}

		// No separating axis found, the boxes intersect
		return true;
	}


};

}

