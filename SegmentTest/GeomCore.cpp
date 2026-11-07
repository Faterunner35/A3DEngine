#include "GeomCore.h"

// Тела пересечений из GeomCore.h и нестатических методов GeomCoreIntersector,
// объявленных в GeomIntersector.h.

namespace mathx
{

std::optional<Ray> intersectionPlanePlane(const Plane& pl_1, const Plane& pl_2)
{
	// Направление прямой: перпендикулярно обеим нормалям
	const Vector<3, double> v = Cross(pl_1.GetNormal(), pl_2.GetNormal());
	const double vv = DotProduct(v, v);

	// Параллельны — прямой пересечения нет
	if (vv <= 1e-12) return std::nullopt;

	// ВНИМАНИЕ: порядок в cross РАЗНЫЙ.
	// Это следует из формулы трёх плоскостей, куда третьей подставлена [v|0]:
	//   p = [d1·(n3 × n2) + d2·(n1 × n3)] / [n1,n2,n3],  где n3 = v, d3 = 0
	const Vector<3, double> t1 = Cross(v, pl_2.GetNormal());   // v × n2
	const Vector<3, double> t2 = Cross(pl_1.GetNormal(), v);   // n1 × v
	const double inv = 1.0 / vv;

	// Скобки: деление применяется ко ВСЕМУ числителю
	const Vector<3, double> point{
		(pl_1.GetDistance() * t1[0] + pl_2.GetDistance() * t2[0]) * inv,
		(pl_1.GetDistance() * t1[1] + pl_2.GetDistance() * t2[1]) * inv,
		(pl_1.GetDistance() * t1[2] + pl_2.GetDistance() * t2[2]) * inv };

	return Ray{ point, v };
}


    double GeomCoreIntersector::RayTriangleIntersect(const Ray& ray, const Triangle& tri)
    {
		
	   auto vec_ray_dist = tri.triangle_data[0] - ray.getOrigin();
	   auto v01_edge = tri.GetArray()[1] - tri.GetArray()[0];
	   auto v02_edge = tri.GetArray()[2] - tri.GetArray()[0];

	   auto dist_ray_edge_02=Cross(ray.getDirection(), v01_edge);
	   auto cross_ray_edge01=Cross(vec_ray_dist,v02_edge);
       
	   auto scalar_ray_intersect=DotProduct(dist_ray_edge_02,v01_edge);

	   auto t1=DotProduct(dist_ray_edge_02,vec_ray_dist)/scalar_ray_intersect;
	   auto t2=DotProduct(ray.getDirection(),cross_ray_edge01)/scalar_ray_intersect;

	   auto t3=1-t1+t2;

	   return t3;

   }

	std::pair<bool, Vector<3,double>> GeomCoreIntersector::RaySphereIntersection(const Ray&ray,const Sphere& s)
	{
		const Vector<3,double>& raydir = ray.getDirection();
		// Adjust ray origin relative to sphere center
		const Vector<3,double>& rayorig = ray.getOrigin() - s.getCenter();
		double radius = s.getRadius();

		// Check origin inside first. Было: Length(rayorig) <= r² — длину
		// сравнивали с квадратом радиуса.
		if (DotProduct(rayorig, rayorig) <= radius * radius)
		{
			return std::pair<bool, Vector<3,double>>(true,Vector<3,double>::Zero());
		}

		// Mmm, quadratics
		// Build coeffs which can be used with std quadratic solver
		// ie t = (-b +/- sqrt(b*b + 4ac)) / 2a
		double a = DotProduct(raydir,raydir);
		double b = 2 * DotProduct(rayorig,raydir);
		double c = DotProduct(rayorig,rayorig) - radius * radius;

		// Calc determinant
		double d = (b * b) - (4 * a * c);
		if (d < 0)
		{
			// No intersection
			Vector<3, double>v{0, 0, 0};
			return std::pair<bool, Vector<3,double>>(false,v);
		}
		else
		{
			// BTW, if d=0 there is one intersection, if d > 0 there are 2
			// But we only want the closest one, so that's ok, just use the
			// '-' version of the solver
			double t = (-b - std::sqrt(d)) / (2 * a);
			if (t < 0)
				t = (-b + std::sqrt(d)) / (2 * a);
			return std::pair<bool, Vector<3,double>>(true, Vector<3,double>(ray.getOrigin()+t*ray.getDirection()));
		}

	}


	std::optional<RayTriangleHit> RayTriIntersectionMollerTrambo(const Ray& ray, const Triangle& tri)
	{
		const Vector<3, double>& d = ray.getDirection();
		const Vector<3, double> e1 = tri.GetVertex(1) - tri.GetVertex(0);
		const Vector<3, double> e2 = tri.GetVertex(2) - tri.GetVertex(0);

		const Vector<3, double> p = Cross(d, e2);
		const double det = DotProduct(e1, p);

		// Параллелен (или треугольник вырожден): относительный порог, см. GeomCore.h
		const double kParallel = 1e-12;
		if (det * det <= kParallel * kParallel * DotProduct(d, d) * DotProduct(e1, e1) * DotProduct(e2, e2))
		{
			return std::nullopt;
		}
		const double invDet = 1.0 / det;

		const Vector<3, double> s = ray.getOrigin() - tri.GetVertex(0);
		const double u = DotProduct(s, p) * invDet;
		if (u < 0.0 || u > 1.0)
		{
			return std::nullopt;
		}

		const Vector<3, double> q = Cross(s, e1);
		const double v = DotProduct(d, q) * invDet;
		if (v < 0.0 || u + v > 1.0)
		{
			return std::nullopt;
		}

		const double t = DotProduct(e2, q) * invDet;
		if (t < 0.0)
		{
			return std::nullopt;                      // позади начала луча
		}

		RayTriangleHit hit;
		hit.t = t;
		hit.u = u;
		hit.v = v;
		hit.point = ray.getPoint(t);
		return hit;
	}



double RayTriangleintersectionFast(const Triangle& tri, const Ray& ray, float weight[3]) {

	Vector<3, double> e1 = tri.GetArray()[1] - tri.GetArray()[0];
	Vector<3, double> e2 = tri.GetArray()[2] - tri.GetArray()[0];
	Vector<3, double> origin_ray_tri = ray.getOrigin() - tri.GetArray()[0];

	auto cross_1 = Cross(ray.getDirection(), e2);

	auto inv_det = 1.0 / DotProduct(e1, cross_1);

	// scalar value of u coordinate (using cross_1!)
	weight[1] = DotProduct(origin_ray_tri, cross_1) * inv_det;

	if (weight[1] < 0.0f || weight[1] > 1.0f)
		return INFINITY;

	auto cross_2 = Cross(origin_ray_tri, e1);

	// scalar value of v coordinate
	weight[2] = DotProduct(ray.getDirection(), cross_2) * inv_det;

	// Checking u + v <= 1.0f
	if (weight[2] < 0.0f || (weight[1] + weight[2]) > 1.0f)
		return INFINITY;

	double t = DotProduct(e2, cross_2) * inv_det;

	// weight array is already float, so just calculate the final barycentric weight
	weight[0] = 1.0f - weight[1] - weight[2];

	return t;
}



std::optional<double> RayRayIntersection(const Ray& ray_1, const Ray& ray_2) {
	auto epsilon = 1e-6;
	Vector<3, double> edge_1 = ray_1.getOrigin() - ray_2.getOrigin(); // w0 vector

	// Gramian Matrix Components
	double d00 = DotProduct(ray_1.getDirection(), ray_1.getDirection());
	double d01 = DotProduct(ray_1.getDirection(), ray_2.getDirection());
	double d11 = DotProduct(ray_2.getDirection(), ray_2.getDirection());


	const Vector<3, double> n = Cross(ray_1.getDirection(), ray_2.getDirection());
	const double denom = DotProduct(n, n);

	

	// A machine epsilon (1e-16) is too strict for geometric stability here.
	// If the denominator is too close to zero, the rays are parallel.
	if (std::abs(denom) < 1e-8) {
		return std::nullopt;
	}

	double w0 = DotProduct(edge_1, ray_1.getDirection());
	double w1 = DotProduct(edge_1, ray_2.getDirection());

	// Fixed parentheses so division applies to the entire numerator!
	double t = (w1 * d01 - w0 * d11) / denom;
	double u = (d00 * w1 - d01 * w0) / denom;

	// Ray boundary check: Ensure intersection isn't behind either origin
	if (t < -epsilon || u < -epsilon) {
		return std::nullopt;
	}

	// Calculate closest points to verify intersection distance
	Vector<3, double> vec_1 = ray_1.getOrigin() + t * ray_1.getDirection();
	Vector<3, double> vec_2 = ray_2.getOrigin() + u * ray_2.getDirection();

	Vector<3, double> dist = vec_1 - vec_2;
	double dist_squared = DotProduct(dist, dist);

	if (dist_squared < epsilon) {
		return t;
	}

	return std::nullopt;
}


//void GeomCoreIntersector::convex_hull_creation(std::vector<Vector<3,double>> points_vector, Triangle triangle, std::vector<Triangle>&_array_triangles)
		//{
		//   int maxdistance=0;
		//   int farthestPointIndex = -1;
		//   for (int i = 0; i < points_vector.size(); i++)
		//   {    
		//       
		//       auto dist_temp=triangle.PointToTriangleDisctance(points_vector[i]);

		//	   if(dist_temp>maxdistance)
		//	   {
  //  	     	   maxdistance=dist_temp;
	 //    		   farthestPointIndex=i;

		//       }
		//	   

		//	   if (farthestPointIndex == -1) {
		//		   hull.push_back(face);
		//		   return;
		//	   }
		//	   
		//	  
		//	   Vector<3,double> farthestPoint = points_vector[farthestPointIndex];

		//	   points_vector.erase(points_vector.begin() + farthestPointIndex);

		//	   Triangle tri1(triangle.GetVertex(0), triangle.GetVertex(1), farthestPoint);
		//	   Triangle tri2(triangle.GetVertex(1), triangle.GetVertex(2), farthestPoint);
		//	   Triangle tri3(triangle.GetVertex(2), triangle.GetVertex(0), farthestPoint);

		//	   // Create three new faces using the farthest point
		//	 GeomCoreIntersector::convex_hull_creation(points_vector, tri1,_array_triangles);
		//	 GeomCoreIntersector::convex_hull_creation(points_vector, tri2,_array_triangles);
		//	 GeomCoreIntersector::convex_hull_creation(points_vector, tri3, _array_triangles);


		//}


		//}

		//  std:: vector<Geom_Exp::Triangle> computeConvexHull3D(vector<Geom_Exp::Vector<3,double>>&points) {
		//	   vector<Geom_Exp::Triangle> hull;

		//	   if (points.size() < 4) {
		//		   cout << "Convex hull not possible with less than 4 points." << endl;
		//		   return hull;
		//	   }

		//	   // Find the extreme points along each axis
		//	   int minX = 0, maxX = 0, minY = 0, maxY = 0, minZ = 0, maxZ = 0;
		//	   for (int i = 1; i < points.size(); i++) {
		//		   if (points[i].x < points[minX].x) minX = i;
		//		   if (points[i].x > points[maxX].x) maxX = i;
		//		   if (points[i].y < points[minY].y) minY = i;
		//		   if (points[i].y > points[maxY].y) maxY = i;
		//		   if (points[i].z < points[minZ].z) minZ = i;
		//		   if (points[i].z > points[maxZ].z) maxZ = i;
		//	   }

		//	   // Create an initial tetrahedron using the extreme points
		//	   Vector<3,double> p1 = points[minX];
		//	   Vector<3,double> p2 = points[maxX];
		//	   Vector<3,double> p3 = points[minY];
		//	   Vector<3,double> p4 = points[maxZ];

		//	   // Recursively find the convex hull
		//	   convex_hull_creation(points, { p1, p2, p3 }, hull);
		//	   convex_hull_creation(points, { p1, p2, p4 }, hull);
		//	   convex_hull_creation(points, { p1, p3, p4 }, hull);
		//	   convex_hull_creation(points, { p2, p3, p4 }, hull);

		//	   // Remove duplicate faces
		//	   sort(hull.begin(), hull.end(), [](const Triangle& f1, const Triangle& f2) {
		//		   return std::tie(f1.GetVertex(0)[0], f1.GetVertex(0)[1], f1.GetVertex(0)[2], f1.b.x, f1.b.y, f1.b.z, f1.c.x, f1.c.y, f1.c.z) <
		//			   std::tie(f2.a.x, f2.a.y, f2.a.z, f2.b.x, f2.b.y, f2.b.z, f2.c.x, f2.c.y, f2.c.z);
		//		   });
		//	   hull.erase(unique(hull.begin(), hull.end()), hull.end());

		//	   return hull;
		//   }
}


     
	


    
