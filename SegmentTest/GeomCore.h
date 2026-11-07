#pragma once
     #include"GeomIntersector.h"
	 #include"Primitives.h"
	 #include"Ray.h"
	 #include<optional>
	 #include<limits>
	 #include<cassert>
	 
namespace mathx
{ 
//  Проекция, отбрасывание и отражение.
//
//  Все три строятся на одной операции — проекции на направление:
//      proj_a(v) = (v·a)·a          <- результат лежит ВДОЛЬ a
//      rej_a(v)  = v − proj_a(v)    <- то, что осталось, ⊥ a
//      refl(v,n) = v − 2·proj_n(v)  <- параллельная часть перевёрнута
//
//  ЧАСТАЯ ОШИБКА: умножить скаляр (v·a) на v вместо a. Результат тогда
//  остаётся вдоль v, то есть направление не меняется вовсе. Ошибка НЕ
//  видна, когда v уже коллинеарен a — там обе формы совпадают.
//
//  ПРЕДУСЛОВИЕ: направление единичное. Иначе (v·a) содержит лишний
//  множитель |a|, и его надо делить на |a|² — см. версии *Unnormalized.
// ============================================================================

// Проекция v на единичное направление a. Результат КОЛЛИНЕАРЕН a.
inline Vector<3, double> VectorProjection(const Vector<3, double>& v,
	const Vector<3, double>& axis)
{
	assert(std::abs(DotProduct(axis, axis) - 1.0) < 1e-9 &&
		"VectorProjection: axis должен быть единичным");
	const double s = DotProduct(v, axis);
	return Vector<3, double>{ axis[0] * s, axis[1] * s, axis[2] * s };
}

// Отбрасывание: то, что осталось от v после снятия проекции. ⊥ axis.
inline Vector<3, double> VectorRejection(const Vector<3, double>& v,
	const Vector<3, double>& axis)
{
	assert(std::abs(DotProduct(axis, axis) - 1.0) < 1e-9 &&
		"VectorRejection: axis должен быть единичным");
	const double s = DotProduct(v, axis);
	return Vector<3, double>{ v[0] - axis[0] * s, v[1] - axis[1] * s, v[2] - axis[2] * s };
}

// Отражение v через плоскость с единичной нормалью n.
// Параллельная нормали часть переворачивается, перпендикулярная остаётся.
// СВОЙСТВА: длина сохраняется, применить дважды = вернуться.
inline Vector<3, double> VectorReflection(const Vector<3, double>& v,
	const Vector<3, double>& normal)
{
	assert(std::abs(DotProduct(normal, normal) - 1.0) < 1e-9 &&
		"VectorReflection: normal должна быть единичной");
	const double s = 2.0 * DotProduct(v, normal);
	return Vector<3, double>{ v[0] - normal[0] * s, v[1] - normal[1] * s, v[2] - normal[2] * s };
}

// Отражение ЛУЧА от поверхности (нормаль направлена наружу).
// Отличается знаком: падающий луч идёт К поверхности, отражённый ОТ неё.
// Формула та же, разница только в том, что здесь v — направление движения.
inline Vector<3, double> ReflectRay(const Vector<3, double>& incident,
	const Vector<3, double>& normal)
{
	return VectorReflection(incident, normal);
}

// --- Версии для НЕединичных направлений: делят на |axis|² ---

inline Vector<3, double> VectorProjectionUnnormalized(const Vector<3, double>& v,
	const Vector<3, double>& axis)
{
	const double aa = DotProduct(axis, axis);
	assert(aa > 1e-24 && "VectorProjectionUnnormalized: нулевая ось");
	const double s = DotProduct(v, axis) / aa;
	return Vector<3, double>{ axis[0] * s, axis[1] * s, axis[2] * s };
}

inline Vector<3, double> VectorRejectionUnnormalized(const Vector<3, double>& v,
	const Vector<3, double>& axis)
{
	const double aa = DotProduct(axis, axis);
	assert(aa > 1e-24 && "VectorRejectionUnnormalized: нулевая ось");
	const double s = DotProduct(v, axis) / aa;
	return Vector<3, double>{ v[0] - axis[0] * s, v[1] - axis[1] * s, v[2] - axis[2] * s };
}

inline double DistancePointToLine(const Vector<3, double>& q, const Ray& line)
{
	const Vector<3, double>& d = line.getDirection();
	const double dd = DotProduct(d, d);

	const Vector<3, double> w{ q[0] - line.getOrigin()[0],
							  q[1] - line.getOrigin()[1],
							  q[2] - line.getOrigin()[2] };

	if (dd < 1e-24) return std::sqrt(DotProduct(w, w));   // вырождена в точку

	const Vector<3, double> area = Cross(w, d);
	return std::sqrt(DotProduct(area, area) / dd);
}


//расстояние от точки до прямой — ЧЕРЕЗ ПЛОЩАДЬ.
//
//  |(q − p) × v| — площадь параллелограмма на этих двух векторах.
//  Площадь = основание × высота, основание это |v|, высота — искомое
//  расстояние. Отсюда:
//
//      d = |(q − p) × v| / |v|
//
//  ПОЧЕМУ НЕ ЧЕРЕЗ ПИФАГОРА. Формула d² = |w|² − (w·v)²/|v|² алгебраически
//  верна, но вычитает два БОЛЬШИХ близких числа. Если точка далеко от p,
//  но близко к самой прямой, оба слагаемых порядка |w|², а разность
//  порядка d². Замерено: при |w| = 1e8 и истинном d = 1 пифагорова форма
//  даёт РОВНО НОЛЬ — оба слагаемых равны 1e16, а в double всего 16
//  значащих цифр.
//
//  Cross тоже вычитает большие числа, но ВНУТРИ компонент, до возведения
//  в квадрат — то есть на величинах порядка координат, а не их квадратов.
//  Скорость та же (замерено 3.80 против 3.82 нс).
// ---------------------------------------------------------------------------


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
	   auto v12_edge = tri.GetArray()[2] - tri.GetArray()[1];
	   auto dist_v0_edge=ray.getDirection() - tri.GetArray()[0];

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

		// Check origin inside first
		if (mathx::Length(rayorig) <= radius * radius)
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


	std::optional<Vector<3, double>> RayTriIntersectionMollerTrambo(const Ray& ray, const Triangle& tri)
	{

		auto epsilon_error = std::numeric_limits<double>::epsilon();

		auto edge_1 = tri.GetVertex(1) - tri.GetVertex(0);
		auto edge_2 = tri.GetVertex(2) - tri.GetVertex(0);

		auto cross_raydest_edge2 = Cross(ray.getDirection(), edge_2);

		auto scalar_1 = DotProduct(cross_raydest_edge2, edge_1);

		if (std::abs(scalar_1) < epsilon_error)
		{
			return {};
		}

		auto invDet = 1 / scalar_1;

		auto rayorigin_triV0 = ray.getOrigin() - tri.GetVertex(0);

		auto cross_rayorigin_edge1 = Cross(rayorigin_triV0, edge_1);

		auto scalar_2 = invDet * DotProduct(cross_rayorigin_edge1, edge_1);

		if (std::abs(scalar_2) < -epsilon_error || scalar_2 - 1 > epsilon_error)
		{
			return {};
		}


		auto scalar_3 = invDet * DotProduct(cross_rayorigin_edge1, edge_1);

		if (scalar_3<-epsilon_error || scalar_3 + scalar_2 - 1>epsilon_error)
		{
			return {};
		}

		auto scalar_4 = invDet * DotProduct(cross_rayorigin_edge1, edge_2);

		if (std::abs(scalar_4) > epsilon_error)
		{
			return Vector<3, double>{ray.getOrigin() + scalar_3 * ray.getDirection()};
		}

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
	double d10 = d01; // Symmetric matrix, d10 = d01
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


inline std::optional<double> RayRayIntersectionRobust(const Ray& ray1, const Ray& ray2,
	double tolerance = 1e-6)
{
	const Vector<3, double>& d1 = ray1.getDirection();
	const Vector<3, double>& d2 = ray2.getDirection();

	const double d00 = DotProduct(d1, d1);
	const double d01 = DotProduct(d1, d2);
	const double d11 = DotProduct(d2, d2);

	const double denom = d00 * d11 - d01 * d01;

	// Порог ОТНОСИТЕЛЬНЫЙ: denom растёт как |d1|²·|d2|², и абсолютный
	// отверг бы перпендикулярные лучи с короткими направляющими.
	if (denom <= 1e-12 * d00 * d11) return std::nullopt;

	const Vector<3, double> w{ ray1.getOrigin()[0] - ray2.getOrigin()[0],
							  ray1.getOrigin()[1] - ray2.getOrigin()[1],
							  ray1.getOrigin()[2] - ray2.getOrigin()[2] };
	const double w0 = DotProduct(w, d1);
	const double w1 = DotProduct(w, d2);

	const double t = (d01 * w1 - d11 * w0) / denom;
	const double u = (d00 * w1 - d01 * w0) / denom;

	if (t < 0.0 || u < 0.0) return std::nullopt;      // позади начала луча

	const Vector<3, double> p1{ ray1.getOrigin()[0] + t * d1[0],
							   ray1.getOrigin()[1] + t * d1[1],
							   ray1.getOrigin()[2] + t * d1[2] };
	const Vector<3, double> p2{ ray2.getOrigin()[0] + u * d2[0],
							   ray2.getOrigin()[1] + u * d2[1],
							   ray2.getOrigin()[2] + u * d2[2] };
	const Vector<3, double> gap{ p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2] };

	if (DotProduct(gap, gap) > tolerance * tolerance) return std::nullopt;
	return t;
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


     
	


    












