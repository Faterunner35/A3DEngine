#pragma once
#include"Geomerty.h"
#include<array>
#include<cstddef>
#include<ostream>
#include<vector>


namespace mathx
{


class Triangle
{

public:
	Triangle(const mathx::Vector<3, double> v0,const Vector<3, double> v1, const Vector<3, double> v2)
	{
		triangle_data[0]=  v0;
		triangle_data[1] = v1;
		triangle_data[2] = v2;
		UpdateNormal();
	}

	Triangle(Vector<3, double> vec_triangle[3])
	{
		triangle_data[0] = vec_triangle[0];
		triangle_data[1] = vec_triangle[1];
		triangle_data[2] = vec_triangle[2];
		UpdateNormal();
	}

	Triangle(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2)
	{
		triangle_data[0] = v0;
		triangle_data[1] = v1;
		triangle_data[2] = v2;
		UpdateNormal();
	}


	Triangle(const Triangle& poly)
	{
		triangle_data[0] = poly.GetVertex(0);
		triangle_data[1] = poly.GetVertex(1);
		triangle_data[2] = poly.GetVertex(2);
		normal = poly.GetNormal();
	}

	Triangle& operator=(const Triangle& poly)
	{
		triangle_data[0] = poly.GetVertex(0);
		triangle_data[1] = poly.GetVertex(1);
		triangle_data[2] = poly.GetVertex(2);
		normal = poly.GetNormal();
		return *this;
	}

	Triangle() = default;

	inline Vector<3, double> GetNormal() const
	{
		return normal;
	}
	inline std::array<Vector<3, double>,3> GetArray() const
	{

		return triangle_data;
	}

	inline Vector<3, double> GetVertex(int vert) const
	{
		return triangle_data[vert];
	}

	inline void SetVertex(int pos, const Vector<3, double>& point)
	{
		triangle_data[pos] = point;
		UpdateNormal();
	}

	// Ближайшая к point точка треугольника (включая границу).
	// Классический разбор по областям Вороного: Ericson, Real-Time Collision
	// Detection, 5.1.5 (ClosestPtPointTriangle). Нормаль не используется,
	// поэтому годится и для неединичной, и для нулевой нормали.
	Vector<3, double> ClosestPoint(const Vector<3, double>& point) const;

	// Расстояние (без знака) от точки до треугольника как множества точек:
	// |point − ClosestPoint(point)|. Не путать с расстоянием до ПЛОСКОСТИ
	// треугольника — оно совпадает только когда проекция попала внутрь.
	// Раньше функция считала расстояние до плоскости и ничего не возвращала.
	double PointToTriangleDisctance(const Vector<3, double>& point) const;


	std::array<Vector<3,double>,3>triangle_data;
	Vector<3, double> normal;

private:
	// Единичная нормаль (v1 − v0) × (v2 − v0); для вырожденного треугольника
	// нулевая. Раньше её вычислял только конструктор из трёх Vector, а
	// остальные оставляли нулевой.
	void UpdateNormal()
	{
		normal = Cross(triangle_data[1] - triangle_data[0], triangle_data[2] - triangle_data[0]);
		Normalize(normal);
	}

};

class Edge
{
public:
	Edge(const Vector<3, double>& p1, const Vector<3, double>& p2)
	{
	   edge_[0]=p1;
	   edge_[1]=p2;

   }

	// Лексикографический порядок по (edge_[0], edge_[1]) — только для
	// упорядоченных контейнеров, геометрического смысла у него нет (как
	// operator<=> у Vector). Раньше результат сравнения отбрасывался.
	inline bool operator >(const Edge& edg_2) const
	{
		return edge_ > edg_2.edge_;
	}

std::array<Vector<3,double>,2>edge_;

};


// Выпуклый плоский многоугольник — список вершин против часовой стрелки
// относительно getNormal(). Семантика методов — как у Ogre::Polygon, откуда
// класс взят; тела в Primitives.cpp.
//
// Удалены объявления без определений, смысл которых из кода не
// восстанавливается и которыми никто не пользовался: список треугольников
// (TriangleList/VertexTriangleList из кортежей <size_t, size_t, Vector>,
// конструктор от него, insertTriangle, GetTriangle, getVertexTriangleList,
// CalculateNormalsInTriangleArray) и HasIntersection с EIntersectionType.
class Polygon
{

    public:
        typedef std::vector<Vector<3,double>> VertexList;

        typedef std::pair<Vector<3,double>, Vector<3,double>> Edge;

    protected:
        VertexList      mVertexList;
        mutable  Vector<3,double> mNormal;
        mutable bool    mIsNormalSet;
        /** Updates the normal.
         */
        void updateNormal() const;

    public:
        Polygon();
        ~Polygon() = default;
        Polygon(const Polygon& cpy) = default;
        Polygon& operator=(const Polygon& cpy) = default;

        /** Inserts a vertex at a specific position.
        @note Vertices must be coplanar.
        */
        void insertVertex(const Vector<3,double>& vdata, size_t vertex);

        /** Inserts a vertex at the end of the polygon.
        @note Vertices must be coplanar.
        */
        void insertVertex(const Vector<3,double>& vdata);

        // Единичная нормаль треугольника (t1 − t0) × (t2 − t0); нулевая для
        // вырожденного.
        static Vector<3,double> CalculateNormalToTriangle(const std::array<Vector<3,double>,3>& triangle);

        /** Returns a vertex.
         */
        const Vector<3,double>& getVertex(size_t vertex) const;

        //sets vertex_polygon list size

        void setVertexListSize(size_t number_of_vertex);


        /** Sets a specific vertex of a polygon.
        @note Vertices must be coplanar.
        */
        void setVertex(const Vector<3,double>& vdata, size_t vertexIndex);

        /** Removes duplicate vertices from a polygon.
         */
        void removeDuplicates();

        /** Vertex count.
         */
        size_t getVertexCount() const;

        /** Returns the polygon normal.
         */
        const Vector<3,double>& getNormal() const;

        /** Deletes a specific vertex.
         */
        void deleteVertex(size_t vertex);

        /** Determines if a point is inside the polygon.
        @remarks
            A point is inside a polygon if it is both on the polygon's plane,
            and within the polygon's bounds. Polygons are assumed to be convex
            and planar.
        */
        bool isPointInside(const Vector<3,double>& point) const;

        /** Stores the edges of the polygon in ccw order.
            The vertices are copied so the user has to take the
            deletion into account.
        */
        //void storeEdges(EdgeMap* edgeMap) const;

        /** Resets the object.
         */
        void reset();

        /** Determines if the current object is equal to the compared one.
         */
        bool operator==(const Polygon& rhs) const;

        /** Determines if the current object is not equal to the compared one.
         */
        bool operator!=(const Polygon& rhs) const { return !(*this == rhs); }

        /** Prints out the polygon data.
         */
         friend std::ostream& operator<<(std::ostream& strm, const Polygon& poly);
    };


    class OOB
    {
      public:
       OOB()=default;
       ~OOB()=default;

    private:
   Vector<3,double>center_oob;
    std::array<Vector<3, double>, 3>oob_axis;
   Vector<3, double> half_extents;

    public:
		    // Раньше проецировался центр ЭТОГО ящика, а не box; полуразмеры брались
		    // из box. Совпадало только при вызове box.projectOntoAxis(box, ...).
		    void projectOntoAxis(const OOB& box, const Vector<3,double>& axis, double& min, double& max) const
            {
			// Project the center of the box ont the axis
			double centerProjection =DotProduct(box.center_oob,axis);

			// Compute the projection of the half-extents onto the axis
			double extent = std::abs(DotProduct(box.oob_axis[0],(axis)) * box.half_extents[0]) +
                std::abs(DotProduct(box.oob_axis[1], (axis)) * box.half_extents[1]) +
                std::abs(DotProduct(box.oob_axis[2], (axis)) * box.half_extents[2]);

			// The min and max projections
			min = centerProjection - extent;
			max = centerProjection + extent;
		    }


       inline std::array<Vector<3, double>,3> GetOOBAxis() const
        {
            return oob_axis;
        }

   };

 }
