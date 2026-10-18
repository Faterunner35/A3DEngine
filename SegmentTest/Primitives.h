#pragma once
#include"Geomerty.h"
#include<array>
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
       

		 Vector<3, double> edge1 = v1 - v0;
		 Vector<3, double> edge2 = v2 - v0;

		normal = Cross(edge1, edge2);
		Normalize(normal);
	}

	Triangle(Vector<3, double> vec_triangle[3])
	{
		triangle_data[0] = vec_triangle[0];
		triangle_data[1] = vec_triangle[1];
		triangle_data[2] = vec_triangle[2];
	}

	Triangle(const std::array<double, 3>& v0, const std::array<double, 3>& v1, const std::array<double, 3>& v2)
	{
		triangle_data[0] = v0;
		triangle_data[1] = v1;
		triangle_data[2] = v2;
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

	inline void SetVertex(int pos, Vector<3, double>& point)
	{
		triangle_data[pos] = point;
	}

    inline double PointToTriangleDisctance(const Vector<3, double>& point)
    { 
      Vector<3,double> Normal=GetNormal();
      Normalize(Normal);
      Vector<3,double>point_vert_dist=point-GetVertex(0);
      auto distance=DotProduct(point_vert_dist,Normal);

    }

    
	std::array<Vector<3,double>,3>triangle_data;
	Vector<3, double> normal;

};

class Edge
{
	Edge(const Vector<3, double>& p1, Vector<3, double>& p2)
	{
	   edge_[0]=p1;
	   edge_[1]=p2;
	   
   }

	inline bool operator >(const Edge& edg_2)
	{
        edge_[0]> edg_2.edge_[0];
		

   }

std::array<Vector<3,double>,2>edge_;

};


class Polygon
{
       
    public:
        typedef std::vector<Vector<3,double>> VertexList;
        typedef std::array<std::tuple<size_t, size_t,Vector<3,double>>, 3> TriangleList;
        typedef std::vector<std::array<std::tuple<size_t, size_t, Vector<3,double>>, 3>> VertexTriangleList;
        
        typedef std::pair<Vector<3,double>, Vector<3,double>> Edge;
        
         typedef std::vector<Vector<3,double>> array_of_normals;
        
        enum EIntersectionType
        {
            eIT_None,
            eIT_Intersection,
            eIT_JustTouch,
        };

    protected:
        VertexList      mVertexList;
        TriangleList        triangle_polygon;
        mutable  VertexTriangleList mTriangleList;
        mutable  Vector<3,double> mNormal;
        mutable bool    mIsNormalSet;
        /** Updates the normal.
         */
        void updateNormal() const;

    public:
        Polygon();
        Polygon(std::array<std::tuple<size_t, size_t, Vector<3,double>>, 3>& triangle);
        ~Polygon();
        Polygon(const Polygon& cpy);

        /** Inserts a vertex at a specific position.
        @note Vertices must be coplanar.
        */
        void insertVertex(const Vector<3,double>& vdata, size_t vertex);

        static EIntersectionType HasIntersection(std::unique_ptr<Polygon>pPolygon0, std::unique_ptr<Polygon>pPolygon1);
        /** Inserts a vertex at the end of the polygon.
        @note Vertices must be coplanar.
        */
        void insertVertex(const Vector<3,double>& vdata);

        std::vector<Vector<3,double>>& CalculateNormalsInTriangleArray();

        inline  Vector<3,double> CalculateNormalToTriangle(std::array<Vector<3,double>,3>& triangle);

        inline TriangleList GetTriangle() const {

            return this->triangle_polygon;

        }

        void insertTriangle(const std::array<std::tuple<size_t, size_t, Vector<3,double>>, 3>& triangle);

        //returns triangles list

        inline  VertexTriangleList* getVertexTriangleList() const
        {
            return &mTriangleList;

        }

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
		    void projectOntoAxis(OOB& box, const Vector<3,double>& axis, float& min, float& max)
            {
			// Project the center of the box ont the axis
			float centerProjection =DotProduct(center_oob,axis);
           
			// Compute the projection of the half-extents onto the axis
			float extent = std::abs(DotProduct(box.oob_axis[0],(axis)) * box.half_extents[0]) +
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


