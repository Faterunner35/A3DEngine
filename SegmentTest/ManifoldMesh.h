#pragma once
// ETManifoldMesh — порт gte::ETManifoldMesh из Geometric Tools Engine
// (Mathematics/ETManifoldMesh.h, версия 6.0.2023.01.17).
//   David Eberly, Geometric Tools, Redmond WA 98052
//   Copyright (c) 1998-2023
//   Distributed under the Boost Software License, Version 1.0.
//   https://www.boost.org/LICENSE_1_0.txt
//   https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// Отличия от оригинала:
//  - LogAssert/LogError (gte/Logger.h) заменены на Check/Fail ниже: они, как и
//    в GTE, бросают std::runtime_error;
//  - Insert проверяет неманифолдность ДО изменения структуры. В GTE отказ
//    посреди цикла по рёбрам оставлял в mEMap уже созданные рёбра со ссылкой
//    на удалённый треугольник (висячий указатель);
//  - operator= защищён от самоприсваивания;
//  - GetBoundaryPolygon проверяет повторный заход в треугольник (как в
//    более новых версиях GTE);
//  - тела методов — в ManifoldMesh.cpp.
#include"TriangleKey.h"
#include"EdgeKey.h"
#include<array>
#include<cstddef>
#include<cstdint>
#include<limits>
#include<map>
#include<memory>
#include<unordered_map>
#include<utility>
#include<vector>

namespace Geom_Exp
{

	class ETManifoldMesh
	{
	public:
		// Edge data types.
		class Edge;
		using ECreator = std::unique_ptr<Edge>(*)(int32_t, int32_t);
		using EMap = std::unordered_map<EdgeKey<false>, std::unique_ptr<Edge>,
			EdgeKey<false>, EdgeKey<false>>;

		// Triangle data types.
		class Triangle;
		using TCreator = std::unique_ptr<Triangle>(*)(int32_t, int32_t, int32_t);
		using TMap = std::unordered_map<TriangleKey<true>, std::unique_ptr<Triangle>,
			TriangleKey<true>, TriangleKey<true>>;

		// Edge object.
		class Edge
		{
		public:
			virtual ~Edge() = default;

			Edge(int32_t v0, int32_t v1)
				:
				V{ v0, v1 },
				T{ nullptr, nullptr }
			{
			}

			// Vertices of the edge.
			std::array<int32_t, 2> V;

			// Triangles sharing the edge.
			std::array<Triangle*, 2> T;
		};

		// Triangle object.
		class Triangle
		{
		public:
			virtual ~Triangle() = default;

			Triangle(int32_t v0, int32_t v1, int32_t v2)
				:
				V{ v0, v1, v2 },
				E{ nullptr, nullptr, nullptr },
				T{ nullptr, nullptr, nullptr }
			{
			}

			// The edge <u0,u1> is directed. Determine whether the triangle
			// has an edge <V[i],V[(i+1)%3]> = <u0,u1> (return +1) or an edge
			// <V[i],V[(i+1)%3]> = <u1,u0> (return -1) or does not have an
			// edge meeting either condition (return 0).
			int32_t WhichSideOfEdge(int32_t u0, int32_t u1) const
			{
				for (size_t i0 = 2, i1 = 0; i1 < 3; i0 = i1++)
				{
					if (V[i0] == u0 && V[i1] == u1)
					{
						return +1;
					}
					if (V[i0] == u1 && V[i1] == u0)
					{
						return -1;
					}
				}
				return 0;
			}

			Triangle* GetAdjacentOfEdge(int32_t u0, int32_t u1)
			{
				for (size_t i0 = 2, i1 = 0; i1 < 3; i0 = i1++)
				{
					if ((V[i0] == u0 && V[i1] == u1) || (V[i0] == u1 && V[i1] == u0))
					{
						return T[i0];
					}
				}
				return nullptr;
			}

			bool GetOppositeVertexOfEdge(int32_t u0, int32_t u1, int32_t& uOpposite)
			{
				for (size_t i0 = 2, i1 = 0; i1 < 3; i0 = i1++)
				{
					if ((V[i0] == u0 && V[i1] == u1) || (V[i0] == u1 && V[i1] == u0))
					{
						uOpposite = V[(i1 + 1) % 3];
						return true;
					}
				}
				return false;
			}

			// Vertices, listed in counterclockwise order (V[0],V[1],V[2]).
			std::array<int32_t, 3> V;

			// Adjacent edges. E[i] points to edge (V[i],V[(i+1)%3]).
			std::array<Edge*, 3> E;

			// Adjacent triangles. T[i] points to the adjacent triangle
			// sharing edge E[i].
			std::array<Triangle*, 3> T;
		};


		// Construction and destruction.
		virtual ~ETManifoldMesh() = default;

		ETManifoldMesh(ECreator eCreator = nullptr, TCreator tCreator = nullptr)
			:
			mECreator(eCreator ? eCreator : CreateEdge),
			mEMap{},
			mTCreator(tCreator ? tCreator : CreateTriangle),
			mTMap{},
			mThrowOnNonmanifoldInsertion(true)
		{
		}

		// Support for a deep copy of the mesh. The mEMap and mTMap objects
		// have dynamically allocated memory for edges and triangles. A
		// shallow copy of the pointers isn't possible with unique_ptr.
		ETManifoldMesh(ETManifoldMesh const& mesh)
			:
			ETManifoldMesh()
		{
			*this = mesh;
		}

		ETManifoldMesh& operator=(ETManifoldMesh const& mesh);

		// Member access.
		inline EMap const& GetEdges() const
		{
			return mEMap;
		}

		inline TMap const& GetTriangles() const
		{
			return mTMap;
		}

		// If the insertion of a triangle fails because the mesh would become
		// nonmanifold, the default behavior is to throw an exception. You
		// can disable this behavior and continue gracefully without an
		// exception. The return value is the previous value of the internal
		// state mAssertOnNonmanifoldInsertion.
		bool ThrowOnNonmanifoldInsertion(bool doException)
		{
			std::swap(doException, mThrowOnNonmanifoldInsertion);
			return doException;  // return the previous state
		}

		// If <v0,v1,v2> is not in the mesh, a Triangle object is created and
		// returned; otherwise, <v0,v1,v2> is in the mesh and nullptr is
		// returned. If the insertion leads to a nonmanifold mesh, the call
		// fails: std::runtime_error is thrown (default) or, after
		// ThrowOnNonmanifoldInsertion(false), nullptr is returned. In both
		// cases the mesh is left unchanged.
		virtual Triangle* Insert(int32_t v0, int32_t v1, int32_t v2);

		// If <v0,v1,v2> is in the mesh, it is removed and 'true' is
		// returned; otherwise, <v0,v1,v2> is not in the mesh and 'false' is
		// returned.
		virtual bool Remove(int32_t v0, int32_t v1, int32_t v2);

		// Destroy the edges and triangles to obtain an empty mesh.
		virtual void Clear()
		{
			mEMap.clear();
			mTMap.clear();
		}

		// A manifold mesh is closed if each edge is shared twice. A closed
		// mesh is not necessarily oriented. For example, you could have a
		// mesh with spherical topology. The upper hemisphere has outer
		// facing normals and the lower hemisphere has inner-facing normals.
		// The discontinuity in orientation occurs on the circle shared by the
		// hemispheres.
		bool IsClosed() const;

		// Test whether all triangles in the mesh are oriented consistently
		// and that no two triangles are coincident. The latter means that
		// you cannot have both triangles <v0,v1,v2> and <v0,v2,v1> in the
		// mesh to be considered oriented.
		bool IsOriented() const;

		// Compute the connected components of the edge-triangle graph that
		// the mesh represents. The first function returns pointers into
		// 'this' object's containers, so you must consume the components
		// before clearing or destroying 'this'. The second function returns
		// triangle keys, which requires three times as much storage as the
		// pointers but allows you to clear or destroy 'this' before consuming
		// the components.
		void GetComponents(std::vector<std::vector<Triangle*>>& components) const;
		void GetComponents(std::vector<std::vector<TriangleKey<true>>>& components) const;

		// Create a compact edge-triangle graph. The vertex indices are those
		// integers passed to an Insert(...) call. These have no meaning to
		// the semantics of maintaining an edge-triangle manifold mesh, so
		// this class makes no assumption about them. The vertex indices do
		// not necessarily start at 0 and they are not necessarily contiguous
		// numbers. The triangles are represented by triples of vertex
		// indices. The compact graph stores these in an array of N triples,
		// say,
		//   T[0] = (v0,v1,v2), T[1] = (v3,v4,v5), ...
		// Each triangle has up to 3 adjacent triangles. The compact graph
		// stores the adjacency information in an array of N triples, say,
		// and the indices of adjacent triangle are stored in an array of Nt
		//   A[0] = (t0,t1,t2), A[1] = (t3,t4,t5), ...
		// where the ti are indices into the array of triangles. For example,
		// the triangle T[0] has edges E[0] = (v0,v1), E[1] = (v1,v2) and
		// E[2] = (v2,v0). The edge E[0] has adjacent triangle T[0]. If E[0]
		// has another adjacent triangle, it is T[A[0][0]. If it does not have
		// another adjacent triangle, then A[0][0] = -1 (represented by
		// std::numeric_limits<size_t>::max()). Similar assignments are made
		// for the other two edges which produces A[0][1] for E[1] and
		// A[0][2] for E[2].
		void CreateCompactGraph(
			std::vector<std::array<size_t, 3>>& triangles,
			std::vector<std::array<size_t, 3>>& adjacents) const;

		// The output of CreateCompactGraph can be used to compute the
		// connected components of the graph, each component having
		// triangles with the same chirality (winding order). Using only
		// the mesh topology, it is not possible to ensure that the
		// chirality is the same for all the components. Additional
		// application-specific geometric information is required.
		//
		// The 'components' contain indices into the 'triangles' array
		// and is partitioned into C subarrays, each representing a
		// connected component. The lengths of the subarrays are
		// stored in 'numComponentTriangles[]'. The number of elements
		// of this array is C. It is the case that the number of triangles
		// in the mesh is sum_{i=0}^{C-1} numComponentTriangles[i].
		//
		// On return, the 'triangles' and 'adjacents' have been modified
		// and have the correct chirality.
		static void GetComponentsConsistentChirality(
			std::vector<std::array<size_t, 3>>& triangles,
			std::vector<std::array<size_t, 3>>& adjacents,
			std::vector<size_t>& components,
			std::vector<size_t>& numComponentTriangles);

		// This is a simple wrapper around CreateCompactGraph(...) and
		// GetComponentsConsistentChirality(...), in particular when you
		// do not need to work directly with the connected components.
		// The mesh is reconstructed, because the bookkeeping details of
		// trying to modify the mesh in-place are horrendous. NOTE: If
		// your mesh has more than 1 connected component, you should read
		// the comments for GetComponentsConsistentChirality(...) about
		// the potential for different chiralities between components.
		virtual void MakeConsistentChirality();

		// Compute the boundary-edge components of the mesh. These are
		// polygons that are simple for the strict definition of manifold
		// mesh that disallows bow-tie configurations. The GTE mesh
		// implementations do allow bow-tie configurations, in which case
		// some polygons might not be simple. If you select duplicateEndpoints
		// to be false, a component has consecutive vertices
		// (v[0], v[1], ..., v[n-1]) and the polygon has edges
		// (v[0],v[1]), (v[1],v[2]), ..., (v[n-2],v[n-1]), (v[n-1],v[0]).
		// If duplicateEndpoints is set to true, a component has consecutive
		// vertices (v[0], v[1], ..., v[n-1], v[0]), emphasizing that the
		// component is closed.
		void GetBoundaryPolygons(std::vector<std::vector<int32_t>>& polygons,
			bool duplicateEndpoints) const;

	protected:
		using TrianglePtrIntMap = std::unordered_map<Triangle*, int32_t>;
		using TrianglePtrSizeTMap = std::unordered_map<Triangle*, size_t>;

		struct BoundaryEdge
		{
			BoundaryEdge()
				:
				triangle(nullptr),
				index(std::numeric_limits<size_t>::max()),
				visited(false)
			{
			}

			BoundaryEdge(Triangle const* inTriangle, size_t inIndex, bool inVisited)
				:
				triangle(inTriangle),
				index(inIndex),
				visited(inVisited)
			{
			}

			Triangle const* triangle;
			size_t index;
			bool visited;
		};

		using BoundaryEdgeMap = std::map<std::array<int32_t, 2>, BoundaryEdge>;

		// Замена LogAssert/LogError из gte/Logger.h. Семантика та же, что в GTE:
		// нарушенное условие — исключение std::runtime_error с сообщением.
		static void Check(bool condition, char const* message)
		{
			if (!condition)
			{
				Fail(message);
			}
		}

		[[noreturn]] static void Fail(char const* message);

		// The edge data and default edge creation.
		static std::unique_ptr<Edge> CreateEdge(int32_t v0, int32_t v1)
		{
			return std::make_unique<Edge>(v0, v1);
		}

		// The triangle data and default triangle creation.
		static std::unique_ptr<Triangle> CreateTriangle(int32_t v0, int32_t v1, int32_t v2)
		{
			return std::make_unique<Triangle>(v0, v1, v2);
		}

		// Support for computing connected components. This is a
		// straightforward depth-first search of the graph but uses a
		// preallocated stack rather than a recursive function that could
		// possibly overflow the call stack.
		void DepthFirstSearch(
			Triangle* tInitial,
			TrianglePtrIntMap& visited,
			std::vector<Triangle*>& component) const;

		void GetBoundaryPolygon(
			Triangle const* initialTriangle,
			size_t initialIndex,
			BoundaryEdgeMap& boundaryEdges,
			std::vector<int32_t>& polygon) const;

		ECreator mECreator;
		EMap mEMap;
		TCreator mTCreator;
		TMap mTMap;
		bool mThrowOnNonmanifoldInsertion;
	};

}
