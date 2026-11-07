#include "ManifoldMesh.h"
// Тела методов ETManifoldMesh — порт gte::ETManifoldMesh (Geometric Tools
// Engine, David Eberly, Boost Software License 1.0; см. шапку ManifoldMesh.h).
// Отличия от оригинала перечислены там же.

#include <queue>
#include <set>
#include <stdexcept>
#include <string>

namespace Geom_Exp
{

	void ETManifoldMesh::Fail(char const* message)
	{
		throw std::runtime_error(std::string("ETManifoldMesh: ") + message);
	}

	ETManifoldMesh& ETManifoldMesh::operator=(ETManifoldMesh const& mesh)
	{
		// В GTE самоприсваивание очищало обе сетки (это один объект) и
		// возвращало пустую.
		if (this == &mesh)
		{
			return *this;
		}

		Clear();

		mECreator = mesh.mECreator;
		mTCreator = mesh.mTCreator;
		mThrowOnNonmanifoldInsertion = mesh.mThrowOnNonmanifoldInsertion;
		for (auto const& element : mesh.mTMap)
		{
			// The typecast avoids warnings about not storing the return
			// value in a named variable. The return value is discarded.
			(void)Insert(element.first.V[0], element.first.V[1], element.first.V[2]);
		}

		return *this;
	}

	ETManifoldMesh::Triangle* ETManifoldMesh::Insert(int32_t v0, int32_t v1, int32_t v2)
	{
		TriangleKey<true> tkey(v0, v1, v2);
		if (mTMap.find(tkey) != mTMap.end())
		{
			// The triangle already exists. Return a null pointer as a
			// signal to the caller that the insertion failed.
			return nullptr;
		}

		// Проверка на неманифолдность ДО каких-либо изменений. В GTE она
		// стояла внутри цикла ниже: если первое ребро уже было создано, а на
		// втором выяснялось, что у него два треугольника, функция выходила,
		// оставляя в mEMap ребро с T[0] = удалённый треугольник. Здесь сетка
		// при отказе остаётся ровно такой, какой была.
		int32_t const vertices[3] = { v0, v1, v2 };
		for (int32_t i0 = 2, i1 = 0; i1 < 3; i0 = i1++)
		{
			auto eiter = mEMap.find(EdgeKey<false>(vertices[i0], vertices[i1]));
			if (eiter != mEMap.end() && eiter->second->T[1] != nullptr)
			{
				if (mThrowOnNonmanifoldInsertion)
				{
					Fail("Attempt to create nonmanifold mesh.");
				}
				return nullptr;
			}
		}

		// Create the new triangle. It will be added to mTMap at the end
		// of the function so that if an assertion is triggered and the
		// function returns early, the (bad) triangle will not be part of
		// the mesh.
		std::unique_ptr<Triangle> newTri = mTCreator(v0, v1, v2);
		Triangle* tri = newTri.get();

		// Add the edges to the mesh if they do not already exist.
		for (int32_t i0 = 2, i1 = 0; i1 < 3; i0 = i1++)
		{
			EdgeKey<false> ekey(tri->V[i0], tri->V[i1]);
			Edge* edge = nullptr;
			auto eiter = mEMap.find(ekey);
			if (eiter == mEMap.end())
			{
				// This is the first time the edge is encountered.
				std::unique_ptr<Edge> newEdge = mECreator(tri->V[i0], tri->V[i1]);
				edge = newEdge.get();
				mEMap[ekey] = std::move(newEdge);

				// Update the edge and triangle.
				edge->T[0] = tri;
				tri->E[i0] = edge;
			}
			else
			{
				// This is the second time the edge is encountered.
				edge = eiter->second.get();
				Check(edge != nullptr, "Unexpected condition.");

				// Update the edge. Третьего треугольника на ребре быть не
				// может: это отсечено проверкой выше.
				Check(edge->T[1] == nullptr, "Unexpected condition.");
				edge->T[1] = tri;

				// Update the adjacent triangles.
				auto adjacent = edge->T[0];
				Check(adjacent != nullptr, "Unexpected condition.");
				for (int32_t j = 0; j < 3; ++j)
				{
					if (adjacent->E[j] == edge)
					{
						adjacent->T[j] = tri;
						break;
					}
				}

				// Update the triangle.
				tri->E[i0] = edge;
				tri->T[i0] = adjacent;
			}
		}

		mTMap[tkey] = std::move(newTri);
		return tri;
	}

	bool ETManifoldMesh::Remove(int32_t v0, int32_t v1, int32_t v2)
	{
		TriangleKey<true> tkey(v0, v1, v2);
		auto titer = mTMap.find(tkey);
		if (titer == mTMap.end())
		{
			// The triangle does not exist.
			return false;
		}

		// Get the triangle.
		Triangle* tri = titer->second.get();

		// Remove the edges and update adjacent triangles if necessary.
		for (int32_t i = 0; i < 3; ++i)
		{
			// Inform the edges the triangle is being deleted.
			auto edge = tri->E[i];
			Check(edge != nullptr, "Unexpected condition.");

			if (edge->T[0] == tri)
			{
				// One-triangle edges always have pointer at index zero.
				edge->T[0] = edge->T[1];
				edge->T[1] = nullptr;
			}
			else if (edge->T[1] == tri)
			{
				edge->T[1] = nullptr;
			}
			else
			{
				Fail("Unexpected condition.");
			}

			// Remove the edge if you have the last reference to it.
			if (!edge->T[0] && !edge->T[1])
			{
				EdgeKey<false> ekey(edge->V[0], edge->V[1]);
				mEMap.erase(ekey);
			}

			// Inform adjacent triangles the triangle is being deleted.
			auto adjacent = tri->T[i];
			if (adjacent)
			{
				for (int32_t j = 0; j < 3; ++j)
				{
					if (adjacent->T[j] == tri)
					{
						adjacent->T[j] = nullptr;
						break;
					}
				}
			}
		}

		mTMap.erase(tkey);
		return true;
	}

	bool ETManifoldMesh::IsClosed() const
	{
		for (auto const& element : mEMap)
		{
			Edge* edge = element.second.get();
			if (!edge->T[0] || !edge->T[1])
			{
				return false;
			}
		}
		return true;
	}

	bool ETManifoldMesh::IsOriented() const
	{
		for (auto const& element : mEMap)
		{
			Edge* edge = element.second.get();
			if (edge->T[0] && edge->T[1])
			{
				// In each triangle, find the ordered edge that
				// corresponds to the unordered edge element.first. Also
				// find the vertex opposite that edge.
				bool edgePositive[2] = { false, false };
				int32_t vOpposite[2] = { -1, -1 };
				for (int32_t j = 0; j < 2; ++j)
				{
					auto tri = edge->T[j];
					for (int32_t i = 0; i < 3; ++i)
					{
						size_t szI = static_cast<size_t>(i);
						if (tri->V[i] == element.first.V[0])
						{
							int32_t vNext = tri->V[(szI + 1) % 3];
							if (vNext == element.first.V[1])
							{
								edgePositive[j] = true;
								vOpposite[j] = tri->V[(szI + 2) % 3];
							}
							else
							{
								edgePositive[j] = false;
								vOpposite[j] = vNext;
							}
							break;
						}
					}
				}

				// To be oriented consistently, the edges must have
				// reversed ordering and the oppositive vertices cannot
				// match.
				if (edgePositive[0] == edgePositive[1] || vOpposite[0] == vOpposite[1])
				{
					return false;
				}
			}
		}
		return true;
	}

	void ETManifoldMesh::GetComponents(std::vector<std::vector<Triangle*>>& components) const
	{
		components.clear();

		// visited: 0 (unvisited), 1 (discovered), 2 (finished)
		TrianglePtrIntMap visited{};
		for (auto const& element : mTMap)
		{
			visited.insert(std::make_pair(element.second.get(), 0));
		}

		for (auto& element : mTMap)
		{
			Triangle* tri = element.second.get();
			if (visited[tri] == 0)
			{
				std::vector<Triangle*> component{};
				DepthFirstSearch(tri, visited, component);
				components.push_back(std::move(component));
			}
		}
	}

	void ETManifoldMesh::GetComponents(std::vector<std::vector<TriangleKey<true>>>& components) const
	{
		components.clear();

		// visited: 0 (unvisited), 1 (discovered), 2 (finished)
		TrianglePtrIntMap visited{};
		for (auto const& element : mTMap)
		{
			visited.insert(std::make_pair(element.second.get(), 0));
		}

		for (auto& element : mTMap)
		{
			Triangle* tri = element.second.get();
			if (visited[tri] == 0)
			{
				std::vector<Triangle*> component{};
				DepthFirstSearch(tri, visited, component);

				std::vector<TriangleKey<true>> keyComponent{};
				keyComponent.reserve(component.size());
				for (auto const& t : component)
				{
					keyComponent.push_back(TriangleKey<true>(t->V[0], t->V[1], t->V[2]));
				}
				components.push_back(std::move(keyComponent));
			}
		}
	}

	void ETManifoldMesh::CreateCompactGraph(
		std::vector<std::array<size_t, 3>>& triangles,
		std::vector<std::array<size_t, 3>>& adjacents) const
	{
		size_t const numTriangles = mTMap.size();
		Check(numTriangles > 0, "Invalid input.");

		triangles.resize(numTriangles);
		adjacents.resize(numTriangles);

		TrianglePtrSizeTMap triIndexMap;
		triIndexMap.insert(std::make_pair(nullptr, std::numeric_limits<size_t>::max()));
		size_t index = 0;
		for (auto const& element : mTMap)
		{
			triIndexMap.insert(std::make_pair(element.second.get(), index++));
		}

		index = 0;
		for (auto const& element : mTMap)
		{
			auto const& triPtr = element.second;
			auto& tri = triangles[index];
			auto& adj = adjacents[index];
			for (size_t j = 0; j < 3; ++j)
			{
				tri[j] = triPtr->V[j];
				adj[j] = triIndexMap[triPtr->T[j]];
			}
			++index;
		}
	}

	void ETManifoldMesh::GetComponentsConsistentChirality(
		std::vector<std::array<size_t, 3>>& triangles,
		std::vector<std::array<size_t, 3>>& adjacents,
		std::vector<size_t>& components,
		std::vector<size_t>& numComponentTriangles)
	{
		Check(
			triangles.size() > 0 && triangles.size() == adjacents.size(),
			"Invalid inputs.");

		// Use a breadth-first search to process the chirality of the
		// triangles. Keep track of the connected components.
		size_t const numTriangles = triangles.size();
		std::vector<bool> visited(numTriangles, false);
		components.reserve(numTriangles);
		components.clear();

		// The 'firstUnvisited' index is the that of the first triangle to
		// process in a connected component of the mesh.
		for (;;)
		{
			// Let n[i] be the number of elements of the i-th connected
			// component. Let C be the number of components. During the
			// execution of this loop, the array numComponentTriangles
			// stores
			//   {0, n[0], n[0]+n[1], ..., n[0]+...+n[C-1]=numTriangles}
			// At the end of this function, the array is modified to
			//   {n[0], n[1], ..., n[C-1]}
			numComponentTriangles.push_back(components.size());

			// Find the starting index of a connected component.
			size_t firstUnvisited = numTriangles;
			for (size_t i = 0; i < numTriangles; ++i)
			{
				if (!visited[i])
				{
					firstUnvisited = i;
					break;
				}
			}
			if (firstUnvisited == numTriangles)
			{
				// All connected components have been found.
				break;
			}

			// Initialize the queue to start at the first unvisited
			// triangle of a connected component.
			std::queue<size_t> triQueue;
			triQueue.push(firstUnvisited);
			visited[firstUnvisited] = true;
			components.push_back(firstUnvisited);

			// Perform the breadth-first search.
			while (!triQueue.empty())
			{
				size_t curIndex = triQueue.front();
				triQueue.pop();

				auto const& curTriangle = triangles[curIndex];
				for (size_t i0 = 0; i0 < 3; ++i0)
				{
					size_t adjIndex = adjacents[curIndex][i0];
					if (adjIndex != std::numeric_limits<size_t>::max() && !visited[adjIndex])
					{
						// The current-triangle has a directed edge
						// <curTriangle[i0],curTriangle[i1]> and there is
						// a triangle adjacent to it.
						size_t i1 = (i0 + 1) % 3;
						auto& adjTriangle = triangles[adjIndex];
						size_t tv0 = curTriangle[i0];
						size_t tv1 = curTriangle[i1];

						// To have the same chirality, it is required
						// that the adjacent triangle have the directed
						// edge <curTriangle[i1],curTriangle[i0]>
						bool sameChirality = true;
						size_t j0{}, j1{};
						for (j0 = 0; j0 < 3; ++j0)
						{
							j1 = (j0 + 1) % 3;
							if (adjTriangle[j0] == tv0)
							{
								if (adjTriangle[j1] == tv1)
								{
									// The adjacent triangle has the same
									// directed edge as the current
									// triangle, so the chiralities do
									// not match.
									sameChirality = false;
								}
								break;
							}
						}
						Check(j0 < 3, "Unexpected condition.");

						if (!sameChirality)
						{
							// Swap the vertices of the adjacent triangle
							// that form the shared directed edge of the
							// current triangle. This requires that the
							// adjacency information for the other two
							// edges of the adjacent triangle be swapped.
							auto& adjAdjacent = adjacents[adjIndex];
							size_t j2 = (j1 + 1) % 3;
							std::swap(adjTriangle[j0], adjTriangle[j1]);
							std::swap(adjAdjacent[j1], adjAdjacent[j2]);
						}

						// The adjacent triangle has been processed, but
						// it might have neighbors that need to be
						// processed. Push the adjacent triangle into the
						// queue to ensure this happens. Insert the
						// adjacent triangle into the active connected
						// component.
						triQueue.push(adjIndex);
						visited[adjIndex] = true;
						components.push_back(adjIndex);
					}
				}
			}
		}

		// Read the comments at the beginning of this function.
		size_t const numSizes = numComponentTriangles.size();
		Check(numSizes > 1, "Unexpected condition.");
		for (size_t i0 = 0, i1 = 1; i1 < numComponentTriangles.size(); i0 = i1++)
		{
			numComponentTriangles[i0] = numComponentTriangles[i1] - numComponentTriangles[i0];
		}
		numComponentTriangles.resize(numSizes - 1);
	}

	void ETManifoldMesh::MakeConsistentChirality()
	{
		std::vector<std::array<size_t, 3>> triangles{};
		std::vector<std::array<size_t, 3>> adjacents{};
		CreateCompactGraph(triangles, adjacents);

		std::vector<size_t> components{}, numComponentTriangles{};
		GetComponentsConsistentChirality(triangles, adjacents,
			components, numComponentTriangles);

		// Only the 'triangles' array is needed to reconstruct the
		// mesh. The other arrays are discarded.
		Clear();

		for (auto const& triangle : triangles)
		{
			int32_t v0 = static_cast<int32_t>(triangle[0]);
			int32_t v1 = static_cast<int32_t>(triangle[1]);
			int32_t v2 = static_cast<int32_t>(triangle[2]);

			// The typecast avoids warnings about not storing the return
			// value in a named variable. The return value is discarded.
			(void)Insert(v0, v1, v2);
		}
	}

	void ETManifoldMesh::GetBoundaryPolygons(std::vector<std::vector<int32_t>>& polygons,
		bool duplicateEndpoints) const
	{
		polygons.clear();

		// Get the boundary edges. The index into the Triangle::T[]
		// adjacency array is also stored to help with the traversal of
		// polygons.
		BoundaryEdgeMap boundaryEdges{};
		for (auto const& triangleElement : mTMap)
		{
			auto const* tri = triangleElement.second.get();
			for (size_t i = 0; i < 3; ++i)
			{
				if (tri->T[i] == nullptr)
				{
					std::array<int32_t, 2> directed{ tri->V[i], tri->V[(i + 1) % 3] };
					BoundaryEdge edge(tri, i, false);
					boundaryEdges.insert(std::make_pair(directed, edge));
				}
			}
		}

		for (auto const& initialEdge : boundaryEdges)
		{
			if (!initialEdge.second.visited)
			{
				std::vector<int32_t> polygon{};
				GetBoundaryPolygon(initialEdge.second.triangle, initialEdge.second.index,
					boundaryEdges, polygon);
				polygons.push_back(std::move(polygon));
			}
		}

		if (!duplicateEndpoints)
		{
			for (auto& polygon : polygons)
			{
				polygon.resize(polygon.size() - 1);
			}
		}
	}

	void ETManifoldMesh::DepthFirstSearch(
		Triangle* tInitial,
		TrianglePtrIntMap& visited,
		std::vector<Triangle*>& component) const
	{
		// Allocate the maximum-size stack that can occur in the
		// depth-first search. The stack is empty when the index top
		// is -1.
		std::vector<Triangle*> tStack(mTMap.size());
		int32_t top = -1;
		tStack[++top] = tInitial;
		while (top >= 0)
		{
			Triangle* tri = tStack[top];
			visited[tri] = 1;
			int32_t i;
			for (i = 0; i < 3; ++i)
			{
				Triangle* adj = tri->T[i];
				if (adj && visited[adj] == 0)
				{
					tStack[++top] = adj;
					break;
				}
			}
			if (i == 3)
			{
				visited[tri] = 2;
				component.push_back(tri);
				--top;
			}
		}
	}

	void ETManifoldMesh::GetBoundaryPolygon(
		Triangle const* initialTriangle,
		size_t initialIndex,
		BoundaryEdgeMap& boundaryEdges,
		std::vector<int32_t>& polygon) const
	{
		auto const* tri = initialTriangle;
		size_t i0 = initialIndex;
		size_t i1 = (i0 + 1) % 3;
		std::array<int32_t, 2> vEdge = { tri->V[i0], tri->V[i1] };
		polygon.push_back(vEdge[0]);
		while (!boundaryEdges[vEdge].visited)
		{
			polygon.push_back(vEdge[1]);
			boundaryEdges[vEdge].visited = true;

			// Traverse the triangle strip with vertex at vEdge[1] until
			// the last triangle is encountered. The final edge of the
			// last triangle is the next boundary edge and starts at
			// vEdge[1].
			std::set<Triangle const*> visited{};
			visited.insert(tri);
			while (tri->T[i1] != nullptr)
			{
				tri = tri->T[i1];
				auto result = visited.insert(tri);

				// It this assertion is trigerred, try calling
				// <mesh>.IsOriented() before calling GetBoundaryPolygons.
				// If <mesh>.IsOriented() returns 'false', the call to
				// GetBoundaryPolygons() will fail.
				Check(result.second, "Triangle already visited. Is the mesh orientable?");

				size_t j{};
				for (j = 0; j < 3; ++j)
				{
					if (vEdge[1] == tri->V[j])
					{
						i1 = j;
						break;
					}
				}
				Check(j < 3, "Unexpected condition.");
			}

			size_t i2 = (i1 + 1) % 3;
			vEdge[0] = vEdge[1];
			vEdge[1] = tri->V[i2];
			i0 = i1;
			i1 = i2;
		}
	}

}
