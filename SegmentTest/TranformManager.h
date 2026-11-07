#pragma once
#include"Geomerty.h"
#include"Matrix3x3.h"
#include"Matrxi4x4.h"
#include"Entity.h"
#include"TransfromMatrix.h"
#include"Affine4X3.h"
#include"Quaternion.h"
#include<vector>

class TransformComponent
{
public:
	static constexpr uint32_t kNone = std::numeric_limits<uint32_t>::max();

	// ---------------- компонент ----------------

	void Create(Entity e,
		const mathx::Vector<3, double>& position = { 0, 0, 0 },
		const mathx::Quaternion& orientation = mathx::Quaternion::Identity(),
		const mathx::Vector<3, double>& scale = { 1, 1, 1 });

	void Destroy(Entity e);

	void GarbageCollect(const EntityPool& pool);

	bool Has(Entity e) const { return IndexOf(e) != kNone; }
	size_t Size() const { return entity_.size(); }

	// ---------------- родитель ----------------

	// Пустой parent — сделать корнем. false — получился бы цикл.
	bool SetParent(Entity child, Entity parent);

	// Родитель с трансформой; удалённый родитель — как отсутствующий.
	Entity Parent(Entity e) const
	{
		const Entity p = parent_[Index(e)];
		return Has(p) ? p : Entity{};
	}

	// ---------------- локальная трансформа ----------------

	void SetPosition(Entity e, const mathx::Vector<3, double>& p) { position_[Index(e)] = p;    worldValid_ = false; }
	void SetOrientation(Entity e, const mathx::Quaternion& q) { orientation_[Index(e)] = q; worldValid_ = false; }
	void SetScale(Entity e, const mathx::Vector<3, double>& s) { scale_[Index(e)] = s;       worldValid_ = false; }

	const mathx::Vector<3, double>& Position(Entity e)    const { return position_[Index(e)]; }
	const mathx::Quaternion& Orientation(Entity e) const { return orientation_[Index(e)]; }
	const mathx::Vector<3, double>& Scale(Entity e)       const { return scale_[Index(e)]; }

	// ---------------- мировая трансформа ----------------

	void Update();

	const mathx::Affine4x3& World(Entity e) const
	{
		assert(worldValid_ && "World() до Update()");
		return world_[Index(e)];
	}

private:
	uint32_t IndexOf(Entity e) const
	{
		if (e.IsNull() || e.index >= sparse_.size()) return kNone;
		const uint32_t i = sparse_[e.index];
		return (i != kNone && entity_[i] == e) ? i : kNone;
	}

	uint32_t Index(Entity e) const
	{
		const uint32_t i = IndexOf(e);
		assert(i != kNone && "у сущности нет трансформы");
		return i;
	}

	// Сортировка по глубине (подсчётом, устойчивая): корни, затем уровень 1, 2, …
	void RebuildOrder();

	template <class T>
	static void Permute(std::vector<T>& v, const std::vector<uint32_t>& order)
	{
		std::vector<T> out;
		out.reserve(v.size());
		for (uint32_t i : order) out.push_back(v[i]);
		v.swap(out);
	}

	std::vector<uint32_t>         sparse_;       // entity.index → индекс в массивах ниже
	std::vector<Entity>           entity_;
	std::vector<Entity>           parent_;
	std::vector<uint32_t>         parentIdx_;    // индекс родителя; kNone — корень
	std::vector<mathx::Vector<3, double>> position_;
	std::vector<mathx::Quaternion>       orientation_;
	std::vector<mathx::Vector<3, double>> scale_;
	std::vector<mathx::Affine4x3>        world_;

	bool orderValid_ = true;
	bool worldValid_ = false;
};