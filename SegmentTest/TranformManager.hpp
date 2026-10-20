#pragma once
#include"Geomerty.h"
#include"Matrix3x3.hpp"
#include"Matrxi4x4.hpp"
#include"Entity.hpp"
#include"TransfromMatrix.hpp"
#include"Affine4X3.hpp"
#include"Quaternion.hpp"
#include<vector>

class TransformComponent
{
public:
	static constexpr uint32_t kNone = std::numeric_limits<uint32_t>::max();

	// ---------------- компонент ----------------

	void Create(Entity e,
		const mathx::Vector<3, double>& position = { 0, 0, 0 },
		const mathx::Quaternion& orientation = mathx::Quaternion::Identity(),
		const mathx::Vector<3, double>& scale = { 1, 1, 1 })
	{
		assert(!e.IsNull() && !Has(e));

		// Индекс сущности переиспользован раньше GarbageCollect: старая запись мертва.
		if (e.index < sparse_.size() && sparse_[e.index] != kNone)
			Destroy(entity_[sparse_[e.index]]);

		if (sparse_.size() <= e.index) sparse_.resize(e.index + 1, kNone);
		sparse_[e.index] = uint32_t(entity_.size());

		entity_.push_back(e);
		parent_.push_back({});
		parentIdx_.push_back(kNone);
		position_.push_back(position);
		orientation_.push_back(orientation);
		scale_.push_back(scale);
		world_.push_back(mathx::MakeTransform(position, orientation, scale));

		worldValid_ = false;
	}

	void Destroy(Entity e)
	{
		const uint32_t i = IndexOf(e);
		if (i == kNone) return;

		const uint32_t last = uint32_t(entity_.size() - 1);
		if (i != last) {
			entity_[i] = entity_[last];
			parent_[i] = parent_[last];
			parentIdx_[i] = parentIdx_[last];
			position_[i] = position_[last];
			orientation_[i] = orientation_[last];
			scale_[i] = scale_[last];
			world_[i] = world_[last];
			sparse_[entity_[i].index] = i;
		}
		entity_.pop_back();
		parent_.pop_back();
		parentIdx_.pop_back();
		position_.pop_back();
		orientation_.pop_back();
		scale_.pop_back();
		world_.pop_back();
		sparse_[e.index] = kNone;

		orderValid_ = false;
		worldValid_ = false;
	}

	void GarbageCollect(const EntityPool& pool)
	{
		for (size_t i = entity_.size(); i-- > 0; )
			if (!pool.IsAlive(entity_[i])) Destroy(entity_[i]);
	}

	bool Has(Entity e) const { return IndexOf(e) != kNone; }
	size_t Size() const { return entity_.size(); }

	// ---------------- родитель ----------------

	// Пустой parent — сделать корнем. false — получился бы цикл.
	bool SetParent(Entity child, Entity parent)
	{
		const uint32_t c = Index(child);
		if (!parent.IsNull()) {
			assert(Has(parent));
			for (Entity a = parent; !a.IsNull(); a = Parent(a))
				if (a == child) return false;
		}

		parent_[c] = parent;
		worldValid_ = false;

		if (orderValid_) {
			const uint32_t p = parent.IsNull() ? kNone : Index(parent);
			if (p == kNone || p < c) parentIdx_[c] = p;   // порядок не нарушен
			else orderValid_ = false;
		}
		return true;
	}

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

	void Update()
	{
		if (!orderValid_) RebuildOrder();

		const size_t n = entity_.size();
		for (size_t i = 0; i < n; ++i) {
			 mathx::Affine4x3 local = mathx::MakeTransform(position_[i], orientation_[i], scale_[i]);
			const uint32_t p = parentIdx_[i];
			world_[i] = (p == kNone) ? local : world_[p] * local;
		}
		worldValid_ = true;
	}

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
	void RebuildOrder()
	{
		const uint32_t n = uint32_t(entity_.size());

		std::vector<uint32_t> parentOf(n);
		for (uint32_t i = 0; i < n; ++i) parentOf[i] = IndexOf(parent_[i]);

		std::vector<uint32_t> depth(n, kNone);
		std::vector<uint32_t> chain;
		uint32_t maxDepth = 0;
		for (uint32_t i = 0; i < n; ++i) {
			uint32_t j = i;
			while (depth[j] == kNone && parentOf[j] != kNone) { chain.push_back(j); j = parentOf[j]; }
			if (depth[j] == kNone) depth[j] = 0;
			for (auto it = chain.rbegin(); it != chain.rend(); ++it) depth[*it] = depth[parentOf[*it]] + 1;
			chain.clear();
			maxDepth = std::max(maxDepth, depth[i]);
		}

		std::vector<uint32_t> start(maxDepth + 2, 0);
		for (uint32_t i = 0; i < n; ++i) ++start[depth[i] + 1];
		for (uint32_t d = 1; d < start.size(); ++d) start[d] += start[d - 1];
		std::vector<uint32_t> order(n);
		for (uint32_t i = 0; i < n; ++i) order[start[depth[i]]++] = i;

		Permute(entity_, order);
		Permute(parent_, order);
		Permute(position_, order);
		Permute(orientation_, order);
		Permute(scale_, order);
		Permute(world_, order);
		for (uint32_t i = 0; i < n; ++i) sparse_[entity_[i].index] = i;

		for (uint32_t i = 0; i < n; ++i) {
			parentIdx_[i] = IndexOf(parent_[i]);
			assert(parentIdx_[i] == kNone || parentIdx_[i] < i);
		}
		orderValid_ = true;
	}

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