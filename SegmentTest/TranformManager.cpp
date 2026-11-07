#include "TranformManager.h"

void TransformComponent::Create(Entity e,
	const mathx::Vector<3, double>& position,
	const mathx::Quaternion& orientation,
	const mathx::Vector<3, double>& scale)
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

void TransformComponent::Destroy(Entity e)
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

void TransformComponent::GarbageCollect(const EntityPool& pool)
{
	for (size_t i = entity_.size(); i-- > 0; )
		if (!pool.IsAlive(entity_[i])) Destroy(entity_[i]);
}

bool TransformComponent::SetParent(Entity child, Entity parent)
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

void TransformComponent::Update()
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

void TransformComponent::RebuildOrder()
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
