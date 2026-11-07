#include "Entity.h"

EntityPool::EntityPool(uint32_t capacity)
	: slots_(capacity)
{
	for (uint32_t i = 0; i < capacity; ++i)
		slots_[i].nextFree = (i + 1 < capacity) ? i + 1 : kEnd;
	freeHead_ = capacity ? 0 : kEnd;
	freeTail_ = capacity ? capacity - 1 : kEnd;
}

Entity EntityPool::Create()
{
	if (freeHead_ == kEnd) return {};

	const uint32_t index = freeHead_;
	Slot& slot = slots_[index];
	freeHead_ = slot.nextFree;
	if (freeHead_ == kEnd) freeTail_ = kEnd;

	slot.nextFree = kEnd;
	slot.alive = true;
	++aliveCount_;
	return { index, slot.generation };
}

void EntityPool::Destroy(Entity e)
{
	if (!IsAlive(e)) return;

	Slot& slot = slots_[e.index];
	slot.alive = false;
	++slot.generation;
	--aliveCount_;

	if (freeTail_ == kEnd) freeHead_ = e.index;
	else slots_[freeTail_].nextFree = e.index;
	freeTail_ = e.index;
}
