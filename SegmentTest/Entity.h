#pragma once
#include <cstdint>
#include <vector>
#include <limits>
#include <cassert>


	struct Entity
	{
		static constexpr uint32_t kInvalid = std::numeric_limits<uint32_t>::max();

		uint32_t index = kInvalid;
		uint32_t generation = 0;

		bool IsNull() const { return index == kInvalid; }
		friend bool operator==(const Entity&, const Entity&) = default;
	};

	// Пул фиксированной ёмкости: память выделяется один раз, Create/Destroy — O(1).
	// Свободные слоты — очередь FIFO внутри массива: освобождённый индекс
	// возвращается в оборот последним, поэтому поколение растёт медленнее.
	class EntityPool
	{
	public:
		explicit EntityPool(uint32_t capacity)
			: slots_(capacity)
		{
			for (uint32_t i = 0; i < capacity; ++i)
				slots_[i].nextFree = (i + 1 < capacity) ? i + 1 : kEnd;
			freeHead_ = capacity ? 0 : kEnd;
			freeTail_ = capacity ? capacity - 1 : kEnd;
		}

		// Пустой Entity, если пул исчерпан.
		Entity Create()
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

		void Destroy(Entity e)
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

		bool IsAlive(Entity e) const
		{
			return e.index < slots_.size()
				&& slots_[e.index].alive
				&& slots_[e.index].generation == e.generation;
		}

		uint32_t Capacity()   const { return uint32_t(slots_.size()); }
		uint32_t AliveCount() const { return aliveCount_; }

	private:
		static constexpr uint32_t kEnd = std::numeric_limits<uint32_t>::max();

		struct Slot
		{
			uint32_t generation = 0;
			uint32_t nextFree = kEnd;
			bool     alive = false;
		};

		std::vector<Slot> slots_;
		uint32_t freeHead_ = kEnd;
		uint32_t freeTail_ = kEnd;
		uint32_t aliveCount_ = 0;
	};
