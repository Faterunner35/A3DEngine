#pragma once
#include"EdgeKey.h"


namespace Geom_Exp
{
	template <bool Ordered>
	class TriangleKey : public FeatureKey<3, Ordered>
	{
	public:
		TriangleKey()
		{
			this->V = { -1, -1, -1 };
		}

		// This constructor is specialized based on Ordered.
		explicit TriangleKey(int32_t v0, int32_t v1, int32_t v2)
		{
			Initialize(v0, v1, v2);
		}

	private:
		template <bool Dummy = Ordered>
		typename std::enable_if<Dummy, void>::type
			Initialize(int32_t v0, int32_t v1, int32_t v2)
		{
			if (v0 < v1)
			{
				if (v0 < v2)
				{
					// v0 is minimum
					this->V[0] = v0;
					this->V[1] = v1;
					this->V[2] = v2;
				}
				else
				{
					// v2 is minimum
					this->V[0] = v2;
					this->V[1] = v0;
					this->V[2] = v1;
				}
			}
			else
			{
				if (v1 < v2)
				{
					// v1 is minimum
					this->V[0] = v1;
					this->V[1] = v2;
					this->V[2] = v0;
				}
				else
				{
					// v2 is minimum
					this->V[0] = v2;
					this->V[1] = v0;
					this->V[2] = v1;
				}
			}
		}

		template <bool Dummy = Ordered>
		typename std::enable_if<!Dummy, void>::type
			Initialize(int32_t v0, int32_t v1, int32_t v2)
		{
			if (v0 < v1)
			{
				if (v0 < v2)
				{
					// v0 is minimum
					this->V[0] = v0;
					this->V[1] = std::min(v1, v2);
					this->V[2] = std::max(v1, v2);
				}
				else
				{
					// v2 is minimum
					this->V[0] = v2;
					this->V[1] = std::min(v0, v1);
					this->V[2] = std::max(v0, v1);
				}
			}
			else
			{
				if (v1 < v2)
				{
					// v1 is minimum
					this->V[0] = v1;
					this->V[1] = std::min(v2, v0);
					this->V[2] = std::max(v2, v0);
				}
				else
				{
					// v2 is minimum
					this->V[0] = v2;
					this->V[1] = std::min(v0, v1);
					this->V[2] = std::max(v0, v1);
				}
			}
		}
	};
}

