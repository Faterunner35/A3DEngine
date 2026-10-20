#pragma once
#include"Geomerty.h"
#include"Quaternion.hpp"
#include"Matrix.h"
#include"Affine4X3.hpp"

namespace mathx
{ 
	struct TransformData
	{
		// ВХОД: локальная трансформа относительно родителя (то, что задаёт пользователь)
		Vector<3, double> position;      // 24 байта
		Quaternion       orientation;   // 32
		Vector<3, double> scale;         // 24
		
		// ВЫХОД: результат композиции с родителем
		//   W = W_parent · T·R·S        (у корня W = T·R·S)
		Affine4x3        world;         // 96

		// СВЯЗЬ: где взять W_parent — индекс в тех же массивах (none — корень)
		uint32_t         parent;        // 4
	};


}