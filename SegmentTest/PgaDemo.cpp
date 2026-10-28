// PgaDemo.cpp — демонстрация PGA Cl(3,0,1): плоскость через три точки, геометрическое
// произведение, дуализация и meet. Отдельная программа (цель pga_demo).
#include <iostream>
#include<iomanip>
#include "Geomerty.h"
#include"WedgeProduct.h"


int main() {
	std::cout << std::fixed << std::setprecision(3);

	// 1) Плоскость через три точки: R ^ S ^ T  (методы-фабрики + оператор ^)
	Multivector R = Multivector::point(0, 0, 0);
	Multivector S = Multivector::point(1, 0, 0);
	Multivector T = Multivector::point(0, 1, 0);
	Multivector plane = Multivector::plane(R, S, T);        // или: R ^ S ^ T
	std::cout << "== 1. plane R^S^T (expect z=0) ==\n";
	std::cout << "  " << plane << "\n";
	float a = plane[0b1110], b = -plane[0b1101],
		c = plane[0b1011], d = -plane[0b0111];
	std::cout << "  (a,b,c,d)=(" << a << "," << b << "," << c << "," << d
		<< ")  [expect (0,0,1,0)]\n";
	std::cout << "  eval(3,5,0)=" << a * 3 + b * 5 + c * 0 + d << "  eval(3,5,2)="
		<< a * 3 + b * 5 + c * 2 + d << "\n\n";

	// 2) geometric product через оператор *  (метрика PGA {+1,+1,+1,0})
	std::cout << "== 2. geometric product a*b, PGA metric ==\n";
	Multivector e1 = Multivector::e(1), e2 = Multivector::e(2), e4 = Multivector::e(4);
	std::cout << "  e1*e1 = " << (e1 * e1) << "  [1]\n";
	std::cout << "  e4*e4 = " << (e4 * e4) << "  [0  (e4^2=0)]\n";
	std::cout << "  e1*e2 = " << (e1 * e2) << "  [e1e2]\n";
	std::cout << "  e2*e1 = " << (e2 * e1) << "  [-e1e2]\n";
	Multivector::set_metric({ +1,+1,+1,+1 });
	std::cout << "  e4*e4 (euclid) = " << (e4 * e4) << "  [1]\n";
	Multivector::set_metric({ +1,+1,+1,0 });
	std::cout << "\n";

	// 3) дуализация .dual() и meet через оператор &
	std::cout << "== 3. dual + meet (operator &) ==\n";
	std::cout << "  e1.dual() = " << e1.dual() << "\n";
	std::cout << "  e1.dual().dual() = " << e1.dual().dual() << "  [+-e1]\n";
	Multivector planeX = Multivector::plane(
		Multivector::point(0, 0, 0), Multivector::point(0, 1, 0), Multivector::point(0, 0, 1));
	std::cout << "  (z=0) & (x=0) = " << (plane & planeX) << "  [line = ось y]\n\n";

	// 4) цепочка операторов: линия ^ точка = плоскость, всё инфиксно
	std::cout << "== 4. chained: (R^S) is line, ^T gives plane ==\n";
	Multivector L = R ^ S;                    // линия
	std::cout << "  line R^S = " << L << "  (grade " << L.top_grade() << ")\n";
	std::cout << "  (R^S)^T  = " << (L ^ T) << "  (grade " << (L ^ T).top_grade() << ")\n";

	std::cout << "\nOK\n";
	return 0;
}
