#pragma once

#include<array>
#include<cmath>
#include<cstdint>
#include<limits>
#include<type_traits>
#include<iostream>

#include<algorithm>

#include<fstream>
#include<string>
#include<sstream>
#include<cassert>



namespace mathx {

	// Тип N — int32_t везде: и в классе, и во всех свободных функциях.
	// При разных типах (uint32_t в классе, int32_t в функциях) вывод шаблона
	// не срабатывает, и operator+, Normalize и т.д. просто не находятся.
	template <int32_t N, typename T>
	class Vector
	{
		static_assert(N > 0, "Vector: N must be positive");

	public:
		// Компоненты всегда инициализированы нулями (см. поле внизу).
		Vector() = default;

		Vector(std::array<T, N> const& values)
			: Vector_Array(values)
		{
		}

		// Меньше N значений — остаток нули. Больше N — ошибка, а не тихое обрезание.
		Vector(std::initializer_list<T> data)
		{
			assert(data.size() <= size_t(N) && "Vector: too many initializers");
			std::copy_n(data.begin(), std::min(data.size(), size_t(N)), Vector_Array.begin());
		}

		// Копирование и присваивание — сгенерированные компилятором (правило нуля).

		static constexpr int32_t GetSize() { return N; }

		T const& operator[](int32_t i) const { return Vector_Array[i]; }
		T& operator[](int32_t i) { return Vector_Array[i]; }

		void SetToVector(int32_t pos, T data) { Vector_Array[pos] = data; }

		std::array<T, N> const& GetVectorArray() const { return Vector_Array; }

		// Особые векторы.
		void MakeZero() { Vector_Array.fill(T(0)); }
		void MakeOnes() { Vector_Array.fill(T(1)); }

		// Компонента d равна 1, остальные 0.
		void MakeUnit(int32_t d)
		{
			MakeZero();
			if (0 <= d && d < N) Vector_Array[d] = T(1);
		}

		// Записать первые data.size() компонент, остальные не трогать.
		void SetNumbers(std::initializer_list<T> data)
		{
			assert(data.size() <= size_t(N) && "SetNumbers: too many values");
			std::copy_n(data.begin(), std::min(data.size(), size_t(N)), Vector_Array.begin());
		}

		static Vector Zero() { return Vector(); }
		static Vector Ones() { Vector v; v.MakeOnes(); return v; }
		static Vector Unit(int32_t d) { Vector v; v.MakeUnit(d); return v; }

		// Лексикографический порядок — только для упорядоченных контейнеров,
		// геометрического смысла у него нет.
		friend bool operator==(Vector const&, Vector const&) = default;
		friend auto operator<=>(Vector const&, Vector const&) = default;

	private:
		std::array<T, N> Vector_Array{};
	};

	// ---------------------------------------------------------------------------
	//  Арифметика. Скаляр — std::type_identity_t<T>, чтобы v * 2 работало для
	//  Vector<3,double>: иначе T выводится и из вектора (double), и из 2 (int),
	//  и вывод шаблона проваливается.
	// ---------------------------------------------------------------------------
	template <int32_t N, typename T>
	Vector<N, T> operator+(Vector<N, T> const& v) { return v; }

	template <int32_t N, typename T>
	Vector<N, T> operator-(Vector<N, T> const& v)
	{
		Vector<N, T> result;
		for (int32_t i = 0; i < N; ++i) result[i] = -v[i];
		return result;
	}

	template <int32_t N, typename T>
	Vector<N, T>& operator+=(Vector<N, T>& v0, Vector<N, T> const& v1)
	{
		for (int32_t i = 0; i < N; ++i) v0[i] += v1[i];
		return v0;
	}

	template <int32_t N, typename T>
	Vector<N, T>& operator-=(Vector<N, T>& v0, Vector<N, T> const& v1)
	{
		for (int32_t i = 0; i < N; ++i) v0[i] -= v1[i];
		return v0;
	}

	template <int32_t N, typename T>
	Vector<N, T>& operator*=(Vector<N, T>& v, std::type_identity_t<T> scalar)
	{
		for (int32_t i = 0; i < N; ++i) v[i] *= scalar;
		return v;
	}

	// Деление на ноль — ошибка вызывающего. Раньше здесь тихо получался нулевой
	// вектор, и ошибка уходила дальше по коду незамеченной.
	template <int32_t N, typename T>
	Vector<N, T>& operator/=(Vector<N, T>& v, std::type_identity_t<T> scalar)
	{
		assert(scalar != T(0) && "Vector: division by zero");
		for (int32_t i = 0; i < N; ++i) v[i] /= scalar;
		return v;
	}

	template <int32_t N, typename T>
	Vector<N, T> operator+(Vector<N, T> v0, Vector<N, T> const& v1) { return v0 += v1; }

	template <int32_t N, typename T>
	Vector<N, T> operator-(Vector<N, T> v0, Vector<N, T> const& v1) { return v0 -= v1; }

	template <int32_t N, typename T>
	Vector<N, T> operator*(Vector<N, T> v, std::type_identity_t<T> scalar) { return v *= scalar; }

	template <int32_t N, typename T>
	Vector<N, T> operator*(std::type_identity_t<T> scalar, Vector<N, T> v) { return v *= scalar; }

	template <int32_t N, typename T>
	Vector<N, T> operator/(Vector<N, T> v, std::type_identity_t<T> scalar) { return v /= scalar; }

	// Покомпонентные операции.
	template <int32_t N, typename T>
	Vector<N, T>& operator*=(Vector<N, T>& v0, Vector<N, T> const& v1)
	{
		for (int32_t i = 0; i < N; ++i) v0[i] *= v1[i];
		return v0;
	}

	// Раньше отсутствовал: operator/(Vector, Vector) вызывал его и не компилировался.
	template <int32_t N, typename T>
	Vector<N, T>& operator/=(Vector<N, T>& v0, Vector<N, T> const& v1)
	{
		for (int32_t i = 0; i < N; ++i) v0[i] /= v1[i];
		return v0;
	}

	template <int32_t N, typename T>
	Vector<N, T> operator*(Vector<N, T> v0, Vector<N, T> const& v1) { return v0 *= v1; }

	template <int32_t N, typename T>
	Vector<N, T> operator/(Vector<N, T> v0, Vector<N, T> const& v1) { return v0 /= v1; }

	// ---------------------------------------------------------------------------
	//  Скалярное произведение, длина, нормирование
	// ---------------------------------------------------------------------------

	// Было: dot = v0[0]*v1[0], потом цикл с i = 0 — первая компонента считалась дважды.
	template <int32_t N, typename T>
	T DotProduct(Vector<N, T> const& v0, Vector<N, T> const& v1)
	{
		T dot = T(0);
		for (int32_t i = 0; i < N; ++i) dot += v0[i] * v1[i];
		return dot;
	}

	// Наибольший модуль компоненты.
	template <int32_t N, typename T>
	T MaxAbsComponent(Vector<N, T> const& v)
	{
		T m = std::fabs(v[0]);
		for (int32_t i = 1; i < N; ++i) m = std::max(m, std::fabs(v[i]));
		return m;
	}

	// robust: делим на наибольшую компоненту перед возведением в квадрат,
	// чтобы квадраты не переполнились и не ушли в денормализованные числа.
	// Было (Lenght): return внутри цикла, множитель — текущая компонента вместо
	// максимальной, для N = 1 функция ничего не возвращала.
	template <int32_t N, typename T>
	T Length(Vector<N, T> const& v, bool robust = false)
	{
		if (!robust) return std::sqrt(DotProduct(v, v));

		const T maxAbs = MaxAbsComponent(v);
		if (maxAbs == T(0)) return T(0);
		const Vector<N, T> scaled = v / maxAbs;
		return maxAbs * std::sqrt(DotProduct(scaled, scaled));
	}

	// Возвращает исходную длину. Нулевой вектор остаётся нулевым.
	template <int32_t N, typename T>
	T Normalize(Vector<N, T>& v, bool robust = false)
	{
		if (robust) {
			const T maxAbs = MaxAbsComponent(v);
			if (maxAbs == T(0)) { v.MakeZero(); return T(0); }
			v /= maxAbs;
			const T length = std::sqrt(DotProduct(v, v));
			v /= length;
			return length * maxAbs;
		}
		const T length = std::sqrt(DotProduct(v, v));
		if (length > T(0)) v /= length;
		else               v.MakeZero();
		return length;
	}

	// Грама–Шмидта (модифицированный: каждая проекция вычитается из уже
	// исправленного v[i]). Возвращает наименьшую из длин до нормирования:
	// близкая к нулю длина значит, что входные векторы почти зависимы.
	// Было: Normalize(v, ...) от указателя — не компилировалось; порядок
	// параметров шаблона <T, N> отличался от остальных функций.
	template <int32_t N, typename T>
	T Orthonormalize(int32_t numInputs, Vector<N, T>* v, bool robust = false)
	{
		if (!v || numInputs < 1 || numInputs > N) return T(0);

		T minLength = Normalize(v[0], robust);
		for (int32_t i = 1; i < numInputs; ++i) {
			for (int32_t j = 0; j < i; ++j) {
				const T dot = DotProduct(v[i], v[j]);
				v[i] -= v[j] * dot;
			}
			minLength = std::min(minLength, Normalize(v[i], robust));
		}
		return minLength;
	}

	// Вектор, перпендикулярный v: берём наибольшую по модулю компоненту imax и
	// следующую за ней inext, результат (…, v[inext], …, −v[imax], …).
	// Скалярное произведение v[imax]·v[inext] − v[inext]·v[imax] = 0.
	// Было: в цикле всегда читалось v[0], поэтому imax всегда был 0, и для
	// векторов вида (0, 0, 1) получался нулевой результат.
	template <int32_t N, typename T>
	Vector<N, T> GetOrthogonal(Vector<N, T> const& v, bool unitLength)
	{
		static_assert(N >= 2, "GetOrthogonal: N must be at least 2");

		int32_t imax = 0;
		T maxAbs = std::fabs(v[0]);
		for (int32_t i = 1; i < N; ++i) {
			const T a = std::fabs(v[i]);
			if (a > maxAbs) { maxAbs = a; imax = i; }
		}
		const int32_t inext = (imax + 1 == N) ? 0 : imax + 1;

		Vector<N, T> result;
		result[imax] = v[inext];
		result[inext] = -v[imax];

		if (unitLength) {
			const T sqrLength = result[imax] * result[imax] + result[inext] * result[inext];
			if (sqrLength > T(0)) {
				const T invLength = T(1) / std::sqrt(sqrLength);
				result[imax] *= invLength;
				result[inext] *= invLength;
			}
		}
		return result;
	}

	// Для N = 4 четвёртая компонента — 0 (вектор, а не точка).
	template <int32_t N, typename T>
	Vector<N, T> Cross(Vector<N, T> const& v0, Vector<N, T> const& v1)
	{
		static_assert(N == 3 || N == 4, "Cross: dimension must be 3 or 4");
		Vector<N, T> result;
		result[0] = v0[1] * v1[2] - v0[2] * v1[1];
		result[1] = v0[2] * v1[0] - v0[0] * v1[2];
		result[2] = v0[0] * v1[1] - v0[1] * v1[0];
		return result;
	}

	template <int32_t N, typename T>
	Vector<N, T> UnitCross(Vector<N, T> const& v0, Vector<N, T> const& v1, bool robust = false)
	{
		Vector<N, T> unitCross = Cross(v0, v1);
		Normalize(unitCross, robust);
		return unitCross;
	}

	// Смешанное произведение v0 · (v1 × v2).
	template <int32_t N, typename T>
	T DotCross(Vector<N, T> const& v0, Vector<N, T> const& v1, Vector<N, T> const& v2)
	{
		return DotProduct(v0, Cross(v1, v2));
	}

	// Покомпонентные минимум и максимум набора векторов (AABB точек).
	template <int32_t N, typename T>
	bool ComputeExtremes(int32_t numVectors, Vector<N, T> const* v,
		Vector<N, T>& vmin, Vector<N, T>& vmax)
	{
		if (!v || numVectors <= 0) return false;

		vmin = v[0];
		vmax = v[0];
		for (int32_t j = 1; j < numVectors; ++j)
			for (int32_t i = 0; i < N; ++i) {
				if (v[j][i] < vmin[i]) vmin[i] = v[j][i];
				else if (v[j][i] > vmax[i]) vmax[i] = v[j][i];
			}
		return true;
	}

	// ---------------------------------------------------------------------------
	//  Подъём и проекция размерности
	// ---------------------------------------------------------------------------

	// (v) → (v, last): однородные координаты.
	template <int32_t N, typename T>
	Vector<N + 1, T> HLift(Vector<N, T> const& v, T last)
	{
		Vector<N + 1, T> result;
		for (int32_t i = 0; i < N; ++i) result[i] = v[i];
		result[N] = last;
		return result;
	}

	// (u, w) → u: отбросить последнюю компоненту (без деления на w).
	template <int32_t N, typename T>
	Vector<N - 1, T> HProject(Vector<N, T> const& v)
	{
		static_assert(N >= 2, "HProject: N must be at least 2");
		Vector<N - 1, T> result;
		for (int32_t i = 0; i < N - 1; ++i) result[i] = v[i];
		return result;
	}

	// Вставить value на место inject: (w0, w1) → (w0, value, w1).
	template <int32_t N, typename T>
	Vector<N + 1, T> Lift(Vector<N, T> const& v, int32_t inject, T value)
	{
		assert(0 <= inject && inject <= N && "Lift: inject out of range");
		Vector<N + 1, T> result;
		for (int32_t i = 0; i < inject; ++i) result[i] = v[i];
		result[inject] = value;
		for (int32_t i = inject; i < N; ++i) result[i + 1] = v[i];
		return result;
	}

	// Убрать компоненту reject: (w0, x, w1) → (w0, w1).
	template <int32_t N, typename T>
	Vector<N - 1, T> Project(Vector<N, T> const& v, int32_t reject)
	{
		static_assert(N >= 2, "Project: N must be at least 2");
		assert(0 <= reject && reject < N && "Project: reject out of range");
		Vector<N - 1, T> result;
		for (int32_t i = 0, j = 0; i < N - 1; ++i, ++j) {
			if (j == reject) ++j;
			result[i] = v[j];
		}
		return result;
	}

	// ---------------------------------------------------------------------------
	//  Отрезок p[0]–p[1]
	// ---------------------------------------------------------------------------
	template <int32_t N, typename T>
	class Segment
	{
	public:
		// По умолчанию: от −e0 до +e0.
		Segment()
		{
			p[1].MakeUnit(0);
			p[0] = -p[1];
		}

		Segment(Vector<N, T> const& p0, Vector<N, T> const& p1)
			: p{ p0, p1 }
		{
		}

		Segment(std::array<Vector<N, T>, 2> const& inP)
			: p(inP)
		{
		}

		Segment(Vector<N, T> const& center, Vector<N, T> const& direction, T extent)
		{
			SetCenteredForm(center, direction, extent);
		}

		// Центр, единичное направление, полудлина.
		void SetCenteredForm(Vector<N, T> const& center, Vector<N, T> const& direction, T extent)
		{
			p[0] = center - extent * direction;
			p[1] = center + extent * direction;
		}

		void GetCenteredForm(Vector<N, T>& center, Vector<N, T>& direction, T& extent) const
		{
			center = T(0.5) * (p[0] + p[1]);
			direction = p[1] - p[0];
			extent = T(0.5) * Normalize(direction);
		}

		std::array<Vector<N, T>, 2> p;

		// Порядок — только для упорядоченных контейнеров.
		friend bool operator==(Segment const&, Segment const&) = default;
		friend auto operator<=>(Segment const&, Segment const&) = default;
	};

} // namespace mathx

	

	//template <typename T>
	//auto readTrianglesFromFile(const std::string& filename)
	//{
	//	std::vector<Triangle> triangles; // To store the parsed triangles
	//	std::ifstream file(filename);   // Open the file for reading

	//	if (!file.is_open()) {
	//		std::cerr << "Error: Could not open file " << filename << "\n";
	//		return triangles; // Return an empty vector if the file cannot be opened
	//	}

	//	std::string line;
	//	while (std::getline(file, line)) {
	//		// Skip empty lines or lines starting with '#'
	//		if (line.empty() || line[0] == '#') {
	//			continue;
	//		}

	//		std::istringstream iss(line);
	//		
	//		std::array<T,3>local_arr{0,0,0};
	//		std::array<T,3>local_arr_2{0,0,0};
	//		std::array<T,3>local_arr_3{0,0,0};

	//		// Read 9 values for the 3 vertices of the triangle
	//		if (!(iss >> local_arr[0] >>local_arr[1] >> local_arr[2] >>
	//			local_arr_2[0] >> local_arr_2[1] >> local_arr_2[2] >>
	//			local_arr_3[0] >> local_arr[1] >>local_arr_3[2])) {
	//			std::cerr << "Error: Invalid line format: " << line << "\n";
	//			continue; // Skip invalid lines
	//		}
 //        

	//	Triangle triangle(local_arr,local_arr_2,local_arr_3);
	//		triangles.push_back(triangle); // Add the parsed triangle to the list
	//	}

	//	file.close(); // Close the file after reading
	//	return triangles;
	//}


	









