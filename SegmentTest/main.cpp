// =============================================================================
//  main.cpp — тесты и бенчмарк ECS-ядра: EntityPool + TransformComponent.
//
//  Часть 1. Корректность
//    1.1  EntityPool: создание, удаление, поколения, устаревшие хэндлы, FIFO
//    1.2  Неравномерный масштаб родителя: точная композиция со сдвигом
//    1.3  20 000 случайных операций против рекурсивного эталона
//    1.4  Удалённый родитель, GarbageCollect, переиспользованный индекс
//    1.5  Matrix3x3: умножение и индексация (регрессия бага m[(r, c)])
//
//  Часть 2. Производительность (десятки тысяч трансформов)
//    2.1  Update(): плоский мир, лес, скелеты (цепочки)
//    2.2  Стоимость перестройки порядка (RebuildOrder)
//    2.3  Запись локальных TRS через Entity (sparse-поиск)
//    2.4  Сравнение с классическим деревом узлов на указателях
//
//  Часть 3. Стресс-нагрузка на иерархию
//    3.1  Глубокая цепочка 1k..20k уровней: точность против аналитического ответа
//    3.2  Широкое дерево: 100k детей у одного корня, гибель корня
//    3.3  Живой мир: рождение, смерть, перевешивание, анимация каждый кадр
//         (10k, 50k; 3.3b — 200k и 500k со сверкой с эталоном и формой иерархии)
//    3.4  Стоимость отдельных операций: SetParent, перестройка порядка, GC
//
//  Режимы: без аргументов — всё; --tests-only; --bench; --stress
//
//  Сборка:  cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
//  Запуск:  build/engine_tests   (Windows: build\Release\engine_tests.exe)
// =============================================================================
#include "TranformManager.hpp"
#include "Projection.hpp"
#include "Frustum.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#ifdef _WIN32
extern "C" __declspec(dllimport) int __stdcall SetConsoleOutputCP(unsigned int codePage);
#endif

using namespace mathx;
using V3 = Vector<3, double>;
static constexpr double kPi = 3.14159265358979323846;

// -----------------------------------------------------------------------------
//  Мини-фреймворк проверок
// -----------------------------------------------------------------------------
static int g_failed = 0;
static int g_checks = 0;

#define CHECK(cond)                                                                 \
    do {                                                                            \
        ++g_checks;                                                                 \
        if (!(cond)) { ++g_failed; std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond); } \
    } while (0)

static double MaxDiff(const Affine4x3& a, const Affine4x3& b)
{
    double e = 0;
    for (int i = 0; i < 12; ++i) e = std::max(e, std::fabs(a[i] - b[i]));
    return e;
}

// Кватернион из оси и угла (ось нормируется)
static Quaternion AxisAngle(double x, double y, double z, double radians)
{
    const double l = std::sqrt(x * x + y * y + z * z);
    const double s = std::sin(radians * 0.5) / l;
    return Quaternion(x * s, y * s, z * s, std::cos(radians * 0.5));
}

static Quaternion RandomQuat(std::mt19937_64& rng)
{
    std::normal_distribution<double> g(0, 1);
    double x = g(rng), y = g(rng), z = g(rng), w = g(rng);
    const double l = std::sqrt(x * x + y * y + z * z + w * w);
    return Quaternion(x / l, y / l, z / l, w / l);
}

// Эталон: W(e) = W(parent) · T·R·S, рекурсивно через публичный API.
static Affine4x3 ReferenceWorld(const TransformComponent& tc, Entity e)
{
    const Affine4x3 local = MakeTransform(tc.Position(e), tc.Orientation(e), tc.Scale(e));
    const Entity p = tc.Parent(e);
    return p.IsNull() ? local : ReferenceWorld(tc, p) * local;
}

// =============================================================================
//  ЧАСТЬ 1. КОРРЕКТНОСТЬ
// =============================================================================
static void TestEntityPool()
{
    std::printf("[1.1] EntityPool\n");
    EntityPool pool(4);

    Entity a = pool.Create(), b = pool.Create();
    CHECK(pool.IsAlive(a) && pool.IsAlive(b));
    CHECK(pool.AliveCount() == 2);

    pool.Destroy(a);
    CHECK(!pool.IsAlive(a));                       // хэндл устарел
    CHECK(pool.AliveCount() == 1);
    pool.Destroy(a);                               // повторное удаление — без эффекта
    CHECK(pool.AliveCount() == 1);

    // FIFO: освобождённый индекс 0 вернётся только после свободных 2 и 3
    Entity c = pool.Create(), d = pool.Create(), e = pool.Create();
    CHECK(c.index == 2 && d.index == 3);
    CHECK(e.index == a.index && e.generation == a.generation + 1);
    CHECK(!pool.IsAlive(a) && pool.IsAlive(e));    // тот же слот, другое поколение

    CHECK(pool.Create().IsNull());                 // пул исчерпан
    CHECK(!pool.IsAlive(Entity{}));                // пустой хэндл
}

static void TestNonUniformScale()
{
    std::printf("[1.2] Неравномерный масштаб родителя: diag(2,1,1), ребёнок повёрнут на 45°\n");
    EntityPool pool(8);
    TransformComponent tc;

    Entity parent = pool.Create(), child = pool.Create();
    tc.Create(parent, V3{ 10, 0, 0 }, Quaternion::Identity(), V3{ 2, 1, 1 });
    tc.Create(child, V3{ 1, 0, 0 }, AxisAngle(0, 0, 1, kPi / 4), V3{ 1, 1, 1 });
    CHECK(tc.SetParent(child, parent));
    tc.Update();

    const Affine4x3& W = tc.World(child);
    const V3 c0 = W.Column(0), c1 = W.Column(1);
    const double dot = c0[0] * c1[0] + c0[1] * c1[1] + c0[2] * c1[2];
    const double l0 = std::sqrt(c0[0] * c0[0] + c0[1] * c0[1] + c0[2] * c0[2]);
    const double l1 = std::sqrt(c1[0] * c1[0] + c1[1] * c1[1] + c1[2] * c1[2]);
    const double angle = std::acos(dot / (l0 * l1)) * 180.0 / kPi;
    const V3 t = GetTranslation(W);

    std::printf("      перенос ребёнка = (%.3f, %.3f, %.3f)   ожидаем (12, 0, 0)\n", t[0], t[1], t[2]);
    std::printf("      угол между осями X и Y ребёнка = %.2f°  (не 90°: сдвиг, TRS его не выражает)\n", angle);
    CHECK(std::fabs(t[0] - 12) < 1e-12 && std::fabs(t[1]) < 1e-12);
    CHECK(std::fabs(dot + 1.5) < 1e-12);
    CHECK(std::fabs(angle - 126.8698976) < 1e-6);
}

static void TestRandomAgainstReference()
{
    std::printf("[1.3] 20 000 случайных операций против рекурсивного эталона\n");
    constexpr uint32_t kCap = 512;
    EntityPool pool(kCap);
    TransformComponent tc;
    std::vector<Entity> alive;
    std::mt19937_64 rng(2026);
    std::uniform_real_distribution<double> u(-3, 3), us(0.3, 2.5);

    int cycles = 0, ops[6] = {};
    double worst = 0;

    for (int step = 0; step < 20000; ++step) {
        const int op = int(rng() % 6);
        ++ops[op];
        if ((op == 0 || alive.empty()) && pool.AliveCount() < kCap) {             // создать
            Entity e = pool.Create();
            tc.Create(e, V3{ u(rng), u(rng), u(rng) }, RandomQuat(rng), V3{ us(rng), us(rng), us(rng) });
            alive.push_back(e);
        }
        else if (op == 1 && !alive.empty()) {                                     // удалить
            const size_t k = rng() % alive.size();
            tc.Destroy(alive[k]);
            pool.Destroy(alive[k]);
            alive[k] = alive.back(); alive.pop_back();
        }
        else if (op == 2 && alive.size() > 1) {                                   // перевесить
            Entity c = alive[rng() % alive.size()];
            Entity p = (rng() % 5 == 0) ? Entity{} : alive[rng() % alive.size()];
            if (!tc.SetParent(c, p)) ++cycles;
        }
        else if (op == 3 && !alive.empty()) {
            tc.SetPosition(alive[rng() % alive.size()], V3{ u(rng), u(rng), u(rng) });
        }
        else if (op == 4 && !alive.empty()) {
            tc.SetOrientation(alive[rng() % alive.size()], RandomQuat(rng));
        }
        else if (op == 5 && !alive.empty()) {
            tc.SetScale(alive[rng() % alive.size()], V3{ us(rng), us(rng), us(rng) });
        }

        if (step % 97 == 0 || step == 19999) {
            tc.Update();
            for (Entity e : alive) worst = std::max(worst, MaxDiff(tc.World(e), ReferenceWorld(tc, e)));
            CHECK(tc.Size() == alive.size());
        }
    }
    std::printf("      живых в конце: %zu, отклонено циклов: %d, max |W - эталон| = %.2e\n",
                alive.size(), cycles, worst);
    CHECK(worst < 1e-9);
    CHECK(cycles > 0);
}

static void TestParentDestroyAndGC()
{
    std::printf("[1.4] Удалённый родитель, GarbageCollect, переиспользованный индекс\n");
    EntityPool pool(3);
    TransformComponent tc;

    Entity p = pool.Create(), c = pool.Create();
    tc.Create(p, V3{ 5, 0, 0 });
    tc.Create(c, V3{ 1, 0, 0 });
    tc.SetParent(c, p);
    tc.Update();
    CHECK(std::fabs(GetTranslation(tc.World(c))[0] - 6) < 1e-12);

    // Родитель удалён из пула, но компонент ещё держит запись → GC её снимает
    pool.Destroy(p);
    tc.GarbageCollect(pool);
    CHECK(!tc.Has(p));
    CHECK(tc.Parent(c).IsNull());                  // ребёнок стал корнем
    tc.Update();
    CHECK(std::fabs(GetTranslation(tc.World(c))[0] - 1) < 1e-12);

    // Индекс удалённой сущности переиспользован ДО GarbageCollect
    pool.Destroy(c);
    Entity x = pool.Create();                      // слот 2
    Entity y = pool.Create();                      // слот 0 (бывший p)
    Entity z = pool.Create();                      // слот 1 (бывший c), поколение +1
    CHECK(z.index == c.index && z.generation == c.generation + 1);
    tc.Create(z, V3{ 7, 0, 0 });                   // старая запись c удаляется внутри Create
    CHECK(!tc.Has(c) && tc.Has(z));
    CHECK(tc.Size() == 1);
    tc.Update();
    CHECK(std::fabs(GetTranslation(tc.World(z))[0] - 7) < 1e-12);
    (void)x; (void)y;
}

// Регрессия: раньше Matrix3x3 индексировалась как m[(r, c)] — это оператор
// «запятая», читается m[c]. Умножение, GetRow/GetColumn, HasScale, IsOrthonormal
// молча считали мусор.
static void TestMatrix3x3Indexing()
{
    std::printf("[1.5] Matrix3x3: умножение и доступ к строкам/столбцам\n");
    const Matrix3x3 A(1, 2, 3, 4, 5, 6, 7, 8, 10);
    const Matrix3x3 B(2, 0, 1, 1, 3, 0, 0, 1, 4);
    const Matrix3x3 C = A * B;
    const Matrix3x3 expect(4, 9, 13, 13, 21, 28, 22, 34, 47);
    double e = 0;
    for (int i = 0; i < 9; ++i) e = std::max(e, std::fabs(C[i] - expect[i]));
    CHECK(e == 0.0);
    const V3 r1 = Matrix3x3::GetRow(A, 1), c2 = Matrix3x3::GetColumn(A, 2);
    CHECK(r1[0] == 4 && r1[1] == 5 && r1[2] == 6);
    CHECK(c2[0] == 3 && c2[1] == 6 && c2[2] == 10);
    CHECK(!HasScale(Matrix3x3(0, -1, 0, 1, 0, 0, 0, 0, 1)));       // поворот на 90°
    CHECK(HasScale(Matrix3x3(2, 0, 0, 0, 1, 0, 0, 0, 1)));
    CHECK(IsOrthonormal(Matrix3x3(0, -1, 0, 1, 0, 0, 0, 0, 1)));
    CHECK(!IsOrthonormal(Matrix3x3(1, 1, 0, 0, 1, 0, 0, 0, 1)));    // сдвиг
}

// [1.6] Проекции FGED (6.36) и обратная бесконечная. Камера: x вправо, y вниз,
// z вперёд. Числа — разбор FGED 6.3: fovY = 60°, s = 16/9, n = 1, f = 100.
static void TestProjection()
{
    std::printf("[1.6] Проекции FGED: (6.36), обратная бесконечная, ClipY\n");
    const double fovY = 3.14159265358979323846 / 3.0, s = 16.0 / 9.0, n = 1.0, f = 100.0;
    const double g = 1.0 / std::tan(fovY * 0.5);                       // 1.7320508
    auto near_ = [](double a, double b, double eps = 1e-12) { return std::fabs(a - b) <= eps; };
    using V4 = Vector<4, double>;

    // (6.36) как в листинге 6.1 (y клипа вниз, Vulkan): элементы матрицы
    const Matrix4x4 Pv = MakeFrustumProjection(fovY, s, n, f, ClipY::Down);
    CHECK(near_(Pv(0, 0), g / s) && near_(Pv(1, 1), g));
    CHECK(near_(Pv(2, 2), f / (f - n)) && near_(Pv(2, 3), -n * f / (f - n)));
    CHECK(Pv(3, 2) == 1.0 && Pv(3, 3) == 0.0);

    // Точка (2,-1,4): NDC (0.487, -0.433, 0.758), как в разборе
    const V3 ndc = PerspectiveDivide(Pv * V4{ 2.0, -1.0, 4.0, 1.0 });
    CHECK(near_(ndc[0], (g / s) * 2.0 / 4.0) && near_(ndc[1], -g / 4.0));
    CHECK(near_(ndc[2], 0.7575757575757576));

    // ближняя -> 0, дальняя -> 1, и обратно через (6.38)
    CHECK(near_(PerspectiveDivide(Pv * V4{ 0.0, 0.0, n, 1.0 })[2], 0.0));
    CHECK(near_(PerspectiveDivide(Pv * V4{ 0.0, 0.0, f, 1.0 })[2], 1.0));
    CHECK(near_(CameraDepthFromFrustum(ndc[2], n, f), 4.0, 1e-10));

    // D3D12: y клипа вверх. Точка выше оси (y камеры < 0) уходит в y_ndc > 0.
    const Matrix4x4 Pd = MakeFrustumProjection(fovY, s, n, f, ClipY::Up);
    CHECK(near_(Pd(1, 1), -g));
    CHECK(PerspectiveDivide(Pd * V4{ 0.0, -1.0, 4.0, 1.0 })[1] > 0.0);
    CHECK(PerspectiveDivide(Pv * V4{ 0.0, -1.0, 4.0, 1.0 })[1] < 0.0);

    // Обратная бесконечная: ближняя -> 1, d = 4 -> n/4, направление -> z_clip = 0 ровно
    const Matrix4x4 Pr = MakeReversedInfiniteProjection(fovY, s, n, ClipY::Up);
    CHECK(near_(PerspectiveDivide(Pr * V4{ 0.0, 0.0, n, 1.0 })[2], 1.0));
    const double zr = PerspectiveDivide(Pr * V4{ 0.0, 0.0, 4.0, 1.0 })[2];
    CHECK(near_(zr, n / 4.0) && near_(CameraDepthFromReversedInfinite(zr, n), 4.0));
    const V4 dir = Pr * V4{ 0.3, 0.2, 0.9, 0.0 };
    CHECK(dir[2] == 0.0 && dir[3] > 0.0);
    CHECK(IsInsideFrustum(Pr * V4{ 0.5, 0.2, 50.0, 1.0 }, DepthRange::ZeroToOne));
    CHECK(!IsInsideFrustum(Pr * V4{ 0.0, 0.0, 0.5 * n, 1.0 }, DepthRange::ZeroToOne)); // ближе near
    CHECK(!IsInsideFrustum(Pr * V4{ 0.0, 0.0, -5.0, 1.0 }, DepthRange::ZeroToOne));    // за камерой

    // Обход: треугольник лицом к камере, CCW в правой тройке ((b-a)x(c-a) = (0,0,-1)),
    // в окне D3D12 (y окна вниз) остаётся CCW -> нужен FrontCounterClockwise = TRUE.
    const V4 tri[3] = { { 0, 0, 1, 1 }, { 1, 0, 1, 1 }, { 0, -1, 1, 1 } };
    double wx[3], wy[3];
    for (int i = 0; i < 3; ++i) { const V3 q = PerspectiveDivide(Pd * tri[i]); wx[i] = q[0]; wy[i] = -q[1]; }
    const double area = (wx[1] - wx[0]) * (wy[2] - wy[0]) - (wy[1] - wy[0]) * (wx[2] - wx[0]);
    CHECK(area < 0.0);   // в координатах с y вниз отрицательная площадь = против часовой
}

// =============================================================================
//  ЧАСТЬ 2. ПРОИЗВОДИТЕЛЬНОСТЬ
// =============================================================================
using Clock = std::chrono::steady_clock;

// printf считает байты, а кириллица в UTF-8 — 2 байта на букву: выравниваем по символам.
static std::string Pad(const char* text, size_t width)
{
    std::string r(text);
    size_t chars = 0;
    for (unsigned char c : r) chars += (c & 0xC0) != 0x80;
    if (chars < width) r.append(width - chars, ' ');
    return r;
}

static double Median(std::vector<double> v)
{
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

// Сумма переносов — чтобы компилятор не выбросил вычисления
static double Checksum(const TransformComponent& tc, const std::vector<Entity>& es)
{
    double s = 0;
    for (size_t i = 0; i < es.size(); i += 61) s += tc.World(es[i])[3];
    return s;
}

enum class Shape { Flat, Forest, Skeletons };

static const char* ShapeName(Shape s)
{
    switch (s) {
    case Shape::Flat:      return "плоский (все корни)";
    case Shape::Forest:    return "лес (случ. родитель)";
    case Shape::Skeletons: return "скелеты (цепочки 32)";
    }
    return "";
}

struct Scene {
    EntityPool                pool;
    TransformComponent        tc;
    std::vector<Entity>       entities;
    explicit Scene(uint32_t n) : pool(n) {}
};

// Сущности создаются в перемешанном порядке, родители назначаются так,
// что родитель часто оказывается в массиве ПОСЛЕ ребёнка — первый Update
// обязан перестроить порядок.
static std::unique_ptr<Scene> BuildScene(uint32_t n, Shape shape, std::mt19937_64& rng)
{
    auto s = std::make_unique<Scene>(n);
    std::uniform_real_distribution<double> u(-1, 1), us(0.5, 1.5);
    s->entities.reserve(n);
    for (uint32_t i = 0; i < n; ++i) {
        Entity e = s->pool.Create();
        s->tc.Create(e, V3{ u(rng), u(rng), u(rng) }, RandomQuat(rng), V3{ us(rng), us(rng), us(rng) });
        s->entities.push_back(e);
    }
    // логический порядок иерархии — случайная перестановка массива
    std::vector<uint32_t> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);

    for (uint32_t k = 1; k < n; ++k) {
        const Entity child = s->entities[perm[k]];
        if (shape == Shape::Forest && k % 8 != 0)                    // 1/8 — корни
            s->tc.SetParent(child, s->entities[perm[rng() % k]]);
        else if (shape == Shape::Skeletons && k % 32 != 0)           // цепочка из 32 костей
            s->tc.SetParent(child, s->entities[perm[k - 1]]);
    }
    return s;
}

// Классическое дерево узлов: каждый узел — отдельная аллокация, дети по указателям.
struct PtrNode {
    V3 position, scale; Quaternion orientation;
    Affine4x3 world;
    std::vector<PtrNode*> children;
};

static void UpdatePtr(PtrNode* n, const Affine4x3& parentWorld)
{
    n->world = parentWorld * MakeTransform(n->position, n->orientation, n->scale);
    for (PtrNode* c : n->children) UpdatePtr(c, n->world);
}

static void BenchUpdate(uint32_t n, Shape shape, int frames)
{
    std::mt19937_64 rng(n * 7 + int(shape));
    auto s = BuildScene(n, shape, rng);

    // первый Update — с перестройкой порядка
    auto t0 = Clock::now();
    s->tc.Update();
    const double firstMs = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();

    std::vector<double> ms;
    double sink = 0;
    for (int f = 0; f < frames; ++f) {
        // анимация: каждый кадр меняется позиция у 1/16 сущностей
        for (size_t i = f % 16; i < s->entities.size(); i += 16)
            s->tc.SetPosition(s->entities[i], V3{ double(f), 0, 0 });
        auto a = Clock::now();
        s->tc.Update();
        ms.push_back(std::chrono::duration<double, std::milli>(Clock::now() - a).count());
        sink += Checksum(s->tc, s->entities);
    }
    const double med = Median(ms);
    std::printf("  %8u  %s %9.3f %9.3f %8.1f %9.1f %8.2f%%\n",
                n, Pad(ShapeName(shape), 22).c_str(), firstMs, med, med * 1e6 / n, n / med / 1e3,
                med / 16.667 * 100.0);
    if (sink == 12345.678) std::printf("%f", sink);
}

static void BenchWrites(uint32_t n)
{
    std::mt19937_64 rng(99);
    auto s = BuildScene(n, Shape::Forest, rng);
    s->tc.Update();

    // доступ в случайном порядке сущностей (как пишет физика или анимация)
    std::vector<Entity> order = s->entities;
    std::shuffle(order.begin(), order.end(), rng);

    std::vector<double> ms;
    for (int f = 0; f < 20; ++f) {
        auto a = Clock::now();
        for (Entity e : order) s->tc.SetPosition(e, V3{ double(f), 1, 2 });
        ms.push_back(std::chrono::duration<double, std::milli>(Clock::now() - a).count());
    }
    const double med = Median(ms);
    std::printf("  %8u  SetPosition через Entity (случ. порядок): %7.3f мс, %5.1f нс/запись\n",
                n, med, med * 1e6 / n);
}

static void BenchVsPointerTree(uint32_t n)
{
    std::mt19937_64 rng(5);
    auto s = BuildScene(n, Shape::Forest, rng);
    s->tc.Update();

    // то же дерево на указателях; узлы аллоцированы вразброс
    std::vector<std::unique_ptr<PtrNode>> owner(n);
    std::vector<uint32_t> allocOrder(n);
    std::iota(allocOrder.begin(), allocOrder.end(), 0);
    std::shuffle(allocOrder.begin(), allocOrder.end(), rng);
    std::vector<std::unique_ptr<char[]>> noise;                     // разносит узлы по куче
    for (uint32_t i : allocOrder) {
        owner[i] = std::make_unique<PtrNode>();
        noise.emplace_back(new char[64 + rng() % 512]);
    }
    std::vector<PtrNode*> roots;
    for (uint32_t i = 0; i < n; ++i) {
        const Entity e = s->entities[i];
        PtrNode* node = owner[i].get();
        node->position = s->tc.Position(e);
        node->orientation = s->tc.Orientation(e);
        node->scale = s->tc.Scale(e);
        const Entity p = s->tc.Parent(e);
        if (p.IsNull()) roots.push_back(node);
        else owner[p.index]->children.push_back(node);       // index совпадает с позицией в entities
    }

    std::vector<double> msEcs, msPtr;
    for (int f = 0; f < 15; ++f) {
        auto a = Clock::now();
        s->tc.Update();
        msEcs.push_back(std::chrono::duration<double, std::milli>(Clock::now() - a).count());
        a = Clock::now();
        for (PtrNode* r : roots) UpdatePtr(r, Affine4x3::Identity());
        msPtr.push_back(std::chrono::duration<double, std::milli>(Clock::now() - a).count());
    }

    // результаты совпадают
    double worst = 0;
    for (uint32_t i = 0; i < n; i += 37)
        worst = std::max(worst, MaxDiff(s->tc.World(s->entities[i]), owner[i]->world));
    CHECK(worst < 1e-12);

    const double e = Median(msEcs), p = Median(msPtr);
    std::printf("  %8u  TransformComponent %8.3f мс | дерево указателей %8.3f мс | x%.1f   (|разница| %.0e)\n",
                n, e, p, p / e, worst);
}

// =============================================================================
//  ЧАСТЬ 3. СТРЕСС-НАГРУЗКА НА ИЕРАРХИЮ
// =============================================================================

// Эталон без рекурсии и с мемоизацией: на 50 000 сущностей — O(n),
// на цепочке глубиной 20 000 — без переполнения стека.
// Возвращает max |W − эталон| по всем живым сущностям.
static double CheckAllAgainstReference(const TransformComponent& tc, const std::vector<Entity>& alive,
                                       uint32_t capacity)
{
    std::vector<Affine4x3> memo(capacity);
    std::vector<uint32_t>  memoGen(capacity, ~0u);
    std::vector<Entity>    chain;
    double worst = 0;

    for (Entity e : alive) {
        // подняться до первого посчитанного предка или до корня
        Entity a = e;
        while (!a.IsNull() && memoGen[a.index] != a.generation) { chain.push_back(a); a = tc.Parent(a); }
        Affine4x3 W = a.IsNull() ? Affine4x3::Identity() : memo[a.index];
        bool haveParent = !a.IsNull();
        for (auto it = chain.rbegin(); it != chain.rend(); ++it) {
            const Affine4x3 L = MakeTransform(tc.Position(*it), tc.Orientation(*it), tc.Scale(*it));
            W = haveParent ? W * L : L;
            haveParent = true;
            memo[it->index] = W;
            memoGen[it->index] = it->generation;
        }
        chain.clear();
        worst = std::max(worst, MaxDiff(tc.World(e), memo[e.index]));
    }
    return worst;
}

// Форма иерархии: число корней, максимальная и средняя глубина. O(n) с мемоизацией.
struct HierarchyShape { size_t roots; uint32_t maxDepth; double avgDepth; };
static HierarchyShape MeasureHierarchy(const TransformComponent& tc, const std::vector<Entity>& alive,
                                       uint32_t capacity)
{
    std::vector<uint32_t> depth(capacity, ~0u), depthGen(capacity, ~0u);
    std::vector<Entity> chain;
    HierarchyShape h{ 0, 0, 0.0 };
    for (Entity e : alive) {
        Entity a = e;
        while (!a.IsNull() && depthGen[a.index] != a.generation) { chain.push_back(a); a = tc.Parent(a); }
        uint32_t d = a.IsNull() ? ~0u : depth[a.index];          // ~0u + 1 == 0 для корня
        for (auto it = chain.rbegin(); it != chain.rend(); ++it) {
            d = d + 1;
            depth[it->index] = d;
            depthGen[it->index] = it->generation;
        }
        chain.clear();
        const uint32_t de = depth[e.index];
        h.roots += (de == 0);
        h.maxDepth = std::max(h.maxDepth, de);
        h.avgDepth += de;
    }
    if (!alive.empty()) h.avgDepth /= double(alive.size());
    return h;
}

struct FrameStats { double p50, p99, max; };
static FrameStats Stats(std::vector<double> v)
{
    std::sort(v.begin(), v.end());
    return { v[v.size() / 2], v[std::min(v.size() - 1, v.size() * 99 / 100)], v.back() };
}

static double Ms(Clock::time_point a, Clock::time_point b)
{
    return std::chrono::duration<double, std::milli>(b - a).count();
}

// -----------------------------------------------------------------------------
// 3.1  Глубокая цепочка: каждый узел сдвинут на (1,0,0) и повёрнут на θ = 2π/1000.
//      Точный ответ известен: узел k стоит в вершине правильного 1000-угольника,
//          A_k = R((k+1)θ),   t_k = Σ_{j=0..k} R(jθ)·(1,0,0)
//      Ошибка показывает накопление округления по глубине иерархии.
// -----------------------------------------------------------------------------
static void StressDeepChain(uint32_t depth)
{
    const double theta = 2 * kPi / 1000.0;
    EntityPool pool(depth);
    TransformComponent tc;
    std::vector<Entity> chain;
    chain.reserve(depth);

    auto t0 = Clock::now();
    for (uint32_t k = 0; k < depth; ++k) {
        Entity e = pool.Create();
        tc.Create(e, V3{ 1, 0, 0 }, AxisAngle(0, 0, 1, theta), V3{ 1, 1, 1 });
        chain.push_back(e);
    }
    auto t1 = Clock::now();
    for (uint32_t k = 1; k < depth; ++k) tc.SetParent(chain[k], chain[k - 1]);
    auto t2 = Clock::now();
    tc.Update();                                    // с перестройкой порядка
    auto t3 = Clock::now();
    tc.Update();
    auto t4 = Clock::now();

    // точный ответ в long double
    double worstA = 0, worstT = 0, worstOrtho = 0;
    long double sx = 0, sy = 0;
    for (uint32_t k = 0; k < depth; ++k) {
        sx += std::cos((long double)k * theta);
        sy += std::sin((long double)k * theta);
        if (k % 97 != 0 && k != depth - 1) continue;
        const Affine4x3& W = tc.World(chain[k]);
        const long double a = (long double)(k + 1) * theta;
        const double c = double(std::cos(a)), s = double(std::sin(a));
        worstA = std::max({ worstA, std::fabs(W(0, 0) - c), std::fabs(W(0, 1) + s),
                            std::fabs(W(1, 0) - s), std::fabs(W(1, 1) - c) });
        worstT = std::max({ worstT, std::fabs(W(0, 3) - double(sx)), std::fabs(W(1, 3) - double(sy)) });
        const V3 x = W.Column(0), y = W.Column(1);
        worstOrtho = std::max({ worstOrtho, std::fabs(x[0] * y[0] + x[1] * y[1] + x[2] * y[2]),
                                std::fabs(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 1) });
    }
    std::printf("  глубина %6u: SetParent всей цепочки %8.1f мс, 1-й Update %6.2f мс, Update %6.2f мс\n",
                depth, Ms(t1, t2), Ms(t2, t3), Ms(t3, t4));
    std::printf("                  ошибка поворота %.1e, переноса %.1e, неортогональность осей %.1e\n",
                worstA, worstT, worstOrtho);
    CHECK(worstA < 1e-9 && worstT < 1e-6);
    (void)t0;
}

// -----------------------------------------------------------------------------
// 3.2  Широкое дерево: один корень и N детей. Затем корень удаляется —
//      все дети одновременно становятся корнями.
// -----------------------------------------------------------------------------
static void StressWideTree(uint32_t children)
{
    EntityPool pool(children + 1);
    TransformComponent tc;
    std::mt19937_64 rng(31);
    std::uniform_real_distribution<double> u(-10, 10);

    Entity root = pool.Create();
    tc.Create(root, V3{ 100, 0, 0 }, AxisAngle(0, 1, 0, 0.7), V3{ 2, 2, 2 });
    std::vector<Entity> kids;
    kids.reserve(children);
    for (uint32_t i = 0; i < children; ++i) {
        Entity e = pool.Create();
        tc.Create(e, V3{ u(rng), u(rng), u(rng) }, RandomQuat(rng));
        tc.SetParent(e, root);
        kids.push_back(e);
    }
    tc.Update();
    auto a = Clock::now();
    tc.Update();
    const double tUpd = Ms(a, Clock::now());

    std::vector<Entity> all = kids;
    all.push_back(root);
    CHECK(CheckAllAgainstReference(tc, all, pool.Capacity()) < 1e-9);

    // корень изменил позицию — сдвинулись все дети
    tc.SetPosition(root, V3{ -100, 0, 0 });
    tc.Update();
    CHECK(CheckAllAgainstReference(tc, all, pool.Capacity()) < 1e-9);

    a = Clock::now();
    pool.Destroy(root);
    tc.GarbageCollect(pool);
    tc.Update();
    const double tKill = Ms(a, Clock::now());

    size_t roots = 0;
    for (Entity e : kids) roots += tc.Parent(e).IsNull();
    CHECK(roots == kids.size());
    CHECK(CheckAllAgainstReference(tc, kids, pool.Capacity()) < 1e-12);
    std::printf("  %u детей у одного корня: Update %.2f мс; удаление корня + GC + Update %.2f мс; "
                "корнями стали %zu/%zu\n", children, tUpd, tKill, roots, kids.size());
}

// -----------------------------------------------------------------------------
// 3.3  Живой мир: каждый кадр часть сущностей умирает, рождается, перевешивается
//      и анимируется. Каждые 50 кадров — полная сверка с эталоном.
// -----------------------------------------------------------------------------
struct ChurnConfig {
    uint32_t population;
    int      frames;
    double   destroyPct, createPct, reparentPct, animatePct;
};

static void StressChurn(const ChurnConfig& cfg)
{
    const uint32_t capacity = cfg.population * 2;
    EntityPool pool(capacity);
    TransformComponent tc;
    std::vector<Entity> alive;
    std::mt19937_64 rng(cfg.population + cfg.frames);
    std::uniform_real_distribution<double> u(-5, 5), us(0.5, 2.0);

    auto spawn = [&]() {
        Entity e = pool.Create();
        if (e.IsNull()) return;
        tc.Create(e, V3{ u(rng), u(rng), u(rng) }, RandomQuat(rng), V3{ us(rng), us(rng), us(rng) });
        if (!alive.empty() && rng() % 8 != 0) tc.SetParent(e, alive[rng() % alive.size()]);
        alive.push_back(e);
    };
    for (uint32_t i = 0; i < cfg.population; ++i) spawn();
    tc.Update();

    const auto count = [&](double pct) { return size_t(alive.size() * pct / 100.0); };
    constexpr size_t kAnimPool = 4096;                     // степень двойки
    std::vector<V3> animPos(kAnimPool);
    std::vector<Quaternion> animRot(kAnimPool);
    for (size_t i = 0; i < kAnimPool; ++i) { animPos[i] = V3{ u(rng), u(rng), u(rng) }; animRot[i] = RandomQuat(rng); }
    size_t poolIdx = 0;
    std::vector<double> frameMs, opsMs, updMs;
    double worst = 0;
    size_t created = 0, destroyed = 0, reparented = 0, rejected = 0;

    for (int f = 0; f < cfg.frames; ++f) {
        auto a = Clock::now();
        // смерть: сущность уходит из пула; компонент чистится сборщиком
        for (size_t k = count(cfg.destroyPct); k > 0 && !alive.empty(); --k) {
            const size_t i = rng() % alive.size();
            pool.Destroy(alive[i]);
            alive[i] = alive.back(); alive.pop_back();
            ++destroyed;
        }
        tc.GarbageCollect(pool);
        // рождение
        for (size_t k = count(cfg.createPct); k > 0; --k) { spawn(); ++created; }
        // перевешивание (циклы отклоняются)
        for (size_t k = count(cfg.reparentPct); k > 0; --k) {
            Entity c = alive[rng() % alive.size()];
            Entity p = (rng() % 10 == 0) ? Entity{} : alive[rng() % alive.size()];
            if (tc.SetParent(c, p)) ++reparented; else ++rejected;
        }
        // анимация: значения берутся из заранее сгенерированного пула, чтобы
        // в «операции» попадала стоимость записи в компонент, а не генератора
        for (size_t k = count(cfg.animatePct); k > 0; --k) {
            Entity e = alive[rng() % alive.size()];
            tc.SetPosition(e, animPos[poolIdx]);
            tc.SetOrientation(e, animRot[poolIdx]);
            poolIdx = (poolIdx + 1) & (kAnimPool - 1);
        }
        auto b = Clock::now();
        tc.Update();
        auto c = Clock::now();
        opsMs.push_back(Ms(a, b));
        updMs.push_back(Ms(b, c));
        frameMs.push_back(Ms(a, c));

        const int checkEvery = std::max(1, std::min(50, cfg.frames / 4));
        if (f % checkEvery == checkEvery - 1 || f == cfg.frames - 1) {
            CHECK(tc.Size() == alive.size());
            worst = std::max(worst, CheckAllAgainstReference(tc, alive, capacity));
        }
    }
    const FrameStats fs = Stats(frameMs), us_ = Stats(updMs), os = Stats(opsMs);
    std::printf("  %6u сущ., %d кадров, за кадр: −%.1f%% +%.1f%% ↻%.1f%% ~%.0f%%\n",
                cfg.population, cfg.frames, cfg.destroyPct, cfg.createPct, cfg.reparentPct, cfg.animatePct);
    std::printf("         кадр p50 %6.2f мс  p99 %6.2f мс  max %6.2f мс   (операции p50 %.2f, Update p50 %.2f)\n",
                fs.p50, fs.p99, fs.max, os.p50, us_.p50);
    std::printf("         родилось %zu, умерло %zu, перевешано %zu, отклонено циклов %zu, |W − эталон| = %.1e\n",
                created, destroyed, reparented, rejected, worst);
    const HierarchyShape hs = MeasureHierarchy(tc, alive, capacity);
    std::printf("         иерархия в конце: корней %zu (%.1f%%), глубина max %u, средняя %.1f; память ≈ %.0f МБ\n",
                hs.roots, 100.0 * hs.roots / std::max<size_t>(1, alive.size()), hs.maxDepth, hs.avgDepth,
                alive.size() * 200.0 / (1024 * 1024));
    CHECK(worst < 1e-9);
}

// -----------------------------------------------------------------------------
// 3.4  Стоимость отдельных операций на большом мире
// -----------------------------------------------------------------------------
static void StressOperationCosts(uint32_t n)
{
    std::mt19937_64 rng(77);
    auto s = BuildScene(n, Shape::Forest, rng);
    s->tc.Update();

    // SetParent без нарушения порядка (родитель левее) и с нарушением
    std::vector<Entity>& es = s->entities;
    const int kOps = 2000;
    auto a = Clock::now();
    for (int i = 0; i < kOps; ++i) s->tc.SetParent(es[i % es.size()], Entity{});
    const double tDetach = Ms(a, Clock::now()) * 1e6 / kOps;

    s->tc.Update();
    a = Clock::now();
    s->tc.Update();
    const double tUpdate = Ms(a, Clock::now());

    // Одна смена родителя, гарантированно ломающая порядок: массивы отсортированы
    // по глубине, поэтому корень стоит левее самого глубокого узла. Корень,
    // подвешенный к нему, получает родителя ПРАВЕЕ себя → Update перестраивает всё.
    auto depthOf = [&](Entity e) { int d = 0; for (Entity a = s->tc.Parent(e); !a.IsNull(); a = s->tc.Parent(a)) ++d; return d; };
    Entity deepest = es[0];
    int maxDepth = 0;
    for (Entity e : es) { const int d = depthOf(e); if (d > maxDepth) { maxDepth = d; deepest = e; } }
    bool broken = false;
    for (Entity e : es)
        if (s->tc.Parent(e).IsNull() && s->tc.SetParent(e, deepest)) { broken = true; break; }
    CHECK(broken);
    a = Clock::now();
    s->tc.Update();
    const double tRebuild = Ms(a, Clock::now());

    // массовая гибель: 30% сущностей уходит из пула, их дети становятся корнями
    for (size_t i = 0; i < es.size(); i += 3) s->pool.Destroy(es[i]);
    a = Clock::now();
    s->tc.GarbageCollect(s->pool);
    const double tGC = Ms(a, Clock::now());
    a = Clock::now();
    s->tc.Update();
    const double tAfter = Ms(a, Clock::now());

    std::vector<Entity> alive;
    for (Entity e : es) if (s->pool.IsAlive(e)) alive.push_back(e);
    CHECK(s->tc.Size() == alive.size());
    CHECK(CheckAllAgainstReference(s->tc, alive, s->pool.Capacity()) < 1e-9);

    std::printf("  N = %u (лес)\n", n);
    std::printf("    SetParent(корень)                         %8.0f нс\n", tDetach);
    std::printf("    Update, порядок цел                       %8.2f мс\n", tUpdate);
    std::printf("    Update после ОДНОГО SetParent «вправо»    %8.2f мс   (перестройка всего порядка)\n", tRebuild);
    std::printf("    GarbageCollect 1/3 мёртвых                %8.2f мс\n", tGC);
    std::printf("    Update после GC                           %8.2f мс   (перестройка + сироты стали корнями)\n", tAfter);
}

// Эталон Ogre: calcProjectionParameters + updateFrustumImpl (камера Ogre: y вверх,
// смотрит в -z, глубина [-1,1]) и перевод в [0,1] рендер-системой: row2 = (row2+row3)/2.
struct OgreRef {
    double fovY = 0, aspect = 1, n = 1, f = 100;
    double offX = 0, offY = 0, focal = 1;          // offY — вверх, как у Ogre
    bool ortho = false; double orthoH = 1;
};
static Matrix4x4 OgreProjectionRS(const OgreRef& o)
{
    double left, right, bottom, top;
    if (!o.ortho) {
        const double tanY = std::tan(o.fovY * 0.5), tanX = tanY * o.aspect;
        const double nearFocal = o.n / o.focal, hw = tanX * o.n, hh = tanY * o.n;
        left = -hw + o.offX * nearFocal; right = hw + o.offX * nearFocal;
        bottom = -hh + o.offY * nearFocal; top = hh + o.offY * nearFocal;
    } else {
        const double hh = o.orthoH * 0.5, hw = hh * o.aspect;
        left = -hw; right = hw; bottom = -hh; top = hh;
    }
    const double iw = 1 / (right - left), ih = 1 / (top - bottom), id = 1 / (o.f - o.n);
    Matrix4x4 P; for (int i = 0; i < 16; ++i) P[i] = 0;
    if (!o.ortho) {
        P(0, 0) = 2 * o.n * iw; P(0, 2) = (right + left) * iw;
        P(1, 1) = 2 * o.n * ih; P(1, 2) = (top + bottom) * ih;
        P(2, 2) = -(o.f + o.n) * id; P(2, 3) = -2 * (o.f * o.n) * id; P(3, 2) = -1;
    } else {
        P(0, 0) = 2 * iw; P(0, 3) = -(right + left) * iw;
        P(1, 1) = 2 * ih; P(1, 3) = -(top + bottom) * ih;
        P(2, 2) = -2 * id; P(2, 3) = -(o.f + o.n) * id; P(3, 3) = 1;
    }
    for (int k = 0; k < 4; ++k) P(2, k) = (P(2, k) + P(3, k)) * 0.5;
    return P;
}

// [1.7] Frustum: статические функции, матрица из окна, обратная, глубина,
// плоскости FGED, видимость, стерео, ортография, эталон Ogre.
static void TestFrustum()
{
    std::printf("[1.7] Frustum: матрица, обратная, плоскости, видимость, эталон Ogre\n");
    using V4 = Vector<4, double>;
    using DM = Frustum::DepthMode;
    const Radian fov(3.14159265358979323846 / 3.0);
    const double s = 16.0 / 9.0, n = 1.0, f = 100.0, inf = std::numeric_limits<double>::infinity();
    auto near_ = [](double a, double b, double eps = 1e-12) { return std::fabs(a - b) <= eps; };
    auto maxDiff = [](const Matrix4x4& A, const Matrix4x4& B) {
        double m = 0; for (int i = 0; i < 16; ++i) m = std::fmax(m, std::fabs(A[i] - B[i])); return m;
    };
    auto isIdentity = [](const Matrix4x4& M) {
        for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j)
            if (std::fabs(M(i, j) - (i == j ? 1.0 : 0.0)) > 1e-11) return false;
        return true;
    };
    auto make = [&](ClipY cy, DM dm, double far_) {
        Frustum F(cy, dm);
        F.SetFOVy(fov); F.SetAspectRatio(s); F.SetNearClipDistance(n); F.SetFarClipDistance(far_);
        return F;
    };

    // 1) матрица класса = статическая CreatePerspectiveMatrix = MakeFrustumProjection
    for (ClipY cy : { ClipY::Down, ClipY::Up }) {
        const Frustum F = make(cy, DM::Standard, f);
        CHECK(maxDiff(F.GetProjectionMatrix(), Frustum::CreatePerspectiveMatrix(fov, s, n, f, cy)) < 1e-14);
        CHECK(maxDiff(F.GetProjectionMatrix(), MakeFrustumProjection(fov.ValueRadians(), s, n, f, cy)) < 1e-14);
        CHECK(maxDiff(F.GetInverseProjectionMatrix(), Frustum::MakeFrustumProjectionInverse(fov, s, n, f, cy)) < 1e-14);
    }
    {
        const Frustum R = make(ClipY::Up, DM::Reversed, inf);
        CHECK(maxDiff(R.GetProjectionMatrix(), MakeReversedInfiniteProjection(fov.ValueRadians(), s, n, ClipY::Up)) < 1e-14);
        CHECK(maxDiff(R.GetInverseProjectionMatrix(), Frustum::MakeReversedInfiniteProjectionInverse(fov, s, n)) < 1e-14);
    }

    // 2) все соглашения: обратная, глубина туда-обратно, видимость
    for (ClipY cy : { ClipY::Down, ClipY::Up })
    for (DM dm : { DM::Standard, DM::Reversed })
    for (double far_ : { f, inf }) {
        const Frustum F = make(cy, dm, far_);
        CHECK(isIdentity(F.GetProjectionMatrix() * F.GetInverseProjectionMatrix()));
        const V4 c = F.GetProjectionMatrix() * V4{ 2.0, -1.0, 4.0, 1.0 };
        const Frustum::DepthParams dp = F.GetProjectionParamsAB();
        CHECK(near_(c[2] / c[3], dp.A + dp.B / 4.0) && near_(dp.B / (c[2] / c[3] - dp.A), 4.0, 1e-9));
        CHECK(near_(c[1] / c[3], cy == ClipY::Up ? 0.25 * std::tan(3.14159265358979323846 / 3.0) : -0.25 * std::tan(3.14159265358979323846 / 3.0)));
        CHECK(F.IsVisible(V3{ 0.0, 0.0, 50.0 }) && !F.IsVisible(V3{ 0.0, 0.0, 0.5 }) && !F.IsVisible(V3{ 0.0, 0.0, -5.0 }));
        CHECK(F.IsVisible(V3{ 0.0, 0.0, 0.5 }, 1.0));                     // сфера задевает near
        CHECK(F.IsVisible(V3{ 0.0, 0.0, 200.0 }) == std::isinf(far_));
    }

    // 3) плоскости = таблица 6.1 FGED
    {
        const Frustum F = make(ClipY::Up, DM::Reversed, f);
        const auto& pl = F.GetFrustumPlanes();
        const double g = 1.0 / std::tan(fov.ValueRadians() * 0.5), L = std::sqrt(g * g + s * s), T = std::sqrt(g * g + 1.0);
        CHECK(near_(pl[Frustum::kLeft].GetNormal()[0], g / L) && near_(pl[Frustum::kLeft].GetNormal()[2], s / L));
        CHECK(near_(pl[Frustum::kTop].GetNormal()[1], g / T) && near_(pl[Frustum::kTop].GetNormal()[2], 1.0 / T));
        CHECK(near_(pl[Frustum::kNear].GetDistance(), -n) && near_(pl[Frustum::kFar].GetDistance(), f));
    }

    // 4) эталон Ogre: P_ogre_RS * diag(1,-1,-1,1) = наша (ClipY::Up, Standard)
    Matrix4x4 S; for (int i = 0; i < 16; ++i) S[i] = 0;
    S(0, 0) = 1; S(1, 1) = -1; S(2, 2) = -1; S(3, 3) = 1;
    {
        OgreRef o; o.fovY = fov.ValueRadians(); o.aspect = s; o.n = n; o.f = f;
        Frustum F = make(ClipY::Up, DM::Standard, f);
        CHECK(maxDiff(OgreProjectionRS(o) * S, F.GetProjectionMatrix()) < 1e-12);
        o.offX = -0.032; o.offY = 0.05; o.focal = 2.0;                    // стерео
        F.SetFrustumOffset(Vector<2, double>{ -0.032, -0.05 }); F.SetFocalLength(2.0);   // y вниз
        CHECK(maxDiff(OgreProjectionRS(o) * S, F.GetProjectionMatrix()) < 1e-12);
        const V4 q = F.GetProjectionMatrix() * V4{ -0.032, -0.05, 2.0, 1.0 };   // центр окна на фокусе
        CHECK(near_(q[0] / q[3], 0.0) && near_(q[1] / q[3], 0.0));
        CHECK(isIdentity(F.GetProjectionMatrix() * F.GetInverseProjectionMatrix()));

        OgreRef oo; oo.ortho = true; oo.orthoH = 10.0; oo.aspect = s; oo.n = n; oo.f = f;   // ортография
        Frustum O(ClipY::Up, DM::Standard);
        O.SetProjectionType(Frustum::ProjectionType::Orthographic);
        O.SetOrthoWindowHeight(10.0); O.SetAspectRatio(s); O.SetNearClipDistance(n); O.SetFarClipDistance(f);
        CHECK(maxDiff(OgreProjectionRS(oo) * S, O.GetProjectionMatrix()) < 1e-12);
        CHECK(isIdentity(O.GetProjectionMatrix() * O.GetInverseProjectionMatrix()));
        const Frustum::DepthParams od = O.GetProjectionParamsAB();
        CHECK(near_((od.A + od.B * 40.0 - od.A) / od.B, 40.0) && near_(od.A + od.B * n, 0.0));
    }
}

int main(int argc, char** argv)
{
#ifdef _WIN32
    SetConsoleOutputCP(65001);   // кириллица в UTF-8 в консоли Windows
#endif
    // без аргументов — всё; --tests-only — только часть 1 (для CTest);
    // --stress — части 1 и 3; --bench — части 1 и 2
    const std::string mode = argc > 1 ? argv[1] : "";
    const bool runBench  = mode.empty() || mode == "--bench";
    const bool runStress = mode.empty() || mode == "--stress";

    std::printf("==================== ЧАСТЬ 1. КОРРЕКТНОСТЬ ====================\n");
    TestEntityPool();
    TestNonUniformScale();
    TestRandomAgainstReference();
    TestParentDestroyAndGC();
    TestMatrix3x3Indexing();
    TestProjection();
    TestFrustum();
    std::printf("  проверок: %d, провалено: %d\n\n", g_checks, g_failed);
    if (mode == "--tests-only") return g_failed ? 1 : 0;

#ifndef NDEBUG
    std::printf("  ВНИМАНИЕ: Debug-сборка, время не показательно. Соберите Release.\n\n");
#endif

    if (runBench) {
        std::printf("================= ЧАСТЬ 2. ПРОИЗВОДИТЕЛЬНОСТЬ =================\n");
        const size_t bytesPerEntity = sizeof(uint32_t) + 2 * sizeof(Entity) + sizeof(uint32_t)
                                    + 2 * sizeof(V3) + sizeof(Quaternion) + sizeof(Affine4x3);
        std::printf("  память на трансформ: %zu байт (вход TRS %zu, выход W %zu, связь и индексы %zu)\n\n",
                    bytesPerEntity, 2 * sizeof(V3) + sizeof(Quaternion), sizeof(Affine4x3),
                    bytesPerEntity - 2 * sizeof(V3) - sizeof(Quaternion) - sizeof(Affine4x3));

        std::printf("[2.1/2.2] Update(): первый кадр с перестройкой порядка, затем медиана по кадрам\n");
        std::printf("  %8s  %s %s %s %s %s %s\n", "N", Pad("иерархия", 22).c_str(),
                    Pad(" 1-й, мс", 9).c_str(), Pad(" кадр, мс", 9).c_str(), Pad("   нс/шт", 8).c_str(),
                    Pad("   тыс/мс", 9).c_str(), Pad(" от 60 FPS", 9).c_str());
        for (uint32_t n : { 10'000u, 50'000u, 100'000u, 1'000'000u })
            for (Shape sh : { Shape::Flat, Shape::Forest, Shape::Skeletons })
                BenchUpdate(n, sh, n >= 1'000'000u ? 10 : 40);

        std::printf("\n[2.3] Запись локальных трансформ\n");
        for (uint32_t n : { 10'000u, 100'000u })
            BenchWrites(n);

        std::printf("\n[2.4] Против дерева узлов на указателях (лес, одинаковые данные)\n");
        for (uint32_t n : { 10'000u, 100'000u, 1'000'000u })
            BenchVsPointerTree(n);
        std::printf("\n");
    }

    if (runStress) {
        std::printf("================ ЧАСТЬ 3. СТРЕСС-НАГРУЗКА ИЕРАРХИИ ================\n");
        std::printf("[3.1] Глубокая цепочка (точный ответ — правильный 1000-угольник)\n");
        for (uint32_t d : { 1'000u, 5'000u, 20'000u })
            StressDeepChain(d);

        std::printf("\n[3.2] Широкое дерево и гибель корня\n");
        StressWideTree(100'000);

        std::printf("\n[3.3] Живой мир: рождение, смерть, перевешивание, анимация каждый кадр\n");
        StressChurn({ 10'000, 600, 0.5, 0.5, 0.5, 10 });
        StressChurn({ 50'000, 300, 0.5, 0.5, 0.5, 10 });
        StressChurn({ 50'000, 300, 0.0, 0.0, 0.0, 30 });     // только анимация
        StressChurn({ 50'000, 300, 2.0, 2.0, 5.0, 30 });     // тяжёлый режим

        std::printf("\n[3.3b] Живой мир на сотнях тысяч сущностей\n");
        StressChurn({ 200'000, 120, 0.5, 0.5, 0.5, 10 });
        StressChurn({ 200'000, 120, 0.0, 0.0, 0.0, 10 });    // без смены структуры — только анимация
        StressChurn({ 500'000,  60, 0.5, 0.5, 0.5, 10 });
        StressChurn({ 500'000,  60, 2.0, 2.0, 5.0, 30 });    // тяжёлый режим

        std::printf("\n[3.4] Стоимость операций\n");
        StressOperationCosts(100'000);
        std::printf("\n");
    }

    std::printf("%s: проверок %d, провалено %d\n", g_failed ? "ЕСТЬ ОШИБКИ" : "ВСЁ ПРОШЛО", g_checks, g_failed);
    return g_failed ? 1 : 0;
}
