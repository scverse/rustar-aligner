#include <random>
#include <cstdio>
int main() {
    std::mt19937 g(19760110u);
    printf("mt19937_first10:");
    for (int i = 0; i < 10; i++) printf(" %u", g());
    printf("\n");

    std::mt19937 g2(19760110u);
    printf("generate_canonical53_first5:");
    for (int i = 0; i < 5; i++)
        printf(" %.17g", std::generate_canonical<double, 53>(g2));
    printf("\n");

    std::mt19937 g3(19760110u);
    std::discrete_distribution<int> d({1.0, 2.0, 3.0, 4.0});
    printf("discrete_1234_first20:");
    for (int i = 0; i < 20; i++) printf(" %d", d(g3));
    printf("\n");
    return 0;
}
