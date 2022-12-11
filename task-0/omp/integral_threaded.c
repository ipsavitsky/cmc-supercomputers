#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

typedef double (*math_func)(double);

double calculate_integral(double from, double to, int n, math_func func) {
    int i;
    double h = (to - from) / (double)n;
    double sum = 0.0, x;
#pragma omp parallel for reduction(+ : sum) private(x)
    for (i = 1; i <= n; ++i) {
        x = h * ((double)i - 0.5);
        sum += func(x);
    }
    return sum * h;
}

int main(int argc, char *argv[]) {
    double start, stop;
    start = omp_get_wtime();
    int n = strtol(argv[1], NULL, 10);
    double res = calculate_integral(0, 10000, n, sqrt);
    stop = omp_get_wtime();
    printf("%d,%lf,%d,%lf\n", omp_get_max_threads(), res, n, stop-start);
    return 0;
}
