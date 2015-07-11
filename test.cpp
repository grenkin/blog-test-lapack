#include <lapacke.h>
#include <iostream>
#include <math.h>
#include <time.h>

const double L = 100;
const double sigma = 0.2;
const double K = 0.01;
const double beta = 1, alpha = 0;

const int M = 2000;
const double h = 2*L/M;

double T[M+1], Told[M+1];

const lapack_int Neq = M + 1;
double d[Neq], dl[Neq-1], du[Neq-1], b[Neq];

double xi(int i)
{
    return -L + i*h;
}

double w(double beta, double T)
{
    return (T - sigma) * (1 - T) * exp(-beta/T);
}

double v(double alpha, double x)
{
  return alpha*(x + L);
}

int main() {
	clock_t start, finish;
	double duration;
	start = clock();

    for (int i = 0; i <= M; ++i)
        Told[i] = 1 + (sigma - 1) * (xi(i) + L) / (2*L);

    d[0] = 1; du[0] = 0;
    d[M] = 1; dl[M-1] = 0;
    b[0] = 1; b[M] = sigma;
    for (int i = 1; i < M; ++i) {
        dl[i-1] = 1/(h*h) - v(alpha, xi(i))/(2*h);
        du[i] = 1/(h*h) + v(alpha, xi(i))/(2*h);
        d[i] = -2/(h*h);
        b[i] = -K*w(beta, Told[i]);
    }

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    std::cout << "matrix construction time: " << duration << " sec\n";

	start = clock();

    lapack_int info = LAPACKE_dgtsv(LAPACK_COL_MAJOR, Neq, 1,
        dl, d, du, b, Neq);
    if (info != 0) {
        std::cout << "Error when solving a linear system, info = " << info << std::endl;
        exit(1);
    }
    for (int i = 0; i <= M; ++i)
        T[i] = b[i];

    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    std::cout << "system solving time: " << duration << " sec\n";

    return 0;
}
