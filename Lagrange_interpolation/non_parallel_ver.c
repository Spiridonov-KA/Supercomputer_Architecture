#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

double dabs(double x) {
	return (x < 0) ? -x : x;
}

double *arange(double start, double stop, double step, int *n) {
	*n = (int) ceil( (stop - start) / step);
	double *res = (double *) malloc(sizeof(double) * (*n));
	res[0] = start;
	for (int i = 1; i < *n; ++i) {
		res[i] = res[i-1] + step;
	}

	return res;
}

double *linspace(double start, double stop, int *n) {
	double *res = (double *) malloc(sizeof(double) * (*n));
	double step = (stop - start) / (*n - 1);
	res[0] = start;
	for (int i = 1; i < *n; ++i) {
		res[i] = res[i-1] + step;
	}

	return res;
}

double lagrange_polinome(double *knownx, double *knowny, int n, double x) {
	double res = 0;
	double polinome;
	for (int i = 0; i < n; ++i) {
		polinome = 1;
		for (int j = 0; j < n; ++j) {
			if (i == j)
				continue;
			polinome *= (x - knownx[j]) / (knownx[i] - knownx[j]);
		}
		res += knowny[i] * polinome;
	}

	return res;
}

double *lagrange_interpolation(double *knownx, double *knowny, int n, double *x, int size_x) {
	double *res = (double *) malloc(sizeof(double) * size_x);

	for (int i = 0; i < size_x; ++i) {
		res[i] = lagrange_polinome(knownx, knowny, n, x[i]);
	}

	return res;
}

void read_input(int *pn, double **knownx, double **knowny, double **x, int *size_x) {
	int way;
	double start, stop, step;
	printf("Input number of known points:\n");
	scanf("%d", pn);
	*knownx = (double *) malloc(sizeof(double) * (*pn));
	*knowny = (double *) malloc(sizeof(double) * (*pn));
	printf("Input known x:\n");
	for (int i = 0; i < *pn; ++i) {
		scanf("%lf", &((*knownx)[i]));
	}

	printf("Input known y:\n");
	for (int i = 0; i < *pn; ++i) {
		scanf("%lf", &((*knowny)[i]));
	}
	printf("Way to create array in which you want to calculate value of function:\n"
		"1. arange\n"
		"2. linspace\n");

	scanf("%d", &way);
	if (way == 1) {
		printf("input start, stop and step:\n");
		scanf("%lf%lf%lf", &start, &stop, &step);
		*x = arange(start, stop, step, size_x);
	}
	else if (way == 2) {
		printf("input start, stop and number of points:\n");
		scanf("%lf%lf%d", &start, &stop, size_x);
		*x = linspace(start, stop, size_x);
	}
	else {
		printf("Error: unknow way");
		abort();
	}
}

int main(int argc, char *argv[]) {
	int n = 0, size_x = 0;
	double *knownx;
	double *knowny;
	double *x;
	double *res;
	double *total_res;
	clock_t start_time, end_time;

	read_input(&n, &knownx, &knowny, &x, &size_x);

		#if DEBUG
		printf("n = %d\n", n);
		printf("My input\n");
		for (int i = 0; i < n; ++i) {
			printf("%lf ", knownx[i]);
		}
		printf("\nsize_x = %d\n", size_x);
		for (int i = 0; i < size_x; ++i) {
			printf("%lf ", x[i]);
		}
		printf("\nMy input\n");
		#endif

	start_time = clock();

	total_res = lagrange_interpolation(knownx, knowny, n, x, size_x);

	end_time = clock();

	printf("Time: %.10lf sec\n", ((double)end_time - start_time) / CLOCKS_PER_SEC);

	#if DEBUG
	printf("TOTAL_RES:\n");
	for (int i = 0; i < size_x; ++i) {
		printf("%lf ", total_res[i]);
	}
	printf("\n");
	#endif

	
	free(knownx);
	free(knowny);
	free(x);

	return 0;
}
