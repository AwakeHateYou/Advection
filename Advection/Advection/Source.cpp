#include <fstream>
#include <iostream>
#include <complex>
#include <math.h>
#include <omp.h>

#define M_PI 3.14159265358979323846
#define N 20
#define D 1.0
#define D1 1.0
#define D2 -3.0
#define W 1.0
#define C 1.0
#define M 1.0
#define P 1.0
using namespace std;


complex<double> matrixInitialization(double x, double y) {
	return cos(P * M_PI * x) + cos(M * M_PI * y);
}
complex<double> GinzburgLandauFunction(
	complex<double> cur,
	complex<double> east, 
	complex<double> south, 
	complex<double> west,
	complex<double> north,
	double DX, double DY, double DT) {
	return cur + DT * (-D * complex<double>(1, W) * (((north - 2.0 * cur + south) / (DX * DX) + (east - 2.0 * cur + west) / (DY * DY))) + 
		(D1 * ((west - cur) / DX)) + (D2 * ((north - cur) / DY))
		+ cur + complex<double>(-1, C) * abs(pow(cur, 2.0))* cur);
}
complex<double> funcWithoutY(
	complex<double> cur,
	complex<double> east,
	complex<double> south,
	complex<double> west,
	complex<double> north,
	double DX, double DY, double DT) {
	return cur + DT * (-D * complex<double>(1, W) * ((north - 2.0 * cur + south) / (DX * DX)) +
		(D1 * ((west - cur) / DX)) + D2
		+ cur + complex<double>(-1, C) * abs(pow(cur, 2.0))* cur);
}
void printMatrix(complex<double> *mass) {
	ofstream fout;
	fout.open("output.txt");
	//fout.precision(3);
	for (int i = 0; i < N*N; i++) {
		if ((i % N == 0) && (i != 0))
			fout << std::endl;
		fout << '{' << real(mass[i]) << ',' << imag(mass[i]) << "} ";
	}
	fout.close();
}
void printGNU(complex<double> *mass, double dx, double dy, char* filename) {
	ofstream fout;
	fout.open(filename);
	//fout.precision(3);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			fout << dx*i << ' ' << dy*j << ' ' << real(mass[i*N + j]) << endl;
			//fout << dx*i << ' ' << dy*j << " {" << real(mass[i*N + j]) << ',' << imag(mass[i*N + j]) << '}' << endl;
		}
		fout << endl;
	}
	fout.close();
}

int main(int argc, char** argv) {
	int i, j, K;
	complex<double> *currentMatrix, *nextMatrix, *swap;
	double t1, t2, DT, DX, DY;
	currentMatrix = new complex<double>[N*N];
	nextMatrix = new complex<double>[N*N];

	DX = 1.0 / N;
	DT = 0.025 * DX * DX;
	DY = DX;
	K = 1 / DT;

	cout << "K = " << K << endl;
	cout << "DT = " << DT << endl;
	cout << "DX = " << DX << endl;
	cout << "DY = " << DY << endl;

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++) {
			currentMatrix[i*N + j] = complex< double >(0.0, 0);
			nextMatrix[i*N + j] = complex< double >(0.0, 0);
		}
	for (i = 1; i < N-1; i++)
		for (j = 1; j < N-1; j++) {
			currentMatrix[i*N + j] = matrixInitialization(i*DX, j*DY);
			nextMatrix[i*N + j] = matrixInitialization(i*DX, j*DY);
		}

	t1 = omp_get_wtime();

	for (int k = 0; k < K; k++) {

		for (i = 1; i < N - 1; i++)
#pragma omp parallel for
			for (j = 1; j < N - 1; j++) {
				nextMatrix[i*N + j] = GinzburgLandauFunction(
					currentMatrix[i*N + j], 
					currentMatrix[i*N + j - 1], 
					currentMatrix[i*N + N],
					currentMatrix[i*N + j + 1],
					currentMatrix[i*N + j - N],
					DX, DY, DT);
			}
		swap = currentMatrix;
		currentMatrix = nextMatrix;
		nextMatrix = swap;
	}
	t2 = omp_get_wtime();
	cout << "Time: " << t2 - t1 << endl;
	printMatrix(currentMatrix);
	printGNU(currentMatrix, DX, DY, "outputGNU.txt");

	return 0;
}
