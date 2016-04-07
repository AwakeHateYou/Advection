#include <fstream>
#include <iostream>
#include <complex>
#include <math.h>
#include <omp.h>

#define N 10
#define D 1.0
#define D1 1.0
#define D2 1.0
#define W 1.0
#define C 1.0
using namespace std;

complex<double> GinzburgLandauFunction(complex<double> cur, complex<double> nextI, complex<double> nextJ, double DX, double DY, double DT) {
	return (-D * complex<double>(1, W) * cur + (D1 * (nextJ - cur) / DX) + (D2 * (nextI - cur) / DY) 
		+ cur + complex<double>(-1, C) * abs(pow(cur, 2.0))* cur) * (DT / (- cur * (1.0 + D * complex<double>(1, W) * DT)));
}
void print(complex<double> **mass) {
	ofstream fout;
	fout.open("output.txt");
	//fout.precision(3);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			fout << mass[i][j] << ' ';
		}
		fout << endl;
	}
	fout.close();
}
void printGNU(complex<double> **mass, double dx, char* filename) {
	ofstream fout;
	fout.open(filename);
	//fout.precision(3);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			fout << dx*j << ' ' << dx*i << ' ' << mass[i][j] << endl; //i - Oy; j - Ox
		}
		fout << endl;
	}
	fout.close();
}

int main(int argc, char** argv) {
	int i, j, K;
	complex<double> **currentMatrix, **nextMatrix, **swap;
	double t1, t2, DT, DX, DY;
	currentMatrix = new complex<double>*[N];
	nextMatrix = new complex<double>*[N];

	DX = 1.0 / N;
	DT = 0.25 * DX * DX;
	DY = DX;
	K = 1 / DT;
	//    cout << "K = " << K << endl;
	//    cout << "DT = " << DT << endl;
	//    cout << "DX = " << DX << endl;
	//    cout << "DY = " << DY << endl;
	for (i = 0; i < N; i++) {
		currentMatrix[i] = new complex<double>[N];
		nextMatrix[i] = new complex<double>[N];
	}

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++) {
			currentMatrix[i][j] = complex< double >(0.0, 0);
			nextMatrix[i][j] = complex< double >(0.0, 0);
		}
	//for (i = 0; i < N; i++) {
	//	M[i][0] = 100.0;
	//	V[i][0] = 100.0;
	//}

	t1 = omp_get_wtime();

	for (int k = 0; k < 2; k++) {
		//#pragma omp parallel for
		for (i = 1; i < N - 1; i++)
			for (j = 1; j < N - 1; j++) {
				nextMatrix[i][j] = GinzburgLandauFunction(currentMatrix[i][j], currentMatrix[i + 1][j], currentMatrix[i][j + 1], DX, DY, DT);
			}
		swap = currentMatrix;
		currentMatrix = nextMatrix;
		nextMatrix = swap;
	}
	t2 = omp_get_wtime();
	cout << "Time: " << t2 - t1 << endl;
	print(nextMatrix);

	return 0;
}
