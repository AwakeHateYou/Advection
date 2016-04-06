#include <fstream>
#include <iostream>
#include <complex>
#include <omp.h>

#define N 100

using namespace std;

complex<double> Ginzburg–LandauFunction(complex<double> cur, complex<double> next) {
	return 1;
}
void print(complex<double> **mass, double dx) {
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
			fout << dx*i << ' ' << dx*j << ' ' << mass[i][j] << endl;
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

	for (int k = 0; k < K; k++) {
		//#pragma omp parallel for
		for (i = 1; i < N - 1; i++)
			for (j = 1; j < N - 1; j++) {
				nextMatrix[i][j] = (DT, DX, DY, M[i][j], M[i - 1][j], M[i][j - 1], M[i][j + 1], M[i + 1][j]);
			}
		SW = M;
		M = V;
		V = SW;
	}
	t2 = omp_get_wtime();
	cout << "Time: " << t2 - t1 << endl;

	FILE *f = fopen("output.txt", "w");
	if (f == NULL)
		return 1;
	else {
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++)
				fprintf(f, "%d %d %3.2f\n", i, j, V[i][j]);
			fprintf(f, "\n");
		}
	}

	fclose(f);

	return 0;
}
