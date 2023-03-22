#include <mpi.h>
#include<iostream>
#include<iomanip>
#include<random>
#include<time.h>
#include<stdio.h>
#include<semaphore.h>
#include <immintrin.h>

using namespace std;
const int ORDER = 3000;
float matrix[ORDER][ORDER + 1];
float temp[ORDER][ORDER + 1];
int order = 100;
int SAMPLE = 0;
int NUM_THREADS = 8;

void tradition(int order, int times) {
	struct timespec sts, ets;
	double time_sum = 0.0;
	for (int time = 0; time < times; time++) {

		//初始化/恢复矩阵 不计入时间
		for (int i = 0; i < order; i++)
			for (int j = 0; j <= order; j++)
				matrix[i][j] = temp[i][j];

		timespec_get(&sts, TIME_UTC);
		for (int k = 0; k < order; k++) {
			for (int j = k + 1; j <= order; j++)
				matrix[k][j] = matrix[k][j] / matrix[k][k];
			matrix[k][k] = 1.0;
			for (int i = k + 1; i < order; i++) {
				for (int j = k + 1; j <= order; j++) {
					matrix[i][j] -= matrix[i][k] * matrix[k][j];
				}
				matrix[i][k] = 0.0;
			}
		}
		timespec_get(&ets, TIME_UTC);
		time_sum += (ets.tv_nsec - sts.tv_nsec) * 0.000001 + (ets.tv_sec - sts.tv_sec) * 1000;
	}

	for (int i = 0; i < order; i++) {
		for (int j = 0; j <= order; j++) {
			cout << fixed << setprecision(2) << matrix[i][j] << " ";
		}
		cout << endl;
	}/**/
	cout << "Traditon order: " << order << "   times:" << times << endl;
	cout << fixed << setprecision(2) << time_sum / times << endl;
	cout << endl;
	return;
}

void MPI_version1(int order, int times) {
	double sts, ets;
	double time_sum = 0.0;


	for (int time = 0; time < times; time++) {

		//初始化/恢复矩阵 不计入时间
		for (int i = 0; i < order; i++)
			for (int j = 0; j <= order; j++)
				matrix[i][j] = temp[i][j];


		int rank;
		int size;
		sts = MPI_Wtime();
		MPI_Status status;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		//
		for (int k = 0; k < order; k++) {
			if (rank == 0) {
				__m256 vector_k_k;
				vector_k_k = _mm256_set1_ps(matrix[k][k]);
				int j;
				for (j = k + 1; j + 7 <= ORDER; j += 8) {
					__m256 vector_k_j = _mm256_loadu_ps(&matrix[k][j]);
					vector_k_j = _mm256_div_ps(vector_k_j, vector_k_k);
					_mm256_store_ps(&matrix[k][j], vector_k_j);
				}
				for (; j <= ORDER; j++) {
					matrix[k][j] = matrix[k][j] / matrix[k][k];
				}
				matrix[k][k] = 1.0;
				for (int j = 1; j < size; j++)
					MPI_Send(&matrix[k][0], order + 1, MPI_FLOAT, j, k + 1, MPI_COMM_WORLD);
			}
			else {
				MPI_Recv(&matrix[k][0], order + 1, MPI_FLOAT, 0, k + 1, MPI_COMM_WORLD, &status);
			}
			for (int i = rank; i < order; i += size) {
				if (i < k + 1) continue;
				else if (i == k + 1) {
					int j;
					__m256 vector_i_k = _mm256_set1_ps(matrix[i][k]);
					for (j = k + 1; j + 7 <= ORDER; j += 8) {
						__m256 vector_k_j = _mm256_loadu_ps(&matrix[k][j]);
						__m256 vector_i_j = _mm256_loadu_ps(&matrix[i][j]);
						__m256 vector_temp = _mm256_mul_ps(vector_i_k, vector_k_j);
						vector_i_j = _mm256_sub_ps(vector_i_j, vector_temp);
						_mm256_store_ps(&matrix[i][j], vector_i_j);
					}
					for (; j <= ORDER; j++) {
						matrix[i][j] -= matrix[i][k] * matrix[k][j];
					}
					matrix[i][k] = 0.0;
					if (rank != 0) //将k+1行的结果更新给0号
						MPI_Send(&matrix[i][0], order + 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
				}
				else {
					for (int j = k + 1; j <= order; j++)
						matrix[i][j] -= matrix[i][k] * matrix[k][j];
					matrix[i][k] = 0;
				}
			}
			if (rank == 0 && ((k + 1) % size != 0)) { //如果k+1行不是由0号处理
				cout << "waiting" << endl;
				MPI_Recv(&matrix[k + 1][0], order + 1, MPI_FLOAT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status); cout << "waitingend" << endl;
			}
		}MPI_Barrier(MPI_COMM_WORLD);
		ets = MPI_Wtime();
		if (rank == 0 && time == 0) {
			cout << "MPI order: " << order << "   times:" << times << endl;
			cout << fixed << setprecision(2) << (ets - sts) * 1000 << endl;
			cout << endl;
			for (int i = 0; i < order; i++) {
				for (int j = 0; j <= order; j++) {
					cout << fixed << setprecision(2) << matrix[i][j] << " ";
				}
				cout << endl;
			}/**/
		}

	}
	/*
	*/
	return;
}

int main() {
	unsigned seed = 47;
	srand(seed);
	for (int i = 0; i < ORDER; i++)
		for (int j = 0; j <= ORDER; j++)
			temp[i][j] = rand() % 10 + rand() % 10000 * 0.0001;
	int sample = 1;
	int samples[8][2] = { {10,1},{300,1},{500,1},{700,1},{1000,1},{1500,1},{2000,1},{3000,1} };
	order = 3000;

	for (int i = 0; i < sample; i++) {
		order = samples[i][0];
		tradition(samples[i][0], samples[i][1]);
	}/**/
	MPI_Init(NULL, NULL);
	for (int i = 0; i < sample; i++) {
		order = samples[i][0];
		MPI_version1(samples[i][0], samples[i][1]);
	}
	MPI_Finalize();
	return 0;
}
