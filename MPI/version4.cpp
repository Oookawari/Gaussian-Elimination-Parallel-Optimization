#include <mpi.h>
#include<iostream>
#include<iomanip>
#include<random>
#include<time.h>
#include<stdio.h>
#include<semaphore.h>
using namespace std;
const int ORDER = 3000;
float matrix[ORDER][ORDER + 1];
float temp[ORDER][ORDER + 1];
int order = 100;
int SAMPLE = 0;
int NUM_THREADS = 8;

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
		int nextline = rank; // 该进程下一次要除法的行
		for (int k = 0; k < order; k++) {
			//本进程需要进行除法
			if (k == nextline) {
				//进行除法
				for (int j = k + 1; j <= order; j++)
					matrix[k][j] = matrix[k][j] / matrix[k][k];
				matrix[k][k] = 1.0;

				//除法结束，广播给下一个进程
				int next_rank = rank + 1;
				if (next_rank == size) { next_rank = 0; }
				if (k + 1 != order) {
					//cout << "sendline: " << k << " from: " << rank << " to: " << next_rank << endl;
					MPI_Send(&matrix[k][0], order + 1, MPI_FLOAT, next_rank, k, MPI_COMM_WORLD);
					//cout << "sendline: " << k << " from: " << rank << " to: " << next_rank << " completed" << endl;
				}
				nextline += size;
				//进行减法消去
				for (int i = nextline; i < order; i += size) {
					for (int j = k + 1; j <= order; j++)
						matrix[i][j] -= matrix[i][k] * matrix[k][j];
					matrix[i][k] = 0;
				}
			}
			else {
				//一定接收，但是不一定发送
				int last_rank = rank - 1;
				int next_rank = rank + 1;
				if (last_rank == -1) { last_rank = size - 1; }
				if (next_rank == size) { next_rank = 0; }
				if (k + 1 != order) {
					//cout << "recvline: " << k << " rank: " << rank << " from: " << last_rank << endl;
					MPI_Recv(&matrix[k][0], order + 1, MPI_FLOAT, last_rank, k, MPI_COMM_WORLD, &status);
					//cout << "recvline: " << k << " rank: " << rank << " from: " << last_rank << " completed" << endl;
				}
				if (k + 1 != order && next_rank % size != k) {
					//cout << "sendline: " << k << " from: " << rank << " to: " << next_rank << endl;
					MPI_Send(&matrix[k][0], order + 1, MPI_FLOAT, next_rank, k, MPI_COMM_WORLD);
					//cout << "sendline: " << k << " from: " << rank << " to: " << next_rank << " completed" << endl;
				}
				for (int i = nextline; i < order; i += size) {
					for (int j = k + 1; j <= order; j++)
						matrix[i][j] -= matrix[i][k] * matrix[k][j];
					matrix[i][k] = 0;
				}
				//进行减法消去
			}
		}
		//0号进程统一接收所有的处理结果
		if (rank == 0) {
			for (int i = 0; i < order; i++) {
				if (i % size != 0) {
					MPI_Recv(&matrix[i][0], order + 1, MPI_FLOAT, MPI_ANY_SOURCE, i, MPI_COMM_WORLD, &status);
				}
			}
		}
		else {
			for (int i = rank; i < order; i += size) {
				MPI_Send(&matrix[i][0], order + 1, MPI_FLOAT, 0, i, MPI_COMM_WORLD);
			}
		}
		ets = MPI_Wtime();
		if (rank == 0 && time == 0) {
			cout << "MPI order: " << order << "   times:" << times << endl;
			cout << fixed << setprecision(2) << (ets - sts) * 1000 << endl;
			cout << endl;
			/*
			for (int i = 0; i < order; i++) {
				for (int j = 0; j <= order; j++) {
					cout << fixed << setprecision(2) << matrix[i][j] << " ";
				}
				cout << endl;
			}*/
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
	int sample = 8;
	int samples[8][2] = { {100,1},{300,1},{500,1},{700,1},{1000,1},{1500,1},{2000,1},{3000,1} };
	order = 3000;

	MPI_Init(NULL, NULL);
	for (int i = 0; i < sample; i++) {
		order = samples[i][0];
		MPI_version1(samples[i][0], samples[i][1]);
	}
	MPI_Finalize();
	return 0;
}
