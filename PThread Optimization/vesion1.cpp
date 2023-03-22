#include<pthread.h>
#include<iostream>
#include<iomanip>
#include<random>
#include<time.h>
#include<stdio.h>
#include<semaphore.h>
#include <immintrin.h>
#pragma comment(lib, "pthreadVC2.lib")

using namespace std;
const int ORDER = 3000;
float matrix[ORDER][ORDER + 1];
float temp[ORDER][ORDER + 1];
const int NUM_THREADS = 7;
int order = 100;

//传统算法
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
	//cout << "Traditon order: " << order << "   times:" << times << endl;
	cout << fixed << setprecision(2) << time_sum / times << endl;
	//cout << endl;
	return;
}

//传统算法
void AVX(int order, int times) {
	struct timespec sts, ets;
	double time_sum = 0.0;
	for (int time = 0; time < times; time++) {

		//初始化/恢复矩阵 不计入时间
		for (int i = 0; i < order; i++)
			for (int j = 0; j <= order; j++)
				matrix[i][j] = temp[i][j];

		timespec_get(&sts, TIME_UTC);

		for (int k = 0; k < order; k++) {
			__m256 vector_k_k;
			vector_k_k = _mm256_set1_ps(matrix[k][k]);
			int j;
			for (j = k + 1; j + 7 <= order; j += 8) {
				__m256 vector_k_j = _mm256_loadu_ps(&matrix[k][j]);
				vector_k_j = _mm256_div_ps(vector_k_j, vector_k_k);
				_mm256_store_ps(&matrix[k][j], vector_k_j);
			}
			for (; j <= order; j++) {
				matrix[k][j] = matrix[k][j] / matrix[k][k];
			}
			matrix[k][k] = 1.0;
			for (int i = k + 1; i < order; i++) {
				__m256 vector_i_k = _mm256_set1_ps(matrix[i][k]);
				for (j = k + 1; j + 7 <= order; j += 8) {
					__m256 vector_k_j = _mm256_loadu_ps(&matrix[k][j]);
					__m256 vector_i_j = _mm256_loadu_ps(&matrix[i][j]);
					__m256 vector_temp = _mm256_mul_ps(vector_i_k, vector_k_j);
					vector_i_j = _mm256_sub_ps(vector_i_j, vector_temp);
					_mm256_store_ps(&matrix[i][j], vector_i_j);
				}
				for (; j <= order; j++) {
					matrix[i][j] -= matrix[i][k] * matrix[k][j];
				}
				matrix[i][k] = 0.0;
			}
		}
		timespec_get(&ets, TIME_UTC);
		time_sum += (ets.tv_nsec - sts.tv_nsec) * 0.000001 + (ets.tv_sec - sts.tv_sec) * 1000;
	}
	cout << "AVX order: " << order << "   times:" << times << endl;
	cout << fixed << setprecision(4) << time_sum << endl;
	cout << "averaged_time:  " << fixed << setprecision(4) << time_sum / times << endl;
	cout << endl;
	return;
}

typedef struct {
	int t_id;   //线程id
}threadparam_t;

sem_t sem_div;  //信号量定义
sem_t sem_elimination;
pthread_barrier_t barrier;

void* threadfunc(void* param) {
	threadparam_t* p = (threadparam_t*)param;
	int t_id = p->t_id;
	for (int k = 0; k < order; k++) {
		if (t_id == 0) {
			for (int j = k + 1; j <= order; j++)
				matrix[k][j] = matrix[k][j] / matrix[k][k];
			matrix[k][k] = 1.0;
			for (int i = 1; i < NUM_THREADS; i++)
				sem_post(&sem_div);
		}
		else sem_wait(&sem_div);
		if (t_id == 0) {

		}
		for (int i = k + 1 + t_id; i < order; i += NUM_THREADS) {
			for (int j = k + 1; j <= order; j++)
				matrix[i][j] -= matrix[i][k] * matrix[k][j];
			matrix[i][k] = 0;
		}
		pthread_barrier_wait(&barrier);
	}
	pthread_exit(NULL);
	return NULL;
}

void pthread(int order, int times) {
	struct timespec sts, ets;
	double time_sum = 0.0;
	for (int time = 0; time < times; time++) {
		//初始化/恢复矩阵 不计入时间
		for (int i = 0; i < order; i++)
			for (int j = 0; j <= order; j++)
				matrix[i][j] = temp[i][j];
		timespec_get(&sts, TIME_UTC);
		sem_init(&sem_div, 0, 0);
		sem_init(&sem_elimination, 0, 0);
		pthread_barrier_init(&barrier, NULL, NUM_THREADS);
		//创建线程
		pthread_t handles[NUM_THREADS];
		threadparam_t param[NUM_THREADS];
		for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
			param[t_id].t_id = t_id;
			pthread_create(&handles[t_id], NULL, threadfunc, (void*)&param[t_id]);
		}
		for (int i = 0; i < NUM_THREADS; i++)
			pthread_join(handles[i], NULL);
		sem_destroy(&sem_div);
		sem_destroy(&sem_elimination);
		pthread_barrier_destroy(&barrier);
		timespec_get(&ets, TIME_UTC);
		time_sum += (ets.tv_nsec - sts.tv_nsec) * 0.000001 + (ets.tv_sec - sts.tv_sec) * 1000;
	}
	//cout << "PThread order: " << order << "   times:" << times << endl;
	//cout << fixed << setprecision(4) << time_sum << endl;
	cout << fixed << setprecision(2) << time_sum / times << endl;
	//cout << endl;
}

int main() {
	unsigned seed = 47;
	srand(seed);
	for (int i = 0; i < ORDER; i++)
		for (int j = 0; j <= ORDER; j++)
			temp[i][j] = rand() % 10 + rand() % 10000 * 0.0001;
	int sample = 8;
	int samples[8][2] = { {100,100},{300,10},{500,5},{700,2},{1000,1},{1500,1},{2000,1},{5000,1} };
	for (int i = 0; i < sample; i++)
		tradition(samples[i][0], samples[i][1]);
	for (int i = 0; i < sample; i++)
		AVX(samples[i][0], samples[i][1]);
	for (int i = 0; i < sample; i++) {
		order = samples[i][0];
		pthread(samples[i][0], samples[i][1]);
	}
	return 0;
}
