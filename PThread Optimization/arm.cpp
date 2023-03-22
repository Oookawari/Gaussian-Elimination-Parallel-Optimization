#include<pthread.h>
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
const int NUM_THREADS = 7;
int order = 100;

typedef struct {
	int t_id;   //线程id
}threadparam_t;

sem_t sem_div;  //信号量定义
sem_t sem_elimination;
sem_t sem_divison[NUM_THREADS - 1];
sem_t sem_leader;
sem_t sem_worker[NUM_THREADS - 1];
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

void* threadfunc2(void* param) {
	threadparam_t* p = (threadparam_t*)param;
	int t_id = p->t_id;
	for (int k = 0; k < order; k++) {
		//for (int i = k + 1 + t_id; i <= order; i += NUM_THREADS) {
			//matrix[k][i] = matrix[k][i] / matrix[k][k];
		//}
		if (t_id == 0) {
			for (int j = k + 1; j <= order; j++)
				matrix[k][j] = matrix[k][j] / matrix[k][k];
			matrix[k][k] = 1.0;
		}
		pthread_barrier_wait(&barrier);
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

void* threadfunc3(void* param) {
	threadparam_t* p = (threadparam_t*)param;
	int t_id = p->t_id;
	for (int k = 0; k < order; k++) {
		if (t_id == 0) {
			for (int j = k + 1; j <= order; j++)
				matrix[k][j] = matrix[k][j] / matrix[k][k];
			matrix[k][k] = 1.0;
			for (int i = 0; i < NUM_THREADS - 1; i++)
				sem_post(&sem_divison[i]);
		}
		else sem_wait(&sem_divison[t_id - 1]);
		for (int i = k + 1 + t_id; i < order; i += NUM_THREADS) {
			for (int j = k + 1; j <= order; j++)
				matrix[i][j] -= matrix[i][k] * matrix[k][j];
			matrix[i][k] = 0;
		}
		if (t_id == 0) {
			for (int i = 1; i < NUM_THREADS; i++)
				sem_wait(&sem_leader);
			for (int i = 0; i < NUM_THREADS - 1; i++) {
				sem_post(&sem_worker[i]); // 通知其它 worker 进入下一轮
			}
		}
		else {
			sem_post(&sem_leader);// 通知 leader, 已完成消去任务
			sem_wait(&sem_worker[t_id - 1]);
		}
	}
	pthread_exit(NULL);
	return NULL;
}

void* threadfunc4(void* param) {
	threadparam_t* p = (threadparam_t*)param;
	int t_id = p->t_id;
	for (int k = 0; k < order; k++) {

		for (int i = k + 1 + t_id; i <= order; i += NUM_THREADS) {
			matrix[k][i] = matrix[k][i] / matrix[k][k];
		}
		if (t_id == 0) {
			for (int i = 0; i < NUM_THREADS - 1; i++)
				sem_post(&sem_divison[i]);
		}
		else {
			sem_wait(&sem_divison[t_id - 1]);
		}
		if (t_id == 0) {
			matrix[k][k] = 1.0;
		}
		for (int i = k + 1 + t_id; i < order; i += NUM_THREADS) {
			for (int j = k + 1; j <= order; j++)
				matrix[i][j] -= matrix[i][k] * matrix[k][j];
			matrix[i][k] = 0;
		}
		if (t_id == 0) {
			for (int i = 1; i < NUM_THREADS; i++)
				sem_wait(&sem_leader);
			for (int i = 0; i < NUM_THREADS - 1; i++) {
				sem_post(&sem_worker[i]); // 通知其它 worker 进入下一轮
			}
		}
		else {
			sem_post(&sem_leader);// 通知 leader, 已完成消去任务
			sem_wait(&sem_worker[t_id - 1]);
		}
	}
	pthread_exit(NULL);
	return NULL;
}

//信号量同步1 barrier同步2
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
	/*
	for (int i = 0; i < order; i++) {
		for (int j = 0; j <= order; j++) {
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
	*/
	cout << "PThread order: " << order << "   times:" << times << endl;
	cout << fixed << setprecision(2) << time_sum / times << endl;
	cout << endl;
}

//barrier同步1 barrier同步2
void pthread2(int order, int times) {
	struct timespec sts, ets;
	double time_sum = 0.0;
	for (int time = 0; time < times; time++) {
		//初始化/恢复矩阵 不计入时间
		for (int i = 0; i < order; i++)
			for (int j = 0; j <= order; j++)
				matrix[i][j] = temp[i][j];
		timespec_get(&sts, TIME_UTC);
		pthread_barrier_init(&barrier, NULL, NUM_THREADS);
		//创建线程
		pthread_t handles[NUM_THREADS];
		threadparam_t param[NUM_THREADS];
		for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
			param[t_id].t_id = t_id;
			pthread_create(&handles[t_id], NULL, threadfunc2, (void*)&param[t_id]);
		}
		for (int i = 0; i < NUM_THREADS; i++)
			pthread_join(handles[i], NULL);
		pthread_barrier_destroy(&barrier);
		timespec_get(&ets, TIME_UTC);
		time_sum += (ets.tv_nsec - sts.tv_nsec) * 0.000001 + (ets.tv_sec - sts.tv_sec) * 1000;
	}
	/*
	for (int i = 0; i < order; i++) {
		for (int j = 0; j <= order; j++) {
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
	*/
	cout << "PThread2 order: " << order << "   times:" << times << endl;
	cout << fixed << setprecision(2) << time_sum / times << endl;
	cout << endl;
}

//信号量同步1 信号量同步2
void pthread3(int order, int times) {
	struct timespec sts, ets;
	double time_sum = 0.0;
	for (int time = 0; time < times; time++) {
		//初始化/恢复矩阵 不计入时间
		for (int i = 0; i < order; i++)
			for (int j = 0; j <= order; j++)
				matrix[i][j] = temp[i][j];
		timespec_get(&sts, TIME_UTC);
		for (int i = 0; i < NUM_THREADS - 1; ++i) {
			sem_init(&sem_divison[i], 0, 0);
			sem_init(&sem_worker[i], 0, 0);
		}



		sem_init(&sem_leader, 0, 0);

		//创建线程
		pthread_t handles[NUM_THREADS];
		threadparam_t param[NUM_THREADS];
		for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
			param[t_id].t_id = t_id;
			pthread_create(&handles[t_id], NULL, threadfunc3, (void*)&param[t_id]);
		}
		for (int i = 0; i < NUM_THREADS; i++)
			pthread_join(handles[i], NULL);
		for (int i = 0; i < NUM_THREADS - 1; ++i) {
			sem_destroy(&sem_divison[i]);
			sem_destroy(&sem_worker[i]);
		}
		sem_destroy(&sem_leader);
		timespec_get(&ets, TIME_UTC);
		time_sum += (ets.tv_nsec - sts.tv_nsec) * 0.000001 + (ets.tv_sec - sts.tv_sec) * 1000;
	}
	/*
	for (int i = 0; i < order; i++) {
		for (int j = 0; j <= order; j++) {
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
	*/
	cout << "PThread3 order: " << order << "   times:" << times << endl;
	cout << fixed << setprecision(2) << time_sum / times << endl;
	cout << endl;
}

//除法进入线程循环
void pthread4(int order, int times) {
	struct timespec sts, ets;
	double time_sum = 0.0;
	for (int time = 0; time < times; time++) {
		//初始化/恢复矩阵 不计入时间
		for (int i = 0; i < order; i++)
			for (int j = 0; j <= order; j++)
				matrix[i][j] = temp[i][j];
		timespec_get(&sts, TIME_UTC);
		for (int i = 0; i < NUM_THREADS - 1; ++i) {
			sem_init(&sem_divison[i], 0, 0);
			sem_init(&sem_worker[i], 0, 0);
		}

		sem_init(&sem_leader, 0, 0);

		//创建线程
		pthread_t handles[NUM_THREADS];
		threadparam_t param[NUM_THREADS];
		for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
			param[t_id].t_id = t_id;
			pthread_create(&handles[t_id], NULL, threadfunc4, (void*)&param[t_id]);
		}
		for (int i = 0; i < NUM_THREADS; i++)
			pthread_join(handles[i], NULL);
		for (int i = 0; i < NUM_THREADS - 1; ++i) {
			sem_destroy(&sem_divison[i]);
			sem_destroy(&sem_worker[i]);
		}
		sem_destroy(&sem_leader);
		timespec_get(&ets, TIME_UTC);
		time_sum += (ets.tv_nsec - sts.tv_nsec) * 0.000001 + (ets.tv_sec - sts.tv_sec) * 1000;
	}

	/*
	for (int i = 0; i < order; i++) {
		for (int j = 0; j <= order; j++) {
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
	*/
	cout << "PThread4 order: " << order << "   times:" << times << endl;
	cout << fixed << setprecision(2) << time_sum / times << endl;
	cout << endl;
}

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
	/*
	for (int i = 0; i < order; i++) {
		for (int j = 0; j <= order; j++) {
			cout << fixed << setprecision(2) << matrix[i][j] << " ";
		}
		cout << endl;
	}*/
	cout << "Traditon order: " << order << "   times:" << times << endl;
	cout << fixed << setprecision(2) << time_sum / times << endl;
	cout << endl;
	return;
}

int main() {
	unsigned seed = 47;
	srand(seed);
	for (int i = 0; i < ORDER; i++)
		for (int j = 0; j <= ORDER; j++)
			temp[i][j] = rand() % 10 + rand() % 10000 * 0.0001;
	int sample = 8;
	int samples[8][2] = { {100,100},{300,10},{500,10},{700,2},{1000,2},{1500,1},{2000,1},{3000,1} };
	for (int i = 0; i < sample; i++) {
		order = samples[i][0];
		tradition(samples[i][0], samples[i][1]);
	}
	for (int i = 0; i < sample; i++) {
		order = samples[i][0];
		pthread(samples[i][0], samples[i][1]);
	}
	for (int i = 0; i < sample; i++) {
		order = samples[i][0];
		pthread3(samples[i][0], samples[i][1]);
	}
	return 0;
}
