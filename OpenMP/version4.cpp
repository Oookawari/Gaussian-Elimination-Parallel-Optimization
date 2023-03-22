#include<omp.h>
#include<iostream>
#include<iomanip>
#include<random>
#include<time.h>
#include<stdio.h>
#include<pthread.h>
#include<semaphore.h>

using namespace std;
const int ORDER = 3000;
float matrix[ORDER][ORDER + 1];
float temp[ORDER][ORDER + 1];
int order = 100;
int SAMPLE = 0;
const int NUM_THREADS = 7;
typedef struct {
    int t_id;   //线程id
}threadparam_t;

sem_t sem_div;  //信号量定义
sem_t sem_elimination;
sem_t sem_divison[NUM_THREADS - 1];
sem_t sem_leader;
sem_t sem_worker[NUM_THREADS - 1];

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

void OpenMP(int order, int times) {
    struct timespec sts, ets;
    double time_sum = 0.0;
    for (int time = 0; time < times; time++) {

        //初始化/恢复矩阵 不计入时间
        for (int i = 0; i < order; i++)
            for (int j = 0; j <= order; j++)
                matrix[i][j] = temp[i][j];

        timespec_get(&sts, TIME_UTC);
        int i, j, k;
#pragma omp parallel num_threads(NUM_THREADS) shared(order) private(k,i,j)
        for (k = 0; k < order; k++) {
#pragma omp single
            {
                for (j = k + 1; j <= order; j++)
                    matrix[k][j] = matrix[k][j] / matrix[k][k];
                matrix[k][k] = 1.0;
            }
#pragma omp for
            for (i = k + 1; i < order; i++) {
                for (j = k + 1; j <= order; j++) {
                    matrix[i][j] -= matrix[i][k] * matrix[k][j];
                }
                matrix[i][k] = 0.0;
            }
#pragma omp barrier
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
    //cout << "OpenMP order: " << order << "   times:" << times << endl;
    cout << fixed << setprecision(2) << time_sum / times << endl;
    //cout << endl;
    return;
}

void OpenMP2(int order, int times) {
    struct timespec sts, ets;
    double time_sum = 0.0;
    for (int time = 0; time < times; time++) {

        //初始化/恢复矩阵 不计入时间
        for (int i = 0; i < order; i++)
            for (int j = 0; j <= order; j++)
                matrix[i][j] = temp[i][j];

        timespec_get(&sts, TIME_UTC);
        int i, j, k;
#pragma omp parallel num_threads(NUM_THREADS) shared(order) private(k,i,j)
        for (k = 0; k < order; k++) {
#pragma omp single
            {
                for (j = k + 1; j <= order; j++)
                    matrix[k][j] = matrix[k][j] / matrix[k][k];
                matrix[k][k] = 1.0;
            }
#pragma omp for simd
            for (i = k + 1; i < order; i++) {
                for (j = k + 1; j <= order; j++) {
                    matrix[i][j] -= matrix[i][k] * matrix[k][j];
                }
                matrix[i][k] = 0.0;
            }

#pragma omp barrier
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
    cout << "OpenMP2 order: " << order << "   times:" << times << endl;
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
    int samples[8][2] = { {100,1},{300,10},{500,10},{700,2},{1000,2},{1500,1},{2000,1},{3000,1} };
    order = 3000;
    for (int i = 0; i < sample; i++) {
        order = samples[i][0];
        // tradition(samples[i][0], samples[i][1]);
    }
    for (int i = 0; i < sample; i++) {
        order = samples[i][0];
        OpenMP(samples[i][0], samples[i][1]);
    }
    for (int i = 0; i < sample; i++) {
        order = samples[i][0];
        OpenMP2(samples[i][0], samples[i][1]);
    }
    return 0;
}
