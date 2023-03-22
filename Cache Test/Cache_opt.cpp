#include<iostream>
#include <windows.h>
using namespace std;
const int SIZE_ = 20000;
const int RE = 10;
int sum[SIZE_];
int b[SIZE_][SIZE_], a[SIZE_];
int main() {
    cout << "矩阵大小: " << SIZE_ << " 重复次数: " << RE << endl;
    LARGE_INTEGER nFreq;
    LARGE_INTEGER nBeginTime;
    LARGE_INTEGER nEndTime;
    double time, time_op;
    QueryPerformanceFrequency(&nFreq);

    //Initializing matrices and vectors
    for (int i = 0; i < SIZE_; i++) {
        a[i] = i;
        for (int j = 0; j < SIZE_; j++) {
            b[i][j] = i + j;
        }
    }

    QueryPerformanceCounter(&nBeginTime);
    //Common algorithm:
    for (int k = 0; k < RE; k++) {
        for (int i = 0; i < SIZE_; i++) {
            sum[i] = 0;
            for (int j = 0; j < SIZE_; j++) {
                sum[i] += b[j][i] * a[j];
            }
        }
    }
    QueryPerformanceCounter(&nEndTime);
    time = (double)(nEndTime.QuadPart - nBeginTime.QuadPart) / (double)nFreq.QuadPart;
    cout << "General algorithm time consumption:   ";
    printf("%f\n", time);

    //Parallel algorithm(cache optimized)
    QueryPerformanceCounter(&nBeginTime);
    for (int k = 0; k < RE; k++) {
        for (int i = 0; i < SIZE_; i++)
            sum[i] = 0;
        for (int j = 0; j < SIZE_; j++) {
            for (int i = 0; i < SIZE_; i++) {
                sum[i] += b[j][i] * a[j];
            }
        }
    }
    QueryPerformanceCounter(&nEndTime);
    time_op = (double)(nEndTime.QuadPart - nBeginTime.QuadPart) / (double)nFreq.QuadPart;
    cout << "Optimized algorithm time consumption:   ";
    printf("%f\n", time_op);
    return 0;
}
